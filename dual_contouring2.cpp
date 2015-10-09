#include "stdafx.h"
#pragma hdrstop
#include "QEF.h"
#include "dual_contouring2.h"

namespace DualContouring
{
void InitNode( Node2& _node )
{
	mxUINT_LOOP_i(8) {
		_node.kids[i] = NIL_NODE;
	}
}
// empty nodes (i.e. internal nodes with no children) should be collapsed
bool IsBadNode( const Node2& _node )
{
	return
		_node.kids[0] == NIL_NODE && _node.kids[1] == NIL_NODE && _node.kids[2] == NIL_NODE && _node.kids[3] == NIL_NODE &&
		_node.kids[4] == NIL_NODE && _node.kids[5] == NIL_NODE && _node.kids[6] == NIL_NODE && _node.kids[7] == NIL_NODE;
}

Octree_DC2::Octree_DC2()
{
}
Octree_DC2::~Octree_DC2()
{
	m_nodes.Empty();
	m_leaves.Empty();
	m_QEFs.Empty();
}
void Octree_DC2::Build(
					   const AVolume* _volume,
					   QEF_Solver* _qef_solver,
					   OctStats &_octree_stats,
					   const Options& options
					   )
{
	mxPROFILE_SCOPE("Build Octree");

	DBGOUT("Building octree (radius: %f, min.cell: %f, QEF threshold: %f)",
		options.radius, options.minSubdiv, options.qef_threshold);

	const UINT32 startTimeMSec = mxGetTimeInMilliseconds();

	m_nodes.Empty();
	m_leaves.Empty();
	m_QEFs.Empty();

	OctCubeF worldBounds;
	worldBounds.x = 0;
	worldBounds.y = 0;
	worldBounds.z = 0;
	worldBounds.radius = options.radius;

	NodeID rootNodeID = 0;

	rootNodeID = BuildOctreeRecursive( worldBounds, 0, _qef_solver, _volume, options, _octree_stats );

#if USE_QEF_SIMPLIFICATION
	rootNodeID = Simplify( rootNodeID, worldBounds, 0, _qef_solver, _volume, options, _octree_stats );
#endif

	_octree_stats.bytesAllocated = m_nodes.GetDataSize() + m_leaves.GetDataSize();

	const UINT32 currentTimeMSec = mxGetTimeInMilliseconds();
	_octree_stats.construction_time_milliseconds = currentTimeMSec - startTimeMSec;
}

// may return NULL
NodeID Octree_DC2::TryCreateLeaf(
					 const OctCubeF& _bounds,
					 const UINT32 _treeLevel,
					 QEF_Solver* _qef_solver,
					 const AVolume* _volume,
					 float _error_threshold,
					 OctStats &stats
					 )
{
	Float3 corners[8];
	OCT_GetCorners( _bounds, corners );

	// Find corner signs and determine zero crossing edges.
	int cornerSigns = 0;
	for( int i = 0; i < 8; i++ )
	{
		const float distance = _volume->GetDistanceAt(corners[i]);
		if( AVolume::IsSolid( distance ) ) {
			cornerSigns |= (1 << i);
		}
	}
	if( cornerSigns == 0xFF || cornerSigns == 0 )
	{
		// voxel is fully inside or outside the volume
		stats.numEmptyLeaves++;
		return NIL_NODE;
	}

	// Calculate the so-called 'Hermite data':
	// intersection points together with their normals
	// (of the contour with the edges of the cube).

	QEF_Solver::Input solverInput;
	solverInput.bounds.min_point = corners[0];
	solverInput.bounds.max_point = corners[CHILD_MASK_X|CHILD_MASK_Y|CHILD_MASK_Z];

	// Calculate intersection points on zero crossing edges.
	const UINT8* edges = CUBE_GetEdgeIndices();
	for( int iEdge = 0; iEdge < 12; iEdge++ )
	{
		int	iPointA = edges[ iEdge*2 + 0 ];
		int	iPointB = edges[ iEdge*2 + 1 ];

		int signA = (cornerSigns >> iPointA) & 1;
		int signB = (cornerSigns >> iPointB) & 1;
		if( signA == signB ) {
			continue;
		}

		const Float3& pointA = corners[ iPointA ];
		const Float3& pointB = corners[ iPointB ];

		Float3 intersectionPoint;
		Float3 intersectionNormal;

		if( _volume->IntersectsLine( pointA, pointB, intersectionPoint, intersectionNormal ) )
		{
			solverInput.positions[ solverInput.numPoints ] = intersectionPoint;
			solverInput.normals[ solverInput.numPoints ] = intersectionNormal;
			solverInput.numPoints++;
		}
	}

	mxASSERT2(solverInput.numPoints > 0, "No intersection points, but differing signs at corners - should happen!");
	mxASSERT(solverInput.numPoints <= 8);

	stats.maxActiveEdges = largest(stats.maxActiveEdges, solverInput.numPoints);

	QEF_Solver::Output solverOutput;
	_qef_solver->Solve( solverInput, solverOutput );

	stats.max_QEF_error = largest( stats.max_QEF_error, solverOutput.error );

	if( solverOutput.error < _error_threshold ) {
		//DBGOUT("solverOutput.error: %f", solverOutput.error);
		return NIL_NODE;
	}

	Float3 centerP, centerN;
	centerP = solverOutput.position;
	centerN = Float3_Normalized(_volume->SampleAt(centerP).normal);

	const UINT16 leafIndex = m_leaves.Alloc();
	Leaf2& leaf = m_leaves[ leafIndex ];
	QuantizePosition( _bounds, centerP, leaf.xyz );
	QuantizeNormal( centerN, leaf.N );
	leaf.signs = cornerSigns;
#if USE_QEF_SIMPLIFICATION
	mxASSERT(solverOutput.qef.numPoints > 0);
	leaf.qefID = m_QEFs.Alloc();
	QEF_t& qef = m_QEFs[ leaf.qefID ];
	(svd::QefData&)qef = solverOutput.qef;
#else
	leaf.qefID = (UINT16)~0;
#endif
	stats.numLeafNodes++;

	stats.maxTreeDepth = largest(stats.maxTreeDepth, _treeLevel);

	return MAKE_LEAF_ID(leafIndex);
}
// may return NULL
NodeID Octree_DC2::TryCreateLeaf2(
					 const OctCubeF& _bounds,
					 const UINT32 _treeLevel,
					 QEF_Solver* _qef_solver,
					 const AVolume* _volume,
					 float _error_threshold,
					 OctStats &stats
					 )
{
	Float3 corners[8];
	OCT_GetCorners( _bounds, corners );

	// Calculate the so-called 'Hermite data':
	// intersection points together with their normals
	// (of the contour with the edges of the cube).

	QEF_Solver::Input solverInput;
	solverInput.bounds.min_point = corners[0];
	solverInput.bounds.max_point = corners[CHILD_MASK_X|CHILD_MASK_Y|CHILD_MASK_Z];

	// Calculate intersection points on zero crossing edges.
	const UINT8* edges = CUBE_GetEdgeIndices();
	for( int iEdge = 0; iEdge < 12; iEdge++ )
	{
		int	iPointA = edges[ iEdge*2 + 0 ];
		int	iPointB = edges[ iEdge*2 + 1 ];
		const Float3& pointA = corners[ iPointA ];
		const Float3& pointB = corners[ iPointB ];

		Float3 intersectionPoint;
		Float3 intersectionNormal;

		if( _volume->IntersectsLine( pointA, pointB, intersectionPoint, intersectionNormal ) )
		{
			solverInput.positions[ solverInput.numPoints ] = intersectionPoint;
			solverInput.normals[ solverInput.numPoints ] = intersectionNormal;
			solverInput.numPoints++;
		}
	}

	if( solverInput.numPoints == 0 )
	{
		// voxel is fully inside or outside the volume
		stats.numEmptyLeaves++;
		return NIL_NODE;
	}

	stats.maxActiveEdges = largest(stats.maxActiveEdges, solverInput.numPoints);

	QEF_Solver::Output solverOutput;
	_qef_solver->Solve( solverInput, solverOutput );

	stats.max_QEF_error = largest( stats.max_QEF_error, solverOutput.error );

	if( solverOutput.error < _error_threshold ) {
		//DBGOUT("solverOutput.error: %f", solverOutput.error);
		return NIL_NODE;
	}

	Float3 centerP, centerN;
	centerP = solverOutput.position;
	centerN = Float3_Normalized(_volume->SampleAt(centerP).normal);

	const UINT16 leafIndex = m_leaves.Alloc();
	Leaf2& leaf = m_leaves[ leafIndex ];
	QuantizePosition( _bounds, centerP, leaf.xyz );
	QuantizeNormal( centerN, leaf.N );

	int cornerSigns = 0;
	for( int i = 0; i < 8; i++ )
	{
		const float distance = _volume->GetDistanceAt(corners[i]);
		if( AVolume::IsSolid( distance ) ) {
			cornerSigns |= (1 << i);
		}
	}
	leaf.signs = cornerSigns;
#if USE_QEF_SIMPLIFICATION
	mxASSERT(solverOutput.qef.numPoints > 0);
	leaf.qefID = m_QEFs.Alloc();
	QEF_t& qef = m_QEFs[ leaf.qefID ];
	(svd::QefData&)qef = solverOutput.qef;
#else
	leaf.qefID = (UINT16)~0;
#endif
	stats.numLeafNodes++;

	stats.maxTreeDepth = largest(stats.maxTreeDepth, _treeLevel);

	return MAKE_LEAF_ID(leafIndex);
}
NodeID Octree_DC2::BuildOctreeRecursive(
										const OctCubeF& _bounds,
										const UINT32 _treeLevel,
										QEF_Solver* _qef_solver,
										const AVolume* _volume,
										const Options& options,
										OctStats &stats
										)
{
	mxPROFILE_SCOPE("BuildOctreeRecursive");

	const bool isLeaf = (_treeLevel >= options.maxDepth)
				|| (_bounds.radius <= options.minSubdiv);

	if( isLeaf )
	{
		const float fake_threshold = -1.0;	// don't throw away leaves with error less than this
		return TryCreateLeaf( _bounds, _treeLevel, _qef_solver, _volume, fake_threshold, stats );
	}
	else
	{
//#if !USE_QEF_SIMPLIFICATION
//		if( _treeLevel > 0 )// don't simplify octree down to a single node-point
//		{
//			NodeID leafID = TryCreateLeaf2( _bounds, _treeLevel, _qef_solver, _volume, options.build_threshold, stats );
//			if( leafID != NIL_NODE )
//			{
//				return leafID;
//			}
//		}
//#endif

		OctCubeF octants[8];
		GetChildOctants(_bounds,octants);
		const UINT32 nextLevel = _treeLevel + 1;

		NodeID nodeID = m_nodes.Alloc();
		for( int i = 0; i < 8; i++ )
		{
			NodeID childID = BuildOctreeRecursive( octants[i], nextLevel, _qef_solver, _volume, options, stats );
			Node2& node = m_nodes[ nodeID ];
			node.kids[i] = childID;
		}

		// collapse empty nodes
		if( IsBadNode( m_nodes[ nodeID ] ) )
		{
			stats.numBadNodes++;
			m_nodes.Free( nodeID );
			nodeID = NIL_NODE;
		}
		else
		{
			stats.numInternalNodes++;
		}		

		return nodeID;
	}
}

#if USE_QEF_SIMPLIFICATION
NodeID Octree_DC2::Simplify(
							const NodeID _nodeIndex,
							const OctCubeF& _bounds,
							const UINT32 _treeLevel,
							QEF_Solver* _qef_solver,
							const AVolume* _volume,
							const Options& options,
							OctStats &stats
							)
{
	NodeID result = _nodeIndex;

	//// don't simplify octree down to a single node-point
	//if( _treeLevel == 1 ) {
	//	return _nodeIndex;
	//}

	if( _nodeIndex != NIL_NODE && !IS_LEAF_ID( _nodeIndex ) )
	{
		OctCubeF octants[8];
		GetChildOctants( _bounds, octants );

		bool isCollapsible = true;

		int signs[8] = { -1, -1, -1, -1, -1, -1, -1, -1 };
		int midsign = -1;


		Float3 corners[8];
		OCT_GetCorners( _bounds, corners );


		svd::QefSolver	solver;

		int edgeCount = 0;


		Float3 averageNormal = Float3_Zero();

		const int nextLevel = _treeLevel + 1;

		for( int i = 0; i < 8; i++ )
		{
			Node2& node = m_nodes[_nodeIndex];

			if( node.kids[i] != NIL_NODE )
			{
				NodeID simplified = Simplify(
					node.kids[i],
					octants[i],
					nextLevel,
					_qef_solver,
					_volume,
					options,
					stats
				);
				mxASSERT( simplified != NIL_NODE );
				if( node.kids[i] != simplified ) {
					FreeNode( node.kids[i] );
					node.kids[i] = simplified;
					stats.nodes_collapsed++;
				}
				if( IS_LEAF_ID(simplified) )
				{
					const Leaf2& leaf = m_leaves[ GET_ID(simplified) ];
					const QEF_t& qef = m_QEFs[ leaf.qefID ];
					solver.add( qef );
		//midsign = (leaf.signs >> (7 - i)) & 1;
		//midsign = (leaf.signs >> i) & 1;
					signs[i] = (leaf.signs >> i) & 1;
					averageNormal += DequantizeNormal( leaf.N );
					edgeCount++;
				}
				else
				{
					// at least one child is an internal node, can't collapse
					isCollapsible = false;
				}
			}
		}

		midsign = _volume->GetDistanceAt(_bounds.XYZ()) >= 0.0f;

Node2& node = m_nodes[_nodeIndex];

		if( isCollapsible )
		{
			const float QEF_ERROR = 1e-6f;
			const int QEF_SWEEPS = 4;

			svd::Vec3 qefPosition;
			solver.solve(qefPosition, QEF_ERROR, QEF_SWEEPS, QEF_ERROR);
			float error = solver.getError();

			// collapse the node if the residual is less than the given tolerance
			if( error <= options.qef_threshold )
			{
				Float3 position = Float3_Set(qefPosition.x,qefPosition.y,qefPosition.z);

				AABB24 aabb;
				aabb.min_point = corners[0];
				aabb.max_point = corners[CHILD_MASK_X|CHILD_MASK_Y|CHILD_MASK_Z];

				if( !AABB_ContainsPoint( aabb, position ) )
				{
					const svd::Vec3& mp = solver.getMassPoint();
					position = Float3_Set(mp.x, mp.y, mp.z);
				}

				const UINT16 leafIndex = m_leaves.Alloc();
				Leaf2& newLeaf = m_leaves[ leafIndex ];
				result = MAKE_LEAF_ID( leafIndex );

				for (int i = 0; i < 8; i++)
				{
					if (signs[i] == -1)
					{
						// Undetermined, use center sign instead
						newLeaf.signs |= (midsign << i);
					}
					else
					{
						newLeaf.signs |= (signs[i] << i);
					}
				}

				QuantizePosition( _bounds, position, newLeaf.xyz );
				QuantizeNormal( Float3_Normalized(averageNormal), newLeaf.N );

				newLeaf.qefID = m_QEFs.Alloc();
				QEF_t& newQEF = m_QEFs[ newLeaf.qefID ];

				(svd::QefData&)newQEF = solver.getData();
			}
		}
	}
	return result;
}
#endif
void Octree_DC2::FreeNode( NodeID nodeID )
{
	const int nodeIndex = GET_ID( nodeID );
	if( IS_LEAF_ID( nodeID ) )
	{
		Leaf2& leaf = m_leaves[ nodeIndex ];
		if( leaf.qefID != (UINT16)~0 ) {
			m_QEFs.Free( leaf.qefID );
		}
		m_leaves.Free( nodeIndex );
	}
	else
	{
		m_nodes.Free( nodeIndex );
	}
}

static void CollectVertices(
							const Octree_DC2& _tree,
							const OctCubeF& _bounds,
							const UINT32 _treeLevel,
							const Options& _options,
							const NodeID _parentID,
							AMeshBuilder & _mesh,
							MeshCtx & _meshCtx,
							OctStats & _stats
							)
{
	mxPROFILE_SCOPE("CollectVertices");

	const bool isLeaf = IS_LEAF_ID(_parentID);
	const UINT32 nodeIndex = GET_ID(_parentID);

	if(isLeaf)
	{
		const Leaf2& leaf = _tree.m_leaves[nodeIndex];

		DrawVertex	vertex;
		DequantizePosition( _bounds, leaf.xyz, &vertex.xyz );
		vertex.N.v = (int&)leaf.N;

		const int vertexIndex = _mesh.AddVertex( vertex );
		_meshCtx.vertexIDs[nodeIndex] = vertexIndex;
	}
	else//if(!isLeaf)
	{
		OctCubeF octants[8];
		GetChildOctants(_bounds,octants);
		const UINT32 nextLevel = _treeLevel + 1;
		mxASSERT(nextLevel <= _options.maxDepth);

		const Node2& node = _tree.m_nodes[nodeIndex];
		for( int i = 0; i < 8; i++ )
		{
			if( node.kids[i] != NIL_NODE ) {
				CollectVertices( _tree, octants[i], nextLevel, _options, node.kids[i], _mesh, _meshCtx, _stats );
			}
		}
	}
}

ERet Octree_DC2::Triangulate(
							 AMeshBuilder & mesh,
							 const Options& options,
							 OctStats &stats
							 )
{
	mxPROFILE_SCOPE("Triangulate");

	mxDO(mesh.Begin());

	const UINT32 startTimeMSec = mxGetTimeInMilliseconds();

	OctCubeF worldBounds;
	worldBounds.x = 0;
	worldBounds.y = 0;
	worldBounds.z = 0;
	worldBounds.radius = options.radius;

	const UINT32 zeroLevel = 0;
	const NodeID rootNodeID = 0;

	// allocate temporary storage
	const UINT32 maxVertices = m_leaves.Num();
	mxASSERT( maxVertices < MAX_UINT16 );

	MeshCtx	meshCtx;

	ScopedStackAlloc	tempAlloc( gCore.frameAlloc );
	meshCtx.vertexIDs = tempAlloc.AllocMany< UINT16 >( maxVertices );

	// create vertices
	CollectVertices( *this, worldBounds, zeroLevel, options, rootNodeID, mesh, meshCtx, stats );

	// create quadrilaterals
	NodeCtx rootCtx = { 0, 0 };
	ProcessNode( rootCtx, mesh, meshCtx );

	int verts, tris;
	mxDO(mesh.End( verts, tris ));


	const UINT32 currentTimeMSec = mxGetTimeInMilliseconds();
	const UINT32 elapsedTimeMSec = currentTimeMSec - startTimeMSec;

	stats.numPolygons += tris;
	stats.contouring_time_milliseconds = elapsedTimeMSec;

	DBGOUT("Octree_DC2::Triangulate(): %d verts, %d tris in %u msec\n", verts, tris, elapsedTimeMSec);

	return ALL_OK;
}

// recursively enumerates each edge of the octree and emits quads at minimal edges when needed
/*
         .-----.-----.		
        /  6  /  7  /|		  Z
       .-----.-----.7|		  |  /Y
      /  4  /  5  /|/|		  | /
     .-----+-----|5|3.		  |/_____X
     |  4  |  5  |/|/		(0,0,0)
     |-----|-----|1/  
     |  0  |  1  |/   
     .-----+-----.    
*/
void Octree_DC2::ProcessNode(
							 const NodeCtx& n,
							 AMeshBuilder& mesh, MeshCtx & ctx
							 )
{
	mxPROFILE_SCOPE("ProcessNode");
	if( n.id != NIL_NODE && !IS_LEAF_ID(n.id) )
	{
		//const Node2& node = m_nodes[ n.id ];
		const NodeCtx c0 = CreateChildCtx( n, 0 );
		const NodeCtx c1 = CreateChildCtx( n, 1 );
		const NodeCtx c2 = CreateChildCtx( n, 2 );
		const NodeCtx c3 = CreateChildCtx( n, 3 );
		const NodeCtx c4 = CreateChildCtx( n, 4 );
		const NodeCtx c5 = CreateChildCtx( n, 5 );
		const NodeCtx c6 = CreateChildCtx( n, 6 );
		const NodeCtx c7 = CreateChildCtx( n, 7 );
		// Call ProcessNode() for each child (8 cell calls).
		{
			ProcessNode( c0, mesh, ctx );
			ProcessNode( c1, mesh, ctx );
			ProcessNode( c2, mesh, ctx );
			ProcessNode( c3, mesh, ctx );
			ProcessNode( c4, mesh, ctx );
			ProcessNode( c5, mesh, ctx );
			ProcessNode( c6, mesh, ctx );
			ProcessNode( c7, mesh, ctx );
		}
		// Call ProcessFace() for each pair of children that share a face (12 face calls).
		{
			// call ProcessFace_X() for children sharing a face in the ZY plane ('vertical' plane with normal pointing left to right)
			ProcessFaces_X( c0, c1, c2, c3, c4, c5, c6, c7, mesh, ctx );	// 4 ProcessFace_X() calls
			// call ProcessFace_Y() for children sharing a face in the XZ plane ('front' plane with normal pointing front to back)
			ProcessFaces_Y( c0, c1, c2, c3, c4, c5, c6, c7, mesh, ctx );	// 4 ProcessFace_Y() calls
			// call ProcessFace_Z() for children sharing a face in the XY plane ('ground' plane with normal pointing upwards)
			ProcessFaces_Z( c0, c1, c2, c3, c4, c5, c6, c7, mesh, ctx );	// 4 ProcessFace_Z() calls
		}
		// Call ProcessEdge() for each of the 6 interior edges (6 edge calls).
		{
			// call ProcessEdge_X for children sharing an edge along the X-axis
			ProcessEdges_X_Axis( c0, c1, c2, c3, c4, c5, c6, c7, mesh, ctx );	// 2 ProcessEdge_X() calls
			// call ProcessEdge_Y for children sharing an edge along the Y-axis
			ProcessEdges_Y_Axis( c0, c1, c2, c3, c4, c5, c6, c7, mesh, ctx );	// 2 ProcessEdge_Y() calls
			// call ProcessEdge_Z for children sharing an edge along the Z-axis
			ProcessEdges_Z_Axis( c0, c1, c2, c3, c4, c5, c6, c7, mesh, ctx );	// 2 ProcessEdge_Z() calls
		}
	}
}
/*
	Processes the two nodes sharing a (vertical) face along the X axis:

  ..Left node  |  Right node..
	  ---.-----.-----.---		
        /  7  /  6  /|	      Z
     --.-----.-----.6|--      |  /Y
      /  5  /  4  /|/|	      | /
  ---.-----+-----|4|2.---     |/_____X*
     |  5  |  4  |/|/		(0,0,0)
     |-----+-----|0/  
     |  1  |  0  |/   
  ---.-----+-----.---    
*/
void Octree_DC2::ProcessFace_X(
	const NodeCtx& n1,	//<= 'left' node
	const NodeCtx& n2,	//<= 'right' node
	AMeshBuilder& mesh, MeshCtx & ctx
	)
{
	// When at least one of the two nodes sharing a face is subdivided,
	// we get four sub-faces and four edges meeting at a single vertex
	// at the center of the face.
	const bool n1Leaf = IS_LEAF_ID(n1.id);
	const bool n2Leaf = IS_LEAF_ID(n2.id);
	// If all nodes are leaves, or one is empty, bail out.
	if( n1.id != NIL_NODE && n2.id != NIL_NODE && (!n1Leaf || !n2Leaf) )
	{
		// left
		const NodeCtx c0 = n1Leaf ? n1 : CreateChildCtx( n1, 1 );
        const NodeCtx c2 = n1Leaf ? n1 : CreateChildCtx( n1, 3 );
		const NodeCtx c6 = n1Leaf ? n1 : CreateChildCtx( n1, 7 );
		const NodeCtx c4 = n1Leaf ? n1 : CreateChildCtx( n1, 5 );
		// right
		const NodeCtx c1 = n2Leaf ? n2 : CreateChildCtx( n2, 0 );
        const NodeCtx c3 = n2Leaf ? n2 : CreateChildCtx( n2, 2 );
        const NodeCtx c7 = n2Leaf ? n2 : CreateChildCtx( n2, 6 );
		const NodeCtx c5 = n2Leaf ? n2 : CreateChildCtx( n2, 4 );        

		// call ProcessFace_YZ for children sharing a face in the ZY plane
		ProcessFaces_X( c0, c1, c2, c3, c4, c5, c6, c7, mesh, ctx );	// 4 ProcessFace_X() calls
		// call ProcessEdge_Y for nodes sharing an edge along the Y-axis
		ProcessEdges_Y_Axis( c0, c1, c2, c3, c4, c5, c6, c7, mesh, ctx );	// 2 ProcessEdge_Y() calls
		// call ProcessEdge_Z for children sharing an edge along the Z-axis
		ProcessEdges_Z_Axis( c0, c1, c2, c3, c4, c5, c6, c7, mesh, ctx );	// 2 ProcessEdge_Z() calls
	}
}
/*
	Processes the two nodes sharing a (vertical) face along the Y axis:

            Back node
	     .-----.-----.   		
        /  4  /  5  /|	      Z
       .-----.-----.5|        |  /Y*
      /  6  /  7  /|/|	      | /
     .-----+-----|7|1.        |/_____X
     |  6  |  7  |/|/		(0,0,0)
     |-----+-----|3/  
     |  2  |  3  |/   
     .-----+-----.          
      Front node
*/
void Octree_DC2::ProcessFace_Y(
	const NodeCtx& n1,	//<= 'front' node
	const NodeCtx& n2,	//<= 'back' node
	AMeshBuilder& mesh, MeshCtx & ctx
	)
{
	const bool n1Leaf = IS_LEAF_ID(n1.id);
	const bool n2Leaf = IS_LEAF_ID(n2.id);
	// If all nodes are leaves, or one is empty, bail out.
	if( n1.id != NIL_NODE && n2.id != NIL_NODE && (!n1Leaf || !n2Leaf) )
	{
		// front
		const NodeCtx c0 = n1Leaf ? n1 : CreateChildCtx( n1, 2 );
        const NodeCtx c4 = n1Leaf ? n1 : CreateChildCtx( n1, 6 );
        const NodeCtx c5 = n1Leaf ? n1 : CreateChildCtx( n1, 7 );
		const NodeCtx c1 = n1Leaf ? n1 : CreateChildCtx( n1, 3 );
		// back
        const NodeCtx c2 = n2Leaf ? n2 : CreateChildCtx( n2, 0 );
        const NodeCtx c6 = n2Leaf ? n2 : CreateChildCtx( n2, 4 );
        const NodeCtx c7 = n2Leaf ? n2 : CreateChildCtx( n2, 5 );
		const NodeCtx c3 = n2Leaf ? n2 : CreateChildCtx( n2, 1 );

		// call ProcessFace_XZ for nodes sharing sub-faces
		ProcessFaces_Y( c0, c1, c2, c3, c4, c5, c6, c7, mesh, ctx );	// 4 ProcessFace_Y() calls
		// call ProcessEdge_X for nodes sharing an edge along the X-axis
		ProcessEdges_X_Axis( c0, c1, c2, c3, c4, c5, c6, c7, mesh, ctx );	// 2 ProcessEdge_X() calls
		// call ProcessEdge_Z for children sharing an edge along the Z-axis
		ProcessEdges_Z_Axis( c0, c1, c2, c3, c4, c5, c6, c7, mesh, ctx );	// 2 ProcessEdge_Z() calls
	}
}
/*
	Processes the two nodes sharing a (horizontal) face along the Z axis:

            Upper node
	     .-----.-----.   		
        /  2  /  3  /|	      Z*
       .-----.-----.3|        |  /Y
      /  0  /  1  /|/|	      | /
     .-----+-----|1|7.        |/_____X
     |  0  |  1  |/|/		(0,0,0)
     |-----+-----|5/  
     |  4  |  5  |/   
     .-----+-----.          
      Lower node
*/
void Octree_DC2::ProcessFace_Z(
	const NodeCtx& n1,	//<= lower node
	const NodeCtx& n2,	//<= upper node
	AMeshBuilder& mesh, MeshCtx & ctx
	)
{
	const bool n1Leaf = IS_LEAF_ID(n1.id);
	const bool n2Leaf = IS_LEAF_ID(n2.id);
	// If all nodes are leaves, or one is empty, bail out.
	if( n1.id != NIL_NODE && n2.id != NIL_NODE && (!n1Leaf || !n2Leaf) )
	{
		// bottom part
		const NodeCtx c0 = n1Leaf ? n1 : CreateChildCtx( n1, 4 );
        const NodeCtx c1 = n1Leaf ? n1 : CreateChildCtx( n1, 5 );
        const NodeCtx c2 = n1Leaf ? n1 : CreateChildCtx( n1, 6 );
        const NodeCtx c3 = n1Leaf ? n1 : CreateChildCtx( n1, 7 );
		// upper part
        const NodeCtx c4 = n2Leaf ? n2 : CreateChildCtx( n2, 0 );
        const NodeCtx c5 = n2Leaf ? n2 : CreateChildCtx( n2, 1 );
        const NodeCtx c6 = n2Leaf ? n2 : CreateChildCtx( n2, 2 );
        const NodeCtx c7 = n2Leaf ? n2 : CreateChildCtx( n2, 3 );

		// call ProcessFace_Z() for nodes sharing sub-faces on the XY plane
		ProcessFaces_Z( c0, c1, c2, c3, c4, c5, c6, c7, mesh, ctx );
		// call ProcessEdge_X for nodes sharing an edge along the X-axis
		ProcessEdges_X_Axis( c0, c1, c2, c3, c4, c5, c6, c7, mesh, ctx );	// 2 ProcessEdge_X() calls
		// call ProcessEdge_Y for nodes sharing an edge along the Y-axis
		ProcessEdges_Y_Axis( c0, c1, c2, c3, c4, c5, c6, c7, mesh, ctx );	// 2 ProcessEdge_Y() calls
	}
}
/*
	Processes the 4 subtrees sharing the edge along the X axis:

         .-----.    Z
        /  2  /|    |  Y
       /-----/2|    | /
      /  3  /|/|    |/
     .-----|3|1|    O-------X*
     |  3  |/|/ 
     |-----|0/  
     |  0  |/   
     .-----.
*/
void Octree_DC2::ProcessEdge_X(
	// 'right-hand rule' order, around X axis, starting from the lowest octant index
	const NodeCtx& n0,const NodeCtx& n1,const NodeCtx& n2,const NodeCtx& n3,
	AMeshBuilder& mesh, MeshCtx & ctx
	)
{
	// If one of the nodes is null, bail out.
	if( n0.id != NIL_NODE && n1.id != NIL_NODE && n2.id != NIL_NODE && n3.id != NIL_NODE )
	{
		const bool n0Leaf = IS_LEAF_ID(n0.id);
		const bool n1Leaf = IS_LEAF_ID(n1.id);
		const bool n2Leaf = IS_LEAF_ID(n2.id);
		const bool n3Leaf = IS_LEAF_ID(n3.id);
		// If all cubes are leaves, stop recursion and emit a quad.
		if( n0Leaf && n1Leaf && n2Leaf && n3Leaf )
		{
			// All nodes are surface leaves - determine the minimal cell/edge (containing no smaller edge).
			const Leaf2& leaf0 = m_leaves[ GET_ID(n0.id) ];
			const Leaf2& leaf1 = m_leaves[ GET_ID(n1.id) ];
			const Leaf2& leaf2 = m_leaves[ GET_ID(n2.id) ];
			const Leaf2& leaf3 = m_leaves[ GET_ID(n3.id) ];

			// Determine the cell with the greatest subdivision
			// and test the corner signs of its edge to determine
			// if a polygon should be generated.

			const NodeCtx* mc;
			mc = (n0.depth >= n1.depth) ? &n0 : &n1;
			mc = (mc->depth >= n2.depth) ? mc : &n2;
			mc = (mc->depth >= n3.depth) ? mc : &n3;

			// left -> right
			int sign1, sign2;	// 0 = outside, 1 = inside surface
			if(mc == &n0) {
				sign1 = (leaf0.signs >> 6) & 1;
				sign2 = (leaf0.signs >> 7) & 1;
			} else if(mc == &n1) {
				sign1 = (leaf1.signs >> 4) & 1;
				sign2 = (leaf1.signs >> 5) & 1;
			} else if(mc == &n2) {
				sign1 = (leaf2.signs >> 0) & 1;
				sign2 = (leaf2.signs >> 1) & 1;
			} else {
				sign1 = (leaf3.signs >> 2) & 1;
				sign2 = (leaf3.signs >> 3) & 1;
			}

			if(sign1 > sign2) {
				/* The first corner is inside the surface, the second is outside */
				EmitQuad( n0, n1, n2, n3, mesh, ctx );
			} else if( sign1 < sign2 ) {
				/* The second corner is inside the surface, the first is outside */
				EmitQuad( n3, n2, n1, n0, mesh, ctx );
			}
		}
		else//if( !n0Leaf || !n1Leaf || !n2Leaf || !n3Leaf )
		{
			const NodeCtx c0 = n0Leaf ? n0 : CreateChildCtx( n0, 6 );
			const NodeCtx c1 = n0Leaf ? n0 : CreateChildCtx( n0, 7 );
			const NodeCtx c2 = n1Leaf ? n1 : CreateChildCtx( n1, 4 );
			const NodeCtx c3 = n1Leaf ? n1 : CreateChildCtx( n1, 5 );

			const NodeCtx c4 = n3Leaf ? n3 : CreateChildCtx( n3, 2 );
			const NodeCtx c5 = n3Leaf ? n3 : CreateChildCtx( n3, 3 );
			const NodeCtx c6 = n2Leaf ? n2 : CreateChildCtx( n2, 0 );
			const NodeCtx c7 = n2Leaf ? n2 : CreateChildCtx( n2, 1 );

			// ...call ProcessEdge_X for nodes sharing smaller edges:
			ProcessEdges_X_Axis( c0, c1, c2, c3, c4, c5, c6, c7, mesh, ctx );	// 2 ProcessEdge_X() calls
		}
	}
}

/*
	Processes the 4 subtrees sharing the edge along the Y axis:

       .-----.-----.     Z
      /  1  /  2  /|     |  Y*
     .-----+-----|2|     | /
     |  1  |  2  |/|     |/
     |-----X-----|3/     O-------X
     |  0  |  3  |/     
     .-----.-----.   
*/
void Octree_DC2::ProcessEdge_Y(
	const NodeCtx& n0,const NodeCtx& n1,const NodeCtx& n2,const NodeCtx& n3,
	AMeshBuilder& mesh, MeshCtx & ctx
	)
{
	// If one of the nodes is null, bail out.
	if( n0.id != NIL_NODE && n1.id != NIL_NODE && n2.id != NIL_NODE && n3.id != NIL_NODE )
	{
		const bool n0Leaf = IS_LEAF_ID(n0.id);
		const bool n1Leaf = IS_LEAF_ID(n1.id);
		const bool n2Leaf = IS_LEAF_ID(n2.id);
		const bool n3Leaf = IS_LEAF_ID(n3.id);
		// If all cubes are leaves, stop recursion and emit a quad.
		if( n0Leaf && n1Leaf && n2Leaf && n3Leaf )
		{
			// All nodes are surface leaves - determine the minimal cell/edge (containing no smaller edge).
			const Leaf2& leaf0 = m_leaves[ GET_ID(n0.id) ];
			const Leaf2& leaf1 = m_leaves[ GET_ID(n1.id) ];
			const Leaf2& leaf2 = m_leaves[ GET_ID(n2.id) ];
			const Leaf2& leaf3 = m_leaves[ GET_ID(n3.id) ];

			const NodeCtx* mc;
			mc = (n0.depth >= n1.depth) ? &n0 : &n1;
			mc = (mc->depth >= n2.depth) ? mc : &n2;
			mc = (mc->depth >= n3.depth) ? mc : &n3;

			// front -> back
			int sign1, sign2;	// 0 = outside, 1 = inside surface
			if(mc == &n0) {
				sign1 = (leaf0.signs >> 5) & 1;
				sign2 = (leaf0.signs >> 7) & 1;
			} else if(mc == &n1) {
				sign1 = (leaf1.signs >> 1) & 1;
				sign2 = (leaf1.signs >> 3) & 1;
			} else if(mc == &n2) {
				sign1 = (leaf2.signs >> 0) & 1;
				sign2 = (leaf2.signs >> 2) & 1;
			} else {
				sign1 = (leaf3.signs >> 4) & 1;
				sign2 = (leaf3.signs >> 6) & 1;
			}

			if( sign1 > sign2 ) {
				// The first corner is inside the surface, the second is outside
				EmitQuad( n0, n1, n2, n3, mesh, ctx );
			} else if( sign1 < sign2 ) {
				// The second corner is inside the surface, the first is outside
				EmitQuad( n3, n2, n1, n0, mesh, ctx );
			}
		}
		else
		{
			const NodeCtx c0 = n0Leaf ? n0 : CreateChildCtx( n0, 5 );
			const NodeCtx c1 = n3Leaf ? n3 : CreateChildCtx( n3, 4 );
			const NodeCtx c2 = n0Leaf ? n0 : CreateChildCtx( n0, 7 );
			const NodeCtx c3 = n3Leaf ? n3 : CreateChildCtx( n3, 6 );

			const NodeCtx c4 = n1Leaf ? n1 : CreateChildCtx( n1, 1 );
			const NodeCtx c5 = n2Leaf ? n2 : CreateChildCtx( n2, 0 );
			const NodeCtx c6 = n1Leaf ? n1 : CreateChildCtx( n1, 3 );
			const NodeCtx c7 = n2Leaf ? n2 : CreateChildCtx( n2, 2 );

			// call ProcessEdge_Y for children sharing an edge along the Y-axis
			ProcessEdges_Y_Axis( c0, c1, c2, c3, c4, c5, c6, c7, mesh, ctx );	// 2 ProcessEdge_Y() calls
		}
	}
}
/*
	Processes the 4 subtrees sharing the edge along the Z axis:

                           Z*
         .-----.-----.     |  Y
        /  3  /  2  /|     | /
       .-----.-----.2|	   |/
      /  0  /  1  /|/	   O-------X
     .-----+-----|1/
     |  0  |  1  |/
     .-----.-----.
*/
void Octree_DC2::ProcessEdge_Z(
	const NodeCtx& n0,const NodeCtx& n1,const NodeCtx& n2,const NodeCtx& n3,
	AMeshBuilder& mesh, MeshCtx & ctx
	)
{
	// If one of the nodes is null, bail out.
	if( n0.id != NIL_NODE && n1.id != NIL_NODE && n2.id != NIL_NODE && n3.id != NIL_NODE )
	{
		const bool n0Leaf = IS_LEAF_ID(n0.id);
		const bool n1Leaf = IS_LEAF_ID(n1.id);
		const bool n2Leaf = IS_LEAF_ID(n2.id);
		const bool n3Leaf = IS_LEAF_ID(n3.id);
		// If all cubes are leaves, stop recursion and emit a quad.
		if( n0Leaf && n1Leaf && n2Leaf && n3Leaf )
		{
			// All nodes are surface leaves - determine the minimal cell/edge (containing no smaller edge).
			const Leaf2& leaf0 = m_leaves[ GET_ID(n0.id) ];
			const Leaf2& leaf1 = m_leaves[ GET_ID(n1.id) ];
			const Leaf2& leaf2 = m_leaves[ GET_ID(n2.id) ];
			const Leaf2& leaf3 = m_leaves[ GET_ID(n3.id) ];

			const NodeCtx* mc;
			mc = (n0.depth >= n1.depth) ? &n0 : &n1;
			mc = (mc->depth >= n2.depth) ? mc : &n2;
			mc = (mc->depth >= n3.depth) ? mc : &n3;

			// bottom -> top
			int sign1, sign2;	// 0 = outside, 1 = inside surface
			if(mc == &n0) {
				sign1 = (leaf0.signs >> 3) & 1;
				sign2 = (leaf0.signs >> 7) & 1;
			} else if(mc == &n1) {
				sign1 = (leaf1.signs >> 2) & 1;
				sign2 = (leaf1.signs >> 6) & 1;
			} else if(mc == &n2) {
				sign1 = (leaf2.signs >> 0) & 1;
				sign2 = (leaf2.signs >> 4) & 1;
			} else {
				sign1 = (leaf3.signs >> 1) & 1;
				sign2 = (leaf3.signs >> 5) & 1;
			}

			if( sign1 > sign2 ) {
				// The first corner is inside the surface, the second is outside
				EmitQuad( n0, n1, n2, n3, mesh, ctx );
			} else if( sign1 < sign2 ) {
				// The second corner is inside the surface, the first is outside
				EmitQuad( n3, n2, n1, n0, mesh, ctx );
			}
		}
		else
		{
			const NodeCtx c0 = n0Leaf ? n0 : CreateChildCtx( n0, 3 );
			const NodeCtx c1 = n1Leaf ? n1 : CreateChildCtx( n1, 2 );
			const NodeCtx c2 = n3Leaf ? n3 : CreateChildCtx( n3, 1 );
			const NodeCtx c3 = n2Leaf ? n2 : CreateChildCtx( n2, 0 );

			const NodeCtx c4 = n0Leaf ? n0 : CreateChildCtx( n0, 7 );
			const NodeCtx c5 = n1Leaf ? n1 : CreateChildCtx( n1, 6 );
			const NodeCtx c6 = n3Leaf ? n3 : CreateChildCtx( n3, 5 );
			const NodeCtx c7 = n2Leaf ? n2 : CreateChildCtx( n2, 4 );

			// call ProcessEdge_Z for children sharing an edge along the Z-axis
			ProcessEdges_Z_Axis( c0, c1, c2, c3, c4, c5, c6, c7, mesh, ctx );	// 2 ProcessEdge_Z() calls
		}
	}
}
// Generates a quad: ( v1, v2, v3, v4 );
void Octree_DC2::EmitQuad(
	const NodeCtx& n0,const NodeCtx& n1,const NodeCtx& n2,const NodeCtx& n3,
	AMeshBuilder& mesh,	MeshCtx & ctx
	)
{
	mxASSERT( IS_LEAF_ID(n0.id) && IS_LEAF_ID(n1.id) && IS_LEAF_ID(n2.id) && IS_LEAF_ID(n3.id) );
	const NodeID n0Index = GET_ID(n0.id);
	const NodeID n1Index = GET_ID(n1.id);
	const NodeID n2Index = GET_ID(n2.id);
	const NodeID n3Index = GET_ID(n3.id);
	//const Leaf2 leaf0 = m_leaves[n0Index];
	//const Leaf2 leaf1 = m_leaves[n1Index];
	//const Leaf2 leaf2 = m_leaves[n2Index];
	//const Leaf2 leaf3 = m_leaves[n3Index];
	mesh.AddTriangle( ctx.vertexIDs[n0Index], ctx.vertexIDs[n1Index], ctx.vertexIDs[n2Index] );
	mesh.AddTriangle( ctx.vertexIDs[n0Index], ctx.vertexIDs[n2Index], ctx.vertexIDs[n3Index] );
}

struct TraceWorks
{
	int	dir_mask;	// direction mask to support negative rays

	Float3 start;
	Float3 direction;

	RayCastResult* result;
};

// Based on "An Efficient Parametric Algorithm for Octree Traversal" [2000]:
// http://wscg.zcu.cz/wscg2000/papers_2000/x31.pdf
// very good explanation:
// http://chiranjivi.tripod.com/octrav.html

// Calculates the index of the first intersected node.
// To determine the entrance and exit planes of the ray,
// we find the maximum of the start t-values (t0's) and the minimum of the existing 3 t1's.
// To determine the first crossed node, we further subdivide the problem using the midplanes.
static inline int GetFirstNode(
							   const Float3& t0,	// times at which the ray enters the node
							   const Float3& tm		// times for the middle planes of the node
							   )
{
	// Remember, the time at which the ray enters the cube is Max( t0.x, t0.y, t0.z ).
	int result = 0;
	// Find the time with maximum value and figure out the plane which is hit first.
	if( t0.x > t0.y )
	{
		if( t0.x > t0.z )
		{
			// t0.x is max - entry plane is YZ,
			// the only possible first-intersected children will be 0, 2, 4, 6.
			if( t0.x > tm.y ) result |= CHILD_MASK_Y;	// node lies above middle-Y (2,6)
			if( t0.x > tm.z ) result |= CHILD_MASK_Z;	// node lies above middle-Z (4,6)
		}
		else
		{
			// t0.z is max - entry plane is XY,
			// the only possible first-intersected children will be 0, 1, 2, 3.
			if( t0.z > tm.x ) result |= CHILD_MASK_X;
			if( t0.z > tm.y ) result |= CHILD_MASK_Y;
		}
	}
	else
	{
		if( t0.y > t0.z )
		{
			// t0.y is max - entry plane is XZ,
			// the only possible first-intersected children will be 0, 1, 4, 5.
			if( t0.y > tm.x ) result |= CHILD_MASK_X;
			if( t0.y > tm.z ) result |= CHILD_MASK_Z;
		}
		else
		{
			// t0.z is max - entry plane is XY,
			// the only possible first-intersected children will be 0, 1, 2, 3.
			if( t0.z > tm.x ) result |= CHILD_MASK_X;
			if( t0.z > tm.y ) result |= CHILD_MASK_Y;
		}
	}
	return result;
}

// Calculates the exit plane and the next node, once the current parent node is exited.
// It simply returns the integer whose corresponding float is the smallest.
static inline int GetNextNode( float x, int nx, float y, int ny, float z, int nz )
{
	if( x < y ) {
		if( x < z ) {
			return nx;	// YZ plane
		} else {
			return nz;	// XY plane
		}
	}
	else
	{
		if( y < z ) {
			return ny;	// XZ plane
		} else {
			return nz;	// XY plane
		}
	}
}

void CastRayRecursive2(
						Octree_DC2& tree,
						NodeID nodeID,
						const Float3& t0,
						const Float3& t1,
						TraceWorks * tr
						)
{
	if( nodeID == NIL_NODE ) {
		return;
	}

	if( t1.x < 0.0f || t1.y < 0.0f || t1.z < 0.0f ) {
		return;	// the node is behind the ray
	}

	int nodeIndex = GET_ID(nodeID);

	if( IS_LEAF_ID(nodeID) )
	{
		tr->result->hitLeaves.Add( nodeIndex );
	}
	else
	{
		// Calculate the time parameters for the middle point of the node.
		const Float3 tm = (t0 + t1) * 0.5f;
		//g_hit_results.Add( tr->start + Float3_Multiply( tr->direction, tm  ));

		const Node2& node = tree.m_nodes[nodeIndex];
		// Because the system only works for positive components,
		// we have to reflect the ray about the midplane of the octree root for that negative axis.
		const int a = tr->dir_mask;

		enum { TERMINAL_NODE = 8 }; // <= this refers to an exit from the parent. 

		// Traverse the children in the order the ray pierces them.
		int currentOctant = GetFirstNode( t0, tm );
		do
		{
			switch( currentOctant )
			{
			case 0:
				CastRayRecursive2(tree, node.kids[0^a], Float3_Set(t0.x, t0.y, t0.z), Float3_Set(tm.x, tm.y, tm.z), tr);
				currentOctant = GetNextNode( tm.x, 1, tm.y, 2, tm.z, 4 );	// from child 0 the ray can go into children 1, 2, 4
				break;
			case 1:
				CastRayRecursive2(tree, node.kids[1^a], Float3_Set(tm.x, t0.y, t0.z), Float3_Set(t1.x, tm.y, tm.z), tr);
				currentOctant = GetNextNode( t1.x, TERMINAL_NODE, tm.y, 3, tm.z, 5 );
				break;
			case 2:
				CastRayRecursive2(tree, node.kids[2^a], Float3_Set(t0.x, tm.y, t0.z), Float3_Set(tm.x, t1.y, tm.z), tr);
				currentOctant = GetNextNode( tm.x, 3, t1.y, TERMINAL_NODE, tm.z, 6 );
				break;
			case 3:
				CastRayRecursive2(tree, node.kids[3^a], Float3_Set(tm.x, tm.y, t0.z), Float3_Set(t1.x, t1.y, tm.z), tr);
				currentOctant = GetNextNode( t1.x, TERMINAL_NODE, t1.y, TERMINAL_NODE, tm.z, 7 );
				break;
			case 4:
				CastRayRecursive2(tree, node.kids[4^a], Float3_Set(t0.x, t0.y, tm.z), Float3_Set(tm.x, tm.y, t1.z), tr);
				currentOctant = GetNextNode( tm.x, 5, tm.y, 6, t1.z, TERMINAL_NODE );
				break;
			case 5:
				CastRayRecursive2(tree, node.kids[5^a], Float3_Set(tm.x, t0.y, tm.z), Float3_Set(t1.x, tm.y, t1.z), tr);
				currentOctant = GetNextNode( t1.x, TERMINAL_NODE, tm.y, 7, t1.z, TERMINAL_NODE );
				break;
			case 6:
				CastRayRecursive2(tree, node.kids[6^a], Float3_Set(t0.x, tm.y, tm.z), Float3_Set(tm.x, t1.y, t1.z), tr);
				currentOctant = GetNextNode( tm.x, 7, t1.y, TERMINAL_NODE, t1.z, TERMINAL_NODE );
				break;
			case 7:
				CastRayRecursive2(tree, node.kids[7^a], Float3_Set(tm.x, tm.y, tm.z), Float3_Set(t1.x, t1.y, t1.z), tr);
				currentOctant = TERMINAL_NODE;	// there are no nodes we can reach from here as our ray is always traveling in a positive direction
				break;
			}
		} while( currentOctant < TERMINAL_NODE );
	}
}

// Based on "An Efficient Parametric Algorithm for Octree Traversal" [2000]
void Octree_DC2::CastRay(
	const Float3& start,
	const Float3& direction,
	RayCastResult &result,
	const Options& options	
)
{
	result.hitAnything = false;
	result.hitLeaves.Empty();

	TraceWorks	tr;
	tr.start = start;
	tr.direction = direction;
	tr.dir_mask = 0;
	tr.result = &result;

	AABB24	aabb;
	aabb.min_point = Float3_Replicate(-options.radius);
	aabb.max_point = Float3_Replicate(+options.radius);
	// the world AABB is centered around the origin
	//Float3 center = Float3_Zero();

	Float3 o = start;
	Float3 dir = Float3_Normalized(direction);

	// Get rid of small ray direction components to avoid division by zero.
	const float epsilon = 1e-4f;

	if( Float_Abs(dir.x) < epsilon ) {
		dir.x = Float_CopySign( epsilon, dir.x );
	}
	if( Float_Abs(dir.y) < epsilon ) {
		dir.y = Float_CopySign( epsilon, dir.y );
	}
	if( Float_Abs(dir.z) < epsilon ) {
		dir.z = Float_CopySign( epsilon, dir.z );
	}

	// 'Normalize' the ray.

	// If one of the components of the direction vector is negative,
	// the ray is mirrored about the corresponding middle plane of the octree
	// and a bit mask is set to recall this fact.
	tr.dir_mask = 0;
	if( dir.x < 0.0f ) {
		o.x = -o.x;
		dir.x = -dir.x;
		tr.dir_mask |= CHILD_MASK_X;
	}
	if( dir.y < 0.0f ) {
		o.y = -o.y;
		dir.y = -dir.y;
		tr.dir_mask |= CHILD_MASK_Y;
	}
	if( dir.z < 0.0f ) {
		o.z = -o.z;
		dir.z = -dir.z;
		tr.dir_mask |= CHILD_MASK_Z;
	}

	mxASSERT(dir.x > 0);
	mxASSERT(dir.y > 0);
	mxASSERT(dir.z > 0);

	const Float3 invDir = Float3_Reciprocal( dir );

	// vector to hold entry time parameters
	Float3 entryTimes = Float3_Multiply( aabb.min_point - o, invDir );
	// entry time parameters
	Float3 exitTimes = Float3_Multiply( aabb.max_point - o, invDir );

	// calculate the time at which we enter and exit the octant
	float tmin = Max3( entryTimes.x, entryTimes.y, entryTimes.z );
	// the smallest of the t-values of the three exit planes
	float tmax = Min3( exitTimes.x, exitTimes.y, exitTimes.z );

	// if the ray hits the root AABB
	if( (tmin < tmax) && ( exitTimes.x > 0 && exitTimes.y > 0 && exitTimes.z > 0 ) )
	{
		Float3 entry = start + direction * tmin;// ray entry point
		Float3 exit = start + direction * tmax;	// ray exit point
		// middle point of the root node
		//Float3 mid = start + Float3_Multiply( direction, (entryTimes + exitTimes) * 0.5f );

		//int entry_octant = 0;
		//if( entry.x > center.x ) {
		//	entry_octant |= CHILD_MASK_X;
		//}
		//if( entry.y > center.y ) {
		//	entry_octant |= CHILD_MASK_Y;
		//}
		//if( entry.z > center.z ) {
		//	entry_octant |= CHILD_MASK_Z;
		//}
		CastRayRecursive2( *this, 0, entryTimes, exitTimes, &tr );

		result.rayEntry = entry;
		result.rayExit = exit;
		result.hitAnything = (result.hitLeaves.Num() > 0);
	}
	else
	{
		//DBGOUT("no intersection");
	}
}

void DebugGatherHitLeaves(
						  const RayCastResult& input,
						  const Octree_DC2& tree,
						  DualContouring::NodeID nodeID,
						  const OctCubeF& bounds,
						  TArray< DebugLeafInfo > &output
						  )
{
	const bool isLeaf = DualContouring::IS_LEAF_ID(nodeID);
	const UINT16 nodeIndex = DualContouring::GET_ID(nodeID);

	AABB24	aabb;
	OctBoundsToAABB(bounds, aabb);

	if( isLeaf )
	{
		const DualContouring::Leaf2& leaf = tree.m_leaves[nodeIndex];

		if( input.hitLeaves.Contains(nodeIndex) )
		{
			DebugLeafInfo& debugLeafInfo = output.Add();
			debugLeafInfo.leafIndex = nodeIndex;
			debugLeafInfo.aabb = aabb;
			DequantizePosition( bounds, leaf.xyz, &debugLeafInfo.xyz );
			debugLeafInfo.N = DequantizeNormal( leaf.N );
		}
	}
	else//if(!isLeaf)
	{
		OctCubeF octants[8];
		GetChildOctants(bounds,octants);

		const DualContouring::Node2& node = tree.m_nodes[nodeIndex];
		for( int i = 0; i < 8; i++ )
		{
			if( node.kids[i] != DualContouring::NIL_NODE ) {
				DebugGatherHitLeaves( input, tree, node.kids[i], octants[i], output );
			}
		}
	}
}

void DBG_CSG_Subtract_R(
						Octree_DC2& tree,
						NodeID nodeID,
						AVolume* vol
						)
{
UNDONE;
}

void DBG_CSG_Subtract( Octree_DC2& tree, const Float3& center, float radius )
{
	SphereSDF vol;
	vol.center = center;
	vol.radius = radius;
	DBG_CSG_Subtract_R(tree,0,&vol);
}

ERet Octree_DC2::Save( AStreamWriter &stream ) const
{
	stream << m_nodes;
	stream << m_leaves;
	stream << m_edges;
	return ALL_OK;
}

ERet Octree_DC2::Load( AStreamReader& stream )
{
	stream >> m_nodes;
	stream >> m_leaves;
	stream >> m_edges;
	return ALL_OK;
}

}//namespace DualContouring
