/*
	Adaptive (octree-based) dual contouring.

	Nodes always contain pointers to 8 children
	(some of which may be internal nodes or leaves).
*/
#pragma once

#include <Base/Memory/FreeList/TPool.h>
#include <Core/VectorMath.h>
#include <Graphics/Device.h>
#include <Graphics/Utils.h>
#include <Meshok/Meshok.h>
#include <Meshok/VoxelEngine.h>
#include <Utility/EngineUtil/EngineUtil.h>

#include "contouring.h"
#include "QEF.h"

// Build octree to a certain subdivision depth/cell size
// and recursively collapse leaves using QEF metric.
#define USE_QEF_SIMPLIFICATION (0)


namespace DualContouring
{
	using namespace Contouring;

	// NOTE: zero material ID is reserved for 'AIR' ('empty space').

	#pragma pack(push,1)
	struct Leaf2
	{
		UINT8 xyz[3];	//3! quantized position of representative vertex
		UINT8 signs;	//1! eight signs for the corner points
		UINT8 N[3];		//3! quantized vertex normal
		UINT8 _pad0;	//1! unused padding
		UINT16 U,V;		//4! texture coordinates
		UINT16 matID;	//2! material id
		UINT16 qefID;	//2! QEF id
		//16!
		UINT16 edges[8];//16! Hermite data: indices to ActiveEdge structs
	};//32!
	#pragma pack(pop)
	mxSTATIC_ASSERT(sizeof(Leaf2) == 32);
	mxDECLARE_POD_TYPE(Leaf2);

	struct QEF_t : svd::QefData//14*4=56!
	{
		UINT8 _pad[8];//8! unused padding
	};//64!
	mxSTATIC_ASSERT(sizeof(QEF_t) == 64);

	// grid points are shared by neighboring nodes
	struct GridPoint
	{
		UINT8	matID;	//2! material id
	};

	// data for a zero-crossing edge
	struct ActiveEdge
	{
		Float3 normal;	//12!
		float distance;	//4! distance from solid to empty grid point, normalized to [0..1]
	};//16!

	enum { NIL_NODE = (UINT16)~0 };

	//NOTE: upper bit = is_leaf flag
	typedef UINT16 NodeID;

	inline NodeID MAKE_LEAF_ID( NodeID nodeID ) {
		return nodeID | (1<<15);
	}
	inline bool IS_LEAF_ID( NodeID nodeID ) {
		return (nodeID & (1<<15));
	}
	inline NodeID GET_ID( NodeID nodeID ) {
		return nodeID & ~(1<<15);
	}

	// Branch cells have eight children and no surface data
	// i tried to make the tree more dynamic
	struct Node2
	{
		NodeID kids[8];	//16! indices to 8 children; upper bit = is_leaf
		//Leaf2 lod;
	};//16!
	mxDECLARE_POD_TYPE(Node2);

	// node context
	struct NodeCtx
	{
		NodeID id;
		UINT32 depth;
	};

	// mesh context
	struct MeshCtx
	{
		UINT16* vertexIDs;	// maps leaf node index to vertex index
	};

	struct Octree_DC2
	{
		Octree_DC2();
		~Octree_DC2();

		// creates an octree from the density source
		void Build(
			const AVolume* _volume,
			QEF_Solver* _qef_solver,
			OctStats &_octree_stats,
			const Options& options = Options()
		);

		void CastRay(
			const Float3& start,
			const Float3& direction,
			RayCastResult &result,
			const Options& options = Options()
		);

		ERet Save( AStreamWriter &stream ) const;
		ERet Load( AStreamReader& stream );

	public:
		TPool< Node2 >	m_nodes;
		TPool< Leaf2 >	m_leaves;
		TPool< QEF_t >	m_QEFs;	// 
		TPool< ActiveEdge >	m_edges;	//@todo: optimize, remove

		NodeCtx CreateChildCtx( const NodeCtx& parent, int child ) const
		{
			mxASSERT(parent.id != NIL_NODE);
			mxASSERT(!IS_LEAF_ID(parent.id));
			NodeCtx result;
			result.id = m_nodes[ parent.id ].kids[ child ];
			result.depth = parent.depth + 1;
			return result;
		}

	public:

		NodeID BuildOctreeRecursive(
			const OctCubeF& _bounds,
			const UINT32 _treeLevel,
			QEF_Solver* _qef_solver,
			const AVolume* _volume,
			const Options& options,
			OctStats &stats
		);
		NodeID TryCreateLeaf(
			const OctCubeF& _bounds,
			const UINT32 _treeLevel,
			QEF_Solver* _qef_solver,
			const AVolume* _volume,
			float _error_threshold,
			OctStats &stats
		);
		NodeID TryCreateLeaf2(
			const OctCubeF& _bounds,
			const UINT32 _treeLevel,
			QEF_Solver* _qef_solver,
			const AVolume* _volume,
			float _error_threshold,
			OctStats &stats
		);
		NodeID Simplify(
			const NodeID _nodeIndex,
			const OctCubeF& _bounds,
			const UINT32 _treeLevel,
			QEF_Solver* _qef_solver,
			const AVolume* _volume,
			const Options& options,
			OctStats &stats
		);

		ERet Triangulate(
			AMeshBuilder & mesh,
			const Options& options,
			OctStats &stats
		);

	private:
		void FreeNode( NodeID nodeID );

		void ProcessNode(
			const NodeCtx& n,
			AMeshBuilder& mesh, MeshCtx & ctx
		);

		// visit faces along the X axis
		void ProcessFaces_X(
			const NodeCtx& n0,const NodeCtx& n1,const NodeCtx& n2,const NodeCtx& n3,
			const NodeCtx& n4,const NodeCtx& n5,const NodeCtx& n6,const NodeCtx& n7,
			AMeshBuilder& mesh, MeshCtx & ctx
			)
		{
			// call ProcessFace_X() for children sharing a face in the ZY plane ('left-right')
			ProcessFace_X( n0, n1, mesh, ctx );
			ProcessFace_X( n2, n3, mesh, ctx );
			ProcessFace_X( n6, n7, mesh, ctx );
			ProcessFace_X( n4, n5, mesh, ctx );
		}
		// visit faces along the Y axis
		void ProcessFaces_Y(
			const NodeCtx& n0,const NodeCtx& n1,const NodeCtx& n2,const NodeCtx& n3,
			const NodeCtx& n4,const NodeCtx& n5,const NodeCtx& n6,const NodeCtx& n7,
			AMeshBuilder& mesh, MeshCtx & ctx
			)
		{
			// call ProcessFace_Y() for children sharing a face in the XZ plane ('front-back')
			ProcessFace_Y( n0, n2, mesh, ctx );
			ProcessFace_Y( n1, n3, mesh, ctx );
			ProcessFace_Y( n5, n7, mesh, ctx );
			ProcessFace_Y( n4, n6, mesh, ctx );
		}
		// visit faces along the Z axis
		void ProcessFaces_Z(
			const NodeCtx& n0,const NodeCtx& n1,const NodeCtx& n2,const NodeCtx& n3,
			const NodeCtx& n4,const NodeCtx& n5,const NodeCtx& n6,const NodeCtx& n7,
			AMeshBuilder& mesh, MeshCtx & ctx
			)
		{
			// call ProcessFace_Z() for children sharing a face in the XY plane ('bottom-top')
			ProcessFace_Z( n0, n4, mesh, ctx );
			ProcessFace_Z( n1, n5, mesh, ctx );
			ProcessFace_Z( n3, n7, mesh, ctx );
			ProcessFace_Z( n2, n6, mesh, ctx );
		}

		void ProcessEdges_X_Axis(
			const NodeCtx& n0,const NodeCtx& n1,const NodeCtx& n2,const NodeCtx& n3,
			const NodeCtx& n4,const NodeCtx& n5,const NodeCtx& n6,const NodeCtx& n7,
			AMeshBuilder& mesh, MeshCtx & ctx
			)
		{
			// call ProcessEdge_X for children sharing an edge along the X-axis
			ProcessEdge_X( n0, n2, n6, n4, mesh, ctx );	// left
			ProcessEdge_X( n1, n3, n7, n5, mesh, ctx );	// right
		}
		void ProcessEdges_Y_Axis(
			const NodeCtx& n0,const NodeCtx& n1,const NodeCtx& n2,const NodeCtx& n3,
			const NodeCtx& n4,const NodeCtx& n5,const NodeCtx& n6,const NodeCtx& n7,
			AMeshBuilder& mesh, MeshCtx & ctx
			)
		{
			// call ProcessEdge_Y for children sharing an edge along the Y-axis
			ProcessEdge_Y( n0, n4, n5, n1, mesh, ctx );	// front
			ProcessEdge_Y( n2, n6, n7, n3, mesh, ctx );	// back
		}
		void ProcessEdges_Z_Axis(
			const NodeCtx& n0,const NodeCtx& n1,const NodeCtx& n2,const NodeCtx& n3,
			const NodeCtx& n4,const NodeCtx& n5,const NodeCtx& n6,const NodeCtx& n7,
			AMeshBuilder& mesh, MeshCtx & ctx
			)
		{
			// call ProcessEdge_Z for children sharing an edge along the Z-axis
			ProcessEdge_Z( n0, n1, n3, n2, mesh, ctx );	// bottom
			ProcessEdge_Z( n4, n5, n7, n6, mesh, ctx );	// top
		}

		void ProcessFace_X(
			const NodeCtx& n1,	//<= 'left' node
			const NodeCtx& n2,	//<= 'right' node
			AMeshBuilder& mesh, MeshCtx & ctx
			);
		void ProcessFace_Y(
			const NodeCtx& n1,	//<= 'front' node
			const NodeCtx& n2,	//<= 'back' node
			AMeshBuilder& mesh, MeshCtx & ctx
			);
		void ProcessFace_Z(
			const NodeCtx& n1,	//<= lower node
			const NodeCtx& n2,	//<= upper node
			AMeshBuilder& mesh, MeshCtx & ctx
			);

		void ProcessEdge_X(
			const NodeCtx& n0,const NodeCtx& n1,const NodeCtx& n2,const NodeCtx& n3,
			AMeshBuilder& mesh, MeshCtx & ctx
			);
		void ProcessEdge_Y(
			const NodeCtx& n0,const NodeCtx& n1,const NodeCtx& n2,const NodeCtx& n3,
			AMeshBuilder& mesh, MeshCtx & ctx
			);
		void ProcessEdge_Z(
			const NodeCtx& n0,const NodeCtx& n1,const NodeCtx& n2,const NodeCtx& n3,
			AMeshBuilder& mesh, MeshCtx & ctx
			);
		void EmitQuad(
			const NodeCtx& n0,const NodeCtx& n1,const NodeCtx& n2,const NodeCtx& n3,
			AMeshBuilder& mesh,	MeshCtx & ctx
			);
	};

	void DBG_CSG_Subtract( Octree_DC2& tree, const Float3& center, float radius );

	void DebugGatherHitLeaves(
		const RayCastResult& input,
		const Octree_DC2& tree,
		DualContouring::NodeID nodeID,
		const OctCubeF& bounds,
		TArray< DebugLeafInfo > &output
	);

}//namespace DualContouring
