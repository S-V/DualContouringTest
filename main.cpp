#include "stdafx.h"
#pragma hdrstop

// for std::multimap
#include <map>
#include <algorithm>

//#include <d3dx11.h>

#include <Base/Util/LogUtil.h>
#include <Base/Util/FPSTracker.h>
#include <Base/Util/Sorting.h>
// for TCallback<>
#include <Base/Template/Delegate/Delegate.h>
#include <Base/Math/NewMath.h>
#include <Core/bitsquid/collection_types.h>
#include <Core/Event.h>
#include <Core/Client.h>
#include <Core/FileSystem.h>
#include <Core/VectorMath.h>
#include <Core/ObjectModel.h>
#include <Core/Serialization.h>

#include <Graphics/Device.h>
#include <Graphics/Utils.h>

#include <Renderer/Mesh.h>
#include <Renderer/Model.h>
#include <Renderer/Light.h>
#include <Renderer/Deferred.h>

#include <Scripting/Scripting.h>
#include <Scripting/FunctionBinding.h>
	#include <Scripting/Lua_Helpers.h>

#include <Utility/TxTSupport/TxTSerializers.h>
#include <Utility/TxTSupport/TxTReader.h>
#include <Utility/TxTSupport/TxTConfig.h>

#include <Utility/Meshok/MeshImporter.h>
#include <Utility/Meshok/BSP.h>
#include <Utility/Meshok/ASDF.h>

#include <Utility/EngineUtil/EngineUtil.h>
#include <Utility/EngineUtil/ImGUI_Renderer.h>
#include <Utility/GameUtil/GameUtil.h>
#include <Utility/DemoFramework/DemoFramework.h>

#include <Utility/EditorSupport/DevAssetFolder.h>

#include <Utility/AssetPipeline/AssetPipeline.h>

// STB_IMAGE_IMPLEMENTATION is defined in Graphics
//#define STB_IMAGE_IMPLEMENTATION
#include <External/stb/stb_image.h>

#include <Meshok/marching_cubes.h>

#include "contouring.h"
#include "dual_contouring.h"
#include "dual_contouring2.h"
#include "dual_marching_cubes.h"
#include "cubical_marching_squares.h"

namespace DC = DualContouring;
namespace DMC = DualMarchingCubes;
namespace CMS = CubicalMarchingSquares;

#define USE_OCTREE_2 (1)

//const float SCENE_RADIUS = 128;
//const float SCENE_RADIUS = 256;
//const float SCENE_RADIUS = 512;
//const float SCENE_RADIUS = 1024;

const char* CONFIG_FILE_NAME = "demo-config-generic.son";
const char* SETTINGS_FILE_NAME = "contouring_demo_settings.son";
const char* OCTREE_SAVE_FILE_NAME = "octree.bin";


void AABB_ScaleSize( AABB24 *aabb, float scale )
{
	Float3 center = AABB_Center(*aabb);
	Float3 radius = AABB_Extent(*aabb);
	radius *= scale;
	aabb->min_point = center - radius;
	aabb->max_point = center + radius;
}

namespace ContouringDemo
{
using namespace Contouring;

static bool g_rebuild_octree = false;
static bool g_regenerate_mesh = false;

struct SolidWireSettings : CStruct
{
	float LineWidth;
	float FadeDistance;
	float PatternPeriod;
	Float4 FillColor;
	Float4 WireColor;
	Float4 PatternColor;
public:
	mxDECLARE_CLASS(SolidWireSettings, CStruct);
	mxDECLARE_REFLECTION;
	SolidWireSettings()
	{
		this->SetDefauls();
	}
	void SetDefauls()
	{
		LineWidth = 0.5;
		FadeDistance = 50;
		PatternPeriod = 1.5;
		FillColor = Float4_Set( 0.1, 0.2, 0.4, 1 );
		WireColor = Float4_Set( 1, 1, 1, 1 );
		PatternColor = Float4_Set( 1, 1, 0.5, 1 );
	}
};
mxDEFINE_CLASS(SolidWireSettings);
mxBEGIN_REFLECTION(SolidWireSettings)
	mxMEMBER_FIELD( LineWidth ),
	mxMEMBER_FIELD( FadeDistance ),
	mxMEMBER_FIELD( PatternPeriod ),
	mxMEMBER_FIELD( FillColor ),
	mxMEMBER_FIELD( WireColor ),
	mxMEMBER_FIELD( PatternColor ),
mxEND_REFLECTION;

struct Settings : CStruct
{
	DualContouring::Options dc;
	SolidWireSettings	solidWire;
	bool visualize_octree;
	bool draw_octree_leaves_only;
	bool drawVertexNormals;
	bool drawWireframeMesh;	
public:
	mxDECLARE_CLASS(Settings, CStruct);
	mxDECLARE_REFLECTION;
	Settings()
	{
		this->SetDefauls();
	}
	void SetDefauls()
	{
		solidWire.SetDefauls();
		visualize_octree = false;
		draw_octree_leaves_only = false;
		drawVertexNormals = false;
		drawWireframeMesh = true;		
	}
};
mxDEFINE_CLASS(Settings);
mxBEGIN_REFLECTION(Settings)
	mxMEMBER_FIELD( dc ),
	mxMEMBER_FIELD( solidWire ),
	mxMEMBER_FIELD( visualize_octree ),
	mxMEMBER_FIELD( draw_octree_leaves_only ),
	mxMEMBER_FIELD( drawVertexNormals ),
	mxMEMBER_FIELD( drawWireframeMesh ),
mxEND_REFLECTION;

static Settings	g_settings;

AABB24 GetSceneBounds()
{
	const float sceneRadius = g_settings.dc.radius;
	AABB24 result;
	result.min_point = Float3_Replicate(-sceneRadius);
	result.max_point = Float3_Replicate(+sceneRadius);
	return result;
}
OctCubeF GetSceneCube()
{
	OctCubeF worldBounds;
	worldBounds.x = 0;	worldBounds.y = 0;	worldBounds.z = 0;
	worldBounds.radius = g_settings.dc.radius;
	return worldBounds;
}

class CsgScene
{
	SphereSDF			sphereSDF;
	SphereSDF			sphereSDF2;
	SphereSDF			sphereSDF3;
	HalfSpaceSDF		planeSDF;
	AxisAlignedBoxSDF	aaboxSDF;
	TorusSDF			torusSDF;
	CSGUnion			union1;
	CSGUnion			union2;
	CSGSubtraction		subtraction;
	CSGSubtraction		subtraction2;

	TransformSDF		boxXform;

public:
	CsgScene()
	{
	}

	void Setup( const Settings& settings )
	{
		{
			sphereSDF.center = Float3_Zero();
			sphereSDF.radius = 96;
		}
		{
			//sphereSDF2.center = Float3_Set(100,20,100);
			sphereSDF2.radius = settings.dc.radius * 0.3f;
		}
		{
			sphereSDF3.center = Float3_Zero();
			sphereSDF3.radius = settings.dc.radius * 0.9f;
		}

		{
			const Float3 planeNormal = 
				//Float3_Normalized(Float3_Set(1,1,1))
				Float4_As_Float3(Float4_Up)
				;
			planeSDF.plane = Float4_Set(planeNormal.x,planeNormal.y,planeNormal.z,-10);
		}
		{
			//aaboxSDF.origin = Float3_Zero();
			//aaboxSDF.extent = Float3_Set(40, 60, 20);
			//aaboxSDF.center = Float3_Set(64, 64, 64);
			//aaboxSDF.size = Float3_Set(128, 128, 128);
			aaboxSDF.size = Float3_Replicate(256);
		}
		{
			//torusSDF.T = Float2_Set(128,128);
			torusSDF.T = Float2_Set(512,128);
		}

		subtraction.opA = &planeSDF;
		subtraction.opB = &sphereSDF;

		//substraction2.opA = &substraction;
		//substraction2.opB = &sphereSDF2;

		subtraction2.opA = &planeSDF;
		subtraction2.opB = &sphereSDF2;

		union1.opA = &subtraction2;
		//union1.opB = &aaboxSDF;
		union1.opB = &torusSDF;

		boxXform.op = &aaboxSDF;
		boxXform.xform = Matrix_Multiply( Matrix_RotationZ(DEG2RAD(45)), Matrix_RotationY(DEG2RAD(45)) );

		union2.opA = &union1;
		//union2.opB = &boxXform;
		union2.opB = &aaboxSDF;
	}
	AVolume* GetRoot() {
		return &planeSDF;
		//return &subtraction2;
		//return &union2;
		//return &aaboxSDF;
		//return &sphereSDF3;
	}
};

class MyMenuState : public MenuState {
public:
	virtual void Update( double deltaSeconds ) override
	{}
	virtual void RenderGUI() override
	{
		if(ImGui::Begin("Debug Visualization"))
		{
			ImGui::Checkbox("Visualize vertex normals",&g_settings.drawVertexNormals);
			ImGui::Checkbox("Draw meshes in wireframe",&g_settings.drawWireframeMesh);
			ImGui::End();
		}

		if(ImGui::Begin("Settings"))
		{
			ImGui_DrawPropertyGrid(&g_settings,mxCLASS_OF(g_settings),"Settings");
			ImGui::End();
		}

		if( ImGui::Begin("Control Panel") )
		{
			if(ImGui::Button("Rebuild"))
			{
				DBGOUT("rebuilding octree!");
				g_rebuild_octree = true;
			}
			if(ImGui::Button("EXIT"))
			{
				WindowsDriver::RequestExit();
			}
			ImGui::End();
		}		

		//ImGui::ShowTestWindow();
	}
};


struct TriangleMeshSDF : AVolume
{
	TPtr< TcMeshData >	meshData;
	TPtr< BspTree >		meshTree;
public:
	Float3 GetClosestPointTo( const Float3& _position ) const
	{
		Float3 closestPointResult;
		float minimumDistanceSquared = BIG_NUMBER;
		for( int iSubMesh = 0; iSubMesh < meshData->sets.Num(); iSubMesh++ )
		{
			const TcTriMesh& subMesh = meshData->sets[ iSubMesh ];

			for( int iTri = 0; iTri < subMesh.indices.Num(); iTri += 3 )
			{
				const UINT32 i0 = subMesh.indices[ iTri + 0 ];
				const UINT32 i1 = subMesh.indices[ iTri + 1 ];
				const UINT32 i2 = subMesh.indices[ iTri + 2 ];
				Float3 triangle[3];
				triangle[0] = subMesh.positions[i0];
				triangle[1] = subMesh.positions[i1];
				triangle[2] = subMesh.positions[i2];

				Float3 closestPoint = ClosestPointOnTriangle( triangle, _position );
				float distanceSquared = Float3_LengthSquared( _position - closestPoint );
				if( distanceSquared < minimumDistanceSquared ) {
					minimumDistanceSquared = distanceSquared;
					closestPointResult = closestPoint;
				}
			}
		}
		return closestPointResult;
	}
	// we compute the distance from a point to M by computing the distance
	// from that point to every triangle in M and picking the smallest distance to any triangle.
	virtual float GetDistanceAt( const Float3& _position ) const override
	{
		Float3 closestPoint = GetClosestPointTo( _position );
		float result = Float3_Length( _position - closestPoint );
		if( meshTree->PointInSolid(_position,1e-6f) ) {
			result *= -1.0f;
		}
		return result;
	}
};


bool PlanePlaneIntersection( const Vector4& a, const Vector4& b, Float3& start, Float3& dir )
{
	const Float3 N1 = Plane_GetNormal( a );
	const Float3 N2 = Plane_GetNormal( b );

	double n00, n01, n11, det, invDet, f0, f1;

	n00 = Float3_LengthSquared( N1 );
	n01 = N1 * N2;
	n11 = Float3_LengthSquared( N2 );
	det = n00 * n11 - n01 * n01;

	if ( Float_Abs(det) < 1e-6f ) {
		return false;
	}

	invDet = 1.0f / det;
	f0 = ( n01 * b.w - n11 * a.w ) * invDet;
	f1 = ( n01 * a.w - n00 * b.w ) * invDet;

	dir = Float3_Cross( N1, N2 );
	start = N1 * f0 + N2 * f1;
	return true;
}

bool PlaneRayIntersection( const Vector4& plane, const Float3& start, const Float3& dir, float &scale )
{
	const Float3 N = Plane_GetNormal( plane );

	float d1, d2;

	d1 = Float3_Dot( N, start ) + plane.w;
	d2 = Float3_Dot( N, dir );
	if ( d2 == 0.0f ) {
		return false;
	}
	scale = -( d1 / d2 );
	return true;
}

//LPStatus
struct Entry2
{
	int plane;
	float dist;
};
float FindClosestPoint( const BspTree& tree, const Entry2* entries, int num_entries, const Float3& _point )
{
	// the optimum lies on the plane
	Vector4 plane1 = tree.m_planes[entries[0].plane];
	int dim = 0;	// 0 - plane, 1 - line, 2 - point

	Float3 lineStart, lineDir;

	Float3 closestPt;

	for( int i = 1; i < num_entries; i++ )
	{
		const Entry2& other = entries[i];
		const Vector4 plane2 = tree.m_planes[entries[i].plane];
		if( dim == 0 ) {
			// plane vs plane
			if( PlanePlaneIntersection(plane1, plane2, lineStart, lineDir) )
			{
				// the optimum point lies on the line
				dim = 1;	// plane -> line
			}
		}
		else if( dim == 1 ) {
			// line vs plane
			// the optimum lies at the line-plane intersection point
			// 

			if( PlaneRayIntersection() )
			{
				dim = 2;
			}
		} else {
			mxASSERT(dim==2);
		}
	}
	//Entry e, e2;
	//e = stack.Pop();
	//if( stack.IsEmpty() ) {
	//	return e.dist;
	//}
	//e2 = stack.Pop();

	UNDONE;
	return 0;
}

struct BSP_SDF : public AVolume
{
	TPtr< BspTree >	tree;

public:
	virtual float GetDistanceAt( const Float3& _point ) const override
	{
		float epsilon = 1e-4f;
	
		TStaticList< Entry2, 64 > stack;

		float neg_dist = 0;

		int nodeIndex = 0;
		float distance = 0.f;
		// If < 0, we are in a leaf node
		while( nodeIndex >= 0 )
		{
			// Find which side of the node we are on
			const BspNode& node = tree->m_nodes[ nodeIndex ];
			const Vector4& plane = tree->m_planes[ node.plane ];
			const float distance = Plane_PointDistance( plane, _point );
			// Go down the appropriate side
			if( distance >= 0.0f ) {
				nodeIndex = (INT16)node.front;
				Entry2& newEntry = stack.Add();
				newEntry.plane = node.plane;
				newEntry.dist = distance;
			} else {
				nodeIndex = (INT16)node.back;
				neg_dist = minf(neg_dist, distance);
			}
		}
		if( nodeIndex == BSP_EMPTY_LEAF )
		{
			// the point is outside the solid
			return FindClosestPoint(*tree, stack.ToPtr(), stack.Num(), _point);
		}
		else
		{
			// the point is inside the solid
			// return the minimum (closest) distance
			mxASSERT(neg_dist < 0.0f);
			return neg_dist;
		}
	}
};

struct App : DemoUtil::DemoApp
	, protected TGlobal< App >	// <= nasty, but needed for static/global functions (CompareDepth())
{
	Clump *		m_sceneData;
	Clump *		m_rendererData;

	FlyingCamera	camera;

	HContext		renderContext;

	DeferredRenderer	m_render_system;
	BatchRenderer		m_batch_renderer;
	GUI_Renderer		m_gui_renderer;

	AxisArrowGeometry	gizmo;

	CsgScene	m_scene;
#if USE_OCTREE_2
	DualContouring::Octree_DC2	m_octree;
#else
	DualContouring::Octree	m_octree;
#endif
	rxMesh *	m_mesh;

	//DualContouring::QEF_Solver_ParticleBased m_qefSolver;
	DualContouring::QEF_Solver_SVD m_qefSolver;


	ASDF::Octree	m_asdf_tree;
	ASDF::Volume	m_asdf_volume;


	TcMeshData		meshData;
	BspTree			meshBspTree;
	BSP_SDF			m_bsp_sdf;
	TriangleMeshSDF	meshSDF;


	MeshBuilder	m_meshBuffer;

	MyMenuState	m_menuState;

	bool m_debugShowNormalsRT;

	Float3			m_lastRayCastStart;
	Float3			m_lastRayCastDirection;
	RayCastResult	m_lastRayCastResult;

	TArray< DebugLeafInfo >	m_hitLeaves;

	TArray< AVolume* >	m_garbage;


	Float3 m_pointPosition;
	Float3 m_closestPointOnMesh;


public:
	App()
	{
		m_debugShowNormalsRT = false;
	}
	ERet Setup()
	{
		SON::LoadFromFile(SETTINGS_FILE_NAME,g_settings);
		//dev_load_state(GameActions::dev_load_state,IS_Pressed,0);

		mxDO(DemoApp::__Initialize());

		this->renderContext = llgl::GetMainContext();

		mxDO(Rendering::RegisterClasses());

//mxDO(LoadClump(MakeAssetID("simple.renderer"), m_rendererData));
mxDO(LoadClump(MakeAssetID("deferred.renderer"), m_rendererData));
		//DemoUtil::DBG_PrintClump(*m_rendererData);

		mxDO(m_render_system.Initialize(m_rendererData));
		mxDO(m_batch_renderer.Initialize());
		mxDO(m_gui_renderer.Initialize(m_rendererData));

		this->gizmo.BuildGeometry();

		return ALL_OK;
	}
	void Shutdown()
	{
		//m_garbage.DeleteContents();
		m_garbage.Empty();

		m_gui_renderer.Shutdown();
		m_batch_renderer.Shutdown();
		m_render_system.Shutdown();

		m_sceneData->Clear();
		delete m_sceneData;
		m_sceneData = NULL;

		DemoApp::__Shutdown();

		//dev_save_state(GameActions::dev_save_state,IS_Pressed,0);
	}

	ERet GenerateMesh( AVolume* vol, rxMesh* mesh )
	{
		const AABB24 worldAabb = GetSceneBounds();
		int gridSize[3] = {32,32,32};
		MarchingCubes::Contour(vol,worldAabb.min_point,worldAabb.max_point,gridSize,m_meshBuffer);
		mxDO(CreateMesh(m_meshBuffer, mesh));
		return ALL_OK;
	}

	ERet ReallocateBuffer( HBuffer* handle, EBufferType type, UINT32 oldSize, UINT32 newSize, const void* data )
	{
		bool reCreateBuffer = handle->IsNull() || (newSize > oldSize);

		if( reCreateBuffer )
		{
			if( handle->IsValid() )
			{
				llgl::DeleteBuffer( *handle );
			}
			//*handle = llgl::CreateBuffer( type, newSize, data );
			*handle = llgl::CreateBuffer( type, newSize );
		}
		//else
		{
			llgl::UpdateBuffer( llgl::GetMainContext(), *handle, newSize, data );
		}

		return ALL_OK;
	}

	ERet CreateMesh( const MeshBuilder& source, rxMesh* mesh )
	{
		{
			const UINT32 oldSize = mesh->m_numVertices * sizeof(source.vertices[0]);
			mxDO(ReallocateBuffer( &mesh->m_vertexBuffer, Buffer_Vertex, oldSize, source.vertices.GetDataSize(), source.vertices.ToPtr() ));
			mesh->m_numVertices = source.vertices.Num();
		}
		{
			const UINT32 oldSize = mesh->m_numIndices * sizeof(source.indices[0]);
			mxDO(ReallocateBuffer( &mesh->m_indexBuffer, Buffer_Index, oldSize, source.indices.GetDataSize(), source.indices.ToPtr() ));
			mesh->m_indexStride = sizeof(source.indices[0]);
			mesh->m_numIndices = source.indices.Num();
		}
		{
			rxSubmesh* submesh0 = NULL;
			if( mesh->m_parts.IsEmpty() ) {
				submesh0 = &mesh->m_parts.Add();
			} else {
				submesh0 = &mesh->m_parts.GetFirst();
			}
			submesh0->startIndex = 0;
			submesh0->indexCount = mesh->m_numIndices;
			submesh0->baseVertex = 0;
			submesh0->vertexCount = mesh->m_numVertices;
		}
		mesh->m_topology = Topology::TriangleList;
		return ALL_OK;
	}

	ERet CreateScene( const Settings& settings )
	{
		mxPROFILE_SCOPE("CreateScene()");

		m_sceneData = new Clump();

		m_scene.Setup( settings );

		{
			rxGlobalLight* globalLight;
			mxDO(m_sceneData->New(globalLight));
			globalLight->m_color = Float3_Set(0.2,0.4,0.3);
		}

		//{
			Float3 testTri[3];
			testTri[0] = Float3_Set(0,0,0);
			testTri[1] = Float3_Set(10,0,0);
			testTri[2] = Float3_Set(10,10,0);

			Float3 testPt = Float3_Set(3,3,-50);
			Float3 closestPt = ClosestPointOnTriangle(testTri, testPt);
			float len = Float3_Length(testPt - closestPt);
		//}


#define USE_MESH_SDF (1)


		{
			const char* path_to_source_assets = gUserINI->GetString("path_to_source_assets");
			String256 pathToMesh;
			Str::CopyS(pathToMesh, path_to_source_assets);
			Str::NormalizePath(pathToMesh);
//Str::AppendS(pathToMesh, "meshes/sm_teapot.mesh");
Str::AppendS(pathToMesh, "meshes/cube.mesh");

			mxDO(SON::LoadFromFile(pathToMesh.c_str(), meshData));

			LogStream(LL_Info) << "Loaded mesh: AABB: " << meshData.bounds;

			// Translate and scale the mesh to fit the scene's AABB.
			{
				AABB24 sceneAABB = GetSceneBounds();
				AABB_ScaleSize( &sceneAABB, 0.6f );

				AABB24 meshAABB = meshData.bounds;
				Float3 center = AABB_Center(meshAABB);
				Float3 radius = AABB_Extent(meshAABB);

				bool preserveRelativeScale = false;

				const Float3 offset = AABB_Center( sceneAABB ) - center;
				const Float3 scale = Float3_Multiply( AABB_Extent( sceneAABB ), Float3_Reciprocal( radius ) );

				for( int iSubMesh = 0; iSubMesh < meshData.sets.Num(); iSubMesh++ )
				{
					TcTriMesh& subMesh = meshData.sets[ iSubMesh ];

					for( int iVertex = 0; iVertex < subMesh.positions.Num(); iVertex++ )
					{
						subMesh.positions[iVertex] = Float3_Multiply( subMesh.positions[iVertex] + offset, scale);
					}
				}

#if !USE_MESH_SDF
				VertexDescription vertexDesc;
				DrawVertex::BuildVertexDescription(vertexDesc);

				RawMeshData rawMeshData;
				mxDO(Meshok::CompileMesh(meshData, vertexDesc, rawMeshData));

				rxMesh* mesh = m_sceneData->New< rxMesh >();
				mxDO(mesh->Create(rawMeshData));

				rxModel* model = rxModel::Create( mesh, m_sceneData );

				rxMaterial*	defaultMaterial = GetAsset< rxMaterial >( MakeAssetID("plain_color.material"), m_sceneData );
				model->m_batches.Add(defaultMaterial);
#endif
			}

			ProcessMeshDataTriangles triMeshCallback(meshData);
			meshBspTree.Build(&triMeshCallback);

			m_bsp_sdf.tree = &meshBspTree;

			meshSDF.meshData = &meshData;
			meshSDF.meshTree = &meshBspTree;
		}

		m_pointPosition = Float3_Zero();
		m_closestPointOnMesh = Float3_Zero();


#if !USE_MESH_SDF
		AVolume* root = m_scene.GetRoot();
		//HalfSpaceSDF* root = new HalfSpaceSDF();
		//root->plane = Float4_Set(Float3_Set(0,0,1),-10);
		m_garbage.Add( root );

		{
			ASDF::Options options;
			options.radius = settings.dc.radius;
			options.max_depth = settings.dc.maxDepth;
			options.min_subdiv = settings.dc.minSubdiv;
			m_asdf_tree.Build(root, options);

			m_asdf_volume.tree = &m_asdf_tree;
			m_asdf_volume.radius = settings.dc.radius;
			root = &m_asdf_volume;
		}
#else
		//AVolume* root = &meshSDF;
		AVolume* root = &m_bsp_sdf;
		m_garbage.Add( root );
#endif


		bool useOctree = 1;

		const UINT32 startTimeMSec = mxGetTimeInMilliseconds();
		DualContouring::OctStats stats;

		if(useOctree)
		{
			m_octree.Build(root, &m_qefSolver, stats, settings.dc);
			stats.Print(mxGetTimeInMilliseconds() - startTimeMSec);
		}

		{
			mxDO(m_sceneData->New( m_mesh ));

			if(useOctree) {
				mxDO(m_octree.Triangulate(m_meshBuffer,settings.dc,stats));
				mxDO(CreateMesh(m_meshBuffer, m_mesh));
			}
			else {
				GenerateMesh(root, m_mesh);
			}

			rxModel* model = rxModel::Create( m_mesh, m_sceneData );

			rxMaterial*	defaultMaterial = GetAsset< rxMaterial >( MakeAssetID("plain_color.material"), m_sceneData );
			model->m_batches.Add(defaultMaterial);
		}

		const UINT32 currentTimeMSec = mxGetTimeInMilliseconds();
		stats.Print(currentTimeMSec - startTimeMSec);

		return ALL_OK;
	}

	ERet Polygonize( const DC::Octree_DC2& tree, rxMesh* mesh, DualContouring::OctStats &stats )
	{
		mxDO(m_octree.Triangulate(m_meshBuffer,g_settings.dc,stats));
		mxDO(CreateMesh(m_meshBuffer, m_mesh));
		return ALL_OK;
	}

	ERet LoadSceneFromCache()
	{
		mxPROFILE_SCOPE("LoadSceneFromCache()");
		UNDONE;
		return ALL_OK;
	}

	void Update(const double deltaSeconds) override
	{
		mxPROFILE_SCOPE("Update()");

		if( g_rebuild_octree )
		{
			DualContouring::OctStats stats;
			m_octree.Triangulate(m_meshBuffer,g_settings.dc,stats);
			CreateMesh(m_meshBuffer, m_mesh);
			g_rebuild_octree = false;
		}
		if( g_regenerate_mesh )
		{
			//
		}

		DemoApp::Update(deltaSeconds);

		camera.Update(deltaSeconds);

		m_gui_renderer.Update();

		//dev_cast_ray(GameActions::dev_cast_ray,EInputState::IS_Pressed,0);
	}

	ERet Render()
	{
		mxPROFILE_SCOPE("Render()");

		OctCubeF worldCube = GetSceneCube();

		const Float4x4 viewMatrix = camera.BuildViewMatrix();

		int screenWidth, screenHeight;
		WindowsDriver::GetWindowSize(&screenWidth, &screenHeight);

		const float aspect_ratio = (float)screenWidth / (float)screenHeight;
		const float FoV_Y_degrees = DEG2RAD(90);
		const float near_clip	= 10.0f;
		const float far_clip	= 5000.0f;

		const Float4x4 projectionMatrix = Matrix_Perspective(FoV_Y_degrees, aspect_ratio, near_clip, far_clip);
		const Float4x4 viewProjectionMatrix = Matrix_Multiply(viewMatrix, projectionMatrix);
		const Float4x4 worldMatrix = Matrix_Identity();
		const Float4x4 worldViewProjectionMatrix = Matrix_Multiply(worldMatrix, viewProjectionMatrix);


		SceneView view;
		view.viewMatrix = viewMatrix;
		view.projectionMatrix = projectionMatrix;
		view.viewProjectionMatrix = viewProjectionMatrix;
		view.worldSpaceCameraPos = camera.m_eyePosition;
		view.viewportWidth = screenWidth;
		view.viewportHeight = screenHeight;
		view.nearClip = near_clip;
		view.farClip = far_clip;

		mxDO(m_render_system.RenderScene(view, *m_sceneData));

		if(g_settings.drawWireframeMesh)
		{
			gfxMARKER(Wireframe);

			llgl::ViewState	viewState;
			{
				viewState.Reset();
				viewState.colorTargets[0].SetDefault();
				viewState.targetCount = 1;
				viewState.depthTarget.SetDefault();
			}
			llgl::SubmitView(renderContext, viewState);

			FxShader* shader;
			mxDO(GetAsset(shader,MakeAssetID("solid_wireframe.shader"),m_sceneData));
			{
				const SolidWireSettings& settings = g_settings.solidWire;
				mxDO(FxSlow_SetUniform(shader,"LineWidth",&settings.LineWidth));
				mxDO(FxSlow_SetUniform(shader,"FadeDistance",&settings.FadeDistance));
				mxDO(FxSlow_SetUniform(shader,"PatternPeriod",&settings.PatternPeriod));
				mxDO(FxSlow_SetUniform(shader,"WireColor",&settings.WireColor));
				mxDO(FxSlow_SetUniform(shader,"PatternColor",&settings.PatternColor));
			}
			mxDO(FxSlow_Commit(renderContext,shader));

			mxDO(FxSlow_SetRenderState(renderContext, *m_rendererData, "SolidWireframe"));
			m_render_system.DBG_Draw_Models_Wireframe(view, *m_sceneData, *shader);
		}
		if(g_settings.drawVertexNormals)
		{
			FxShader* shader;
			mxDO(GetAsset(shader,MakeAssetID("debug_draw_normals.shader"),m_sceneData));
			{
				float lineLength = 100;
				mxDO(FxSlow_SetUniform(shader,"g_lineLength",&lineLength));
			}
			mxDO(FxSlow_Commit(renderContext,shader));

			mxDO(FxSlow_SetRenderState(renderContext, *m_rendererData, "Default"));
			m_render_system.DBG_Draw_Models_With_Custom_Shader(view, *m_sceneData, *shader, Topology::PointList);
		}

		{
			gfxMARKER(DebugDrawOctree);
			mxDO(FxSlow_SetRenderState(renderContext, *m_rendererData, "NoCulling"));

			FxShader* shader;
			mxDO(GetByName(*m_rendererData, "auxilary_rendering", shader));
			m_batch_renderer.SetShader(shader);
			mxDO(FxSlow_SetUniform(shader,"u_modelViewProjMatrix",&worldViewProjectionMatrix));
			mxDO(FxSlow_SetUniform(shader,"u_modelMatrix",&worldMatrix));
			mxDO(FxSlow_SetUniform(shader,"u_color",RGBAf::WHITE.ToPtr()));
			mxDO(FxSlow_Commit(renderContext,shader));



			m_batch_renderer.DrawLine(m_pointPosition, m_closestPointOnMesh, RGBAf::RED, RGBAf::RED);



			//OctCubeF octants[8];
			//GetChildOctants(worldCube,octants);
			//for( int i = 0; i < 8; i++ )
			//{
			//	AABB24 aabb;
			//	OctBoundsToAABB(octants[i], aabb);
			//	m_batch_renderer.DrawAABB(aabb.min_point, aabb.max_point, RGBAf::lightblue);
			//}

			if(g_settings.visualize_octree)
			{
				DebugDrawOctree(0,worldCube,0,0,m_batch_renderer);
			}

			// 'crosshair'
			m_batch_renderer.DrawCircle(
				camera.m_eyePosition + camera.m_lookDirection * 100,
				camera.m_rightDirection,
				camera.m_upDirection,
				RGBAf::RED,
				1, 8
			);

			if( m_lastRayCastResult.hitAnything )
			{
				m_batch_renderer.DrawLine( m_lastRayCastStart, m_lastRayCastResult.rayExit, RGBAf::WHITE, RGBAf::RED );
				m_batch_renderer.DrawCircle(m_lastRayCastResult.rayEntry, camera.m_rightDirection, camera.m_upDirection, RGBAf::YELLOW, 10, 8);
				m_batch_renderer.DrawCircle(m_lastRayCastResult.rayExit, camera.m_rightDirection, camera.m_upDirection, RGBAf::GREEN, 10, 8);
				//DrawHitLeaves( 0, worldBounds, m_batch_renderer );
				for( int i = 0; i < m_hitLeaves.Num(); i++ )
				{
					const DebugLeafInfo& leafInfo = m_hitLeaves[i];
					m_batch_renderer.DrawAABB(leafInfo.aabb.min_point, leafInfo.aabb.max_point, RGBAf::RED);
					m_batch_renderer.DrawCircle( leafInfo.xyz, camera.m_rightDirection, camera.m_upDirection, RGBAf::ORANGE, 20, 4 );
					m_batch_renderer.DrawLine( leafInfo.xyz, leafInfo.xyz + leafInfo.N * 100, RGBAf::WHITE, RGBAf::RED );
				}
			}
		}

		{
			gfxMARKER(Gizmos);
			mxDO(FxSlow_SetRenderState(renderContext, *m_rendererData, "NoCulling"));

			FxShader* shader;
			mxDO(GetByName(*m_rendererData, "auxilary_rendering", shader));
			m_batch_renderer.SetShader(shader);
			mxDO(FxSlow_SetUniform(shader,"u_modelViewProjMatrix",&worldViewProjectionMatrix));
			mxDO(FxSlow_SetUniform(shader,"u_modelMatrix",&worldMatrix));
			mxDO(FxSlow_SetUniform(shader,"u_color",RGBAf::WHITE.ToPtr()));
			mxDO(FxSlow_Commit(renderContext,shader));

			DrawGizmo(gizmo, worldMatrix, camera.m_eyePosition, m_batch_renderer);

			const AABB24 worldAabb = GetSceneBounds();
			m_batch_renderer.DrawAABB(worldAabb.min_point, worldAabb.max_point, RGBAf::LIGHT_YELLOW_GREEN);

			m_batch_renderer.Flush();
		}

		m_batch_renderer.Flush();

		if(m_debugShowNormalsRT)
		{
			{
				mxDO(FxSlow_SetRenderState(renderContext, *m_rendererData, "NoCulling"));
				FxShader* shader;
				mxDO(GetByName(*m_rendererData, "full_screen_triangle", shader));
				mxDO(FxSlow_SetResource(shader, "t_sourceTexture", llgl::AsResource(m_render_system.m_colorRT1->handle), Rendering::g_samplers[PointSampler]));
				mxDO(FxSlow_Commit(renderContext,shader));
				mxDO(m_render_system.DrawFullScreenTriangle(shader));
			}
		}


		m_gui_renderer.Start_Rendering();
		{
			mxDO(FxSlow_SetRenderState(renderContext, *m_rendererData, "NoCulling"));

			if( m_lastRayCastResult.hitAnything && m_hitLeaves.Num() > 0 )
			{
#if 0
				std::sort( m_hitLeaves.ToPtr(), m_hitLeaves.ToPtr() + m_hitLeaves.Num(), CompareDepth );

				for( int i = 0; i < m_hitLeaves.Num(); i++ )
				{
					const DebugLeafInfo& leafInfo = m_hitLeaves[i];

					if(ImGui::Begin("Hit Leaf"))
					{
						ImGui::SetWindowPos(ImGui_GetWindowPosition(leafInfo.xyz,view));
						ImGui_Window_Text_Float3(leafInfo.xyz);
						const float dist = Float3_LengthSquared( leafInfo.xyz - Get().camera.m_eyePosition );
						ImGui::Text("dist=%f",dist);
						ImGui::End();
					}
				}
#endif
			}

			AClientState* currentState = stateMgr.GetCurrentState();

			{
				ImGui::SetNextWindowPos(ImVec2(0, 0), ImGuiSetCond_FirstUseEver);
				static bool opened = true;
				if( ImGui::Begin("Stats",&opened) )
				{
					// Calculate and show frame rate
					static float ms_per_frame[120] = { 0 };
					static int ms_per_frame_idx = 0;
					static float ms_per_frame_accum = 0.0f;
					ms_per_frame_accum -= ms_per_frame[ms_per_frame_idx];
					ms_per_frame[ms_per_frame_idx] = ImGui::GetIO().DeltaTime * 1000.0f;
					ms_per_frame_accum += ms_per_frame[ms_per_frame_idx];
					ms_per_frame_idx = (ms_per_frame_idx + 1) % 120;
					const float ms_per_frame_avg = ms_per_frame_accum / 120;
					ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", ms_per_frame_avg, 1000.0f / ms_per_frame_avg);

					// Octree/Voxels/Memory
					{
						const int num_nodes = m_octree.m_nodes.NumValidItems();
						const int num_leaves = m_octree.m_leaves.NumValidItems();
						const int nodes_used_bytes = num_nodes * sizeof(m_octree.m_nodes[0]);
						const int leaves_used_bytes = num_leaves * sizeof(m_octree.m_leaves[0]);
						const int nodes_alloc_bytes = m_octree.m_nodes.GetDataSize();
						const int leaves_alloc_bytes = m_octree.m_leaves.GetDataSize();
						ImGui::Text("Bytes used: %d (nodes: %d(%d), leaves: %d(%d))",
							nodes_used_bytes+leaves_used_bytes, nodes_used_bytes,num_nodes, leaves_used_bytes,num_leaves);
						ImGui::Text("Bytes allocated: %d (nodes: %d, leaves: %d)",
							nodes_alloc_bytes+leaves_alloc_bytes, nodes_alloc_bytes, leaves_alloc_bytes);
						ImGui::Text("QEFs: %d (%d bytes)", m_octree.m_QEFs.NumValidItems(), m_octree.m_QEFs.NumValidItems() * sizeof(m_octree.m_QEFs[0]));
					}
					// Mesh/Polygons/Memory
					{
						const int num_vertices = m_meshBuffer.vertices.Num();
						const int num_indices = m_meshBuffer.indices.Num();
						const int vertices_used_bytes = num_vertices * sizeof(m_meshBuffer.vertices[0]);
						const int indices_used_bytes = num_indices * sizeof(m_meshBuffer.indices[0]);
						ImGui::Text("Vertices: %d (%d bytes, %d quads)", num_vertices, vertices_used_bytes, num_vertices/4);
						ImGui::Text("Indices: %d (%d bytes)", num_indices, indices_used_bytes);
					}

					String256 temp;
					const int stateMask = currentState->GetStateMask();
					Dbg_FlagsToString(&stateMask,GetTypeOf_GameStateF(),temp);
					ImGui::Text("Current State: %s (mask=%s)", typeid(*currentState).name(), temp.c_str());

					{
						String256 tmp;
						StringStream(tmp) << camera.m_eyePosition;
						ImGui::Text("Camera Pos: %s", tmp.c_str());
					}
					ImGui::End();
				}				
			}

			//if(ImGui::Begin("G-Buffer"))
			//{
			//ImGui::Checkbox("Show normals",&m_debugShowNormalsRT);
			//ImGui::End();
			//}

			currentState->RenderGUI();
		}
		m_gui_renderer.Finish_Rendering();

		return ALL_OK;
	}

	static bool CompareDepth( const DebugLeafInfo& a, const DebugLeafInfo& b )
	{
		const float distA = Float3_LengthSquared( a.xyz - Get().camera.m_eyePosition );
		const float distB = Float3_LengthSquared( b.xyz - Get().camera.m_eyePosition );
		return distA > distB;
	}

	void DrawHitLeaves( DualContouring::NodeID nodeID, const OctCubeF& bounds, BatchRenderer & renderer )
	{
		const bool isLeaf = DualContouring::IS_LEAF_ID(nodeID);
		const UINT16 nodeIndex = DualContouring::GET_ID(nodeID);
	
		AABB24	aabb;
		OctBoundsToAABB(bounds, aabb);

		if( isLeaf )
		{
			const DualContouring::Leaf2& leaf = m_octree.m_leaves[nodeIndex];

			if( m_lastRayCastResult.hitLeaves.Contains(nodeIndex) )
			{
				renderer.DrawAABB(aabb.min_point, aabb.max_point, RGBAf::RED);

				Float3 xyz;
				DequantizePosition( bounds, leaf.xyz, &xyz );
				m_batch_renderer.DrawCircle( xyz, camera.m_rightDirection, camera.m_upDirection, RGBAf::ORANGE, 20, 4 );

				Float3 N = DequantizeNormal( leaf.N );
				m_batch_renderer.DrawLine( xyz, xyz + N * 100, RGBAf::WHITE, RGBAf::RED );
			}
		}
		else//if(!isLeaf)
		{
			OctCubeF octants[8];
			GetChildOctants(bounds,octants);

			const DualContouring::Node2& node = m_octree.m_nodes[nodeIndex];
			for( int i = 0; i < 8; i++ )
			{
				if( node.kids[i] != DualContouring::NIL_NODE ) {
 					DrawHitLeaves( node.kids[i], octants[i], renderer );
				}
			}
		}
	}

	void DebugDrawOctree( DualContouring::NodeID nodeID, const OctCubeF& bounds, const UINT32 code, int depth, BatchRenderer & renderer )
	{
		static const UINT32 s_octTreeLevelColors[8] = {
			SColor::RED,
			SColor::ORANGE,
			SColor::YELLOW,
			SColor::GREEN,
			SColor::LIGHTBLUE,
			SColor::BLUE,
			SColor::VIOLET,
			SColor::WHITE,
		};

		int colorIndex = Clamp<int>(depth,0,mxCOUNT_OF(s_octTreeLevelColors));
		const SColor argb(s_octTreeLevelColors[colorIndex]);
		const RGBAf color(argb);


		const bool isLeaf = DualContouring::IS_LEAF_ID(nodeID);
		const UINT32 nodeIndex = DualContouring::GET_ID(nodeID);


		AABB24	aabb;
		OctBoundsToAABB(bounds, aabb);


#if 1
		if( g_settings.draw_octree_leaves_only )
		{
			if( isLeaf ) {
				renderer.DrawAABB(aabb.min_point, aabb.max_point, color);
			}
		}
		else
		{
			renderer.DrawAABB(aabb.min_point, aabb.max_point, color);
		}
#endif


		if(isLeaf)
		{
			//const DualContouring::Leaf2& leaf = m_octree.m_leaves[nodeIndex];
		}
		else//if(!isLeaf)
		{
			OctCubeF octants[8];
			GetChildOctants(bounds,octants);

			const DualContouring::Node2& node = m_octree.m_nodes[nodeIndex];
			for( int i = 0; i < 8; i++ )
			{
				if( node.kids[i] != DualContouring::NIL_NODE ) {
 					DebugDrawOctree( node.kids[i], octants[i], (code<<3)|i, depth+1, renderer );
				}
			}
		}
	}

	static const int CSG_RADIUS = 600;

	void CSG_Subtract( const Float3& position )
	{
		//DC::DBG_CSG_Subtract( m_octree, firstHit.xyz, 100 );

		const UINT32 startTimeMSec = mxGetTimeInMilliseconds();

#if 1
		SphereSDF* sphere = new SphereSDF();
		sphere->center = position;
		sphere->radius = CSG_RADIUS;
#else
		SphereSDF* sphere = new SphereSDF();
		sphere->center = position;
		sphere->radius = 500;
#endif

		CSGSubtraction* subtract = new CSGSubtraction();
		subtract->opA = m_garbage.GetLast();
		subtract->opB = sphere;

		m_garbage.Add( sphere );
		m_garbage.Add( subtract );

		DualContouring::OctStats stats;
		m_octree.Build( subtract, &m_qefSolver, stats, g_settings.dc );
		m_octree.Triangulate(m_meshBuffer,g_settings.dc,stats);
		CreateMesh(m_meshBuffer, m_mesh);

		const UINT32 currentTimeMSec = mxGetTimeInMilliseconds();
		stats.Print(currentTimeMSec - startTimeMSec);
	}
	void CSG_Add( const Float3& position )
	{
		const UINT32 startTimeMSec = mxGetTimeInMilliseconds();

		SphereSDF* sphere = new SphereSDF();
		sphere->center = position;
		sphere->radius = CSG_RADIUS;

		CSGUnion* csgOp = new CSGUnion();
		csgOp->opA = m_garbage.GetLast();
		csgOp->opB = sphere;

		m_garbage.Add( sphere );
		m_garbage.Add( csgOp );

		DualContouring::OctStats stats;
		m_octree.Build( csgOp, &m_qefSolver, stats, g_settings.dc );
		m_octree.Triangulate(m_meshBuffer,g_settings.dc,stats);
		CreateMesh(m_meshBuffer, m_mesh);

		const UINT32 currentTimeMSec = mxGetTimeInMilliseconds();
		stats.Print(currentTimeMSec - startTimeMSec);
	}

	void attack1( GameActionID action, EInputState status, float value )
	{
		m_pointPosition = camera.m_eyePosition;
		m_closestPointOnMesh = meshSDF.GetClosestPointTo(m_pointPosition);
		LogStream(LL_Debug) << "Closest point: " << m_closestPointOnMesh << ", dist: " << Float3_Length(m_pointPosition-m_closestPointOnMesh);
		return;

		RayCastResult result;
		m_octree.CastRay(camera.m_eyePosition, camera.m_lookDirection, result, g_settings.dc);
		if( result.hitLeaves.Num() > 0 )
		{
			TLocalArray< DebugLeafInfo, 32 >	hitLeaves;
			DC::DebugGatherHitLeaves( result, m_octree, 0, GetSceneCube(), hitLeaves );
			const DebugLeafInfo& firstHit = hitLeaves.GetFirst();
			//int leafIndex = result.hitLeaves.GetFirst();
			CSG_Subtract( firstHit.xyz );
		}
	}
	void attack2( GameActionID action, EInputState status, float value )
	{
		RayCastResult result;
		m_octree.CastRay(camera.m_eyePosition, camera.m_lookDirection, result, g_settings.dc);
		if( result.hitLeaves.Num() > 0 )
		{
			TLocalArray< DebugLeafInfo, 32 >	hitLeaves;
			DC::DebugGatherHitLeaves( result, m_octree, 0, GetSceneCube(), hitLeaves );
			const DebugLeafInfo& firstHit = hitLeaves.GetFirst();
			//int leafIndex = result.hitLeaves.GetFirst();
			CSG_Add( firstHit.xyz );
		}
	}
	void dev_cast_ray( GameActionID action, EInputState status, float value )
	{
		//DBGOUT("Cast ray");
		m_lastRayCastStart = camera.m_eyePosition;
		m_lastRayCastDirection = camera.m_lookDirection;
		m_octree.CastRay(camera.m_eyePosition, camera.m_lookDirection, m_lastRayCastResult, g_settings.dc);

		m_hitLeaves.Empty();
		DC::DebugGatherHitLeaves( m_lastRayCastResult, m_octree, 0, GetSceneCube(), m_hitLeaves );

		if( m_hitLeaves.Num() > 0 )
		{
			DBGOUT("%d hits",m_hitLeaves.Num());

			const DebugLeafInfo& firstHit = m_hitLeaves.GetFirst();
			LogStream(LL_Debug) << "firstHit.pos=" << firstHit.xyz;

			for( int i = 0; i < m_hitLeaves.Num(); i++ )
			{
				const DebugLeafInfo& leafInfo = m_hitLeaves[i];
				const DC::Leaf2& leaf = m_octree.m_leaves[ leafInfo.leafIndex ];
				LogStream(LL_Debug) << "Leaf[" << leafInfo.leafIndex << "]: pos=" << leafInfo.xyz << "; quantized: " << (int)leaf.xyz[0] << "," << (int)leaf.xyz[1] << "," << (int)leaf.xyz[2];
			}

			CSG_Subtract( firstHit.xyz );
		}
		else
		{
			DBGOUT("no intersection");
		}
	}
	void dev_toggle_menu( GameActionID action, EInputState status, float value ) override
	{
		ChangeGameState( &m_menuState );
	}
	void dev_reset_state( GameActionID action, EInputState status, float value ) override
	{
		camera.m_eyePosition = Float3_Set(300,300,300);
	}
	void dev_save_state( GameActionID action, EInputState status, float value ) override
	{
		//SON::SaveToFile(g_settings,SETTINGS_FILE_NAME);
		FileWriter file;
		if(mxSUCCEDED(file.Open(OCTREE_SAVE_FILE_NAME))) {
			m_octree.Save(file);
		}		
	}
	void dev_load_state( GameActionID action, EInputState status, float value ) override
	{
		//SON::LoadFromFile(SETTINGS_FILE_NAME,g_settings);
		FileReader file;
		if(mxSUCCEDED(file.Open(OCTREE_SAVE_FILE_NAME))) {
			m_octree.Load(file);
			DualContouring::OctStats stats;
			m_octree.Triangulate(m_meshBuffer,g_settings.dc,stats);
			CreateMesh(m_meshBuffer, m_mesh);
		}	
	}
	void take_screenshot( GameActionID action, EInputState status, float value )
	{
		CalendarTime	localTime( CalendarTime::GetCurrentLocalTime() );

		String64	timeOfDay;
		GetTimeOfDayString( timeOfDay, localTime.hour, localTime.minute, localTime.second );

		String128	timeStamp;
		Str::SPrintF(timeStamp, "%d.%d.%d - %s",
			localTime.year, localTime.month, localTime.day, timeOfDay.c_str());

		String256	fileName("screenshots/");
		Str::Append(fileName, timeStamp);
		Str::AppendS(fileName, ".jpg");

		DBGOUT("Saving screenshot to file '%s'", fileName.c_str());
		llgl::SaveScreenshot( fileName.c_str() );
	}
};

}//namespace ContouringDemo

int mxAppMain()
{
	SetupBaseUtil	setupBase;
	FileLogUtil		fileLog;

	SON::TextConfig	engineINI;
	engineINI.LoadFromFile("EngineConfig.son");
	gINI = &engineINI;
	
	SON::TextConfig	demoINI;
	demoINI.LoadFromFile(CONFIG_FILE_NAME);
	gUserINI = &demoINI;

	SetupCoreUtil	setupCore;


	ContouringDemo::App	app;

	FlyingCamera& camera = app.camera;
	camera.m_eyePosition = Float3_Set(0, -100, 20);
	camera.m_movementSpeed = 100;
	camera.m_strafingSpeed = 100;
//camera.m_rotationSpeed = 0.03;
camera.m_rotationSpeed = 2;
	camera.m_invertYaw = true;
	camera.m_invertPitch = true;

	int screen_width = 800, screen_height = 600;
	gUserINI->GetInteger("screen_width",screen_width,800,4096);
	gUserINI->GetInteger("screen_height",screen_height,600,4096);

	WindowsDriver::Settings cInfo;
	cInfo.screenWidth = screen_width;
	cInfo.screenHeight = screen_height;
	mxDO(WindowsDriver::Initialize(cInfo));

	WindowsDriver::SetRelativeMouseMode(true);

	llgl::Settings	graphicsOptions;
	graphicsOptions.window = WindowsDriver::GetNativeWindowHandle();
	mxDO(llgl::Initialize(graphicsOptions));

	mxDO(app.Setup());

	{
		ClientCallbacks::Data& input = ClientCallbacks::GetData();

		input.handlers.Add( ActionBindingT(GameActions::move_forward, mxBIND_MEMBER_FUNCTION(FlyingCamera,move_forward,camera)) );
		input.handlers.Add( ActionBindingT(GameActions::move_back, mxBIND_MEMBER_FUNCTION(FlyingCamera,move_back,camera)) );
		input.handlers.Add( ActionBindingT(GameActions::move_left, mxBIND_MEMBER_FUNCTION(FlyingCamera,strafe_left,camera)) );
		input.handlers.Add( ActionBindingT(GameActions::move_right, mxBIND_MEMBER_FUNCTION(FlyingCamera,strafe_right,camera)) );
		input.handlers.Add( ActionBindingT(GameActions::move_up, mxBIND_MEMBER_FUNCTION(FlyingCamera,move_up,camera)) );
		input.handlers.Add( ActionBindingT(GameActions::move_down, mxBIND_MEMBER_FUNCTION(FlyingCamera,move_down,camera)) );
		input.handlers.Add( ActionBindingT(GameActions::rotate_pitch, mxBIND_MEMBER_FUNCTION(FlyingCamera,rotate_pitch,camera)) );
		input.handlers.Add( ActionBindingT(GameActions::rotate_yaw, mxBIND_MEMBER_FUNCTION(FlyingCamera,rotate_yaw,camera)) );

		input.handlers.Add( ActionBindingT(GameActions::dev_cast_ray, mxBIND_MEMBER_FUNCTION(ContouringDemo::App,dev_cast_ray,app)) );
	}

	mxDO(app.CreateScene( ContouringDemo::g_settings ));

	SON::LoadFromFile(gUserINI->GetString("path_to_camera_state"),camera);


	UINT64	prevTimeMicroseconds = mxGetTimeInMicroseconds();

	const double updateFrequency = 60;	// Hz
	const double timeStepMs = 1000.f / updateFrequency;
	double timeAccumulatedMs = 0;

	while( !WindowsDriver::IsAboutToExit() )
	{
		mxPROFILE_SCOPE("MainLoop");

		UINT64	currTimeMicroseconds = mxGetTimeInMicroseconds();
		UINT64	deltaTimeMicroseconds = currTimeMicroseconds - prevTimeMicroseconds;
		prevTimeMicroseconds = currTimeMicroseconds;

		timeAccumulatedMs += deltaTimeMicroseconds * 1e-3;



		while( timeAccumulatedMs >= timeStepMs )
		{
			timeAccumulatedMs -= timeStepMs;

			//const bool bHasFocus = WindowsDriver::HasFocus();
			AClientState* currentState = app.stateMgr.GetCurrentState();
			{
				int stateMask = currentState->GetStateMask();
				WindowsDriver::GenerateSystemEvents(stateMask);
				WindowsDriver::ProcessSystemEvents(stateMask);
			}
			//const double deltaTimeSeconds = WindowsDriver::GetDeltaSeconds();
			//const double deltaTimeSeconds = deltaTimeMicroseconds * 1e-6;
			const double deltaTimeSeconds = timeStepMs * 1e-3;

			app.Update(deltaTimeSeconds);

			app.Render();

			llgl::NextFrame();

			mxSleepMilliseconds(5);

			mxPROFILE_INCREMENT_FRAME_COUNTER;
		}
	}

	SON::SaveToFile(camera,gUserINI->GetString("path_to_camera_state"));

	app.Shutdown();

	llgl::Shutdown();

	WindowsDriver::Shutdown();

	return 0;
}
mxAPPLICATION_ENTRY_POINT;
