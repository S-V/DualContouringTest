#pragma once

#include <Core/VectorMath.h>
#include <Renderer/Vertex.h>
#include <QEF.h>

namespace Contouring
{

struct QEF_Solver
{
	enum { MAX_POINTS = 16 };

	struct Input
	{
		AABB24 bounds;	// can be used for clamping the solution
		Float3 positions[MAX_POINTS];	// intersection points
		Float3 normals[MAX_POINTS];		// intersection normals
		int numPoints;		// number of intersections
		float threshold;	// error threshold
		int maxIterations;
	public:
		Input();
	};
	struct Output
	{
		Float3 position;// position with least error
		float error;	// always positive (least-squares error)
		svd::QefData qef;
	};
	virtual void Solve( const Input& input, Output &output ) = 0;

	//virtual void const svd::QefData& GetQefData() { return svd::QefData(); }
};

struct Options : CStruct
{
	float	radius;	// half-size of the entire world cube
	UINT32	maxDepth;
	float	minSubdiv;
	float	threshold;	// error threshold
	float	qef_threshold;	// threshold for QEF-based simplification (collapse nodes if error is less than this value)
	float	build_threshold;// error threshold to use for building
public:
	mxDECLARE_CLASS(Options, CStruct);
	mxDECLARE_REFLECTION;
	Options();
	void SetDefauls();
};

class OctStats {
public:
	UINT32		numPolygons;	// number of resulting polygons
	UINT32		numLeafNodes;	// number of leaf ('gray', surface) nodes
	UINT32		numEmptyLeaves;	// number of leaves fully inside/outside the surface
	UINT32		numInternalNodes;// number of internal ('glue') nodes
	UINT32		numBadNodes;	// number of collapsed empty nodes
	UINT32		maxTreeDepth;
	UINT32		bytesAllocated;

	// maximum number of zero-crossing edges in a cell
	UINT32		maxActiveEdges;

	// maximum QEF error
	float		max_QEF_error;

	UINT32		construction_time_milliseconds;
	UINT32		contouring_time_milliseconds;

	// nodes merged during QEF-based simplification
	UINT32		nodes_collapsed;

public:
	OctStats();
	void Print( UINT32 elapsedTimeMSec );
};

// Uses Leonardo Augusto Schmitz's easy-to-understand method
// with exact normals at intersection points to reduce complexity.
//
struct QEF_Solver_ParticleBased : QEF_Solver
{
	virtual void Solve( const Input& input, Output &output ) override;

	static const int MAX_ITERATIONS = 100;
};

struct QEF_Solver_SVD : QEF_Solver
{
	virtual void Solve( const Input& input, Output &output ) override;
};

struct RayCastResult
{
	Float3	rayEntry;	//for debugging
	Float3	rayExit;	//for debugging
	bool hitAnything;
	TLocalArray< UINT16, 32 > hitLeaves;	// in the same order the ray passes through them

	RayCastResult()
	{
		hitAnything = false;
	}
	PREVENT_COPY(RayCastResult);
};

inline void QuantizePosition( const OctCubeF& _bounds, const Float3& xyz, UINT8 output[3] )
{
	AABB24 aabb;
	OctBoundsToAABB(_bounds, aabb);
	AABB_EncodePosition(aabb, xyz, output);
}
inline void DequantizePosition( const OctCubeF& _bounds, const UINT8 input[3], Float3 *xyz )
{
	AABB24 aabb;
	OctBoundsToAABB(_bounds, aabb);
	AABB_DecodePosition(aabb, input, xyz);
}
inline void QuantizeNormal( const Float3& normal, UINT8 output[3] )
{
	UByte4 packed = PackNormal( normal );
	output[0] = packed.x;
	output[1] = packed.y;
	output[2] = packed.z;
}
inline Float3 DequantizeNormal( const UINT8 input[3] )
{
	UByte4 packed = { input[0], input[1], input[2], 0 };
	return UnpackNormal( packed );
}


struct DebugLeafInfo
{
	int leafIndex;
	AABB24 aabb;
	Float3 xyz;
	Float3 N;
};

}//namespace Contouring
