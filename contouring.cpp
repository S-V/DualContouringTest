#include "stdafx.h"
#pragma hdrstop
#include "QEF.h"
#include "contouring.h"

namespace Contouring
{

QEF_Solver::Input::Input()
{
	AABB24_Clear(&bounds);
	numPoints = 0;
	threshold = 1e-6f;
	maxIterations = ~0;
}

mxDEFINE_CLASS(Options);
mxBEGIN_REFLECTION(Options)
	mxMEMBER_FIELD( radius ),
	mxMEMBER_FIELD( maxDepth ),
	mxMEMBER_FIELD( minSubdiv ),
	mxMEMBER_FIELD( threshold ),
	mxMEMBER_FIELD( qef_threshold ),
	mxMEMBER_FIELD( build_threshold ),
mxEND_REFLECTION;
Options::Options()
{
	this->SetDefauls();

}
void Options::SetDefauls()
{
	radius = 1.0f;
	maxDepth = 16;
	minSubdiv = 1.0f;
	threshold = 1.0f;
	qef_threshold = -1.0f;	// don't simplify at all
	build_threshold = -1.0f;// build to the very last level
}

OctStats::OctStats()
{
	mxZERO_OUT( *this );
}
void OctStats::Print( UINT32 elapsedTimeMSec )
{
	DBGOUT( "\n=== Octree statistics ========" );
	DBGOUT( "Polygons:           %u", numPolygons	 );
	DBGOUT( "Bad Nodes:          %u", numBadNodes );
	DBGOUT( "Leaf Nodes:         %u", numLeafNodes );	
	DBGOUT( "Empty Leaves:       %u", numEmptyLeaves );
	DBGOUT( "Internal Nodes:     %u", numInternalNodes );
	DBGOUT( "Max. Tree Depth:    %u", maxTreeDepth );
	DBGOUT( "Max. Intersections: %u", maxActiveEdges );
	DBGOUT( "Nodes Collapsed:    %u", nodes_collapsed );
	DBGOUT( "Max. QEF error:     %f", max_QEF_error );
	DBGOUT( "Construction Time: %u msec", construction_time_milliseconds );
	DBGOUT( "Contouring Time:   %u msec", contouring_time_milliseconds );
	DBGOUT( "Total Memory Used:	%u bytes", bytesAllocated );
	//DBGOUT( "Memory used:		%u bytes (leaves: %u, nodes %u)", bytesAllocated );
	DBGOUT( "Time elapsed:      %u msec", elapsedTimeMSec		);
	DBGOUT( "==== End ====================" );
}

// Basically, you start the vertex in the centre of the cell.
// Then you average all the vectors taken from the vertex to each plane
// and move the vertex along that resultant, and repeat this step a fixed number of times.
// http://gamedev.stackexchange.com/a/83757
//
void QEF_Solver_ParticleBased::Solve( const Input& input, Output &output )
{
	// Center the particle on the masspoint.
	Float3 masspoint = Float3_Zero();
	for( int i = 0; i < input.numPoints; i++ )
	{
		masspoint += input.positions[ i ];
	}
	masspoint /= input.numPoints;

	Float3 particlePosition = masspoint;

	// find the force that starts moving the particle from the masspoint towards the iso-surface
	Float3 force;

	Float4 planes[MAX_POINTS];
	for( int i = 0; i < input.numPoints; i++ )
	{
		const Float3& planePoint = input.positions[ i ];
		const Float3& planeNormal = input.normals[ i ];
		planes[ i ] = Plane_FromPointNormal( planePoint, planeNormal );
	}

    //const float FORCE_TRESHOLD = 0.00001f;
    const float forceRatio = 0.75f;

	int iteration = 0;
	while( iteration++ < MAX_ITERATIONS )
	{
		force = Float3_Zero();

		// For each intersection point:
		for (int i = 0; i < input.numPoints; i++)
		{
			const Float3& point = input.positions[ i ];
			force += Plane_GetNormal( planes[i] ) * Plane_PointDistance( planes[i], point );
		}

		// Average the force over all the intersection points, and multiply 
		// with a ratio and some damping to avoid instabilities.
		float damping = 1.0f - ((float) iteration) / MAX_ITERATIONS;

		force *= forceRatio * damping / input.numPoints;

		// Apply the force.
		particlePosition += force;

		// If the force was almost null, break.
		if( Float3_LengthSquared(force) < squaref(input.threshold) )
		{
			break;
		}
	}

	output.position = particlePosition;
	output.error = Float3_Length(force);
}

void QEF_Solver_SVD::Solve( const Input& input, Output &output )
{
	svd::QefSolver	qef;

	for( int i = 0; i < input.numPoints; i++ )
	{
		const Float3& planePoint = input.positions[ i ];
		const Float3& planeNormal = input.normals[ i ];
		qef.add(
			planePoint.x, planePoint.y, planePoint.z,
			planeNormal.x, planeNormal.y, planeNormal.z
		);
	}

	const int QEF_SWEEPS = 4;

	svd::Vec3 qefPosition;
	const float error = qef.solve( qefPosition, input.threshold, QEF_SWEEPS, input.threshold );

	output.position = Float3_Set(qefPosition.x, qefPosition.y, qefPosition.z);
	output.error = error;
	output.qef = qef.getData();

	//NOTE: this is important!
	if( !AABB_ContainsPoint( input.bounds, output.position ) )
	{
		const svd::Vec3& mp = qef.getMassPoint();
		output.position = Float3_Set(mp.x, mp.y, mp.z);
	}
}

}//namespace Contouring
