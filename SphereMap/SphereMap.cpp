/*
Copyright (c) 2018, Michael Kazhdan, Alex Baden, and Keenan Crane
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution.

Neither the name of the Johns Hopkins University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <map>
#include "Misha/CmdLineParser.h"
#include "Misha/Ply.h"
#include "Misha/SparseMatrix.h"
#include "Misha/HalfEdge.h"
#include "Misha/Solver.h"
#include "Misha/Timer.h"
#include "Misha/MemoryUsage.h"
#include "Misha/MeshStuff.h"
#include "Misha/SphericalGeometry.h"

cmdLineParameter< char* > In( "in" ) , Out( "out" ) , OutGrid( "outG" ) , OutTessellation( "outT" );
cmdLineParameter< int > Iterations( "iters" , 100 ) , Threads( "threads" , omp_get_num_procs() ) , Resolution( "res" , 256 );
cmdLineParameter< float > StepSize( "stepSize" , 0.1f ) , CutOff( "cutOff" , 1e-10f ) , Smooth( "smooth" , 5e-4f );
cmdLineReadable Verbose( "verbose" ) , ASCII( "ascii" ) , Randomize( "random" ) , NoCenter( "noCenter" ) , Collapse( "collapse" );

void Usage( const char* ex )
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input mesh>\n" , In.name );
	printf( "\t[--%s <output mesh>]\n" , Out.name );
	printf( "\t[--%s <output spherical grid>]\n" , OutGrid.name );
	printf( "\t[--%s <output spherical tessellation>]\n" , OutTessellation.name );
	printf( "\t[--%s <CMCF iterations>=%d]\n" , Iterations.name , Iterations.value );
	printf( "\t[--%s <CMCF step size>=%f]\n" , StepSize.name , StepSize.value );
	printf( "\t[--%s <Moebius centering cut-off>=%g]\n" , CutOff.name , CutOff.value );
	printf( "\t[--%s <spherical resolution>=%d]\n" , Resolution.name , Resolution.value );
	printf( "\t[--%s <spherical diffusion time>=%g]\n" , Smooth.name , Smooth.value );
	printf( "\t[--%s]\n" , Randomize.name );
	printf( "\t[--%s]\n" , NoCenter.name );
	printf( "\t[--%s]\n" , ASCII.name );
	printf( "\t[--%s]\n" , Collapse.name );
	printf( "\t[--%s]\n" , Verbose.name );
}
cmdLineReadable* params[] = { &In , &Out , &OutGrid , &OutTessellation , &Iterations , &StepSize , &Threads , & Verbose , &Randomize , &ASCII , &NoCenter , &Resolution , &Smooth , &Collapse , NULL };


template< class Real >
Point3D< Real > NormalColor( Point3D< Real > n )
{
	Point3D< Real > c = ( -n + Point3D< Real >( 1 , 1 , 1 ) ) * Real(128);
	for( int d=0 ; d<3 ; d++ )
	{
		if( c[d]>Real(255) ) c[d] = Real(255);
		if( c[d]<Real(  0) ) c[d] = Real(  0);
	}
	return c;
}

template< class Real >
void PoseMesh( const EmptyHEMesh& mesh , std::vector< Point3D< Real > >& vertices , Point3D< Real >& center , Real& scale )
{
	center = AreaCenter< Real , Point3D< Real > , EmptyHEMesh >( &vertices[0] , mesh );
	for( int i=0 ; i<(int)vertices.size() ; i++ ) vertices[i] -= center;
	scale = 0;
	for( int i=0 ; i<(int)vertices.size() ; i++ ) for( int j=0 ; j<3 ; j++ ) scale = std::max< Real >( scale , (Real)fabs(vertices[i][j]) );
	scale *= 2;
	for( int i=0 ; i<(int)vertices.size() ; i++ ) vertices[i] /= scale;
}
template< class Real >
void PoseMesh( const EmptyHEMesh& mesh , std::vector< Point3D< Real > >& vertices )
{
	Point3D< Real > center;
	Real scale;
	PoseMesh( mesh , vertices , center , scale );
}

template< class Real >
void CMCF( const EmptyHEMesh& mesh , std::vector< Point3D< Real > >& vertices , int iters )
{
	double stepSize = StepSize.value;
	std::vector< Point3D< Real > > b( vertices.size() ) , x( vertices.size() ) , oldX( vertices.size() ) , n( vertices.size() );
	Vector< Real > _b( b.size() ) , _x( x.size() );
	std::vector< MetricData< Real > > cData;
	for( int i=0 ; i<(int)vertices.size() ; i++ ) x[i] = Point3D< Real >( vertices[i] );

	InitMetricData( vertices , mesh , cData );

	Point3D< Real > center;
	Real scale = 1.;
	PoseMesh( mesh , x , center , scale );

	// Check if the mesh has any boundary vertices.
	for( int i=0 ; i<(int)mesh.vertex_size() ; i++ )
	{
		bool isBoundary = false;
		EmptyHEMesh::Halfedge_around_vertex_const_circulator iter = mesh.halfedge_around_vertex_begin( i );
		EmptyHEMesh::Halfedge_around_vertex_const_circulator end = iter;
		do
		{
			if( !iter.facet() ) isBoundary = true;
			iter++;
		}
		while( end!=iter );
		if( isBoundary ) fprintf( stderr , "[ERROR] Mesh should be water-tight\n" ) , exit( 0 );
	}

	scale *= MakeUnitArea( x , mesh );

	// Initialize the matrices.
	// M is the current matrix which is initialized with the sparse matrix topology.
	// D stores the mass matrix.
	// L stores the stiffness matrix.
	// The stiffness matrix is unchanged throughout the flow so we can initialize it once.
	SparseMatrix< Real , int > D , L , M;
	GetMatrices< Real , Real >( x , mesh , NULL , &L , &M , true , true );
	if( Randomize.set ) for( int i=0 ; i<x.size() ; i++ ) x[i] = RandomSpherePoint< Real >();

	double radius;

	EigenSolverCholeskyLLt< Real >* solver = NULL;
	for( int i=0 ; i<iters ; i++ )
	{
		double t;
		double tt = Timer::Time();

		t = Timer::Time();
		// Update the matrices that are changing
		GetMatrices< Real , Real >( x , mesh , &D , NULL , NULL , i==0 , true );
		double mTime = Timer::Time()-t;
		// Set the new constraint vector: b = D * x
		D.MultiplyParallel( &x[0] , &b[0] , Threads.value , MULTIPLY_CLEAR );
		// Set the system matrix: M = D + t * L
#pragma omp parallel for num_threads( Threads.value )
		for( int j=0 ; j<D.rows ; j++ )
			for( int k=0  ; k<D.rowSizes[j] ; k++ )
				M[j][k].Value = D[j][k].Value + L[j][k].Value * stepSize;

		t = Timer::Time();

		// If this is the first solve perform both the symbolic and numerical factorization.
		if( !solver ) solver = new EigenSolverCholeskyLLt< Real >( M );
		// Otherwise, if the system matrix has changed, update the numerical factorization.
		else solver->update( M );

		double sTime1 = Timer::Time()-t;

		t = Timer::Time();
		// Solve for each of the x, y, z coefficients independently
		for( int d=0 ; d<3 ; d++ )
		{
#pragma omp parallel for num_threads( Threads.value )
			for( int j=0 ; j<(int)b.size() ; j++ ) _b[j] = b[j][d] , _x[j] = oldX[j][d] = x[j][d];
			solver->solve( &_b[0] , &_x[0] );
#pragma omp parallel for num_threads( Threads.value )
			for( int j=0 ; j<(int)x.size() ; j++ ) x[j][d] = _x[j];
		}

		double sTime2 = Timer::Time()-t;
		double mem = MemoryInfo::UsageMB();
		// Re-scale/center the mash
		TranslateVolumeCenterToOrigin( x , mesh ) , MakeUnitArea( x , mesh );
		double sphericalError = GetSphericalVariation< Real , Point3D< Real > , EmptyHEMesh >( x , mesh , radius , Threads.value );
		double conformalRatio = GetInitializedConformalRatio< Real , Point3D< Real > , EmptyHEMesh >( cData , x , mesh , Threads.value , true );
		double differenceNorm=0 , oldNorm=0;

		for( int j=0 ; j<D.rows ; j++ )
			for( int k=0  ; k<D.rowSizes[j] ; k++ )
			{
				int jj = D[j][k].N;
				differenceNorm += Point3D< Real >::Dot( oldX[j]-x[j] , oldX[jj]-x[jj] ) * D[j][k].Value;
				oldNorm += Point3D< Real >::Dot( oldX[j] , oldX[jj] ) * D[j][k].Value;
			}
		double deformationScale = sqrt( differenceNorm / oldNorm );

		if( Verbose.set )
		{
			tt = Timer::Time() - tt;
			printf( "\rCMCF[%d / %d] %4.2f(s): D-Norm=%6.5f / QC-Ratio=%6.5f / R-Var=%6.5f     " , i+1 , iters , tt , deformationScale , conformalRatio , sphericalError );
			fflush(stdout);
		}
	}
	if( Verbose.set && iters>0 ) printf( "\n" );
	PoseMesh( mesh , x );

	for( int i=0 ; i<(int)x.size() ; i++ ) x[i] = x[i] * scale + center;
	for( int i=0 ; i<(int)vertices.size() ; i++ ) vertices[i] = Point3D< float >( x[i] );
	if( solver ) delete solver;
}


template< class Real >
void Execute( const std::vector< TriangleIndex >& triangles , std::vector< PlyColorVertex< float > >& vertices )
{
	std::vector< Point3D< Real > > sphericalCoordinates( vertices.size() );
	for( int i=0 ; i<vertices.size() ; i++ ) sphericalCoordinates[i] = Point3D< Real >( vertices[i].point );
	EmptyHEMesh heMesh;
	heMesh.SetHalfEdgeData( triangles );
	CMCF( heMesh , sphericalCoordinates , Iterations.value );

	// Normalize the parameterization
	{
		Point3D< double > center;
		double area = 0;
		for( int i=0 ; i<triangles.size() ; i++ )
		{
			Point3D< double > v[] = { sphericalCoordinates[ triangles[i][0] ] , sphericalCoordinates[ triangles[i][1] ] , sphericalCoordinates[ triangles[i][2] ] };
			double a = Length( Point3D< double >::CrossProduct( v[1]-v[0] , v[2]-v[0] ) ) / 2;
			area += a;
			center += ( v[0] + v[1] + v[2] ) / 3 * a;
		}
		center /= area;
#pragma omp parallel for
		for( int i=0 ; i<sphericalCoordinates.size() ; i++ ) sphericalCoordinates[i] -= center;
#pragma omp parallel for
		for( int i=0 ; i<sphericalCoordinates.size() ; i++ ) sphericalCoordinates[i] /= Length( sphericalCoordinates[i] );
	}

	SphericalGeometry::Mesh< Real > mesh;
	{
		std::vector< Point3D< Real > > _vertices( vertices.size() );
		for( int i=0 ; i<vertices.size() ; i++ ) _vertices[i] = vertices[i].point;
		mesh = SphericalGeometry::Mesh< Real >( _vertices , sphericalCoordinates , triangles );
		mesh.makeUnitMass();
	}


	if( !NoCenter.set )
	{
		Timer t;
		static const int MAX_ITERATIONS = 1000;
		if( mesh.normalize( MAX_ITERATIONS , CutOff.value , true , Verbose.set )==MAX_ITERATIONS ) fprintf( stderr , "[WARNING] Failed to meet centering threshold after %d iterations\n" , MAX_ITERATIONS );
		sphericalCoordinates = mesh.vertices;
		if( Verbose.set ) printf( "Centered: %.2f (s)\n" , t.elapsed() );
	}

	if( Out.set )
	{
		std::vector< PlyParametrizedColorVertex< float > > outVertices( vertices.size() );
		for( int i=0 ; i<vertices.size() ; i++ )
		{
			outVertices[i].point = vertices[i].point;
			outVertices[i].color = vertices[i].color;
			outVertices[i].param = sphericalCoordinates[i];
		}
		PlyWriteTriangles( Out.value , outVertices , triangles , PlyParametrizedColorVertex< float >::WriteProperties , PlyParametrizedColorVertex< float >::WriteComponents , ASCII.set ? PLY_ASCII : PLY_BINARY_NATIVE , NULL , 0 );
	}

	if( OutGrid.set || OutTessellation.set )
	{
		struct VertexAndColor
		{
			VertexAndColor( void ){}
			VertexAndColor( Point3D< Real > v , Point3D< Real > c ) : vertex(v) , color(c) {}
			Point3D< Real > vertex , color;
			VertexAndColor operator * ( Real s ) const { return VertexAndColor( vertex*s , color*s ); }
			VertexAndColor operator / ( Real s ) const { return VertexAndColor( vertex/s , color/s ); }
			VertexAndColor operator + ( const VertexAndColor& vc ) const { return VertexAndColor( vertex + vc.vertex , color + vc.color ); }
		};
		Timer t;

		std::vector< VertexAndColor > verticesAndColors( vertices.size() );
		for( int i=0 ; i<vertices.size() ; i++ ) verticesAndColors[i].vertex = Point3D< Real >( vertices[i].point ) , verticesAndColors[i].color = Point3D< Real >( vertices[i].color );
		SphericalGeometry::Tessellation< Real > tessellator( mesh , verticesAndColors , Resolution.value , Verbose.set );
		if( Collapse.set ) tessellator.collapseCells( verticesAndColors );

		if( OutTessellation.set )
		{
			std::vector< Point3D< Real > > _vertices( verticesAndColors.size() ) , _colors( verticesAndColors.size() );
			for( int i=0 ; i<verticesAndColors.size() ; i++ ) _vertices[i] = verticesAndColors[i].vertex , _colors[i] = verticesAndColors[i].color;
			tessellator.write( OutTessellation.value , _vertices , _colors , !ASCII.set );
		}

		SphericalGrid< Real > sGrid;
		sGrid.resize( Resolution.value );
		tessellator.createSGrid( sGrid , Smooth.value , false , false );

		if( OutGrid.set )
		{
			char* ext = GetFileExtension( OutGrid.value );
			if     ( !strcasecmp( ext , "ply" ) ) tessellator.writeGrid( OutGrid.value , Smooth.value , false , false , []( Real v ){ return v<0 ? (Real)-sqrt(-v) : (Real)sqrt(v); } , !ASCII.set );
			else if( !strcasecmp( ext , "sgrid" ) )
			{
				SphericalGrid< float > _sGrid;
				_sGrid.resize( sGrid.resolution() );
				for( int i=0 ; i<sGrid.resolution() ; i++ ) for( int j=0 ; j<sGrid.resolution() ; j++ ) _sGrid(i,j) = (float)sGrid(i,j);
				_sGrid.write( OutGrid.value );
			}
			delete[] ext;
		}
		if( Verbose.set ) printf( "Tessellated: %.2f(s)\n" , t.elapsed() );
	}

}
int main
(
	int argc ,
	char* argv[]
)
{
	cmdLineParse( argc-1 , argv+1 , params );
	if( !In.set )
	{
		Usage( argv[0] );
		return EXIT_FAILURE;
	}
	int fileType;
	std::vector< TriangleIndex > triangles;
	std::vector< PlyColorVertex< float > > vertices;
	bool propertyFlags[ PlyColorVertex< float >::ReadComponents ];
	PlyReadTriangles( In.value , vertices , triangles , PlyColorVertex< float >::ReadProperties , propertyFlags , PlyColorVertex< float >::ReadComponents , fileType );
	bool hasColor = (propertyFlags[3]||propertyFlags[6]) && (propertyFlags[4]||propertyFlags[7]) && (propertyFlags[5]||propertyFlags[8]);

	if( !hasColor )
	{
		std::vector< Point3D< float > > normals( vertices.size() );
		for( int i=0 ; i<triangles.size() ; i++ )
		{
			Point3D< float > v[] = { vertices[ triangles[i][0] ].point , vertices[ triangles[i][1] ].point , vertices[ triangles[i][2] ].point };
			Point3D< float > n = Point3D< float >::CrossProduct( v[1]-v[0] , v[2]-v[0] );
			for( int j=0 ; j<3 ; j++ ) normals[ triangles[i][j] ] += n;
		}
#pragma omp parallel for
		for( int i=0 ; i<vertices.size() ; i++ ) normals[i] /= (float)Length( normals[i] );
		for( int i=0 ; i<(int)normals.size() ; i++ )
		{
			// Negating the normal to comply with earlier version.
			normals[i] = -normals[i];
			normals[i] /= (float)Length( normals[i] );
			vertices[i].color = NormalColor( normals[i] );
		}
	}

	size_t euler = vertices.size() - triangles.size() / 2;
	if( euler!=2 ) fprintf( stderr , "[ERROR] Assuming genus-0, water-tight mesh\n" ) , exit( 0 );

	Execute< double >( triangles , vertices );

	return EXIT_SUCCESS;
}
