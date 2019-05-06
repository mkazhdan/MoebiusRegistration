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
#include <unordered_map>
#include <Eigen/Dense>
#include "Misha/CmdLineParser.h"
#include "Misha/Ply.h"
#include "Misha/SparseMatrix.h"
#include "Misha/Solver.h"
#include "Misha/Timer.h"
#include "Misha/MemoryUsage.h"
#include "Misha/SphericalGeometry.h"

enum
{
	MESH_TYPE_TRIANGLE ,
	MESH_TYPE_POLYGON ,
	MESH_TYPE_TRIANGULATED_POLYGON ,
	MESH_TYPE_COUNT
};
enum
{
	HOLE_FILL_TYPE_NONE ,
	HOLE_FILL_TYPE_CENTER_TRIANGULATION ,
	HOLE_FILL_TYPE_MINIMAL_AREA_TRIANGULATION ,
	HOLE_FILL_TYPE_POLYGON ,
	HOLE_FILL_TYPE_COUNT
};
const char *MeshTypeNames[] = { "triangle mesh" , "polygon mesh" , "implicitly triangulated polygon mesh" };
const char* HoleFillTypeNames[] = { "none" , "center triangulation" , "minimal area triangulation" , "polygon" };
cmdLineParameter< char* > In( "in" ) , Out( "out" ) , OutGrid( "outG" ) , OutTessellation( "outT" );
cmdLineParameter< int > Iterations( "iters" , 100 ) , Threads( "threads" , omp_get_num_procs() ) , Resolution( "res" , 256 ) , SHDegree( "degree" , 1 ) , AdvectionSteps( "aSteps" , 10 ) , MeshType( "mesh" , MESH_TYPE_TRIANGLE+1 ) , HoleFillType( "fill" , HOLE_FILL_TYPE_NONE ) , CenterToInversion( "c2i" , 2 );
cmdLineParameter< float > StepSize( "stepSize" , 0.1f ) , CutOff( "cutOff" , 1e-10f ) , Smooth( "smooth" , 5e-4f ) , AdvectionStepSize( "aStepSize" , 0.05f ) , GSSTolerance( "gssTolerance" , (float)1e-6 ) , PoincareMaxNorm( "poincareMaxNorm" , 2.f );
cmdLineReadable Verbose( "verbose" ) , FullVerbose( "fullVerbose" ) , ASCII( "ascii" ) , Randomize( "random" ) , NoCenter( "noCenter" ) , Collapse( "collapse" ) , NoOrient( "noOrient" ) , NoGridScale( "noGridScale" ) , Lump( "lump" ) , Spherical( "spherical" );

void Usage( const char* ex )
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input mesh>\n" , In.name );
	printf( "\t[--%s <output mesh>]\n" , Out.name );
	printf( "\t[--%s <output spherical grid>]\n" , OutGrid.name );
	printf( "\t[--%s <output spherical tessellation>]\n" , OutTessellation.name );
	printf( "\t[--%s <mesh type>=%d]\n" , MeshType.name , MeshType.value );
	for( int i=0 ; i<MESH_TYPE_COUNT ; i++ ) printf( "\t\t%d] %s\n" , i+1 , MeshTypeNames[i] );
	printf( "\t[--%s <hole fill type>=%d]\n" , HoleFillType.name , HoleFillType.value );
	for( int i=0 ; i<HOLE_FILL_TYPE_COUNT ; i++ ) printf( "\t\t%d] %s\n" , i , HoleFillTypeNames[i] );
	printf( "\t[--%s <CMCF iterations>=%d]\n" , Iterations.name , Iterations.value );
	printf( "\t[--%s <CMCF step size>=%f]\n" , StepSize.name , StepSize.value );
	printf( "\t[--%s <Moebius centering cut-off>=%g]\n" , CutOff.name , CutOff.value );
	printf( "\t[--%s <spherical harmonic degree>=%d]\n" , SHDegree.name , SHDegree.value );
	printf( "\t[--%s <advection steps>=%d]\n" , AdvectionSteps.name , AdvectionSteps.value );
	printf( "\t[--%s <advection step size>=%f]\n" , AdvectionStepSize.name , AdvectionStepSize.value );
	printf( "\t[--%s <spherical resolution>=%d]\n" , Resolution.name , Resolution.value );
	printf( "\t[--%s <spherical diffusion time>=%g]\n" , Smooth.name , Smooth.value );
	printf( "\t[--%s <center to inversion type>=%d]\n" , CenterToInversion.name , CenterToInversion.value );
	printf( "\t\t0] Trivial\n" );
	printf( "\t\t1] Golden section search\n" );
	printf( "\t\t2] Poincare\n" );
	printf( "\t[--%s <golden section search tolerance>=%e]\n" , GSSTolerance.name , GSSTolerance.value );
	printf( "\t[--%s <Poincare max norm>=%e]\n" , PoincareMaxNorm.name , PoincareMaxNorm.value );
	printf( "\t[--%s]\n" , Randomize.name );
	printf( "\t[--%s]\n" , NoCenter.name );
	printf( "\t[--%s]\n" , ASCII.name );
	printf( "\t[--%s]\n" , Collapse.name );
	printf( "\t[--%s]\n" , NoOrient.name );
	printf( "\t[--%s]\n" , NoGridScale.name );
	printf( "\t[--%s]\n" , Spherical.name );
	printf( "\t[--%s]\n" , Lump.name );
	printf( "\t[--%s]\n" , Verbose.name );
	printf( "\t[--%s]\n" , FullVerbose.name );
}
cmdLineReadable* params[] = { &In , &Out , &OutGrid , &OutTessellation , &Iterations , &StepSize , &Threads , &Verbose , &FullVerbose , &Randomize , &ASCII , &NoCenter , &Resolution , &Smooth , &Collapse , &NoOrient , &SHDegree , &AdvectionSteps , &AdvectionStepSize , &CenterToInversion , &GSSTolerance , &PoincareMaxNorm , &Lump , &MeshType , &HoleFillType , &Spherical , &NoGridScale , NULL };

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

template< typename Real >
SquareMatrix< Real , 3 > MassMatrix( const Point3D< Real > vertices[] )
{
	SquareMatrix< Real , 3 > mass;
	Real area = (Real)Length( Point3D< Real >::CrossProduct( vertices[1]-vertices[0] , vertices[2]-vertices[0] ) ) / 2;
	mass(1,0) = mass(0,1) = mass(2,0) = mass(0,1) = mass(0,2) = mass(1,2) = mass(2,1) = area / 12;
	mass(0,0) = mass(1,1) = mass(2,2) = area / 6;
	return mass;
}
template< typename Real >
SquareMatrix< Real , 3 > StiffnessMatrix( const Point3D< Real > vertices[] )
{
	SquareMatrix< Real , 3 > stiffness;
	Real area = (Real)Length( Point3D< Real >::CrossProduct( vertices[1]-vertices[0] , vertices[2]-vertices[0] ) ) / 2;
	for( int i=0 ; i<3 ; i++ )
	{
		int i0 = (i+1)%3 , i1 = (i+2)%3;
		Real dot = Point3D< Real >::Dot( vertices[i0]-vertices[i] , vertices[i1]-vertices[i] );
		stiffness(i0,i1) = stiffness(i1,i0) = - dot / area;
	}
	stiffness(0,0) = - ( stiffness(0,1) + stiffness(0,2) );
	stiffness(1,1) = - ( stiffness(1,2) + stiffness(1,0) );
	stiffness(2,2) = - ( stiffness(2,0) + stiffness(2,1) );
	return stiffness;
}


template< typename Real >
struct Mesh
{
	virtual int faces( void ) const = 0;
	virtual Point3D< Real > normal( int f ) const = 0;
	virtual Point3D< Real > center( int f ) const = 0;
	virtual Real area( int f ) const = 0;
	virtual SparseMatrix< Real , int > massMatrix( bool lump ) const = 0;
	virtual SparseMatrix< Real , int > stiffnessMatrix( void ) const = 0;
	virtual std::vector< Real > vertexAreas( void ) const = 0;
	virtual std::vector< Point3D< Real > > vertexNormals( void ) const = 0;
	virtual void push_back( int v1 , int v2 , int v3 ) = 0;
	virtual void pop_back( void ) = 0;
	virtual void update( void ) {}

	std::vector< Point3D< Real > > vertices;

	Real area( void ) const
	{
		Real a=0;
#pragma omp parallel for reduction( + : a )
		for( int f=0 ; f<faces() ; f++ ) a += area( f );
		return a;
	}
	Point3D< Real > center( void ) const
	{
		Real area=0 , x=0 , y=0 , z=0;
#pragma omp parallel for reduction( + : area , x , y , z )
		for( int f=0 ; f<faces() ; f++ )
		{
			Real a = this->area(f);
			Point3D< Real > c = center(f)*a;
			area += a , x += c[0] , y += c[1] , z += c[2];
		}
		return Point3D< Real >( x , y , z ) / area;
	}
	Real makeUnitArea( void )
	{
		Real scl = (Real)( 1./sqrt( area() ) );
#pragma omp parallel for
		for( int i=0 ; i<vertices.size() ; i++ ) vertices[i] *= scl;
		return (Real)1./scl;
	}
	void pose( Point3D< Real >& center , Real& scale )
	{
		center = this->center();
		for( int i=0 ; i<vertices.size() ; i++ ) vertices[i] -= center;
		scale = 0;
		for( int i=0 ; i<vertices.size() ; i++ ) for( int j=0 ; j<3 ; j++ ) scale = std::max< Real >( scale , (Real)fabs(vertices[i][j]) );
		scale *= 2;
		for( int i=0 ; i<vertices.size() ; i++ ) vertices[i] /= scale;
	}
	void pose( void )
	{
		Point3D< Real > center;
		Real scale;
		pose( center , scale );
	}

	Real radialDeviation( void ) const
	{
		Real area = 0;
		std::vector< Real > vertexAreas = this->vertexAreas();
#pragma omp parallel for reduction( + : area )
		for( int i=0 ; i<vertexAreas.size() ; i++ ) area += vertexAreas[i];
		Point3D< Real > center;
		for( int i=0 ; i<vertices.size() ; i++ ) center += vertices[i] * vertexAreas[i];
		center /= area;

		// set v = \sum_i[ w_i * v_i ]
		// var = \sum_i[ w_i * ( v_i -v )^2 ]
		//     = \sum_i[ w_i * ( v_i^2 + v^2 - 2 * v_i * v ) ]
		//     = \sum_i[ w_i * v_i^2 ] + \sum_i[ w_i * v^2 ] - \sum_i[ 2 * w_i * v_i * v ) ]
		//     = \sum_i[ w_i * v_i^2 ] + v^2 - 2 * v^2
		//     = \sum_i[ w_i * v_i^2 ] - v^2
		Real average = 0 , var = 0;
#pragma omp parallel for reduction( + : average , var )
		for( int i=0 ; i<vertices.size() ; i++ )
		{
			Real r = (Real)sqrt( Point3D< Real >::SquareNorm( Point3D< Real >( vertices[i]-center ) ) );
			average += r * vertexAreas[i];
			var += r * r * vertexAreas[i];
		}
		average /= area , var /= area;
		var -= average * average;
		var = (Real)sqrt( var );
		return var / average;
	}
};

template< typename Real >
struct TriangleMesh : public Mesh< Real >
{
	std::vector< TriangleIndex > triangles;
	using Mesh< Real >::center;
	using Mesh< Real >::vertices;

	int faces( void ) const { return (int)triangles.size(); }
	void push_back( int v1 , int v2 , int v3 ){ triangles.push_back( TriangleIndex( v1 , v2 , v3 ) ); }
	void pop_back( void ){ triangles.pop_back(); }

	Point3D< Real > normal( int t ) const { return Point3D< Real >::CrossProduct( vertices[ triangles[t][1] ] - vertices[ triangles[t][0] ] , vertices[ triangles[t][2] ] - vertices[ triangles[t][0] ] ); }
	Point3D< Real > center( int t ) const { return ( vertices[ triangles[t][0] ] + vertices[ triangles[t][1] ] + vertices[ triangles[t][2] ] ) / 3; }
	Real area( int t ) const { return (Real)Length( normal(t) )/2.; }

#if 0
	template< typename TriangleMatrixFunctor >
	SparseMatrix< Real , int > matrix( TriangleMatrixFunctor F ) const
#else
	SparseMatrix< Real , int > matrix( SquareMatrix< Real , 3 >(*F)( const Point3D< Real > [] ) ) const
#endif
	{
		struct Entry
		{
			Entry( void ) : row(-1) , col(-1) , value(0){}
			Entry( int r , int c , Real v ) : row(r) , col(c) , value(v){}
			int row , col;
			Real value;
		};
		SparseMatrix< Real , int > M;
		M.resize( (int)vertices.size() );
		std::vector< std::vector< Entry > > entries( omp_get_max_threads() );
#pragma omp parallel for
		for( int t=0 ; t<triangles.size() ; t++ )
		{
			int thread = omp_get_thread_num();
			Point3D< Real > v[] = { vertices[ triangles[t][0] ] , vertices[ triangles[t][1] ] , vertices[ triangles[t][2] ] };
			SquareMatrix< Real , 3 > m = F( v );
			for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) entries[thread].push_back( Entry( triangles[t][i] , triangles[t][j] , m(i,j) ) );
		}
		for( int i=0 ; i<entries.size() ; i++ ) for( int j=0 ; j<entries[i].size() ; j++ )	M.rowSizes[ entries[i][j].row ]++;
#pragma omp parallel for
		for( int i=0 ; i<M.rows ; i++ )
		{
			int rowSize = M.rowSizes[i];
			M.rowSizes[i] = 0;
			M.SetRowSize( i , rowSize );
			M.rowSizes[i] = 0;
		}
		for( int i=0 ; i<entries.size() ; i++ ) for( int j=0 ; j<entries[i].size() ; j++ ) M[ entries[i][j].row ][ M.rowSizes[entries[i][j].row]++ ] = MatrixEntry< Real , int >( entries[i][j].col , entries[i][j].value );
#pragma omp parallel for
		for( int i=0 ; i<M.rows ; i++ )
		{
			std::unordered_map< int , Real > row;
			for( int j=0 ; j<M.rowSizes[i] ; j++ ) row[ M[i][j].N ] += M[i][j].Value;
			M.SetRowSize( i , (int)row.size() );
			int j=0;
			for( typename std::unordered_map< int , Real >::const_iterator iter=row.begin() ; iter!=row.end() ; iter++ ) M[i][j++] = MatrixEntry< Real , int >( iter->first , iter->second );
		}
		return M;
	}
	SparseMatrix< Real , int > massMatrix( bool lump ) const
	{
		SparseMatrix< Real , int > mass = matrix( MassMatrix );
		if( lump )
		{
			for( int i=0 ; i<mass.rows ; i++ )
			{
				Real sum = 0;
				for( int j=0 ; j<mass.rowSizes[i] ; j++ ){ sum += mass[i][j].Value ; mass[i][j].Value = 0; }
				for( int j=0 ; j<mass.rowSizes[i] ; j++ ) if( mass[i][j].N==i ) mass[i][j].Value = sum;
			}
		}
		return mass;
	}
	SparseMatrix< Real , int > stiffnessMatrix( void ) const { return matrix( StiffnessMatrix ); }
	
	std::vector< Real > vertexAreas( void ) const
	{
		std::vector< Real > areas( vertices.size() , 0 );
		for( int t=0 ; t<triangles.size() ; t++ )
		{
			Real a = area(t)/3;
			for( int j=0 ; j<3 ; j++ ) areas[ triangles[t][j] ] += a;
		}
		return areas;
	}
	std::vector< Point3D< Real > > vertexNormals( void ) const
	{
		std::vector< Point3D< Real > > normals( vertices.size() );
		for( int t=0 ; t<triangles.size() ; t++ )
		{
			Point3D< Real > n = normal(t);
			Real a = area(t);
			n *= a / (Real)Length(n);
			for( int j=0 ; j<3 ; j++ ) normals[ triangles[t][j] ] += n;
		}
		return normals;
	}
};

template< typename Real >
struct PolygonMesh : public Mesh< Real >
{
protected:
	std::unordered_map< unsigned long long , int > _edges;
public:
	std::vector< std::vector< int > > polygons;
	using Mesh< Real >::center;
	using Mesh< Real >::vertices;

	int faces( void ) const { return (int)polygons.size(); }
	void push_back( int v1 , int v2 , int v3 ){ polygons.push_back( std::vector< int >( { v1 , v2 , v3 } ) ); }
	void pop_back( void ){ polygons.pop_back(); }

	Point3D< Real > normal( int p ) const
	{
		Point3D< Real > n;
		for( int i=0 ; i<polygons[p].size() ; i++ )
		{
			Point3D< Real > b = vertices[ polygons[p][ (i+1)%polygons[p].size() ] ] + vertices[ polygons[p][i] ];
			Point3D< Real > e = vertices[ polygons[p][ (i+1)%polygons[p].size() ] ] - vertices[ polygons[p][i] ];
			n += Point3D< Real >::CrossProduct( b , e ) / 2;
		}
		return n;
	}
	Point3D< Real > center( int p ) const
	{
		Point3D< Real > c;
		for( int i=0 ; i<polygons[p].size() ; i++ ) c += vertices[ polygons[p][i] ];
		return c / (int)polygons[p].size();
	}
	Real area( int p ) const { return (Real)Length( normal( p ) )/2.; }

	SparseMatrix< Real , int > massMatrix( bool lump ) const
	{
		SparseMatrix< Real , int > mass;
		mass.resize( (int)vertices.size() );
		Eigen::DiagonalMatrix< double , Eigen::Dynamic > M0 = mass0();
#pragma omp parallel for
		for( int i=0 ; i<mass.rows ; i++ )
		{
			mass.SetRowSize( i , 1 );
			mass[i][0] = MatrixEntry< Real , int >( i , M0.diagonal()(i) );
		}
		return mass;
	}

	SparseMatrix< Real , int > stiffnessMatrix( void ) const
	{
		SparseMatrix< Real , int > stiffness;
		stiffness.resize( (int)vertices.size() );
		Eigen::SparseMatrix< double > D = coBoundary();
		Eigen::SparseMatrix< double > M = mass1();
		Eigen::SparseMatrix< double > S = D.transpose() * M * D;

		for( int i=0 ; i<S.outerSize() ; i++ ) for( Eigen::SparseMatrix< double >::InnerIterator it(S,i) ; it ; ++it ) stiffness.rowSizes[it.row()]++;
		for( int i=0 ; i<stiffness.rows ; i++ )
		{
			int cols = stiffness.rowSizes[i];
			stiffness.rowSizes[i] = 0;
			stiffness.SetRowSize(i,cols);
			stiffness.rowSizes[i] = 1;
		}
		for( int i=0 ; i<S.outerSize() ; i++ ) for( Eigen::SparseMatrix< double >::InnerIterator it(S,i) ; it ; ++it ) 
		{
			if( it.row()==it.col() ) stiffness[it.row()][0] = MatrixEntry< Real , int >( it.col() , it.value() );
			else
			{
				int idx = stiffness.rowSizes[it.row()]++;
				stiffness[it.row()][idx] = MatrixEntry< Real , int >( it.col() , it.value() );
			}
		}
		return stiffness;
	}
	std::vector< Real > vertexAreas( void ) const
	{
		std::vector< Real > areas( vertices.size() , 0 );
		for( int p=0 ; p<polygons.size() ; p++ )
		{
			Real a = area(p)/(int)polygons[p].size();
			for( int j=0 ; j<polygons[p].size() ; j++ ) areas[ polygons[p][j] ] += a;
		}
		return areas;
	}
	std::vector< Point3D< Real > > vertexNormals( void ) const
	{
		std::vector< Point3D< Real > > normals( vertices.size() );
		for( int p=0 ; p<polygons.size() ; p++ )
		{
			Point3D< Real > n = normal(p);
			Real a = area(p);
			n *= a / (Real)Length(n);
			for( int j=0 ; j<polygons.size() ; j++ ) normals[ polygons[p][j] ] += n;
		}
		return normals;
	}

	void update( void )
	{
		for( int i=0 ; i<polygons.size() ; i++ ) for( int j=0 ; j<polygons[i].size() ; j++ ) for( int k=0 ; k<polygons[i].size() ; k++ ) if( j!=k ) _edges[ EdgeKey( polygons[i][j] , polygons[i][k] ) ] = 0;
		int idx = 0;
		for( auto e=_edges.begin() ; e!=_edges.end() ; e++ ) e->second = idx++;
	}

	Eigen::DiagonalMatrix< double , Eigen::Dynamic > mass0( void ) const
	{
		Eigen::DiagonalMatrix< double , Eigen::Dynamic > M0( vertices.size() );
		std::vector< double > areas( polygons.size() );
#pragma omp parallel for
		for( int i=0 ; i<polygons.size() ; i++ ) areas[i] = area( i ) / polygons[i].size();
#pragma omp parallel for
		for( int i=0 ; i<vertices.size() ; i++ ) M0.diagonal()(i) = 0;
		for( int i=0 ; i<polygons.size() ; i++ ) for( int j=0 ; j<polygons[i].size() ; j++ ) M0.diagonal()( polygons[i][j] ) += areas[i];
		return M0;
	}
	Eigen::SparseMatrix< double > mass1( double lambda=1 ) const
	{
		std::vector< Eigen::Triplet< double > > triplets;
		for( int i=0 ; i<polygons.size() ; i++ )
		{
			Eigen::MatrixXd S = faceMass1( i , lambda );
			const std::vector< int > &p = polygons[i];
			int pSize = (int)p.size();
			for( int j=0 ; j<pSize ; j++ ) for( int k=0 ; k<pSize ; k++ )
			{
				auto it1 = _edges.find( EdgeKey( p[j] , p[(j+1)%pSize] ) );
				auto it2 = _edges.find( EdgeKey( p[k] , p[(k+1)%pSize] ) );
				double sign1 = p[j]<p[(j+1)%pSize] ? 1 : -1;
				double sign2 = p[k]<p[(k+1)%pSize] ? 1 : -1;
				triplets.push_back( Eigen::Triplet< double >( it1->second , it2->second , S(j,k)*sign1*sign2 ) );
			}
		}
		Eigen::SparseMatrix< double > M1( (int)_edges.size() , (int)_edges.size() );
		M1.setFromTriplets( triplets.begin() , triplets.end() );
		return M1;
	}
	Eigen::SparseMatrix< double > coBoundary( void ) const
	{
		std::vector< Eigen::Triplet< double > > triplets;
		triplets.reserve( _edges.size()*2 );

		int idx = 0;
		for( auto e=_edges.begin() ; e!=_edges.end() ; e++ )
		{
			int i , j;
			FactorEdgeKey( e->first , i , j );
			triplets.push_back( Eigen::Triplet< double >( e->second , i , -1 ) );
			triplets.push_back( Eigen::Triplet< double >( e->second , j ,  1 ) );
		}
		Eigen::SparseMatrix< double > D( (int)_edges.size() , (int)vertices.size() );
		D.setFromTriplets( triplets.begin() , triplets.end() );
		return D;
	}

	Eigen::MatrixXd faceMass1( int p , double lambda ) const
	{
		std::vector< Point3D< double > > _vertices( polygons[p].size() );
		for( int i=0 ; i<polygons[p].size() ; i++ ) _vertices[i] = Point3D< double >( vertices[ polygons[p][i] ] );

		Eigen::MatrixXd B( polygons[p].size() , 3 ) , E( polygons[p].size()  , 3 ) , _E( polygons[p].size() , 3 );
		for( int i=0 ; i<polygons[p].size() ; i++ )
		{
			Point3D< double > v[] = { _vertices[i] , _vertices[ (i+1)%polygons[p].size() ] };
			for( int j=0 ; j<3 ; j++ )
			{
				E(i,j) =   v[1][j]-v[0][j];
				B(i,j) = ( v[1][j]+v[0][j] ) / 2.;
			}
		}
		Eigen::MatrixXd A = E.transpose() * B;
		Eigen::MatrixXd _S = B * B.transpose() * sqrt( 2 ) / A.norm();
		Point3D< double > n( -A(1,2) , A(0,2) , - A(0,1) );
		n /= Length( n );
		for( int i=0 ; i<polygons[p].size() ; i++ ) _vertices[i] -= n * Point3D< double >::Dot( n , _vertices[i] );
		for( int i=0 ; i<polygons[p].size() ; i++ )
		{
			Point3D< double > v[] = { _vertices[i] , _vertices[ (i+1)%polygons[p].size() ] };
			for( int j=0 ; j<3 ; j++ ) _E(i,j) = v[1][j] - v[0][j];
		}

		Eigen::MatrixXd C = _E.transpose().fullPivLu().kernel();
		return _S + C * C.transpose() * lambda;
	}
};

// Triangulate the polygons by adding the mid-point but constrain the value at the mid-point to be the average of the values along the vertices of the face
template< typename Real >
struct TriangulatedPolygonMesh : public Mesh< Real >
{
	std::vector< std::vector< int > > polygons;
	using Mesh< Real >::center;
	using Mesh< Real >::vertices;

	int faces( void ) const { return (int)polygons.size(); }
	void push_back( int v1 , int v2 , int v3 ){ polygons.push_back( std::vector< int >( { v1 , v2 , v3 } ) ); }
	void pop_back( void ){ polygons.pop_back(); }

	Point3D< Real > normal( int p ) const
	{
		Point3D< Real > n;
		for( int i=0 ; i<polygons[p].size() ; i++ )
		{
			Point3D< Real > b = vertices[ polygons[p][ (i+1)%polygons[p].size() ] ] + vertices[ polygons[p][i] ];
			Point3D< Real > e = vertices[ polygons[p][ (i+1)%polygons[p].size() ] ] - vertices[ polygons[p][i] ];
			n += Point3D< Real >::CrossProduct( b , e ) / 2;
		}
		return n;
	}
	Point3D< Real > center( int p ) const
	{
		Point3D< Real > c;
		for( int i=0 ; i<polygons[p].size() ; i++ ) c += vertices[ polygons[p][i] ];
		return c / (int)polygons[p].size();
	}
	Real area( int p ) const
	{
		if( polygons[p].size()==3 )
		{
			Point3D< Real > v[] = { vertices[ polygons[p][0] ] , vertices[ polygons[p][1] ] , vertices[ polygons[p][2] ] };
			Real a = (Real)Length( Point3D< Real >::CrossProduct( v[1]-v[0] , v[2]-v[0] ) )/2;
			return a;
		}
		else
		{
			Real a = 0;
			Point3D< Real > c = center(p);
			for( int i0=0 ; i0<polygons[p].size() ; i0++ )
			{
				int i1 = (i0+1)%polygons[p].size();
				Point3D< Real > v[] = { vertices[ polygons[p][i0] ] , vertices[ polygons[p][i1] ] , c };
				a += (Real)Length( Point3D< Real >::CrossProduct( v[1]-v[0] , v[2]-v[0] ) )/2;
			}
			return a;
		}
	}

	std::vector< Real > vertexAreas( void ) const
	{
		std::vector< Real > areas( vertices.size() , 0 );
		for( int p=0 ; p<polygons.size() ; p++ )
		{
			if( polygons[p].size()==3 )
			{
				Point3D< Real > v[] = { vertices[ polygons[p][0] ] , vertices[ polygons[p][1] ] , vertices[ polygons[p][2] ] };
				Real a = (Real)Length( Point3D< Real >::CrossProduct( v[1]-v[0] , v[2]-v[0] ) )/2;
				a /= 3;
				for( int i=0 ; i<3 ; i++ ) areas[ polygons[p][i] ] += a;
			}
			else
			{
				Point3D< Real > c = center(p);
				for( int i0=0 ; i0<polygons[p].size() ; i0++ )
				{
					int i1 = (i0+1)%polygons[p].size();
					Point3D< Real > v[] = { vertices[ polygons[p][i0] ] , vertices[ polygons[p][i1] ] , c };
					Real a = (Real)Length( Point3D< Real >::CrossProduct( v[1]-v[0] , v[2]-v[0] ) )/2;
					a /= 2;
					areas[ polygons[p][i0] ] += a , areas[ polygons[p][i1] ] += a;
				}
			}
		}
		return areas;
	}
	std::vector< Point3D< Real > > vertexNormals( void ) const
	{
		std::vector< Point3D< Real > > normals( vertices.size() );
		for( int p=0 ; p<polygons.size() ; p++ )
		{
			Point3D< Real > n = normal(p);
			Real a = area(p);
			n *= a / (Real)Length(n);
			for( int j=0 ; j<polygons.size() ; j++ ) normals[ polygons[p][j] ] += n;
		}
		return normals;
	}

	SparseMatrix< Real , int > massMatrix( bool lump ) const
	{
		SparseMatrix< Real , int > mass = matrix( MassMatrix );
		if( lump )
		{
#pragma omp parallel for
			for( int i=0 ; i<mass.rows ; i++ )
			{
				Real sum = 0;
				for( int j=0 ; j<mass.rowSizes[i] ; j++ ){ sum += mass[i][j].Value ; mass[i][j].Value = 0; }
				for( int j=0 ; j<mass.rowSizes[i] ; j++ ) if( mass[i][j].N==i ){ mass[i][j].Value = sum ; break; }
			}
		}
		return mass;
	}
	SparseMatrix< Real , int > stiffnessMatrix( void ) const { return matrix( StiffnessMatrix ); }

#if 0
	template< typename TriangleMatrixFunctor >
	SparseMatrix< Real , int > matrix( TriangleMatrixFunctor F ) const
#else
	SparseMatrix< Real , int > matrix( SquareMatrix< Real , 3 >(*F)( const Point3D< Real > [] ) ) const
#endif
	{
		struct Entry
		{
			Entry( void ) : row(-1) , col(-1) , value(0){}
			Entry( int r , int c , Real v ) : row(r) , col(c) , value(v){}
			int row , col;
			Real value;
		};
		SparseMatrix< Real , int > M;
		M.resize( (int)vertices.size() );
		std::vector< std::vector< Entry > > entries( omp_get_max_threads() );
#pragma omp parallel for
		for( int p=0 ; p<polygons.size() ; p++ )
		{
			int thread = omp_get_thread_num();
			int pSize = (int)polygons[p].size();
			Eigen::MatrixXd m = faceMatrix( p , F );
			for( int i=0 ; i<pSize ; i++ ) for( int j=0 ; j<pSize ; j++ ) entries[thread].push_back( Entry( polygons[p][i] , polygons[p][j] , m(i,j) ) );
		}
		for( int i=0 ; i<entries.size() ; i++ ) for( int j=0 ; j<entries[i].size() ; j++ )	M.rowSizes[ entries[i][j].row ]++;
#pragma omp parallel for
		for( int i=0 ; i<M.rows ; i++ )
		{
			int rowSize = M.rowSizes[i];
			M.rowSizes[i] = 0;
			M.SetRowSize( i , rowSize );
			M.rowSizes[i] = 0;
		}
		for( int i=0 ; i<entries.size() ; i++ ) for( int j=0 ; j<entries[i].size() ; j++ ) M[ entries[i][j].row ][ M.rowSizes[entries[i][j].row]++ ] = MatrixEntry< Real , int >( entries[i][j].col , entries[i][j].value );
#pragma omp parallel for
		for( int i=0 ; i<M.rows ; i++ )
		{
			std::unordered_map< int , Real > row;
			for( int j=0 ; j<M.rowSizes[i] ; j++ ) row[ M[i][j].N ] += M[i][j].Value;
			M.SetRowSize( i , (int)row.size() );
			int j=0;
			for( typename std::unordered_map< int , Real >::const_iterator iter=row.begin() ; iter!=row.end() ; iter++ ) M[i][j++] = MatrixEntry< Real , int >( iter->first , iter->second );
		}
		return M;
	}

	Eigen::MatrixXd centerConstraintMatrix( int p ) const
	{
		int pSize = (int)polygons[p].size();
		Eigen::MatrixXd C(pSize+1,pSize);
		for( int i=0 ; i<C.rows() ; i++ ) for( int j=0 ; j<C.cols() ; j++ ) C(i,j) = 0;
		for( int i=0 ; i<pSize ; i++ ) C(i,i) = 1 , C(pSize,i) = 1./pSize;
		return C;
	}

	template< typename TriangleMatrixFunctor >
	Eigen::MatrixXd faceMatrix( int p , TriangleMatrixFunctor F ) const
	{
		int pSize = (int)polygons[p].size();
		if( pSize==3 )
		{
			Eigen::MatrixXd M( pSize , pSize );
			Point3D< Real > v[] = { vertices[ polygons[p][0] ] , vertices[ polygons[p][1] ] , vertices[ polygons[p][2] ] };
			SquareMatrix< Real , 3 > m = F( v );
			for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) M(i,j) = m(i,j);
			return M;
		}
		else
		{
			Eigen::MatrixXd M( pSize+1 , pSize+1 ) , C = centerConstraintMatrix(p);
			for( int i=0 ; i<M.rows() ; i++ ) for( int j=0 ; j<M.cols() ; j++ ) M(i,j) = 0;

			Point3D< Real > c = center(p);
			int idx[] = { 0 , 0 , pSize };
			for( int i0=0 ; i0<pSize ; i0++ )
			{
				int i1 = (i0+1)%pSize;
				idx[0] = i0 , idx[1] = i1;
				Point3D< Real > v[] = { vertices[ polygons[p][i0] ] , vertices[ polygons[p][i1] ] , c };
				SquareMatrix< Real , 3 > m = F( v );
				for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ ) M( idx[j] , idx[k] ) += m(j,k);
			}
			return C.transpose() * M * C;
		}
	}
};

template< typename Real >
struct QuasiConformalRatio
{
	static const Real ConformalCutOff;
	std::vector< TriangleIndex > triangles;
	std::vector< SquareMatrix< Real , 2 > > massInvs;
	std::vector< Real > areas;

	QuasiConformalRatio( const            TriangleMesh< Real > & tMesh , int tCount ) : QuasiConformalRatio( & tMesh.vertices[0] , & tMesh.triangles[0] , tCount ){ }
	QuasiConformalRatio( const             PolygonMesh< Real > & pMesh , int tCount ) : QuasiConformalRatio( & pMesh.vertices[0] , & pMesh.polygons [0] , tCount ){ }
	QuasiConformalRatio( const TriangulatedPolygonMesh< Real > &tpMesh , int tCount ) : QuasiConformalRatio( &tpMesh.vertices[0] , &tpMesh.polygons [0] , tCount ){ }

	QuasiConformalRatio( const Point3D< Real > *v , const TriangleIndex *t , int tCount )
	{
		triangles.resize( tCount );
		for( int i=0 ; i<tCount ; i++ ) triangles[i] = t[i];
		_init( v );
	}
	QuasiConformalRatio( const Point3D< Real > *v , const std::vector< int > *p , int pCount )
	{
		MinimalAreaTriangulation< Real > MAT;
		int count = 0;
		for( int i=0 ; i<pCount ; i++ ) count += (int)p[i].size()-2;
		triangles.reserve( count );
		for( int i=0 ; i<pCount ; i++ )
		{
			const std::vector< int >& poly = p[i];
			std::vector< Point3D< Real > > _vertices( poly.size() );
			std::vector< TriangleIndex > _triangles;
			for( int j=0 ; j<poly.size() ; j++ ) _vertices[j] = v[ poly[j] ];
			MAT.GetTriangulation( _vertices , _triangles );
			for( int j=0 ; j<_triangles.size() ; j++ ) triangles.push_back( TriangleIndex( poly[ _triangles[j][0] ] , poly[ _triangles[j][1] ] , poly[ _triangles[j][2] ] ) );
		}
		_init( v );
	}

	Real operator()( const std::vector< Point3D< Real > > &vertices ) const
	{
		Real error = 0 , area = 0;
#pragma omp parallel for reduction( + : error , area )
		for( int i=0 ; i<triangles.size() ; i++ )
		{
			Point3D< Real > v[] = { vertices[ triangles[i][0] ] , vertices[ triangles[i][1] ] , vertices[ triangles[i][2] ] };
			SquareMatrix< Real , 2 > m = massInvs[i] * _MassMatrix( v );
			Real det = m.determinant() , trace = m.trace();
			trace /= 2;
			if( det/4<=ConformalCutOff*ConformalCutOff ) continue;
			// The eigenvectors of this matrix are the roots of
			// P(x) = [x-d(0,0)]*[x-d(1,1)] - d(1,0)*d(0,1)
			//      = x^2 - x * Tr(d) + Det(d)
			// x = Tr(d)/2 +/- sqrt( Tr^2(d)/4 - Det(d) )
			Real disc  = trace*trace-det;
			if( disc<=Real(0) ) disc = 0;
			else                disc = sqrt( disc );
			Real x1 = trace - disc;
			Real x2 = trace + disc;
			if( x1<=0 ) continue;
			Real tError = sqrt( x2/x1 ) * areas[i];
			error += tError;
			area += areas[i];
		}
		return error / area;
	}

protected:
	static SquareMatrix< Real , 2 > _MassMatrix( const Point3D< Real > v[3] )
	{
		Point3D< Real > t[] = { v[1]-v[0] , v[2]-v[0] };
		SquareMatrix< Real , 2 > mass;
		for( int  j=0 ; j<2 ; j++ ) for( int k=0 ; k<2 ; k++ ) mass(j,k) = Point3D< Real >::Dot( t[j] , t[k] );
		return mass;
	}
	void _init( const Point3D< Real > *vertices )
	{
		massInvs.resize( triangles.size() );
		areas.resize( triangles.size() , 0 );

		Real area = 0;
#pragma omp parallel for reduction( + : area )
		for( int i=0 ; i<triangles.size() ; i++ )
		{
			Point3D< Real > v[] = { vertices[ triangles[i][0] ] , vertices[ triangles[i][1] ] , vertices[ triangles[i][2] ] };
			areas[i] = (Real)Length( Point3D< Real >::CrossProduct( v[1]-v[0] , v[2]-v[0] ) ) / 2;
			area += areas[i];
			if( areas[i]<=ConformalCutOff ) continue;
			massInvs[i] = _MassMatrix( v ).inverse();
		}
#pragma omp parallel for
		for( int i=0 ; i<areas.size() ; i++ ) areas[i] /= area;

	}
};
template< typename Real > const Real QuasiConformalRatio< Real >::ConformalCutOff = (Real)1e-15;

int FaceSize( const TriangleIndex& tri ){ return 3; }
int FaceSize( const std::vector< int > &poly ){ return (int)poly.size(); }

template< typename Face >
std::vector< std::vector< int > > BoundaryLoops( const std::vector< Face > &faces )
{
	auto EdgeKey = []( unsigned int v1 , unsigned int v2 ){ return ( (unsigned long long)v1 )<<32 | ( (unsigned long long)v2 ); };
	auto FactorEdgeKey = []( unsigned long long key , int& v1 , int& v2 ){ v1=(int)(key>>32) , v2 = (int)( (key<<32)>>32 ); };

	std::unordered_set< unsigned long long > edges;
	for( int f=0 ; f<faces.size() ; f++ )
	{
		int fSize = FaceSize( faces[f] );
		for( int v=0 ; v<fSize ; v++ )
		{
			int v1 = faces[f][v] , v2 = faces[f][(v+1)%fSize];
			edges.insert( EdgeKey( v1 , v2 ) );
		}
	}

	std::vector< std::pair< int , int > > boundaryEdges;
	for( auto i=edges.begin() ; i!=edges.end() ; i++ )
	{
		int v1 , v2;
		FactorEdgeKey( *i , v1 , v2 );
		if( edges.find( EdgeKey( v2 , v1 ) )==edges.end() ) boundaryEdges.push_back( std::pair< int , int >( v2 , v1 ) );
	}

	std::vector< std::vector< int > > boundaryLoops;
	while( boundaryEdges.size() )
	{
		std::vector< int > boundaryLoop;
		std::pair< int , int > boundaryEdge = boundaryEdges.back();
		boundaryEdges.pop_back();
		int start = boundaryEdge.first , v = boundaryEdge.second;
		boundaryLoop.push_back( v );
		while( v!=start )
		{
			for( int i=0 ; i<boundaryEdges.size() ; i++ ) if( boundaryEdges[i].first==v )
			{
				v = boundaryEdges[i].second;
				boundaryLoop.push_back( v );
				boundaryEdges[i] = boundaryEdges.back();
				boundaryEdges.pop_back();
				break;
			}
		}
		boundaryLoops.push_back( boundaryLoop );
	}
	return boundaryLoops;
}
template< typename Real >
std::vector< std::vector< int > > BoundaryLoops( const TriangleMesh< Real > &mesh ){ return BoundaryLoops( mesh.triangles ); }
template< typename Real >
std::vector< std::vector< int > > BoundaryLoops( const PolygonMesh< Real > &mesh ){ return BoundaryLoops( mesh.polygons ); }
template< typename Real >
std::vector< std::vector< int > > BoundaryLoops( const TriangulatedPolygonMesh< Real > &mesh ){ return BoundaryLoops( mesh.polygons ); }

template< class Real , typename Mesh >
void CMCF( Mesh &mesh , int iters , bool lump , int fCount )
{
	double stepSize = StepSize.value;
	QuasiConformalRatio< Real > qcRatio( mesh , fCount );
	std::vector< Point3D< Real > > &vertices = mesh.vertices;
	std::vector< Point3D< Real > > b( vertices.size() ) , oldX( vertices.size() ) , n( vertices.size() );
	Vector< Real > _b( b.size() ) , _x( vertices.size() );
	Point3D< Real > center;
	Real scale = 1.;
	mesh.pose( center , scale );
	scale *= mesh.makeUnitArea();

	// Initialize the matrices.
	// M is the current matrix which is initialized with the sparse matrix topology.
	// D stores the mass matrix.
	// L stores the stiffness matrix.
	// The stiffness matrix is unchanged throughout the flow so we can initialize it once.
	SparseMatrix< Real , int > D , L , M;
	M = L = mesh.stiffnessMatrix();
	if( Randomize.set ) for( int i=0 ; i<vertices.size() ; i++ ) vertices[i] = RandomSpherePoint< Real >();
	EigenSolverCholeskyLLt< Real >* solver = NULL;

	auto SetStats = [&]( Real &deformationScale , Real &quasiConformalRatio , Real &radialDeviation )
	{
		Real differenceNorm=0 , oldNorm=0;
#pragma omp parallel for
		for( int j=0 ; j<D.rows ; j++ )
		{
			differenceNorm += Point3D< Real >::Dot( oldX[j]-vertices[j] , oldX[j]-vertices[j] ) * D[j][0].Value;
			oldNorm += Point3D< Real >::Dot( oldX[j] , oldX[j] ) * D[j][0].Value;
		}
		deformationScale = sqrt( differenceNorm / oldNorm );
		quasiConformalRatio = qcRatio( vertices );
		radialDeviation = mesh.radialDeviation();
	};
	Timer timer;
	for( int i=0 ; i<iters ; i++ )
	{
		double t;
		double tt = Timer::Time();

		t = Timer::Time();
		D = mesh.massMatrix( lump );
		double mTime = Timer::Time()-t;
		// Set the new constraint vector: b = D * x
		D.MultiplyParallel( &vertices[0] , &b[0] , Threads.value , MULTIPLY_CLEAR );
		// Set the system matrix: M = D + t * L
#pragma omp parallel for num_threads( Threads.value )
		for( int j=0 ; j<M.rows ; j++ )
		{
			for( int k=0 ; k<M.rowSizes[j] ; k++ ) M[j][k].Value  = L[j][k].Value * stepSize;
			for( int k=0 ; k<D.rowSizes[j] ; k++ ) M[j][k].Value += D[j][k].Value;
		}

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
			for( int j=0 ; j<b.size() ; j++ ) _b[j] = b[j][d] , _x[j] = oldX[j][d] = vertices[j][d];
			solver->solve( &_b[0] , &_x[0] );
#pragma omp parallel for num_threads( Threads.value )
			for( int j=0 ; j<vertices.size() ; j++ ) vertices[j][d] = _x[j];
		}
		double sTime2 = Timer::Time()-t;

		double mem = MemoryInfo::UsageMB();
		// Re-scale/center the mesh
		mesh.pose() , mesh.makeUnitArea();

		if( FullVerbose.set )
		{
			tt = Timer::Time() - tt;
			Real deformationScale , quasiConformalRatio , radialDeviation;
			SetStats( deformationScale , quasiConformalRatio , radialDeviation );
			printf( "CMCF[%d / %d] %4.2f = %4.2f + %4.2f + %4.2f(s): D-Norm=%6.5f / QC-Ratio=%6.5f / R-Deviation=%6.5f\r" , i+1 , iters , tt , mTime , sTime1 , sTime2 , deformationScale , quasiConformalRatio , radialDeviation );
		}
	}
	if( FullVerbose.set && iters>0 ) printf( "\n" );
	if( Verbose.set && iters>0 )
	{
		Real deformationScale , quasiConformalRatio , radialDeviation;
		SetStats( deformationScale , quasiConformalRatio , radialDeviation );
		printf( "CMCF[%d] %4.2f(s): D-Norm=%6.5f / QC-Ratio=%6.5f / R-Deviation=%6.5f\n" , iters , timer.elapsed() , deformationScale , quasiConformalRatio , radialDeviation );
	}
	mesh.pose();

	for( int i=0 ; i<vertices.size() ; i++ ) vertices[i] = vertices[i] * scale + center;
	if( solver ) delete solver;
}

template< typename Real >
void Center( SphericalGeometry::Mesh< Real >& mesh )
{
	static const int MAX_ITERATIONS = 50;
	int iters;
	int verbose = FullVerbose.set ? 2 : ( Verbose.set ? 1 : 0 );
	switch( CenterToInversion.value )
	{
		case 0:  iters = mesh.normalize( MAX_ITERATIONS , CutOff.value , true , typename SphericalGeometry::Mesh< Real >::CenterToInversion() , verbose ) ; break;
		case 1:  iters = mesh.normalize( MAX_ITERATIONS , CutOff.value , true , typename SphericalGeometry::Mesh< Real >::GSSCenterToInversion( (Real)GSSTolerance.value ) , verbose ) ; break;
		default: iters = mesh.normalize( MAX_ITERATIONS , CutOff.value , true , typename SphericalGeometry::Mesh< Real >::PoincareCenterToInversion( (Real)PoincareMaxNorm.value ) , verbose ) ; break;
	}
	if( iters==MAX_ITERATIONS && CenterToInversion.value!=1 )
	{
		fprintf( stderr , "[WARNING] Failed to meet centering threshold after %d iterations, trying golden-section search\n" , MAX_ITERATIONS );
		iters = mesh.normalize( MAX_ITERATIONS , CutOff.value , true , typename SphericalGeometry::Mesh< Real >::GSSCenterToInversion( (Real)GSSTolerance.value ) , verbose );
	}
	if( iters==MAX_ITERATIONS ) fprintf( stderr , "[WARNING] Failed to meet centering threshold after %d iterations\n" , MAX_ITERATIONS );
}
template< unsigned int SHDegree , typename Real >
void SHCenter( SphericalGeometry::Mesh< Real >& mesh )
{
	static const int MAX_ITERATIONS = 50;
	int verbose = FullVerbose.set ? 2 : ( Verbose.set ? 1 : 0 );
	if( mesh.template normalizeSH< SHDegree >( MAX_ITERATIONS , AdvectionSteps.value , AdvectionStepSize.value , CutOff.value , true , verbose )==MAX_ITERATIONS ) fprintf( stderr , "[WARNING] Failed to meet centering threshold after %d iterations\n" , MAX_ITERATIONS );
}

template< class Real >
void Execute( SphericalGeometry::Mesh< Real > &mesh , const std::vector< Point3D< Real > > &vertices , const std::vector< Point3D< Real > > &colors )
{
	// Center the spherical mesh
	if( !NoCenter.set )
	{
		Timer t;
		Center( mesh );
		if( SHDegree.value>=2 ) SHCenter< 2 >( mesh );
		if( SHDegree.value>=3 ) SHCenter< 3 >( mesh );
		if( SHDegree.value>=4 ) SHCenter< 4 >( mesh );
	}

	if( Out.set )
	{
		if( Spherical.set )
		{
			std::vector< PlyColorVertex< float > > outVertices( vertices.size() );
			for( int i=0 ; i<vertices.size() ; i++ )
			{
				outVertices[i].point = Point3D< float >( mesh.vertices[i] );
				outVertices[i].color = Point3D< float >( colors[i] );
			}
			PlyWritePolygons( Out.value , outVertices , mesh.polygons , PlyColorVertex< float >::WriteProperties , PlyColorVertex< float >::WriteComponents , ASCII.set ? PLY_ASCII : PLY_BINARY_NATIVE , NULL , 0 );
		}
		else
		{
			std::vector< PlyParametrizedColorVertex< float , Real > > outVertices( vertices.size() );
			for( int i=0 ; i<vertices.size() ; i++ )
			{
				outVertices[i].point = Point3D< float >( vertices[i] );
				outVertices[i].color = Point3D< float >( colors[i] );
				outVertices[i].param = mesh.vertices[i];
			}
			PlyWritePolygons( Out.value , outVertices , mesh.polygons , PlyParametrizedColorVertex< float , Real >::WriteProperties , PlyParametrizedColorVertex< float , Real >::WriteComponents , ASCII.set ? PLY_ASCII : PLY_BINARY_NATIVE , NULL , 0 );
		}
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
		Timer timer;

		std::vector< VertexAndColor > verticesAndColors( vertices.size() );
		for( int i=0 ; i<vertices.size() ; i++ ) verticesAndColors[i].vertex = Point3D< Real >( vertices[i] ) , verticesAndColors[i].color = Point3D< Real >( colors[i] );
		SphericalGeometry::Tessellation< Real > tessellator( mesh , verticesAndColors , Resolution.value );
		if( Collapse.set ) tessellator.collapseCells( verticesAndColors );

		if( OutTessellation.set )
		{
			std::vector< Point3D< Real > > _vertices( verticesAndColors.size() ) , _colors( verticesAndColors.size() );
			for( int i=0 ; i<verticesAndColors.size() ; i++ ) _vertices[i] = verticesAndColors[i].vertex , _colors[i] = verticesAndColors[i].color;
			tessellator.write( OutTessellation.value , _vertices , _colors , !ASCII.set );
		}

		SphericalGrid< Real > sGrid;
		sGrid.resize( Resolution.value );
		tessellator.createSGrid( sGrid , Smooth.value );

		if( OutGrid.set )
		{
			char* ext = GetFileExtension( OutGrid.value );
			if     ( !strcasecmp( ext , "ply" ) )
			{
				if( NoGridScale.set ) tessellator.writeGrid( OutGrid.value , Smooth.value , []( Real v ){ return (Real)1; } , !ASCII.set );	
				else                  tessellator.writeGrid( OutGrid.value , Smooth.value , []( Real v ){ return v<0 ? (Real)-sqrt(-v) : (Real)sqrt(v); } , !ASCII.set );
			}
			else if( !strcasecmp( ext , "sgrid" ) )
			{
				SphericalGrid< float > _sGrid;
				_sGrid.resize( sGrid.resolution() );
				for( int i=0 ; i<sGrid.resolution() ; i++ ) for( int j=0 ; j<sGrid.resolution() ; j++ ) _sGrid(i,j) = (float)sGrid(i,j);
				_sGrid.write( OutGrid.value );
			}
			delete[] ext;
		}
		if( Verbose.set ) printf( "Tessellated: %.2f(s)\n" , timer.elapsed() );
	}
}
template< class Real >
void Execute( TriangleMesh< Real > &mesh , const std::vector< Point3D< Real > > &vertices , const std::vector< Point3D< Real > > &colors )
{
	SphericalGeometry::Mesh< Real > sMesh;
	{
		// Convert the triangles to polygons
		std::vector< std::vector< int > > polygons( mesh.triangles.size() );
		for( int i=0 ; i<mesh.triangles.size() ; i++ )
		{
			polygons[i].resize(3);
			for( int j=0 ; j<3 ; j++ ) polygons[i][j] = mesh.triangles[i][j];
		}
		sMesh = SphericalGeometry::Mesh< Real >( vertices , mesh.vertices , polygons );
		sMesh.makeUnitMass();
	}
	Execute( sMesh , vertices , colors );
}
template< class Real >
void Execute( PolygonMesh< Real > &mesh , const std::vector< Point3D< Real > > &vertices , const std::vector< Point3D< Real > > &colors )
{
	SphericalGeometry::Mesh< Real > sMesh;
	{
		sMesh = SphericalGeometry::Mesh< Real >( vertices , mesh.vertices , mesh.polygons );
		sMesh.makeUnitMass();
	}
	Execute( sMesh , vertices , colors );
}
template< class Real >
void Execute( TriangulatedPolygonMesh< Real > &mesh , const std::vector< Point3D< Real > > &vertices , const std::vector< Point3D< Real > > &colors )
{
	SphericalGeometry::Mesh< Real > sMesh;
	{
		sMesh = SphericalGeometry::Mesh< Real >( vertices , mesh.vertices , mesh.polygons );
		sMesh.makeUnitMass();
	}
	Execute( sMesh , vertices , colors );
}
template< typename Real >
void HoleFillCenterTriangulation( Mesh< Real > &mesh , const std::vector< int > &loop )
{
	Point3D< Real > c;
	for( int j=0 ; j<loop.size() ; j++ ) mesh.push_back( loop[j] , loop[(j+1)%loop.size()] , (int)mesh.vertices.size() ) , c += mesh.vertices[ loop[j] ];
	c /= (int)loop.size();
	mesh.vertices.push_back( c );
}
template< typename Real >
void HoleFillMinimalAreaTriangulation( Mesh< Real > &mesh , const std::vector< int > &loop )
{
	MinimalAreaTriangulation< Real > MAT;
	std::vector< TriangleIndex > triangles( loop.size()-2 );
	std::vector< Point3D< Real > > vertices( loop.size() );
	for( int j=0 ; j<loop.size() ; j++ ) vertices[j] = mesh.vertices[ loop[j] ];
	MAT.GetTriangulation( vertices , triangles );
	for( int j=0 ; j<triangles.size() ; j++ ) mesh.push_back( loop[ triangles[j][0] ] , loop[ triangles[j][1] ] , loop[ triangles[j][2] ] );;
}
template< typename Real > void HoleFillPolygon(            TriangleMesh< Real >  &tMesh , const std::vector< int > &loop ){ fprintf( stderr , "[WARNING] Polygon hole-filling not supported for triangle mesh\n" ); }
template< typename Real > void HoleFillPolygon(             PolygonMesh< Real >  &pMesh , const std::vector< int > &loop ){  pMesh.polygons.push_back( loop ); }
template< typename Real > void HoleFillPolygon( TriangulatedPolygonMesh< Real > &tpMesh , const std::vector< int > &loop ){ tpMesh.polygons.push_back( loop ); }

template< class Real , typename Mesh >
void Execute( Mesh &mesh , const std::vector< Point3D< Real > > &colors )
{
	int fCount = mesh.faces() , vCount = (int)mesh.vertices.size();
	if( HoleFillType.value!=HOLE_FILL_TYPE_NONE )
	{
		Timer t;
		std::vector< std::vector< int > > boundaryLoops = BoundaryLoops( mesh );
		for( int i=0 ; i<boundaryLoops.size() ; i++ )
		{
			if( FullVerbose.set ) printf( "\t%d] Hole size: %d\n" , i+1 , (int)boundaryLoops[i].size() );
			switch( HoleFillType.value )
			{
				case HOLE_FILL_TYPE_CENTER_TRIANGULATION:       HoleFillCenterTriangulation     ( mesh , boundaryLoops[i] ) ; break;
				case HOLE_FILL_TYPE_MINIMAL_AREA_TRIANGULATION: HoleFillMinimalAreaTriangulation( mesh , boundaryLoops[i] ) ; break;
				case HOLE_FILL_TYPE_POLYGON:                    HoleFillPolygon                 ( mesh , boundaryLoops[i] ) ; break;
				default: fprintf( stderr , "[WARNING] Unrecognizd hole-fill type:%d\n" , HoleFillType.value );
			}
		}
		mesh.update();
		if( Verbose.set ) printf( "Filled in %d holes: %.2f(s)\n" , (int)boundaryLoops.size() , t.elapsed() );
	}
	std::vector< Point3D< Real > > vertices = mesh.vertices;

	CMCF< Real >( mesh , Iterations.value , Lump.set , fCount );

	if( mesh.faces()!=fCount || mesh.vertices.size()!=vCount )
	{
		while( mesh.faces()>fCount ) mesh.pop_back();
		while( mesh.vertices.size()>vCount ) mesh.vertices.pop_back();
		mesh.update();
	}

	// Normalize the parameterization
	{
		Point3D< Real > center = mesh.center();
#pragma omp parallel for
		for( int i=0 ; i<mesh.vertices.size() ; i++ ) mesh.vertices[i] -= center;
#pragma omp parallel for
		for( int i=0 ; i<mesh.vertices.size() ; i++ ) mesh.vertices[i] /= Length( mesh.vertices[i] );
	}

	// Orient the spherical parameterization so that it's outward facing
	if( !NoOrient.set )
	{
		Real pArea = 0 , nArea = 0;
#pragma omp parallel for reduction ( + : pArea , nArea )
		for( int f=0 ; f<mesh.faces() ; f++ )
		{
			Point3D< Real > c = mesh.center(f);
			Point3D< Real > n = mesh.normal(f);
			Real l = (Real)Length(n);
			if( Point3D< Real >::Dot( c , n )>0 ) pArea += l;
			else                                  nArea += l;
		}
		if( nArea>pArea )
#pragma omp parallel for
			for( int i=0 ; i<mesh.vertices.size() ; i++ ) mesh.vertices[i] = - mesh.vertices[i];
	}
	Execute( mesh , vertices , colors );
}

int main( int argc , char* argv[] )
{
	typedef double Real;

	cmdLineParse( argc-1 , argv+1 , params );
	if( !In.set )
	{
		Usage( argv[0] );
		return EXIT_FAILURE;
	}
	if( FullVerbose.set ) Verbose.set = true;
	MeshType.value -= 1;
	if( HoleFillType.value<0 || HoleFillType.value>=HOLE_FILL_TYPE_COUNT ) HoleFillType.value = HOLE_FILL_TYPE_NONE;
	if( MeshType.value==MESH_TYPE_TRIANGLE && HoleFillType.value==HOLE_FILL_TYPE_POLYGON ) fprintf( stderr , "[WARNING] Polygon hole-filling not supported for triangle meshes. Disabling hole-filling\n" ) , HoleFillType.value = HOLE_FILL_TYPE_NONE;
	int fileType;
	TriangleMesh< Real > tMesh;
	TriangulatedPolygonMesh< Real > tpMesh;
	PolygonMesh< Real > pMesh;
	std::vector< Point3D< Real > > colors;
	std::vector< PlyColorVertex< float > > vertices;
	bool propertyFlags[ PlyColorVertex< float >::ReadComponents ];
	switch( MeshType.value )
	{
	case MESH_TYPE_POLYGON:              PlyReadPolygons ( In.value , vertices ,  pMesh.polygons , PlyColorVertex< float >::ReadProperties , propertyFlags , PlyColorVertex< float >::ReadComponents , fileType ) ; break;
	case MESH_TYPE_TRIANGULATED_POLYGON: PlyReadPolygons ( In.value , vertices , tpMesh.polygons , PlyColorVertex< float >::ReadProperties , propertyFlags , PlyColorVertex< float >::ReadComponents , fileType ) ; break;
	default:                             PlyReadTriangles( In.value , vertices , tMesh.triangles , PlyColorVertex< float >::ReadProperties , propertyFlags , PlyColorVertex< float >::ReadComponents , fileType );
	}
	bool hasColor = (propertyFlags[3]||propertyFlags[6]) && (propertyFlags[4]||propertyFlags[7]) && (propertyFlags[5]||propertyFlags[8]);

	colors.resize( vertices.size() );
	if( hasColor ) for( int i=0 ; i<vertices.size() ; i++ ) colors[i] = Point3D< Real >( vertices[i].color );

	if( MeshType.value==MESH_TYPE_TRIANGULATED_POLYGON )
	{
		tpMesh.update();
		tpMesh.vertices.resize( vertices.size() );
		for( int i=0 ; i<vertices.size() ; i++ ) tpMesh.vertices[i] = Point3D< Real >( vertices[i].point );
		if( !hasColor )
		{
			for( int p=0 ; p<tpMesh.polygons.size() ; p++ )
			{
				Point3D< Real > n = tpMesh.normal(p);
				for( int j=0 ; j<tpMesh.polygons[p].size() ; j++ ) colors[ tpMesh.polygons[p][j] ] += n;
			}
#pragma omp parallel for
			for( int i=0 ; i<colors.size() ; i++ ) colors[i] = NormalColor( -colors[i]/(Real)Length( colors[i] ) );
		}
		Execute( tpMesh , colors );
	}
	else if( MeshType.value==MESH_TYPE_POLYGON )
	{
		pMesh.update();
		pMesh.vertices.resize( vertices.size() );
		for( int i=0 ; i<vertices.size() ; i++ ) pMesh.vertices[i] = Point3D< Real >( vertices[i].point );
		if( !hasColor )
		{
			for( int p=0 ; p<pMesh.polygons.size() ; p++ )
			{
				Point3D< Real > n = pMesh.normal(p);
				for( int j=0 ; j<pMesh.polygons[p].size() ; j++ ) colors[ pMesh.polygons[p][j] ] += n;
			}
#pragma omp parallel for
			for( int i=0 ; i<colors.size() ; i++ ) colors[i] = NormalColor( -colors[i]/(Real)Length( colors[i] ) );
		}
		Execute( pMesh , colors );
	}
	else
	{
		tMesh.vertices.resize( vertices.size() );
		for( int i=0 ; i<vertices.size() ; i++ ) tMesh.vertices[i] = Point3D< Real >( vertices[i].point );
		if( !hasColor )
		{
			for( int t=0 ; t<tMesh.triangles.size() ; t++ )
			{
				Point3D< Real > n = tMesh.normal(t);
				for( int j=0 ; j<3 ; j++ ) colors[ tMesh.triangles[t][j] ] += n;
			}
#pragma omp parallel for
			for( int i=0 ; i<colors.size() ; i++ ) colors[i] = NormalColor( -colors[i]/(Real)Length( colors[i] ) );
		}
		Execute( tMesh , colors );
	}

	return EXIT_SUCCESS;
}
