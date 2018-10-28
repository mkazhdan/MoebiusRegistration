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

#include <cstdlib>
#include <vector>
#include <omp.h>

#include "Misha/Geometry.h"
#include "Misha/Algebra.h"
#include "Misha/Ply.h"
#include "Misha/CmdLineParser.h"
#include "Misha/Timer.h"
#include "Misha/SphericalGeometry.h"


enum
{
	CORRELATE_ROTATION   = 1 ,
	CORRELATE_REFLECTION = 2
};

cmdLineParameterArray< char* , 2 > In( "in" );
cmdLineParameter< char* > Out( "out" );
cmdLineParameter< int > Resolution( "res" , 256 ) , CorrelationType( "cType" , CORRELATE_ROTATION );
cmdLineParameter< float > Smooth( "smooth" , 5e-4f );
cmdLineReadable Correspondence( "correspondence" ) , Verbose( "verbose" ) , GradientDomain( "gradientDomain" );

void Usage( const char* ex ) 
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input source/target>\n" , In.name );
	printf( "\t[--%s <output mesh>]\n" , Out.name );
	printf( "\t[--%s <spherical resolution>=%d]\n" , Resolution.name , Resolution.value );
	printf( "\t[--%s <spherical smoothing>=%g]\n" , Smooth.name , Smooth.value );
	printf( "\t[--%s <correlation type>=%d]\n" , CorrelationType.name , CorrelationType.value );
	printf( "\t\t%d] Rotation\n" , CORRELATE_ROTATION );
	printf( "\t\t%d] Reflection\n" , CORRELATE_REFLECTION );
	printf( "\t\t%d] Rotation + Reflection\n" , CORRELATE_ROTATION | CORRELATE_REFLECTION );
	printf( "\t[--%s]\n" , GradientDomain.name );
	printf( "\t[--%s]\n" , Correspondence.name );
	printf( "\t[--%s]\n" , Verbose.name );
}

cmdLineReadable* params[] = { &In , &Out , &Verbose , &Resolution , &Smooth , &GradientDomain , &CorrelationType , &Correspondence , NULL };

template< class Real >
Point3D< Real > NearestPointOnEdge( Point3D< Real > point , const Point3D< Real > edge[2] , Real& b0 , Real& b1 )
{
	Point3D< Real > d = edge[1] - edge[0] , p = point - edge[0];
	Real dot = Point3D< Real >::Dot( p , d );
	if( dot<0 ) 
	{
		b0 = 1.;
		return edge[0];
	}
	else if( dot>Point3D< Real >::SquareNorm( d ) ) { 
		b1 = 1.; 
		return edge[1];
	}
	else
	{
		// Solve for the minimizer of:
		//                            E(t) = || p - t*d ||^2
		//                                 = || p ||^2 - 2 t < p , d > + t^2 || d ||^2
		//            =>             0 = -< p , d > + t || d ||^2
		//            <=>    t = < p , d > / || d ||^2
		Real t = dot / Point3D< Real >::SquareNorm( d );
		b0 = 1.-t , b1 = t;
		return edge[0] + d * t;
	}
}
template< class Real >
Point3D< Real > NearestPointOnTriangle( Point3D< Real > point , const Point3D< Real > triangle[3] , Real* b )
{

	b[0] = b[1] = b[2] = 0;
	Point3D< Real > d[] = { triangle[1]-triangle[0] , triangle[2]-triangle[0] } , p = point - triangle[0] , n = CrossProduct( d[0] , d[1] );

	if( !Length(n) ) return ( triangle[0] + triangle[1] + triangle[2] ) / (Real)3.;


	if     ( Point3D< Real >::Dot( point-triangle[0] , CrossProduct( n , triangle[1]-triangle[0] ) )<0 ){ Point3D< Real > edge[] = { triangle[0] , triangle[1] } ; return NearestPointOnEdge( point , edge , b[0] , b[1] ); }
	else if( Point3D< Real >::Dot( point-triangle[1] , CrossProduct( n , triangle[2]-triangle[1] ) )<0 ){ Point3D< Real > edge[] = { triangle[1] , triangle[2] } ; return NearestPointOnEdge( point , edge , b[1] , b[2] ); }
	else if( Point3D< Real >::Dot( point-triangle[2] , CrossProduct( n , triangle[0]-triangle[2] ) )<0 ){ Point3D< Real > edge[] = { triangle[2] , triangle[0] } ; return NearestPointOnEdge( point , edge , b[2] , b[0] ); }
	else
	{
		// Solve for the minimizer of:
		//                            E(s,t) = || p - s*d[0]-t*d[1] ||^2
		//                                   = || p ||^2 - 2 s < p , d[0] > - 2 t < p , d[1] > + 2 s t < d[0] , d[1] > + s^2 || d[0] ||^2 + t^2 || d[1] ||^2
		//   =>  (0,0) = ( -< p , d[0] > + t < d[0] , d[1] > + s || d[0] ||^2 , -< p , d[1] > + s < d[0] , d[1] > + t || d[1] ||^2
		//            <=> | < p , d[0] > | = | < d[0] , d[0] >   < d[0] , d[1] > | | s |
		//                | < p , d[1] > |   | < d[0] , d[1] >   < d[1] , d[1] > | | t |
		SquareMatrix< Real , 2 > M , M_inverse;
		M(0,0) = Point3D< Real >::SquareNorm( d[0] ) , M(1,0) = M(0,1) = Point3D< Real >::Dot( d[0] , d[1] ) , M(1,1) = Point3D< Real >::SquareNorm( d[1] );
		Real det = M(0,0)*M(1,1) - M(0,1)*M(1,0);
		M_inverse(0,0) = M(1,1) , M_inverse(0,1) = -M(0,1) , M_inverse(1,0) = -M(1,0) , M_inverse(1,1) = M(0,0);
		M_inverse /= det;
		Point2D< Real > st = M_inverse * Point2D< Real >( Point3D< Real >::Dot( p , d[0] ) , Point3D< Real >::Dot( p , d[1] ) );
		// Sanity check
		if( st[0]<0 || st[1]<0 || st[0]+st[1]>1 ) fprintf( stderr , "[ERROR] Bad barycentric coordinates: %g %g %g (%g)\n" , 1.-st[0]-st[1] , st[0] , st[1] , Length(n) ) , exit( 0 );
		b[0] = 1. - st[0] - st[1] , b[1] = st[0] , b[2] = st[1];
		Point3D< Real > ret = triangle[0] * ( Real )( 1. - st[0] - st[1] ) + d[0] * st[0] + d[1] * st[1];

		return triangle[0] * ( Real )( 1. - st[0] - st[1] ) + d[0] * st[0] + d[1] * st[1];
	}
}
template< class Real >
Point3D< Real > NearestPointOnPolygon( Point3D< Real > point , const std::vector< Point3D< Real > > &vertices , TriangleIndex &tri , Real* b )
{
	MinimalAreaTriangulation< Real > MAT;
	std::vector< TriangleIndex > triangles;
	MAT.GetTriangulation( vertices , triangles );

	Real dist=-1 , _b[3];
	for( int i=0 ; i<triangles.size() ; i++ )
	{
		Point3D< Real > triangle[] = { vertices[ triangles[i][0] ] , vertices[ triangles[i][1] ] , vertices[ triangles[i][2] ] };
		Point3D< Real > p = NearestPointOnTriangle( point , triangle , _b );
		Real d2 = Point3D< Real >::SquareNorm( p - point );
		if( dist<0 || d2<dist )
		{
			dist = d2;
			b[0] = _b[0] , b[1] = _b[1] , b[2] = _b[2];
			tri = triangles[i];
		}
	}
	if( dist<0 ) fprintf( stderr , "[ERROR] NearestPointOnTriangle: could not find nearest point\n" ) , exit( 0 );
	return vertices[ tri[0] ] * b[0] + vertices[ tri[1] ] * b[1] + vertices[ tri[2] ] * b[2];
}

template< class Real >
void Execute( void ) 
{
	bool hasGeometry;
	SphericalGrid< Real > sourceGrid , targetGrid;
	SphericalGeometry::Mesh< Real > sourceMesh , targetMesh;
	std::vector< Point3D< Real > > sourceVertices , sourceColors , targetVertices , targetColors;
	std::vector< std::vector< SphericalGeometry::Polygon< Real > > > targetCells( Resolution.value * Resolution.value );

	char* sExt = GetFileExtension( In.values[0] );
	char* tExt = GetFileExtension( In.values[1] );

	// If the inputs are both spherical grids
	if( !strcasecmp( sExt , "sgrid" ) && !strcasecmp( tExt , "sgrid" ) )
	{
		hasGeometry = false;
		SphericalGrid< float > _sourceGrid , _targetGrid;
		_sourceGrid.read( In.values[0] ) , _targetGrid.read( In.values[1] );
		if( _sourceGrid.resolution()!=_targetGrid.resolution() ) fprintf( stderr , "[ERROR] Source and target grid resolutions don't match: %d != %d\n" , _sourceGrid.resolution() , _targetGrid.resolution() ) , exit( 0 );
		sourceGrid.resize( _sourceGrid.resolution() ) , targetGrid.resize( _targetGrid.resolution() );
#pragma omp parallel for
		for( int j=0 ; j<sourceGrid.resolution() ; j++ ) for( int i=0 ; i<sourceGrid.resolution() ; i++ ) sourceGrid(i,j) = (Real)_sourceGrid(i,j) , targetGrid(i,j) = (Real)_targetGrid(i,j);
	}
	// If the inputs are both parameterized meshes
	else if( !strcasecmp( sExt , "ply" ) && !strcasecmp( tExt , "ply" ) )
	{
		hasGeometry = true;
		sourceMesh.read( In.values[0] , sourceVertices , sourceColors , Verbose.set );
		targetMesh.read( In.values[1] , targetVertices , targetColors , Verbose.set );
		sourceMesh.makeUnitMass() , targetMesh.makeUnitMass();

		struct VertexAndColor
		{
			VertexAndColor( void ){}
			VertexAndColor( Point3D< Real > v , Point3D< Real > c ) : vertex(v) , color(c) {}
			Point3D< Real > vertex , color;
			VertexAndColor operator * ( Real s ) const { return VertexAndColor( vertex*s , color*s ); }
			VertexAndColor operator / ( Real s ) const { return VertexAndColor( vertex/s , color/s ); }
			VertexAndColor operator + ( const VertexAndColor& vc ) const { return VertexAndColor( vertex + vc.vertex , color + vc.color ); }
		};
		// Compute the spherical tesselations
		{
			std::vector< VertexAndColor > verticesAndColors;

			Timer t;
			targetGrid.resize( Resolution.value );
			verticesAndColors.resize( targetVertices.size() );
			for( int i=0 ; i<targetVertices.size() ; i++ ) verticesAndColors[i].vertex = targetVertices[i] , verticesAndColors[i].color = targetColors[i];
			SphericalGeometry::Tessellation< Real > targetTessellator( targetMesh , verticesAndColors , Resolution.value );
			targetVertices.resize( verticesAndColors.size() ) , targetColors.resize( verticesAndColors.size() );
			for( int i=0 ; i<verticesAndColors.size() ; i++ ) targetVertices[i] = verticesAndColors[i].vertex , targetColors[i] = verticesAndColors[i].color;
			targetTessellator.createSGrid( targetGrid , Smooth.value );
			const std::vector< SphericalGeometry::Polygon< Real > > p = targetTessellator.polygons();
			for( int i=0 ; i<targetTessellator.polygons().size() ; i++ ) targetCells[ p[i].theta*Resolution.value + p[i].phi ].push_back( p[i] );

			sourceGrid.resize( Resolution.value );
			verticesAndColors.resize( sourceVertices.size() );
			for( int i=0 ; i<sourceVertices.size() ; i++ ) verticesAndColors[i].vertex = sourceVertices[i] , verticesAndColors[i].color = sourceColors[i];
			SphericalGeometry::Tessellation< Real > sourceTessellator( sourceMesh , verticesAndColors , Resolution.value );
			sourceTessellator.createSGrid( sourceGrid , Smooth.value );

			if( Verbose.set ) printf( "Tessellated: %.2f(s)\n" , t.elapsed() );
		}
	}
	else fprintf( stderr , "[ERROR] File extensions must either both be \"sgrid\" or \"ply\": %s , %s\n" , sExt , tExt ) , exit( 0 );
	delete[] sExt;
	delete[] tExt;

	// Correlate the spherical grids and transform the source
	{
		Timer t;
		Real error , rotError , refError;
		SquareMatrix< Real , 3 > R , rotR , refR;
		if( CorrelationType.value & CORRELATE_ROTATION   ) rotR = SphericalGeometry::Correlate( sourceGrid , targetGrid , rotError , GradientDomain.set , false );
		if( CorrelationType.value & CORRELATE_REFLECTION ) refR = SphericalGeometry::Correlate( sourceGrid , targetGrid , refError , GradientDomain.set , true  );
		switch( CorrelationType.value )
		{
		case CORRELATE_ROTATION:
			R = rotR , error = rotError;
			break;
		case CORRELATE_REFLECTION:
			R = refR , error = refError;
			break;
		case CORRELATE_ROTATION | CORRELATE_REFLECTION:
			if( rotError<refError ) error = rotError , R = rotR;
			else                    error = refError , R = refR;
			break;
		default:
			fprintf( stderr , "[ERROR] Invalid correlation type: %d\n" , CorrelationType.value ) , exit( 0 );
		}
#pragma omp parallel for
		for( int i=0 ; i<sourceMesh.vertices.size() ; i++ ) sourceMesh.vertices[i] = R * sourceMesh.vertices[i];
		if( Verbose.set )  printf( "Correlated %.2f(s)\t Error: %g\n" , t.elapsed() , error );
		else printf( "Error: %g\n" , error );
	}


	if( hasGeometry && Out.set )
	{
		if( Correspondence.set )
		{
			// Re-map the source vertex positions
#pragma omp parallel for
			for( int i=0 ; i<sourceMesh.vertices.size() ; i++ )
			{
				Real theta , phi; 

				SphericalGrid< Real >::SetCoordinates( sourceMesh.vertices[i].coords , theta , phi );
				if( theta<0 ) theta += (Real)( 2. * M_PI );

				Real x = (Real)( (theta * Resolution.value) / (2.*M_PI) );
				Real y = (Real)( (phi * Resolution.value ) / M_PI );
				int _x = (int)floor( x ) , _y = (int)floor( y );

				int idx = _x * Resolution.value + _y;

				TriangleIndex tri;
				int tIdx = -1;
				Real dist = -1;
				Point3D< Real > b;
				for( int j=0 ; j<targetCells[idx].size() ; j++ )
				{
					Real _b[3];
					const std::vector< int > poly = targetMesh.polygons[ targetCells[idx][j].sourceID ];
					std::vector< Point3D< Real > > _vertices( poly.size() );
					for( int k=0 ; k<poly.size() ; k++ ) _vertices[k] = targetMesh.vertices[ poly[k] ];

					TriangleIndex _tri;
					Point3D< Real > p = NearestPointOnPolygon( sourceMesh.vertices[i] , _vertices , _tri , _b );

					Real d2 = Point3D< Real >::SquareNorm( p - sourceMesh.vertices[i] );
					if( dist<0 || d2<dist )
					{
						dist = d2;
						tri = TriangleIndex( poly[ _tri[0] ] , poly[ _tri[1] ] , poly[ _tri[2] ] );
						tIdx = targetCells[idx][j].sourceID;
						b = Point3D< Real >( _b[0] , _b[1] , _b[2] );
					}
				}
				if( dist<0 ) fprintf( stderr , "[ERROR] Could not find a polygon in cell: %d %d\n" , _x , _y ) , exit( 0 );
				sourceVertices[i] = targetVertices[ tri[0] ] * b[0] + targetVertices[ tri[1] ] * b[1] + targetVertices[ tri[2] ] * b[2];
			}
		}
		sourceMesh.write( Out.value , sourceVertices , sourceColors );
	}
}

int main( int argc , char* argv[] )
{
	cmdLineParse( argc-1 , argv+1 , params );
	if( !In.set )
	{
		Usage( argv[0] );
		return EXIT_FAILURE;
	}
	Execute< double >();
	return EXIT_SUCCESS;    
}