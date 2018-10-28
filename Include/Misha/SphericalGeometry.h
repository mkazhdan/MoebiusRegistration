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

#ifndef SPHERICAL_GEOMETRY_INCLUDED
#define SPHERICAL_GEOMETRY_INCLUDED

#include <iostream>
#include "Misha/Ply.h"
#include "Misha/SphericalGrid.h"
#include "Misha/SphericalHarmonics.h"

template< class Real > std::ostream& operator << ( std::ostream& os , const Point2D< Real > p ){ return os << "(" << p[0] << "," << p[1] << ")"; }
template< class Real > std::ostream& operator << ( std::ostream& os , const Point3D< Real > p ){ return os << "(" << p[0] << "," << p[1] << "," << p[2] << ")"; }

namespace SphericalGeometry
{
	template< class Real > Point2D< Real > StereographicProjection( Point3D< Real > p ){ return Point2D< Real >( p[0] / ( 1-p[2] ) , p[1] / ( 1-p[2] ) ); }
	template< class Real > Point3D< Real > IStereographicProjection( Point2D< Real > p ){ return Point3D< Real >( 2*p[0] , 2*p[1] , Point2D< Real >::SquareNorm( p )-1 ) / ( Point2D< Real >::SquareNorm( p ) + 1 ) ; }
	template< typename Real > Real Area( const std::vector< Point3D< Real > > &vertices , const std::vector< int > &polygon );
	template< typename Real > Real Area( const std::vector< Point3D< Real > > &vertices );
	unsigned long long Key( int v1 , int v2 );

	template< class Real > SquareMatrix< Real , 3 > Correlate( SphericalGrid< Real >& source , SphericalGrid< Real >& target , Real& error , bool gradientDomain , bool invert );

	template< class Real >
	struct SphericalInversion
	{
		Point3D< Real > center;
		SphericalInversion( Point3D< Real > c=Point3D< Real >() ) : center(c) {}
		Point3D< Real > operator() ( Point3D< Real > p ) const;
	};
	template< class Real > std::ostream& operator << ( std::ostream& os , const SphericalInversion< Real > si ){ return os << si.center; }

	template< class Real >
	struct FractionalLinearTransformation
	{
		SquareMatrix< std::complex< Real > , 2 > matrix;

		FractionalLinearTransformation( void ){ matrix(0,0) = matrix(1,1) = 1 , matrix(0,1) = matrix(1,0) = 0; }
		FractionalLinearTransformation( SphericalInversion< Real > si );
		FractionalLinearTransformation( SquareMatrix< Real , 3 > m );
		FractionalLinearTransformation( const std::pair< Point2D< Real > , Point2D< Real > > pairs[3] );

		Point2D< Real > operator()( Point2D< Real > p ) const;
		Point3D< Real > operator()( Point3D< Real > p ) const;
		FractionalLinearTransformation operator * ( FractionalLinearTransformation flt ) const;
		FractionalLinearTransformation& operator *= ( FractionalLinearTransformation flt );
		FractionalLinearTransformation inverse( void ) const;
	protected:
		// Returns the matrix taking p1 -> 0 , p2 -> infty , p3 -> 1
		static SquareMatrix< std::complex< Real > , 2 > _Transformation( Point2D< Real > p1 , Point2D< Real > p2 , Point2D< Real > p3 );

		// Returns the matrix taking 0 -> p1 , infty -> p2 , 1 -> p3
		static SquareMatrix< std::complex< Real > , 2 > _ITransformation( Point2D< Real > p1 , Point2D< Real > p2 , Point2D< Real > p3 );
		static const Point3D< Real > _ZERO , _INFINITY , _ONE;
	};
	template< class Real > std::ostream& operator << ( std::ostream& os , const FractionalLinearTransformation< Real > fls ){ return os << "[ " << fls.matrix(0,0) << " * z + " << fls.matrix(1,0) << " ] / [ " << fls.matrix(0,1) << " * z + " << fls.matrix(1,1) << " ] "; }

	template< class Real >
	struct Mesh
	{
		std::vector< std::vector< int > > polygons;
		std::vector< Point3D< Real > > vertices;
		std::vector< Real > masses;

		Mesh( void ){}
		Mesh( const std::vector< Point3D< Real > >& vertices , const std::vector< Point3D< Real > >& sVertices , const std::vector< std::vector< int > > &polygons );

		void write( const char* fileName , const std::vector< Point3D< Real > >& vertices ,                                                bool binary=true ) const;
		void write( const char* fileName , const std::vector< Point3D< Real > >& vertices , const std::vector< Point3D< Real > >& colors , bool binary=true ) const;
		void read( const char* fileName , std::vector< Point3D< Real > >& vertices ,                                          bool verbose=false );
		bool read( const char* fileName , std::vector< Point3D< Real > >& vertices , std::vector< Point3D< Real > >& colors , bool verbose=false );

		double area( int p ) const;
		Point3D< Real > center( int p ) const;
		template< typename F > double area( int p , F f ) const;
		template< typename F > Point3D< Real > center( int p , F f ) const;
		Real constant( void ) const;
		Point3D< Real > center( void ) const;
		template< typename F > Point3D< Real > center( F f ) const;
		SquareMatrix< Real , 3 > dCenter( void ) const;
		template< typename F > SquareMatrix< Real , 3 > dCenter( F f ) const;

		template< unsigned int SHDegree > Point< Real , SphericalHarmonics::Dimension< SHDegree >() > centerSH( void ) const;
		template< unsigned int SHDegree > SquareMatrix< Real , SphericalHarmonics::Dimension< SHDegree >() > dCenterSH( void ) const;
		template< typename VF > void advect( VF vf , int steps );

		template< typename F > Mesh &operator *= ( F f );

		void makeUnitMass( void );

		struct CenterToInversion
		{
			virtual SphericalInversion< Real > operator()( const Mesh& mesh , Point3D< Real > c ) const { return SphericalInversion< Real >( c ); }
		};
		struct GSSCenterToInversion : public CenterToInversion
		{
			Real tolerance , maxNorm;
			GSSCenterToInversion( Real t , Real mn=(Real)0.99 ) : tolerance(t) , maxNorm(mn) {}
			SphericalInversion< Real > operator()( const Mesh& mesh , Point3D< Real > c ) const
			{
				Real len = (Real)Length(c);
				if( len>=maxNorm ) c *= maxNorm / len;
				_NormalizationFunctor nf( mesh , c );
				return SphericalInversion< Real >( c * GoldenSectionSearch( nf , (Real)0 , maxNorm , tolerance/len ).second );
			}
		};
		struct PoincareCenterToInversion : public CenterToInversion
		{
			Real maxNorm;
			PoincareCenterToInversion( Real mn=(Real)2 ) : maxNorm(mn) {}
			SphericalInversion< Real > operator()( const Mesh& mesh , Point3D< Real > c ) const
			{
				Real len = (Real)Length(c);
				c /= len;
				len = std::min< Real >( len , maxNorm );
				return SphericalInversion< Real >( c * tanh(len) );
			}
		};

		FractionalLinearTransformation< Real > normalizer( int iters , double cutOff , bool gaussNewton , bool verbose=false ) const;
		int normalize( int iters , double cutOff , bool gaussNewton , const CenterToInversion& c2i=CenterToInversion() , int verbose=0 );

		template< unsigned int SHDegree >
		int normalizeSH( int iters , int advectionSteps , Real advectionStepSize , double cutOff , bool gaussNewton , int verbose=0 );

		static Point3D< Real > SphericalInvert( Point3D< Real > p , Point3D< Real > c );
	protected:
		struct _NormalizationFunctor
		{
			const SphericalGeometry::Mesh< Real >& mesh;
			Point3D< Real > dir;
			_NormalizationFunctor( const SphericalGeometry::Mesh< Real >& m , Point3D< Real > d ) : mesh(m) , dir(d) {}
			Real operator()( Real s ) const { return Point3D< Real >::SquareNorm( mesh.center( SphericalInversion< Real >( dir*s ) ) ); }
		};
	};

	template< class Real >
	class Polygon
	{
	public:
		struct PlanarSources
		{
			int idx[2];
			PlanarSources( void ){ idx[0] = idx[1] = -1; }
			static int SharedIndex( const PlanarSources &ps1 , const PlanarSources &ps2 )
			{
				if( ps1.idx[0]>=0 )
					if( ps1.idx[0]==ps2.idx[0] || ps1.idx[0]==ps2.idx[1] ) return ps1.idx[0];
					else if( ps1.idx[1]>=0 ) if( ps1.idx[1]==ps2.idx[0] || ps1.idx[1]==ps2.idx[1] ) return ps1.idx[1];
					return -1;
			}
			void addIndex( int i )
			{
				if     ( idx[0]==-1 ) idx[0] = i;
				else if( idx[1]==-1 ) idx[1] = i;
				else if( i>=0 ) fprintf( stderr , "[WARNING] PlanarSources::addIndex: both indices assigned: %d %d <- %d\n" , idx[0] , idx[1] , i );
			}
			bool full( void ) const { return idx[1]>=0; }
		};
	protected:
		std::vector< int > _vertices;

		template< class VertexData >
		void _split( Point3D< Real > pNormal , Real pOffset , int pIdx , Polygon &back , Polygon &front , std::vector< Point3D< Real > >& vertices , std::vector< PlanarSources > &pSources , std::vector< VertexData >* vData , std::unordered_map< unsigned long long , int >& vMap ) const;
	public:
		int sourceID , theta , phi;
		Real mass;
		Polygon( int id=-1 , Real m=0 ) : sourceID(id) , mass(m) { theta = phi = -1; }
		Polygon( int id , Real m , const std::vector< int > &vIndices ) : sourceID(id) , mass(m) , _vertices( vIndices ) , theta(-1) , phi(-1) { }
		int& operator[] ( int idx ) { return _vertices[(idx+_vertices.size() )%_vertices.size()]; }
		const int& operator[] ( int idx ) const { return _vertices[(idx+_vertices.size() )%_vertices.size()]; }
		void push_back( int idx ){ _vertices.push_back(idx); }
		size_t size( void ) const { return _vertices.size(); }

		void split( Point3D< Real > pNormal , Real pOffset , int pIdx , Polygon &back , Polygon &front , std::vector< Point3D< Real > >& vertices , std::vector< PlanarSources > &pSources ,                                    std::unordered_map< unsigned long long , int >& vMap ) const { _split( pNormal , pOffset , pIdx , back , front , vertices , pSources , (std::vector< Real >* )NULL , vMap ); }
		template< class VertexData >
		void split( Point3D< Real > pNormal , Real pOffset , int pIdx , Polygon &back , Polygon &front , std::vector< Point3D< Real > >& vertices , std::vector< PlanarSources > &pSources , std::vector< VertexData >& vData , std::unordered_map< unsigned long long , int >& vMap ) const { _split( pNormal , pOffset , pIdx , back , front , vertices , pSources , &vData                      , vMap ); }

		double area( const std::vector< Point3D< Real > >& vertices ) const;
	};

	template< class Real >
	class Tessellation 
	{
	protected:
		int _resolution;
		std::vector< Point3D< Real > > _vertices;
		std::vector< Polygon< Real > > _polygons;
		template< class VertexData > void _splitPolygon   ( Polygon< Real > p ,             std::unordered_map< unsigned long long , int >& vMap , std::vector< typename Polygon< Real >::PlanarSources > &pSources , std::vector< VertexData >* vData ); 
		template< class VertexData > void _splitPolygonPhi( Polygon< Real > p , int theta , std::unordered_map< unsigned long long , int >& vMap , std::vector< typename Polygon< Real >::PlanarSources > &pSources , std::vector< VertexData >* vData );
		template< class VertexData > void _collapseCells( std::vector< VertexData >* vData );
		template< class VertexData > void _init( const Mesh< Real >& mesh , std::vector< VertexData >* vData , int resolution );
	public:
		Tessellation( const Mesh< Real >& mesh , int resolution ){ _init( mesh , (std::vector< char >*)NULL , resolution ); }
		template< class VertexData > Tessellation( const Mesh< Real >& mesh , std::vector< VertexData >& vData , int resolution ){ _init( mesh , &vData , resolution ); }

		void createSGrid( SphericalGrid< Real >& sGrid , Real smoothValue ) const;

		const std::vector< Polygon< Real > >& polygons( void ) const { return _polygons; }
		const std::vector< Point3D< Real > >& vertices( void ) const { return _vertices; }

		void collapseCells( void ){ _collapseCells( (std::vector< char >*)NULL ); }
		template< class VertexData > void collapseCells( std::vector< VertexData >& vData ){ _collapseCells( &vData ); }

		void write( const char* fileName , bool binary=true ) const;
		void write( const char* fileName , const std::vector< Point3D< Real > >& vertices , const std::vector< Point3D< Real > >& colors , bool binary=true ) const;

		template< typename PointScaleFunction >
		void writeGrid( const char* fileName , Real smooth , PointScaleFunction f=[]( Real v ){ return v; } , bool binary=true ) const;

		template< typename PointScaleFunction >
		static void WriteGrid( const char* fileName , const SphericalGrid< Real >& sGrid , PointScaleFunction f=[]( Real v ){ return v; } , bool binary=true );
	};

	template< class Real >
	struct CircumscribedSphere
	{
		static Point3D< Real > Center( const std::vector< Point3D< Real > >& vertices , int iters );
	protected:
		static void _RandomCenter( const std::vector< Point3D< Real > >& vertices , Point3D< Real >& c , Real& a );
	};
}
#include "Misha/SphericalGeometry.inl"
#endif // SPHERICAL_GEOMETRY_INCLUDED