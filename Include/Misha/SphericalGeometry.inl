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

#include <algorithm>
#include <unordered_set>
#include <random>
#include "Misha/Fourier.h"

#define ALIGN_INVERSION_TO_CENTER

template< class Real >
Real SphericalGeometry::Area( const Point3D< Real >& v1, const Point3D< Real >& v2 , const Point3D< Real >& v3 )
{
	Point3D< Real > areaVector = Point3D< Real >::CrossProduct( v2-v1 , v3-v1 );
	return (Real)sqrt( Point3D< Real >::SquareNorm( areaVector ) ) / 2.;
};

template< class Real >
void SphericalGeometry::BarycentricCoordinates( const Point3D< Real >& p , const Point3D< Real >& v1 , const Point3D< Real >& v2, const Point3D< Real >& v3 , Real& a0 , Real& a1 , Real& a2 )
{
	Point3D< Real > p0 =  p - v1;
	Point3D< Real > p1 = v2 - v1;
	Point3D< Real > p2 = v3 - v1;
	Point3D< Real >  n  = Point3D<Real>::CrossProduct( p1 , p2 );
	Point3D< Real > _v1 = Point3D<Real>::CrossProduct( p2 , n  );
	Point3D< Real > _v2 = Point3D<Real>::CrossProduct( p1 , n  );

	_v1 /= Point3D<Real>::Dot( _v1 , p1 );
	_v2 /= Point3D<Real>::Dot( _v2 , p2 );

	a1 = Point3D<Real>::Dot( _v1 , p0 );
	a2 = Point3D<Real>::Dot( _v2 , p0 );
	a0 = Real(1.0) - a1 - a2;
}

inline unsigned long long SphericalGeometry::Key( int v1 , int v2 )
{
	if( v1<v2 ) return ( (unsigned long long)v1 )<<32 | ( (unsigned long long )v2);
	else        return ( (unsigned long long)v2 )<<32 | ( (unsigned long long )v1);
}

template< class Real >
double SphericalGeometry::TriangleArea( Point3D< Real > p1 , Point3D< Real > p2 , Point3D< Real > p3 )
{
	Point3D< double > p[] { Point3D< double >( p1 ) , Point3D< double >( p2 ) , Point3D< double >( p3 ) };
	double theta[3];
	for( int i=0 ; i<3 ; i++ )
	{
		Point3D< double > v[] = { p[(i+2)%3] - p[i] , p[(i+1)%3] - p[i] };
		for( int j=0 ; j<2 ; j++ ){ v[j] -= p[i] * Point3D< double >::Dot( v[j] , p[i] ) ; v[j] /= Length( v[j] ); }
		theta[i] = acos( std::max< Real >( -1 , std::min< Real >( 1, Point3D< double >::Dot( v[0] , v[1] ) ) ) );
	}
	return theta[0] + theta[1] + theta[2] - M_PI;
}

template< class Real >
SquareMatrix< Real , 3 > SphericalGeometry::Correlate( SphericalGrid< Real >& source , SphericalGrid< Real >& target , Real& error , bool gradientDomain , bool invert )
{
	FourierKeySO3< Real> keySO3; 
	HarmonicTransform< Real > hForm; 
	WignerTransform< Real > wForm; 
	FourierKeyS2<Real> sourceKey, targetKey; 
	RotationGrid<Real> euler; 
	Real matrix[3][3];

	hForm.ForwardFourier( source, sourceKey );
	hForm.ForwardFourier( target, targetKey );

	if( gradientDomain )
#pragma omp parallel for
		for( int i=0 ; i<sourceKey.bandWidth() ; i++ )
		{
			Real scl = (Real)sqrt( i*(i+1) );
			for( int j=0 ; j<=i ; j++ ) sourceKey(i,j) *= scl , targetKey(i,j) *= scl;
		}

	// Allocate memory 
	if( !keySO3.resize( source.resolution() ) ) fprintf( stderr , "[ERROR] Correlate: Could not allocate key: %d\n" , source.resolution() ) , exit( 0 ); 

	// Compute the aligning rotation
#pragma omp parallel for
	for( int i=0 ; i<keySO3.bandWidth(); i++ ) for( int j=0 ; j<=i ; j++ ) for( int k=0 ; k<=i ; k++ ) 
	{
		// set each Wigner-D coefficient to be the cross multiplication 
		// within each frequency 
		int sign = 1;
		if( invert && (i%2) )sign = -1;
		keySO3(i,j,k) = sourceKey(i,k) * targetKey(i,j).conjugate() * sign;
		if( k ) keySO3(i,j,-k) = sourceKey(i,k).conjugate() * targetKey(i,j).conjugate() * sign;
	}

	// Take the inverse Wigner D transform 
	wForm.InverseFourier( keySO3 , euler );

	Real max = 0;
	int maxcoords[] = { 0 , 0 , 0 };
	for( int x=0 ; x<euler.resolution() ; x++ ) for( int y=0 ; y<euler.resolution() ; y++ ) for( int z=0 ; z<euler.resolution() ; z++ ) if( euler(x,y,z)>max )
	{
		max = euler(x,y,z);
		maxcoords[0] = x;
		maxcoords[1] = y;
		maxcoords[2] = z; 
	}

	euler.setCoordinates( maxcoords[0] , maxcoords[1] , maxcoords[2] , matrix );
	SquareMatrix< Real , 3 > m;
	for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) m(i,j) = matrix[i][j];
	error = ( sourceKey.squareNorm() + targetKey.squareNorm() - 2.*max ) / ( sourceKey.squareNorm() + targetKey.squareNorm() );
	if( invert ) return -m;
	else         return  m;
}

///////////////////////////////////////////
// SphericalGeometry::SphericalInversion //
///////////////////////////////////////////
template< class Real >
Point3D< Real > SphericalGeometry::SphericalInversion< Real >::operator() ( Point3D< Real > p ) const
{
#ifdef ALIGN_INVERSION_TO_CENTER
	p += center;
#else // !ALIGN_INVERSION_TO_CENTER
	p -= center;
#endif // ALIGN_INVERSION_TO_CENTER
	p /= Point3D< Real >::SquareNorm( p );
	p *= 1 - Point3D< Real >::SquareNorm( center );
#ifdef ALIGN_INVERSION_TO_CENTER
	p += center;
#else // !ALIGN_INVERSION_TO_CENTER
	p -= center;
#endif // ALIGN_INVERSION_TO_CENTER
	return p;
}

///////////////////////////////////////////////////////
// SphericalGeometry::FractionalLinearTransformation //
///////////////////////////////////////////////////////
template< class Real >
SphericalGeometry::FractionalLinearTransformation< Real >::FractionalLinearTransformation( SphericalInversion< Real > si )
{
	matrix = _ITransformation( StereographicProjection( si(_ZERO) ) , StereographicProjection( si(_INFINITY) ) , StereographicProjection( si(_ONE) ) );
}
template< class Real >
SphericalGeometry::FractionalLinearTransformation< Real >::FractionalLinearTransformation( SquareMatrix< Real , 3 > m )
{
	matrix = _ITransformation( StereographicProjection( m * _ZERO ) , StereographicProjection( m * _INFINITY ) , StereographicProjection( m * _ONE ) );
}

template< class Real >
Point2D< Real > SphericalGeometry::FractionalLinearTransformation< Real >::operator()( Point2D< Real > p ) const
{
	std::complex< Real > c( p[0] , p[1] );
	c = ( c * matrix(0,0) + matrix(1,0) ) / ( c * matrix(0,1) + matrix(1,1) );
	return Point2D< Real >( c.real() , c.imag() );
}
template< class Real >
Point3D< Real > SphericalGeometry::FractionalLinearTransformation< Real >::operator()( Point3D< Real > p ) const
{
	return IStereographicProjection( (*this)( StereographicProjection( p ) ) );
}
template< class Real >
SphericalGeometry::FractionalLinearTransformation< Real > SphericalGeometry::FractionalLinearTransformation< Real >::operator * ( SphericalGeometry::FractionalLinearTransformation< Real > flt ) const
{
	FractionalLinearTransformation _flt = *this;
	return _flt *= flt;
}
template< class Real >
SphericalGeometry::FractionalLinearTransformation< Real >& SphericalGeometry::FractionalLinearTransformation< Real >::operator *= ( SphericalGeometry::FractionalLinearTransformation< Real > flt )
{
	matrix *= flt.matrix;
	return *this;
};
template< class Real >
SphericalGeometry::FractionalLinearTransformation< Real > SphericalGeometry::FractionalLinearTransformation< Real >::inverse( void ) const
{
	FractionalLinearTransformation flt;
	std::complex< Real > det = matrix(0,0) * matrix(1,1) - matrix(1,0) * matrix(0,1);
	flt.matrix(0,0) =  matrix(1,1) / det , flt.matrix(1,1) =  matrix(0,0) / det;
	flt.matrix(1,0) = -matrix(1,0) / det , flt.matrix(0,1) = -matrix(0,1) / det;
	return flt;
}

template< class Real >
SquareMatrix< std::complex< Real > , 2 > SphericalGeometry::FractionalLinearTransformation< Real >::_Transformation( Point2D< Real > p1 , Point2D< Real > p2 , Point2D< Real > p3 )
{
	SquareMatrix< std::complex< Real > , 2 > matrix;
	std::complex< Real > _p1( p1[0] , p1[1] ) , _p2( p2[0] , p2[1] ) , _p3( p3[0] , p3[1] );
	if( std::isfinite( p2[0] ) && std::isfinite( p2[1] ) )
	{
		matrix(0,0) =         ( _p3 - _p2 );
		matrix(1,0) = - _p1 * ( _p3 - _p2 );
		matrix(0,1) =         ( _p3 - _p1 );
		matrix(1,1) = - _p2 * ( _p3 - _p1 );
	}
	else
	{
		matrix(0,0) =     1;
		matrix(1,0) = - _p1;
		matrix(0,1) =     0;
		matrix(1,1) = _p3 - _p1;
	}

	return matrix;
}
template< class Real >
SquareMatrix< std::complex< Real > , 2 > SphericalGeometry::FractionalLinearTransformation< Real >::_ITransformation( Point2D< Real > p1 , Point2D< Real > p2 , Point2D< Real > p3 )
{
	SquareMatrix< std::complex< Real > , 2 > matrix , _matrix = _Transformation( p1 , p2 , p3 );
	std::complex< Real > det = _matrix(0,0) * _matrix(1,1) - _matrix(1,0) * _matrix(0,1);
	matrix(0,0) = _matrix(1,1) , matrix(1,1) = _matrix(0,0);
	matrix(0,1) = -_matrix(0,1) , matrix(1,0) = -_matrix(1,0);
	return matrix / det;
}
template< class Real > const Point3D< Real > SphericalGeometry::FractionalLinearTransformation< Real >::_ZERO    (  0 ,  0 , -1 );
template< class Real > const Point3D< Real > SphericalGeometry::FractionalLinearTransformation< Real >::_INFINITY(  0 ,  0 ,  1 );
template< class Real > const Point3D< Real > SphericalGeometry::FractionalLinearTransformation< Real >::_ONE     (  1 ,  0 ,  0 );

/////////////////////////////
// SphericalGeometry::Mesh //
/////////////////////////////
template< class Real >
SphericalGeometry::Mesh< Real >::Mesh( const std::vector< Point3D< Real > >& vertices , const std::vector< Point3D< Real > >& sVertices , const std::vector< TriangleIndex >& triangles )
{
	this->vertices = sVertices;
	this->triangles = triangles;
	masses.resize( triangles.size() );
	for( int i=0 ; i<triangles.size() ; i++ )
	{
		Point3D< Real > v[] = { vertices[ triangles[i][0] ] , vertices[ triangles[i][1] ] , vertices[ triangles[i][2] ] };
		masses[i] = (Real)Point3D< Real >::Length( Point3D< Real >::CrossProduct( v[1]-v[0] , v[2]-v[0] ) ) / (Real)2.;
	}
}
template< class Real >
void SphericalGeometry::Mesh< Real >::write( const char* fileName , const std::vector< Point3D< Real > >& vertices , bool binary ) const
{
	std::vector< PlyParametrizedVertex< float > > _vertices( vertices.size() );
	for( int i=0 ; i<vertices.size() ; i++ ) _vertices[i].point = Point3D< float >( vertices[i] ) , _vertices[i].param = Point3D< float >( this->vertices[i] );
	PlyWriteTriangles( fileName , _vertices , triangles , PlyParametrizedVertex< float >::WriteProperties , PlyParametrizedVertex< float >::WriteComponents , binary ? PLY_BINARY_NATIVE : PLY_ASCII , NULL , 0 );
}
template< class Real >
void SphericalGeometry::Mesh< Real >::write( const char* fileName , const std::vector< Point3D< Real > >& vertices , const std::vector< Point3D< Real > >& colors , bool binary ) const
{
	std::vector< PlyParametrizedColorVertex< float > > _vertices( vertices.size() );
	for( int i=0 ; i<vertices.size() ; i++ ) _vertices[i].point = Point3D< float >( vertices[i] ) , _vertices[i].param = Point3D< float >( this->vertices[i] ) , _vertices[i].color = Point3D< float >( colors[i] );
	PlyWriteTriangles( fileName , _vertices , triangles , PlyParametrizedColorVertex< float >::WriteProperties , PlyParametrizedColorVertex< float >::WriteComponents , binary ? PLY_BINARY_NATIVE : PLY_ASCII , NULL , 0 );
}
template< class Real >
void SphericalGeometry::Mesh< Real >::read( const char* fileName , std::vector< Point3D< Real > >& vertices , bool verbose , bool normalize )
{
	int plyType;
	std::vector< PlyParametrizedVertex< float > > _vertices;
	PlyReadTriangles( fileName , _vertices , triangles , PlyParametrizedVertex< float >::WriteProperties , NULL , PlyParametrizedVertex< float >::WriteComponents , plyType , NULL , 0 );
	vertices.resize( _vertices.size() );
	for( int i=0 ; i<vertices.size() ; i++ ) this->vertices[i] = Point3D< Real >( _vertices[i].param ) , vertices[i] = Point3D< Real >( _vertices[i].point );
	masses.resize( triangles.size() );
	for( int i=0 ; i<triangles.size() ; i++ )
	{
		Point3D< Real > v[] = { vertices[ triangles[i][0] ] , vertices[ triangles[i][1] ] , vertices[ triangles[i][2] ] };
		masses[i] = (Real)Point3D< Real >::Length( Point3D< Real >::CrossProduct( v[1]-v[0] , v[2]-v[0] ) ) / (Real)2.;
	}
	if( normalize ) _normalize( verbose );
}
template< class Real >
void SphericalGeometry::Mesh< Real >::read( const char* fileName , std::vector< Point3D< Real > >& vertices , std::vector< Point3D< Real > >& colors , bool verbose , bool normalize )
{
	int plyType;
	std::vector< PlyParametrizedColorVertex< float > > _vertices;
	PlyReadTriangles( fileName , _vertices , triangles , PlyParametrizedColorVertex< float >::WriteProperties , NULL , PlyParametrizedColorVertex< float >::WriteComponents , plyType , NULL , 0 );
	this->vertices.resize( _vertices.size() ) , vertices.resize( _vertices.size() ) , colors.resize( _vertices.size() );
	for( int i=0 ; i<_vertices.size() ; i++ ) this->vertices[i] = Point3D< Real >( _vertices[i].param ) , vertices[i] = Point3D< Real >( _vertices[i].point ) , colors[i] = Point3D< Real >( _vertices[i].color );
	masses.resize( triangles.size() );
	for( int i=0 ; i<triangles.size() ; i++ )
	{
		Point3D< Real > v[] = { vertices[ triangles[i][0] ] , vertices[ triangles[i][1] ] , vertices[ triangles[i][2] ] };
		masses[i] = (Real)Point3D< Real >::Length( Point3D< Real >::CrossProduct( v[1]-v[0] , v[2]-v[0] ) ) / (Real)2.;
	}
	if( normalize ) _normalize( verbose );
}

template< class Real >
double SphericalGeometry::Mesh< Real >::area( int t ) const
{
	Point3D< Real > v[] = { vertices[ triangles[t][0] ] , vertices[ triangles[t][1] ] , vertices[ triangles[t][2] ] };
	return Length( Point3D< Real >::CrossProduct( v[1]-v[0] , v[2]-v[0] ) ) / 2.;
}
template< class Real >
double SphericalGeometry::Mesh< Real >::area( int t , FractionalLinearTransformation< Real > flt ) const
{
	Point3D< Real > v[3];
	for( int j=0 ; j<3 ; j++ ) v[j] = flt( vertices[ triangles[t][j] ] );
	return Length( Point3D< Real >::CrossProduct( v[1]-v[0] , v[2]-v[0] ) ) / 2.;
}
template< class Real >
Point3D< Real > SphericalGeometry::Mesh< Real >::center( int t ) const
{
	Point3D< Real > v[] = { vertices[ triangles[t][0] ] , vertices[ triangles[t][1] ] , vertices[ triangles[t][2] ] };
	Point3D< Real > c = ( v[0] + v[1] + v[2] ) / (Real)3.;
	return c / Length( c );
}
template< class Real >
Point3D< Real > SphericalGeometry::Mesh< Real >::center( int t , FractionalLinearTransformation< Real > flt ) const
{
	Point3D< Real > v[3];
	for( int j=0 ; j<3 ; j++ ) v[j] = flt( vertices[ triangles[t][j] ] );
	Point3D< Real > c = ( v[0] + v[1] + v[2] ) / (Real)3.;
	return c / Length( c );
}
template< class Real >
Real SphericalGeometry::Mesh< Real >::constant( void ) const
{
	Real c = 0;
	for( int i=0 ; i<masses.size() ; i++ ) c += masses[i];
	return c;
}
template< class Real >
Point3D< Real > SphericalGeometry::Mesh< Real >::center( void ) const
{
	// We are looking for the vector v that minimizes
	// E(v) = \sum_T ( < v , p_T > - f(p_T) )^2 * |T|
	//      = v^t ( \sum_T (p_T) * (p_T)^T * |T| ) v - v^t ( 2 * \sum_T (p_T) * f(p_T) * |T| ) + ...
	// => ( \sum_T (p_T) * (p_T)^T * |T| ) v = \sum_T (p_T) * f(p_T) * |T|
	// Noting that the integral of p * p_T over the sphere is 4*Pi/3 * Id., this gives:
	//    v = 3 / (4*PI) * \sum_T (p_T) * f(p_T) * |T|
	Point3D< Real > c;
	for( int t=0 ; t<triangles.size() ; t++ ) c += center(t) * masses[t]; // = p * ( mA / sA ) * sA;
	return c;
}
template< class Real >
Point3D< Real > SphericalGeometry::Mesh< Real >::center( SphericalGeometry::FractionalLinearTransformation< Real > flt ) const
{
	Point3D< Real > c;
	for( int t=0 ; t<triangles.size() ; t++ ) c += center( t , flt ) * masses[t];
	return c;
}

template< class Real >
template< unsigned int SHDegree >
Point< Real , SphericalHarmonics::Dimension< SHDegree >() > SphericalGeometry::Mesh< Real >::centerSH( void ) const
{
	static_assert( SHDegree<=SphericalHarmonics::MaxDegree , "[ERROR] Degree exceeds maximum spherical harmonic degree" );
	std::vector< Point< double , SphericalHarmonics::Dimension< SHDegree >() > > projections( omp_get_max_threads() );
#pragma omp parallel for
	for( int t=0 ; t<masses.size() ; t++ )
	{
		int thread = omp_get_thread_num();
		double _p[ SphericalHarmonics::Dimension< SHDegree >() ];
		SphericalHarmonics::HarmonicValues< SHDegree >( Point3D< double >( center(t) ) , _p );
		for( unsigned int i=0 ; i<SphericalHarmonics::Dimension< SHDegree >() ; i++ ) projections[thread][i] += _p[i] * masses[t];
	}
	Point< Real , SphericalHarmonics::Dimension< SHDegree >() > p;
	for( int i=0 ; i<projections.size() ; i++ ) for( unsigned int j=0 ; j<SphericalHarmonics::Dimension< SHDegree >() ; j++ ) p[j] += (Real)projections[i][j];
	return p;
}
template< class Real >
template< unsigned int SHDegree >
SquareMatrix< Real , SphericalHarmonics::Dimension< SHDegree >() > SphericalGeometry::Mesh< Real >::dCenterSH( void ) const
{
	static_assert( SHDegree<=SphericalHarmonics::MaxDegree , "[ERROR] Degree exceeds maximum spherical harmonic degree" );
	std::vector< SquareMatrix< double , SphericalHarmonics::Dimension< SHDegree >() > > dProjections( omp_get_max_threads() );
#pragma omp parallel for
	for( int t=0 ; t<masses.size() ; t++ )
	{
		int thread = omp_get_thread_num();
		Point3D< double > _p[ SphericalHarmonics::Dimension< SHDegree >() ];
		SphericalHarmonics::HarmonicGradients< SHDegree >( Point3D< double >( center(t) ) , _p );
		for( unsigned int i=0 ; i<SphericalHarmonics::Dimension< SHDegree >() ; i++ ) for( unsigned int j=0 ; j<=i ; j++ ) dProjections[thread](i,j) += Point3D< Real >::Dot( _p[i] , _p[j] ) * masses[t];
	}
	SquareMatrix< Real , SphericalHarmonics::Dimension< SHDegree >() > d;
	for( int i=0 ; i<dProjections.size() ; i++ ) for( unsigned int j=0 ; j<SphericalHarmonics::Dimension< SHDegree >() ; j++ ) for( unsigned int k=0 ; k<=j ; k++ ) d(j,k) += (Real)dProjections[i](j,k);
	for( unsigned int i=0 ; i<SphericalHarmonics::Dimension< SHDegree >() ; i++ ) for( unsigned int j=0 ; j<=i ; j++ ) d(j,i) = d(i,j);
	return d;
}
template< class Real >
template< typename VF >
void SphericalGeometry::Mesh< Real >::advect( VF vf , int steps )
{
#pragma omp parallel for
	for( int i=0 ; i<vertices.size() ; i++ )
	{
		Point3D< Real >& v = vertices[i];
		for( int j=0 ; j<steps ; j++ )
		{
			Point3D< Real > _vf = Point3D< Real >( vf( v ) );
			Real l = Length( _vf );
			if( l ) v = v * cos( l ) + _vf / l * sin( l );
		}
	}
}

template< class Real >
SquareMatrix< Real , 3 > SphericalGeometry::Mesh< Real >::dCenter( void ) const
{
	SquareMatrix< Real , 3 > D;
	for( int t=0 ; t<triangles.size() ; t++ )
	{
		SquareMatrix< Real , 3 > _D = SquareMatrix< Real , 3 >::Identity();
		Point3D< Real > p = center( t );
		for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) _D(i,j) -= p[i] * p[j];
		D += _D * masses[t];
	}
#ifdef ALIGN_INVERSION_TO_CENTER
	return D * (Real)2;
#else // !ALIGN_INVERSION_TO_CENTER
	return -D * (Real)2;
#endif // ALIGN_INVERSION_TO_CENTER
};
template< class Real >
SquareMatrix< Real , 3 > SphericalGeometry::Mesh< Real >::dCenter( SphericalGeometry::FractionalLinearTransformation< Real > flt ) const
{
	SquareMatrix< Real , 3 > D;
	for( int t=0 ; t<triangles.size() ; t++ )
	{
		SquareMatrix< Real , 3 > _D = SquareMatrix< Real , 3 >::Identity();
		Point3D< Real > p = center( t , flt );
		for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) _D(i,j) -= p[i] * p[j];
		D += _D * masses[t];
	}
#ifdef ALIGN_INVERSION_TO_CENTER
	return D * (Real)2;
#else // !ALIGN_INVERSION_TO_CENTER
	return -D * (Real)2;
#endif // ALIGN_INVERSION_TO_CENTER
};

template< class Real >
SphericalGeometry::Mesh< Real >& SphericalGeometry::Mesh< Real >::operator *= ( SphericalGeometry::FractionalLinearTransformation< Real > flt )
{
#pragma omp parallel for
	for( int i=0 ; i<vertices.size() ; i++ ) vertices[i] = flt( vertices[i] );
	return *this;
}

template< class Real >
void SphericalGeometry::Mesh< Real >::makeUnitMass( void )
{
	Real m = constant();
	for( int i=0 ; i<masses.size() ; i++ ) masses[i] /= m;
}

template< class Real >
SphericalGeometry::FractionalLinearTransformation< Real > SphericalGeometry::Mesh< Real >::normalizer( int iters , double cutOff , bool gaussNewton , bool verbose ) const
{
	FractionalLinearTransformation< Real > flt;
	for( int i=0 ; i<iters ; i++ )
	{
		Point3D< Real > c = center( flt );
		if( verbose ) printf( "Moebius register[%d]: %8.2e\t(%8.1e %8.1e %8.1e)\n" , i , Length( c ) , c[0] , c[1] , c[2] );
		if( Point3D< Real >::SquareNorm( c )<=cutOff*cutOff ) return flt;
		SquareMatrix< Real , 3 > D = dCenter( flt );
		Point3D< Real > center;
		if( gaussNewton )
		{
			SquareMatrix< Real , 3 > D2 = ( D.transpose() * D ).inverse();
			center = - D2 * ( D * c );
		}
		else center = - D * c * 2;
		// [WARNING] Accumulating double inversions through fractional linear transformations can result in a loss of precision
		flt = FractionalLinearTransformation< Real >( SphericalInversion< Real >( center ) ) * flt;
	}
	if( verbose )
	{
		Point3D< Real > c = center( flt );
		printf( "\tMoebius register[%d]: %8.2e\t(%8.1e %8.1e %8.1e)\n" , iters , Length( c ) , c[0] , c[1] , c[2] );
	}
	return flt;
}
template< class Real >
int SphericalGeometry::Mesh< Real >::normalize( int iters , double cutOff , bool gaussNewton , bool verbose )
{
	for( int i=0 ; i<iters ; i++ )
	{
		Point3D< Real > c = center( );
		if( verbose ) printf( "Moebius register[%d]: %8.2e\t(%8.1e %8.1e %8.1e)\n" , i , Length( c ) , c[0] , c[1] , c[2] );
		if( Point3D< Real >::SquareNorm( c )<=cutOff*cutOff ) return i;
		SquareMatrix< Real , 3 > D = dCenter( );
		if( gaussNewton )
		{
			SquareMatrix< Real , 3 > D2 = ( D.transpose() * D ).inverse();
			c = - D2 * ( D * c );
		}
		else c = - D * c * 2;
		SphericalInversion< Real > inv( c );
#pragma omp parallel for
		for( int i=0 ; i<vertices.size() ; i++ ) vertices[i] = inv( vertices[i] );
	}
	if( verbose )
	{
		Point3D< Real > c = center( );
		printf( "Moebius register[%d]: %8.2e\t(%8.1e %8.1e %8.1e)\n" , iters , Length( c ) , c[0] , c[1] , c[2] );
	}
	return iters;
}

template< class Real >
template< unsigned int SHDegree >
int SphericalGeometry::Mesh< Real >::normalizeSH( int iters , int advectionSteps , Real advectionStepSize , double cutOff , bool gaussNewton , bool verbose )
{
	static const unsigned int Dim = SphericalHarmonics::Dimension< SHDegree >() - 1;
	auto SubPoint = [&]( const Point< Real , SphericalHarmonics::Dimension< SHDegree >() >& p )
	{
		Point< Real , Dim > _center , _p;
		for( int i=0 ; i<Dim ; i++ ) _p[i] = p[i+1];
		return _p;
	};
	auto SubMatrix = [&]( const SquareMatrix< Real , SphericalHarmonics::Dimension< SHDegree >() >& M )
	{
		SquareMatrix< Real , Dim > _M;
		for( int i=0 ; i<Dim ; i++ ) for( int j=0 ; j<Dim ; j++ ) _M(i,j) = M(i+1,j+1);
		return _M;
	};

	for( int i=0 ; i<iters ; i++ )
	{
		Point< Real , SphericalHarmonics::Dimension< SHDegree >() > center;
		Point< Real , Dim > c = SubPoint ( centerSH< SHDegree >() );

		if( verbose )
		{
			printf( "Moebius register[%d]: %8.2e\t" , i , sqrt( Point< Real , Dim >::SquareNorm( c ) ) );
			printf( "( " );
			for( int i=0 ; i<Dim ; i++ ) printf( " %8.1e" , c[i] );
			printf( " )\n" );
		}
		if( Point< Real , Dim >::SquareNorm( c )<=cutOff*cutOff ) return i;

		SquareMatrix< Real , Dim > D = SubMatrix( dCenterSH< SHDegree >() );
		if( gaussNewton )
		{
			SquareMatrix< Real , Dim > D2 = ( D.transpose() * D ).inverse();
			c = - D2 * ( D * c );
		}
		else c = - D * c * 2;

		for( int i=0 ; i<Dim ; i++ ) center[i+1] = c[i];

		advect( [&]( Point3D< Real > p ){ return SphericalHarmonics::Gradient< SHDegree >( center*advectionStepSize , p ); } , advectionSteps );
	}
	if( verbose )
	{
		Point< Real , Dim > c = SubPoint ( centerSH< SHDegree >() );
		printf( "Moebius register[%d]: %8.2e\t" , iters , sqrt( Point< Real , Dim >::SquareNorm( c ) ) );
		printf( "( " );
		for( int i=0 ; i<Dim ; i++ ) printf( " %8.1e" , c[i] );
		printf( " )\n" );
	}
	return iters;
}

template< class Real >
void SphericalGeometry::Mesh< Real >::_normalize( bool verbose )
{
	Point3D< Real > center;
	Real area = 0;
	for( int i=0 ; i<triangles.size() ; i++ )
	{
		Point3D< Real > v[] = { vertices[ triangles[i][0] ] , vertices[ triangles[i][1] ] , vertices[ triangles[i][2] ] };
		Real a = (Real)Length( Point3D< Real >::CrossProduct( v[1]-v[0] , v[2]-v[0] ) ) / 2;
		area += a;
		center += ( v[0] + v[1] + v[2] ) / 3 * a;
	}
	center /= area;
#pragma omp parallel for
	for( int i=0 ; i<vertices.size() ; i++ ) vertices[i] -= center;
#pragma omp parallel for
	for( int i=0 ; i<vertices.size() ; i++ ) vertices[i] /= (Real)Length( vertices[i] );
	if( verbose ) std::cout << "Center: " << center << std::endl;
}

////////////////////////////////
// SphericalGeometry::Polygon //
////////////////////////////////
template< class Real >
template< class VertexData , typename F >
void SphericalGeometry::Polygon< Real >::_split( Point3D< Real > pNormal , Real pOffset , SphericalGeometry::Polygon< Real >& back , SphericalGeometry::Polygon< Real >& front , std::vector< Point3D< Real > >& vertices , std::vector< VertexData >* vData , std::unordered_map< unsigned long long , int >& vMap , F VertexOnPlane ) const
{
	std::vector< Real > values( _vertices.size() );
	bool frontSet = false , backSet = false;

	front.sourceID = back.sourceID = sourceID;

	// Evaluate the plane's function at the vertices and mark if front/back vertices have been found
	for( int i=0 ; i<_vertices.size() ; i++ )
	{
		values[i] = Point3D< Real >::Dot( vertices[ _vertices[i] ] , pNormal ) - pOffset;

		// Don't split on edges near the pole
		if( VertexOnPlane( i ) || VertexOnPlane( i ) ) values[i] = 0.; 

		backSet |= ( values[i]<0 ) , frontSet |= ( values[i]>0 );
	}
	if( !frontSet )
	{
		for( int i=0 ; i<_vertices.size() ; i++ ) back._vertices.push_back( _vertices[i] );
		back.mass = mass;
if( back.mass!=back.mass || front.mass!=front.mass )
{
	printf( "1: %g %g %g\n" , back.mass , front.mass , mass );
	exit( 0 );
}
		return;
	}
	if( !backSet )
	{
		for( int i=0 ; i<_vertices.size() ; i++ ) front._vertices.push_back( _vertices[i] );
		front.mass = mass;
if( back.mass!=back.mass || front.mass!=front.mass )
{
	printf( "2: %g %g %g\n" , back.mass , front.mass , mass );
	exit( 0 );
}
		return;
	}
	frontSet = backSet = false;
	for( int i=0 ; i<size() ; i++ )
	{
		int i1 = ( i-1+(int)size() ) % (int)size();	
		if( values[i]*values[i1]<0 )
		{
			int vIdx;
			unsigned long long key = Key( _vertices[i] , _vertices[i1] );
			if( vMap.find( key )==vMap.end() )
			{
				Real t = values[i] / ( values[i] - values[i1] );
				Point3D< Real > v = vertices[ _vertices[i1] ]*t + vertices[ _vertices[i] ]*(Real)(1.-t);
				v /= (Real)Length( v );
				vIdx = (int)vertices.size();
				vertices.push_back( v );
				if( vData ) vData->push_back( (*vData)[ _vertices[i1] ]*t + (*vData)[ _vertices[i] ]*(Real)(1.-t) );
				vMap[key] = vIdx;
			}
			else vIdx = vMap[key];
			back._vertices.push_back( vIdx );
			front._vertices.push_back( vIdx );
			frontSet = backSet = true;
		}
		if( values[i]==0 )
		{
			back._vertices.push_back( _vertices[i] );
			front._vertices.push_back( _vertices[i] );
		}
		else
			if( values[i]<0 ) back._vertices.push_back( _vertices[i] ) , backSet = true;
			else			  front._vertices.push_back( _vertices[i] ) , frontSet = true;
	}
	if( !backSet ) back._vertices.resize( 0 );
	if( !frontSet ) front._vertices.resize( 0 );
	Real a[] = { back.area( vertices ) , front.area( vertices ) };
if( a[0]!=a[0] || a[1]!=a[1] )
{
	printf( "%g %g\n" , a[0] , a[1] );
	printf( "%d %d\n" , (int)back.size() , (int)front.size() );
	for( int i=0 ; i< back.size() ; i++ ) printf( "B[%d] %g %g %g\n" , i , vertices[ back[i]][0] , vertices[ back[i]][1] , vertices[ back[i]][2] );
	for( int i=0 ; i<front.size() ; i++ ) printf( "F[%d] %g %g %g\n" , i , vertices[ front[i]][0] , vertices[ front[i]][1] , vertices[ front[i]][2] );
	exit( 0 );
}
	Real _a = a[0] + a[1];
	if( _a ) a[0] /= _a , a[1] /= _a;
	else a[0] = a[1] = (Real)1./2;
	back.mass = mass * a[0] , front.mass = mass * a[1];
if( back.mass!=back.mass || front.mass!=front.mass )
{
	printf( "%g %g <- %g %g %g\n" , back.mass , front.mass , mass , a[0] , a[1] );
	exit( 0 );
}
}
template< class Real >
double SphericalGeometry::Polygon< Real >::area( const std::vector< Point3D< Real > >& vertices ) const
{
	double area = 0;
	for( int i=2 ; i<size() ; i++ ) area += TriangleArea( vertices[ _vertices[0] ] , vertices[ _vertices[i-1] ] , vertices[ _vertices[i] ] );
	return area;
}

/////////////////////////////////////
// SphericalGeometry::Tessellation //
/////////////////////////////////////
template< class Real >
template< class VertexData >
void SphericalGeometry::Tessellation< Real >::_splitPolygonPhi( Polygon< Real > p , int theta , std::unordered_map< unsigned long long , int >& vMap , std::vector< VertexData >* vData )
{
	auto Phi = [&] ( Point3D< Real > p )
	{
		if     ( p[1]<=-1 ) return M_PI;
		else if( p[1]>= 1 ) return 0.;
		else                return acos( p[1] );
	};
	Real minPhi = Phi( _vertices[ p[0] ] ) , maxPhi = Phi( _vertices[ p[0] ] );
	for( int j=1 ; j<p.size() ; j++ ) minPhi = std::min< Real >( minPhi , Phi( _vertices[ p[j] ] ) ) , maxPhi = std::max< Real >( maxPhi , Phi( _vertices[ p[j] ] ) );
	minPhi *= (Real)_resolution / M_PI , maxPhi *= (Real)_resolution / M_PI;
	int start = (int)floor( minPhi ) , end = (int)ceil( maxPhi );

	for( int i=start+1 ; i<=end-1 ; i++ )
	{
		Real phi = ( (Real)i ) * M_PI / _resolution;
		Point3D< Real > normal( 0 , -1 , 0 );
		Real offset = -cos( phi );
		Polygon< Real > back , front;

		if( vData ) p.split( normal , offset , back , front , _vertices , *vData , vMap , [&]( int ){ return false; } );
		else        p.split( normal , offset , back , front , _vertices ,          vMap , [&]( int ){ return false; } );
		back.theta = theta<0 ? theta+_resolution : theta;
		back.phi = i-1;
		_polygons.push_back( back );
		p = front;
	}
	p.theta = theta<0 ? theta+_resolution : theta;
	p.phi = end-1;
	_polygons.push_back( p );
}

template< class Real >
template< class VertexData >
void SphericalGeometry::Tessellation< Real >::_splitPolygon( SphericalGeometry::Polygon< Real > p , std::unordered_map< unsigned long long , int >& vMap , const int poles[2] , std::vector< VertexData >* vData )
{
	Polygon< Real > _p[2];
	if( vData ) p.split( Point3D< Real >( 0 , 0 , 1 ) , 0. , _p[0] , _p[1] , _vertices , *vData , vMap , [&]( int vIdx ) { return vIdx==poles[0] || vIdx==poles[1]; } );
	else        p.split( Point3D< Real >( 0 , 0 , 1 ) , 0. , _p[0] , _p[1] , _vertices ,          vMap , [&]( int vIdx ) { return vIdx==poles[0] || vIdx==poles[1]; } );

	for( int w=0 ; w<2 ; w++ ) if( _p[w].size() )
	{
		// Recall that the spherical parameterization is given by:
		//		Phi( theta , phi ) = ( cos(theta) * sin(phi) , cos(phi) , sin(theta) * sin(phi) )
		// with theta in [0,2*pi] and phi in [0,pi]
		auto Theta = [&] ( int vIdx )
		{
			Point3D< Real > v = _vertices[ vIdx ];
			if( ( w==0 && v[2]>0 ) || ( w==1 && v[2]<0 ) ) v[2] = 0;
			Real theta = atan2( v[2] , v[0] );
			if( w==0 ) // [-Pi,0]
			{
				if( theta>0 ) theta -= 2 * M_PI;
			}
			else // [0,Pi]
			{
				if( theta<0 ) theta += 2 * M_PI;
			}
			return theta;
		};
		Polygon< Real > p = _p[w];

		Real minTheta = Theta( p[0] ) , maxTheta = Theta( p[0] );
		int vSkip = -1;  
		if( p[0]==poles[0] || p[0]==poles[1] ) minTheta = Theta( p[1] ) , maxTheta = Theta( p[1] ) , vSkip = 0;
		else for( int j=1 ; j<p.size() ; j++ ) if( p[j]==poles[0] || p[j]==poles[1] ) vSkip = j;

		for( int j=1 ; j<p.size() ; j++ ) 
		{
			if( j==vSkip ) continue; 
			minTheta = std::min< Real >( minTheta , Theta( p[j] ) );
			maxTheta = std::max< Real >( maxTheta , Theta( p[j] ) );	
		}

		minTheta *= (Real)_resolution / ( 2.*M_PI ) , maxTheta *= (Real)_resolution / ( 2.*M_PI );
		int start = (int)floor( minTheta ) , end = (int)ceil( maxTheta );
		for( int i=start+1 ; i<=end-1 ; i++ )
		{
			Real theta = ( (Real)i ) * 2. * M_PI / _resolution;
			Polygon< Real > back , front;
			if( vData ) p.split( Point3D< Real >( -sin(theta) , 0 , cos(theta) ) , 0. , back , front, _vertices , *vData , vMap , [&]( int vIdx ){ return vIdx==poles[0] || vIdx==poles[1]; } );
			else        p.split( Point3D< Real >( -sin(theta) , 0 , cos(theta) ) , 0. , back , front, _vertices ,          vMap , [&]( int vIdx ){ return vIdx==poles[0] || vIdx==poles[1]; } );
			_splitPolygonPhi( back , i-1 , vMap , vData );
			p = front;
		}
		_splitPolygonPhi( p , end-1 , vMap , vData );
	}
}

template< class Real >
template< class VertexData >
void SphericalGeometry::Tessellation< Real >::_init( const SphericalGeometry::Mesh< Real >& mesh , std::vector< VertexData >* vData , int resolution , bool verbose )
{
	_resolution = resolution;
	int poles[] = { -1 , -1 };
	_vertices = mesh.vertices;
	_polygons.resize( mesh.triangles.size() );
#pragma omp parallel for
	for( int i=0 ; i<mesh.triangles.size() ; i++ ) _polygons[i] = Polygon< Real >( i , mesh.masses[i] , mesh.triangles[i][0] , mesh.triangles[i][1] , mesh.triangles[i][2] );
#pragma omp parallel for
	for( int i=0 ; i<_vertices.size() ; i++ ) _vertices[i] /= (Real)Length( _vertices[i] );
	auto Normal = [&] ( int t )
	{
		Point3D< Real > v[3] = { _vertices[ _polygons[t][0] ] , _vertices[ _polygons[t][1] ] , _vertices[ _polygons[t][2] ] };
		return Point3D< Real >::CrossProduct( v[1]-v[0] , v[2]-v[0] );
	};

	// Determine if a triangle is covering a pole. 
	// If so, preemptively split and store the pole index. 
	{
		int count = 0; 
		for( int i=(int)_polygons.size()-1 ; i>=0 ; i-- )
		{
			int v0 = _polygons[i][0] , v1 = _polygons[i][1] , v2 = _polygons[i][2];

			// 1. Find the plane containing this triangle 
			Point3D< Real > n = Normal(i);
			// Ignore triangles with very, very, very small area 
			if( n[0]==0 && n[1]==0 && n[2]==0 ) continue;
			n /= (Real)sqrt( Point3D< Real >::SquareNorm( n ) ); 

			// No need to test if the normal is perpendicular to the y-axis 
			if( !Point3D< Real >::Dot( n , Point3D<Real>( 0. , 1. , 0. ) ) ) continue;
			Point3D< Real > v[] = { _vertices[ v0 ] , _vertices[ v1 ] , _vertices[ v2 ] };

			Real d = -n[0]*v[0][0] - n[1]*v[0][1] - n[2]*v[0][2]; 

			// 2. Compute intersection of P with the y-axis 
			Point3D< Real > p; 
			p[1] = - d / n[1]; 

			// 3. Determine whether p is inside or outside of the triangle 
			Real a[3];
			BarycentricCoordinates< Real >( p , v[0] , v[1] , v[2] , a[0] , a[1] , a[2] );
			if( a[0]>=0 && a[1]>=0 && a[2]>=0 ) 
			{
				// This triangle is at one of the poles. Manually split by placing a vertex at p. 
				int vIdx = (int)_vertices.size();
				_vertices.push_back( p/Length(p) );
				if( vData ) vData->push_back( BarycentricInterpolate( a , (*vData)[v0] , (*vData)[v1] , (*vData)[v2] ) );

				if( count>1 ) fprintf( stderr , "Error: Not expecting more than 2 poles!\n" ) , exit( 0 );
				poles[count++] = vIdx; 

				// Create new polygons 
				Polygon< Real > p1( _polygons[i].sourceID , 0 , vIdx , v0 , v1 );
				Polygon< Real > p2( _polygons[i].sourceID , 0 , vIdx , v1 , v2 );
				Polygon< Real > p3( _polygons[i].sourceID , 0 , vIdx , v2 , v0 );
				Real m1 , m2 , m3;
				{
					m1 = p1.area( _vertices ) , m2 = p2.area( _vertices ) , m3 = p3.area( _vertices );
					Real m = m1+m2+m3;
					if( m ) m1 /= m , m2 /= m , m3 /= m;
					else m1 = m2 = m3 = (Real)1./3;
					m1 *= _polygons[i].mass , m2 *= _polygons[i].mass , m3 *= _polygons[i].mass;
				}
				p1.mass = m1 , p2.mass = m2 , p3.mass = m3;
				_polygons[i] = p1 , _polygons.push_back( p2 ) , _polygons.push_back( p3 );
			}
		}
		if( count!=2 ) fprintf( stderr , "[ERROR] Couldn't find both poles\n" ) , exit( 0 );
	}

	std::unordered_map< unsigned long long , int > vMap;
	for( int i=(int)_polygons.size()-1 ; i>=0 ; i-- )
	{
		Polygon< Real > p = _polygons[i];
		_polygons[i] = _polygons.back() , _polygons.pop_back();
		_splitPolygon( p , vMap , poles , vData );
	}
}

template< class Real >
void SphericalGeometry::Tessellation< Real >::createSGrid( SphericalGrid< Real >& sGrid , Real smooth ) const
{
	if( sGrid.resolution()!=_resolution ) 
	{
		fprintf( stderr , "[WARNING] Input sGrid resolution (%d) did not match resolution of tessellation (%d). Resizing.\n" , sGrid.resolution() , _resolution );
		sGrid.resize( _resolution ); 
	}

	auto Theta = [&]( Real i ){ return Real( 2.0*M_PI*i/sGrid.resolution() ); };
	auto Phi   = [&]( Real i ){ return Real( M_PI*(2.0*i+1)/(2.0*sGrid.resolution()) ); };

	for( int x=0 ; x<sGrid.resolution() ; x++ ) for( int y=0 ; y<sGrid.resolution() ; y++ ) sGrid(x,y) = 0;
	for( int i=0 ; i<_polygons.size() ; i++ ) sGrid( _polygons[i].theta , _polygons[i].phi ) += _polygons[i].mass;
	for( int x=0 ; x<sGrid.resolution() ; x++ )  for( int y=0 ; y<sGrid.resolution() ; y++ ) 
	{
		// Normalize the bucket values by area 
		Real theta[] = { Theta(x-0.5) , Theta(x+0.5) };
		Real phi[] = { Phi(y-0.5) , Phi(y+0.5) };
		Real area = ( theta[1] - theta[0] ) * ( -cos(phi[1]) + cos(phi[0]) );
		sGrid(x,y) = sqrt( sGrid(x,y) / area );
	}

	if( smooth>0 )
	{
		HarmonicTransform< Real > hForm; 
		FourierKeyS2< Real > key;

		hForm.ForwardFourier( sGrid , key );
		for( int i=0 ; i<key.bandWidth() ; i++ ) for( int j=0 ; j<=i ; j++ ) key(i,j) *= exp( -smooth * i * ( (Real)(i+1) ) * 1. );
		hForm.InverseFourier( key , sGrid );
	}
}

template< class Real >
template< class VertexData >
void SphericalGeometry::Tessellation< Real >::_collapseCells( std::vector< VertexData >* vData )
{
	std::vector< std::vector< Polygon< Real > > > polygons( _resolution*_resolution );
	for( int i=(int)_polygons.size()-1 ; i>=0 ; i-- ) polygons[ _polygons[i].theta*_resolution + _polygons[i].phi ].push_back( _polygons[i] ) , _polygons.pop_back();
	auto EdgeKey = []( int v1 , int v2 )
	{
		unsigned long long _v1 = v1 , _v2 = v2;
		return ( _v1<<32 ) | _v2;
	};
	auto FactorEdgeKey = []( unsigned long long key , int& v1 , int& v2 )
	{
		v1 = (int)((key>>32)    );
		v2 = (int)((key<<32)>>32);
	};
#pragma omp parallel for
	for( int i=0 ; i<polygons.size() ; i++ )
	{
		Real mass = 0;
		std::unordered_set< unsigned long long > eSet;
		for( int j=0 ; j<polygons[i].size() ; j++ ) mass += polygons[i][j].mass;
		for( int j=0 ; j<polygons[i].size() ; j++ ) for( int k=0 ; k<polygons[i][j].size() ; k++ ) eSet.insert( EdgeKey( polygons[i][j][k] , polygons[i][j][k+1] ) );
		std::vector< std::pair< int , int > > bEdges;
		for( std::unordered_set< unsigned long long >::iterator i=eSet.begin() ; i!=eSet.end() ; i++ )
		{
			int v1 , v2;
			FactorEdgeKey( *i , v1 , v2 );
			if( eSet.find( EdgeKey( v2 , v1 ) )==eSet.end() ) bEdges.push_back( std::pair< int , int >( v1 , v2 ) );
		}
		int pCount = 0;
		while( bEdges.size() )
		{
			std::pair< int , int > e = bEdges.back() ; bEdges.pop_back();
			int first = e.first , last = e.second;
			Polygon< Real > p( -1 , mass );
			p.theta = polygons[i][0].theta;
			p.phi = polygons[i][0].phi;
			p.push_back( e.first );
			while( first!=last )
			{
				bool found = false;
				for( int b=0 ; b<bEdges.size() ; b++ ) if( bEdges[b].first==last )
				{
					p.push_back( last );
					last = bEdges[b].second;
					bEdges[b] = bEdges.back() ; bEdges.pop_back();
					found = true;
					break;
				}
				if( !found ) fprintf( stderr , "[ERROR] Couldn't find next edge\n" ) , exit( 0 );
			}
#pragma omp critical
			{
				_polygons.push_back( p );
			}
			pCount++;
		}
		if( pCount!=1 ) fprintf( stderr , "[WARNING] Found %d polygons in cell\n" , pCount );
	}
	std::vector< int > vMap( _vertices.size() , -1 );
	int idx = 0;
	for( int i=0 ; i<_polygons.size() ; i++ ) for( int j=0 ; j<_polygons[i].size() ; j++ ) if( vMap[ _polygons[i][j] ]==-1 ) vMap[ _polygons[i][j] ] = idx++;

	// Remove valence-2 vertices
	std::vector< int > vCount( _vertices.size() , 0 );
	for( int i=0 ; i<_polygons.size() ; i++ ) for( int j=0 ; j<_polygons[i].size() ; j++ ) vCount[ _polygons[i][j] ]++;
	for( int i=0 ; i<_polygons.size() ; i++ )
	{
		Polygon< Real > p( -1 , _polygons[i].mass );
		p.theta = _polygons[i].theta , p.phi = _polygons[i].phi;
		for( int j=0 ; j<_polygons[i].size() ; j++ ) if( vCount[ _polygons[i][j] ]>2 ) p.push_back( _polygons[i][j] );
		_polygons[i] = p;
	}

	// Remove unused vertices
	std::vector< Point3D< Real > > vertices( idx );
	for( int i=0 ; i<vMap.size() ; i++ ) if( vMap[i]>=0 ) vertices[ vMap[i] ] = _vertices[i];
	for( int i=0 ; i<_polygons.size() ; i++ ) for( int j=0 ; j<_polygons[i].size() ; j++ ) _polygons[i][j] = vMap[ _polygons[i][j] ];
	_vertices = vertices;
	if( vData )
	{
		std::vector< VertexData > _vData( idx );
		for( int i=0 ; i<vMap.size() ; i++ ) if( vMap[i]>=0 ) _vData[ vMap[i] ] = (*vData)[i];
		*vData = _vData;
	}
}
template< class Real >
void SphericalGeometry::Tessellation< Real >::write( const char* fileName , bool binary ) const
{
	std::vector< PlyVertex< float > > v( _vertices.size() );
	std::vector< std::vector< int > > p( _polygons.size() );
	for( int i=0 ; i<_vertices.size() ; i++ ) v[i].point = Point3D< float >( _vertices[i] );
	for( int i=0 ; i<_polygons.size() ; i++ )
	{
		p[i].resize( _polygons[i].size() );
		for( int j=0 ; j<_polygons[i].size() ; j++ ) p[i][j] = _polygons[i][j];
	}
	std::vector< float > m;
	m.resize( _polygons.size() );
	for( int i=0 ; i<_polygons.size() ; i++ ) m[i] = (float)_polygons[i].mass;
	PlyProperty massProperties[] = { { "mass" , PLY_FLOAT , PLY_FLOAT , 0 , 0 , 0 , 0 , 0 } };
	PlyWritePolygons( fileName , v , p , m , PlyVertex< float >::WriteProperties , PlyVertex< float >::WriteComponents , massProperties , 1 , PLY_BINARY_NATIVE , NULL , 0 );
}
template< class Real >
void SphericalGeometry::Tessellation< Real >::write( const char* fileName , const std::vector< Point3D< Real > >& vertices , const std::vector< Point3D< Real > >& colors , bool binary ) const
{
	std::vector< PlyParametrizedColorVertex< float > > v( _vertices.size() );
	std::vector< std::vector< int > > p( _polygons.size() );
	for( int i=0 ; i<_vertices.size() ; i++ )
	{
		v[i].param = Point3D< float >( _vertices[i] );
		v[i].point = Point3D< float >(  vertices[i] );
		for( int j=0 ; j<3 ; j++ ) v[i].color[j] = (float)colors[i][j];
	}
	for( int i=0 ; i<_polygons.size() ; i++ )
	{
		p[i].resize( _polygons[i].size() );
		for( int j=0 ; j<_polygons[i].size() ; j++ ) p[i][j] = _polygons[i][j];
	}
	PlyWritePolygons( fileName , v , p , PlyParametrizedColorVertex< float >::WriteProperties , PlyParametrizedColorVertex< float >::WriteComponents , PLY_BINARY_NATIVE , NULL , 0 );
}

template< class Real >
template< typename F >
void SphericalGeometry::Tessellation< Real >::writeGrid( const char* fileName , Real smooth , F f , bool binary ) const
{
	std::vector< PlyColorVertex< float > > v( _resolution*(_resolution-1)+2 );
	std::vector< std::vector< int > > polygons;

	auto Index = [&]( int i , int j )
	{
		i %= _resolution;
		if( j==0 ) return 0;
		else if( j==_resolution ) return _resolution*(_resolution-1)+1;
		else return 1 + (j-1)*_resolution+i;
	};
	auto Position = [&]( int i , int j )
	{
		Real theta = (2.*M_PI*i)/_resolution;
		Real phi = (PI*j)/_resolution;
		Point3D< Real > p;
		SphericalGrid< Real >::SetCoordinates( theta , phi , &p[0] );
		return p;
	};
	auto Color = []( Real value )
	{
		Point3D< Real > blue( 0 , 0 , 1 ) , red( 1 , 0 , 0 ) , gray( 0.5 , 0.5 , 0.5 );
		value = std::max< Real >( -1 , std::min< Real >( value , 1 ) );
		if( value<0 ) return - blue * value + gray * ( 1 + value );
		else          return    red * value + gray * ( 1 - value );
	};

	for( int i=0 ; i<_resolution ; i++ ) for( int j=0 ; j<=_resolution ; j++ ) v[ Index(i,j) ].point = Position( i , j );
	// North triangles
	{
		int j=0;
		for( int i=0 ; i<_resolution ; i++ )
		{
			std::vector< int > triangle(3);
			triangle[0] = Index(0,j);
			triangle[1] = Index(i+1,j+1);
			triangle[2] = Index(i,j+1);
			polygons.push_back( triangle );
		}
	}
	// South triangles
	{
		int j=_resolution;
		for( int i=0 ; i<_resolution ; i++ )
		{
			std::vector< int > triangle(3);
			triangle[0] = Index(0,j);
			triangle[1] = Index(i,j-1);
			triangle[2] = Index(i+1,j-1);
			polygons.push_back( triangle );
		}
	}
	// Interior quads
	for( int j=1 ; j<_resolution-1 ; j++ ) for( int i=0 ; i<_resolution ; i++ )
	{
		std::vector< int > quad(4);
		quad[0] = Index(i,j);
		quad[1] = Index(i+1,j);
		quad[2] = Index(i+1,j+1);
		quad[3] = Index(i,j+1);
		polygons.push_back( quad );
	}

	SphericalGrid< Real > sGrid( _resolution );
	createSGrid( sGrid , smooth );
	std::vector< Real > values( v.size() );
	for( int i=0 ; i<values.size() ; i++ )
	{
		Point3D< Real > p = v[i].point;
		Real x , y;
		sGrid.setCoordinates( &p[0] , x , y );
		values[i] = sGrid(x,y);
	}
	Real max = 0;
	for( int i=0 ; i<values.size() ; i++ ) max = std::max< Real >( max , (Real)fabs( values[i] ) );
	for( int i=0 ; i<values.size() ; i++ ) values[i] = (Real)f( fabs(values[i]/max) ) * ( values[i]>0 ? 1 : -1 );

	for( int i=0 ; i<v.size() ; i++ )
	{
		v[i].point *= (float)fabs( values[i] );
		v[i].color = Point3D< float >( Color( values[i] ) ) * 255.f;
	}

	PlyWritePolygons( fileName , v , polygons , PlyColorVertex< float >::WriteProperties , PlyColorVertex< float >::WriteComponents , PLY_ASCII , NULL , 0 );
}

////////////////////////////////////////////
// SphericalGeometry::CircumscribedSphere //
////////////////////////////////////////////
template< class Real >
void SphericalGeometry::CircumscribedSphere< Real >::_RandomCenter( const std::vector< Point3D< Real > >& vertices , Point3D< Real >& c , Real& a )
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> dis(0, (int)vertices.size() - 1);

	Point3D< Real > v[4]; 
	for( int i=0 ; i<4 ; i++ ) v[i] = vertices[dis(gen)];

	// Compute determinant matrices 
	SquareMatrix< Real , 4 > A , Dx , Dy , Dz; 
	Real dx , dy , dz;

	for( int i=0 ; i<4 ; i++ ) 
	{
		for( int j=0 ; j<3 ; j++ ) A(i, j) = v[i][j]; 
		A( i , 3 ) = 1; 
	}
	a = A.determinant(); 

	for( int i=0 ; i<4 ; i++ )
	{
		Dx( i , 0 ) = v[i].squareNorm(); 
		Dx( i , 1 ) = v[i][1];
		Dx( i , 2 ) = v[i][2];
		Dx( i , 3 ) = 1; 
	}
	dx = Dx.determinant(); 

	for( int i=0 ; i<4 ; i++ )
	{
		Dy( i , 0 ) = v[i].squareNorm(); 
		Dy( i , 1 ) = v[i][0];
		Dy( i , 2 ) = v[i][2];
		Dy( i , 3 ) = 1; 
	}
	dy = -Dy.determinant(); 

	for( int i=0 ; i<4 ; i++ )
	{
		Dz( i , 0 ) = v[i].squareNorm(); 
		Dz( i , 1 ) = v[i][0];
		Dz( i , 2 ) = v[i][1];
		Dz( i , 3 ) = 1; 
	}
	dz = Dz.determinant(); 

	c[0] = dx/(2.*a); 
	c[1] = dy/(2.*a);
	c[2] = dz/(2.*a);
}

template< class Real >
Point3D< Real > SphericalGeometry::CircumscribedSphere< Real >::Center( const std::vector< Point3D< Real > >& vertices , int iters )
{
	Real w=0 , x=0 , y=0 , z=0;
#pragma omp parallel for reduction( + : w , x , y , z )
	for( int i=0 ; i<iters ; i++ )
	{
		Point3D< Real > c;
		Real _w;
		_RandomCenter( vertices , c , _w );
		if( _w>0 ) x += c[0] *_w , y += c[1] * _w , z += c[2] * _w , w +=_w;
	}
	return Point3D< Real >( x/w , y/w , z/w );
}
