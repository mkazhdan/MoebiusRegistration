#ifndef SPHERICAL_HARMONICS_INCLUDED
#define SPHERICAL_HARMONICS_INCLUDED
#include <vector>
#include <omp.h>
#include "Misha/Geometry.h"

namespace HomogeneousPolynomials
{
	template< unsigned int Degree > constexpr unsigned int Dimension( void ){ return ( (Degree+2)*(Degree+1) ) / 2;}
	unsigned int Index( unsigned int x , unsigned int y , unsigned int z ){ return (x+y+z)*x - x * ( x-1 )/2 + x + y; }

	// Evaluates the different powers of the x-, y-, and z-coefficientss
	template< unsigned int Degree >
	static Point< double , Dimension< Degree >() > Monomials( Point3D< double > p )
	{
		Point< double , Dimension< Degree >() > m;
		double _x = 1;
		for( unsigned int x=0 ; x<=Degree ; x++ )
		{
			double _y = 1;
			for( unsigned int y=0 ; y<=Degree-x ; y++ ) 
			{
				unsigned int z = Degree - x - y;
				{
					double _z = pow( p[2] , z );
					m[ Index( x , y , z ) ] = _x * _y * _z;
				}
				_y *= p[1];
			}
			_x *= p[0];
		}
		return m;
	}
	template< unsigned int Degree >
	struct Function
	{
		static const unsigned int Dim = Dimension< Degree >();
		Point< double , Dim > coefficients;
		Function( void ){}
		Function( unsigned int x , unsigned int y , unsigned int z , double v=1. ){ coefficients[ Index( x , y , z ) ] = v; }
		Function( Point< double , Dim > c ) : coefficients(c) { }
		double value( Point< double , Dim > p ) const{ return Point< double , Dim >::Dot( coefficients , p ); }
		double operator()( Point3D< double > p ) const { return value( Monomials< Degree >( p ) ); }
		Function operator + ( const Function& f ) const { return Function( coefficients + f.coefficients ); }
		Function operator - ( const Function& f ) const { return Function( coefficients - f.coefficients ); }
		Function operator * ( double s ) const { return Function( coefficients*s ); }
		Function operator / ( double s ) const { return Function( coefficients/s ); }

		void print( bool newline=true ) const
		{
			bool first = true;
			for( unsigned int x=0 , idx=0 ; x<=Degree ; x++ ) for( unsigned int y=0 ; y<=Degree-x ; y++ , idx++ )
			{
				unsigned int z = Degree-x-y;
				if( coefficients[idx] )
				{
					if( coefficients[idx]>0 ) 
						if( first ) printf(   " %.2f" ,  coefficients[idx] ); 
						else        printf( " + %.2f" ,  coefficients[idx] );
					else            printf( " - %.2f" , -coefficients[idx] );
					if( x ) printf( " * X^%d" , x );
					if( y ) printf( " * Y^%d" , y );
					if( z ) printf( " * Z^%d" , z );
					first = false;
				}
			}
			if( newline ) printf( "\n" );
		}
	};

	template< unsigned int Degree >
	Function< Degree > operator * ( double scale , const Function< Degree >& f ){ return Function< Degree >( f.coefficients*scale ); }

	template< unsigned int Degree1 , unsigned int Degree2 >
	Function< Degree1+Degree2 > operator * ( const Function< Degree1 >& f1 , const Function< Degree2 >& f2 )
	{
		Function< Degree1+Degree2 > f;
		for( unsigned int idx1=0 , x1=0 ; x1<=Degree1 ; x1++ ) for( unsigned int y1=0 ; y1<=Degree1-x1 ; y1++ , idx1++ )
		{
			unsigned int z1 = Degree1 - x1 - y1;
			for( unsigned int idx2=0 , x2=0 ; x2<=Degree2 ; x2++ ) for( unsigned int y2=0 ; y2<=Degree2-x2 ; y2++ , idx2++ )
			{
				unsigned int z2 = Degree2 - x2 - y2;
				f.coefficients[ Index( x1+x2 , y1+y2 , z1+z2 ) ] += f1.coefficients[idx1] * f2.coefficients[idx2];
			}
		}
		return f;
	}

	template< unsigned int Degree >
	struct VectorField
	{
		static const unsigned int Dim = Dimension< Degree >();
		Point< double , Dim > dx , dy , dz;
		VectorField( void ){}
		// Sets the vector field to be the gradient of the scalar function
		VectorField( const Function< Degree+1 >& f )
		{
			int idx = 0;
			for( unsigned int x=0 ; x<=Degree+1 ; x++ ) for( unsigned int y=0 ; y<=Degree+1-x ; y++ , idx++ ) 
			{
				int z = Degree+1 - x - y;
				if( x ) dx[ Index(x-1,y,z) ] += f.coefficients[idx] * x;
				if( y ) dy[ Index(x,y-1,z) ] += f.coefficients[idx] * y;
				if( z ) dz[ Index(x,y,z-1) ] += f.coefficients[idx] * z;
			}
		}
		Point3D< double > value( Point< double , Dim > p ) const { return Point3D< double >( Point< double , Dim >::Dot( dx , p ) , Point< double , Dim >::Dot( dy , p ) , Point< double , Dim >::Dot( dz , p ) ); }
		Point3D< double > operator()( Point3D< double > p ) const { return value( Monomials< Degree >( p ) ); }
	};
};

namespace SphericalHarmonics
{
	const unsigned int MaxDegree = 4;
	template< unsigned int Degree > constexpr unsigned int Dimension( void ){ return (Degree+1)*(Degree+1); }

	HomogeneousPolynomials::Function< 0 > ONE = HomogeneousPolynomials::Function< 0 >( 0 , 0 , 0 );
	HomogeneousPolynomials::Function< 1 > X   = HomogeneousPolynomials::Function< 1 >( 1 , 0 , 0 );
	HomogeneousPolynomials::Function< 1 > Y   = HomogeneousPolynomials::Function< 1 >( 0 , 1 , 0 );
	HomogeneousPolynomials::Function< 1 > Z   = HomogeneousPolynomials::Function< 1 >( 0 , 0 , 1 );

	HomogeneousPolynomials::Function< 2 > R2 = X*X + Y*Y + Z*Z;

	HomogeneousPolynomials::Function< 0 > Harmonic0[] =
	{
		ONE * sqrt(1./(4.*M_PI) )
	};
	HomogeneousPolynomials::Function< 1 > Harmonic1[] =
	{
		Y * sqrt(3./(4.*M_PI) ) ,
		Z * sqrt(3./(4.*M_PI) ) ,
		X * sqrt(3./(4.*M_PI) ) ,
	};
	HomogeneousPolynomials::Function< 2 > Harmonic2[] =
	{
		( X*Y )					* sqrt(15./( 4.*M_PI) ) ,
		( Y*Z )					* sqrt(15./( 4.*M_PI) ) ,
		( 2*Z*Z - X*X - Y*Y )	* sqrt( 5./(16.*M_PI) ) ,
		( X*Z )					* sqrt(15./( 4.*M_PI) ) ,
		( X*X - Y*Y )			* sqrt(15./(16.*M_PI) ) ,
	};
	HomogeneousPolynomials::Function< 3 > Harmonic3[] =
	{
		( 3*X*X - Y*Y ) * Y				* sqrt( 35./(32.*M_PI) ) ,
		( X*Y*Z )						* sqrt(105./( 4.*M_PI) ) ,
		( 4*Z*Z - X*X - Y*Y ) * Y		* sqrt( 21./(32.*M_PI) ) ,
		( 2*Z*Z - 3*X*X - 3*Y*Y ) * Z	* sqrt(  7./(16.*M_PI) ) ,
		( 4*Z*Z - X*X - Y*Y ) * X		* sqrt( 21./(32.*M_PI) ) ,
		( X*X - Y*Y ) * Z				* sqrt(105./(16.*M_PI) ) ,
		( X*X - 3*Y*Y ) * X				* sqrt( 35./(32.*M_PI) ) ,
	};
	HomogeneousPolynomials::Function< 4 > Harmonic4[] =
	{
		( X*X   - Y*Y           ) * X*Y										* sqrt( 35./(   M_PI) ) * (3./ 4) , 
		( 3*X*X - Y*Y           ) * Y*Z										* sqrt( 35./(2.*M_PI) ) * (3./ 4) ,
		( 6*Z*Z - X*X   - Y*Y   ) * X*Y										* sqrt(  5./(   M_PI) ) * (3./ 4) ,
		( 4*Z*Z - 3*X*X - 3*Y*Y ) * Y*Z										* sqrt(  5./(2.*M_PI) ) * (3./ 4) ,
		( 35*Z*Z*Z*Z - 30*(X*X+Y*Y+Z*Z)*Z*Z + (X*X+Y*Y+Z*Z)*(X*X+Y*Y+Z*Z) )	* sqrt(  1./(   M_PI) ) * (3./16) ,
		( 4*Z*Z - 3*X*X - 3*Y*Y ) * X*Z										* sqrt(  5./(2.*M_PI) ) * (3./ 4) ,
		( 6*Z*Z - X*X - Y*Y ) * ( X*X - Y*Y )								* sqrt(  5./(   M_PI) ) * (3./ 8) ,
		( X*X - 3*Y*Y ) * X*Z												* sqrt( 35./(2.*M_PI) ) * (3./ 4) ,
		( ( X*X - 3*Y*Y ) * X*X - ( 3*X*X - Y*Y ) * Y*Y )					* sqrt( 35./(   M_PI) ) * (3./16) ,
	};

	HomogeneousPolynomials::VectorField< 0 > DHarmonic1[] =
	{
		HomogeneousPolynomials::VectorField< 0 >( Harmonic1[0] ) , HomogeneousPolynomials::VectorField< 0 >( Harmonic1[1] ) , HomogeneousPolynomials::VectorField< 0 >( Harmonic1[2] ) ,
	};
	HomogeneousPolynomials::VectorField< 1 > DHarmonic2[] =
	{
		HomogeneousPolynomials::VectorField< 1 >( Harmonic2[0] ) , HomogeneousPolynomials::VectorField< 1 >( Harmonic2[1] ) , HomogeneousPolynomials::VectorField< 1 >( Harmonic2[2] ) , HomogeneousPolynomials::VectorField< 1 >( Harmonic2[3] ) , HomogeneousPolynomials::VectorField< 1 >( Harmonic2[4] ) ,
	};
	HomogeneousPolynomials::VectorField< 2 > DHarmonic3[] =
	{
		HomogeneousPolynomials::VectorField< 2 >( Harmonic3[0] ) , HomogeneousPolynomials::VectorField< 2 >( Harmonic3[1] ) , HomogeneousPolynomials::VectorField< 2 >( Harmonic3[2] ) , HomogeneousPolynomials::VectorField< 2 >( Harmonic3[3] ) , HomogeneousPolynomials::VectorField< 2 >( Harmonic3[4] ) , HomogeneousPolynomials::VectorField< 2 >( Harmonic3[5] ) , HomogeneousPolynomials::VectorField< 2 >( Harmonic3[6] ) ,
	};
	HomogeneousPolynomials::VectorField< 3 > DHarmonic4[] =
	{
		HomogeneousPolynomials::VectorField< 3 >( Harmonic4[0] ) , HomogeneousPolynomials::VectorField< 3 >( Harmonic4[1] ) , HomogeneousPolynomials::VectorField< 3 >( Harmonic4[2] ) , HomogeneousPolynomials::VectorField< 3 >( Harmonic4[3] ) , HomogeneousPolynomials::VectorField< 3 >( Harmonic4[4] ) , HomogeneousPolynomials::VectorField< 3 >( Harmonic4[5] ) , HomogeneousPolynomials::VectorField< 3 >( Harmonic4[6] ) , HomogeneousPolynomials::VectorField< 3 >( Harmonic4[7] ) , HomogeneousPolynomials::VectorField< 3 >( Harmonic4[8] ) ,
	};

	template< unsigned int Degree >
	void HarmonicValues( Point3D< double > p , double values[ SphericalHarmonics::Dimension< Degree >() ] )
	{
		static_assert( Degree<=SphericalHarmonics::MaxDegree , "[ERROR] Degree exceeds maximum spherical harmonic degree" );
		int idx = 0;

		{
			Point< double ,  1 > m0 = HomogeneousPolynomials::Monomials< 0 >( p );
			for( int i=0 ; i<1 ; i++ ) values[idx++] = Harmonic0[i].value( m0 );
		}
		if( Degree>0 )
		{
			Point< double ,  3 > m1 = HomogeneousPolynomials::Monomials< 1 >( p );
			for( int i=0 ; i<3 ; i++ ) values[idx++] = Harmonic1[i].value( m1 );
		}
		if( Degree>1 )
		{
			Point< double ,  6 > m2 = HomogeneousPolynomials::Monomials< 2 >( p );
			for( int i=0 ; i<5 ; i++ ) values[idx++] = Harmonic2[i].value( m2 );
		}
		if( Degree>2 )
		{
			Point< double , 10 > m3 = HomogeneousPolynomials::Monomials< 3 >( p );
			for( int i=0 ; i<7 ; i++ ) values[idx++] = Harmonic3[i].value( m3 );
		}
		if( Degree>3 )
		{
			Point< double , 15 > m4 = HomogeneousPolynomials::Monomials< 4 >( p );
			for( int i=0 ; i<9 ; i++ ) values[idx++] = Harmonic4[i].value( m4 );
		}
	}
	template< unsigned int Degree >
	void HarmonicGradients( Point3D< double > p , Point3D< double > gradients[ SphericalHarmonics::Dimension< Degree >() ] )
	{
		static_assert( Degree<=SphericalHarmonics::MaxDegree , "[ERROR] Degree exceeds maximum spherical harmonic degree" );
		int idx = 0;

		{
			for( int i=0 ; i<1 ; i++ ) gradients[idx++] = Point3D< double >();
		}
		if( Degree>0 )
		{
			Point< double ,  1 > m1 = HomogeneousPolynomials::Monomials< 0 >( p );
			for( int i=0 ; i<3 ; i++ ) gradients[idx++] = DHarmonic1[i].value( m1 );
		}
		if( Degree>1 )
		{
			Point< double ,  3 > m2 = HomogeneousPolynomials::Monomials< 1 >( p );
			for( int i=0 ; i<5 ; i++ ) gradients[idx++] = DHarmonic2[i].value( m2 );
		}
		if( Degree>2 )
		{
			Point< double ,  6 > m3 = HomogeneousPolynomials::Monomials< 2 >( p );
			for( int i=0 ; i<7 ; i++ ) gradients[idx++] = DHarmonic3[i].value( m3 );
		}
		if( Degree>3 )
		{
			Point< double , 10 > m4 = HomogeneousPolynomials::Monomials< 3 >( p );
			for( int i=0 ; i<9 ; i++ ) gradients[idx++] = DHarmonic4[i].value( m4 );
		}
		for( unsigned int i=0 ; i<SphericalHarmonics::Dimension< Degree >() ; i++ ) gradients[i] -= p * Point3D< double >::Dot( gradients[i] , p );
	}

	template< unsigned int Degree >
	Point3D< double > Gradient( Point< double , SphericalHarmonics::Dimension< Degree >() > c , Point3D< double > p )
	{
		static_assert( Degree<=SphericalHarmonics::MaxDegree , "[ERROR] Degree exceeds maximum spherical harmonic degree" );
		Point3D< double > g , gradients[ SphericalHarmonics::Dimension< Degree >() ];
		SphericalHarmonics::HarmonicGradients< Degree >( p , gradients );
		for( unsigned int i=0 ; i<SphericalHarmonics::Dimension< Degree >() ; i++ ) g += gradients[i] * c[i];
		return g;
	}
}

#endif // SPHERICAL_HARMONICS_INCLUDED