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

#ifndef FOURIER_INCLUDED
#define FOURIER_INCLUDED
#include <fftw3.h>
#include "Array.h"

#ifndef PI
#define PI     3.1415926535897932384
#endif
#ifndef SQRT_2
#define SQRT_2 1.4142135623709504880
#endif

#define FIX_SCALING 0 // Corrects the FFTs so that everything is done w.r.t. an orthonormal basis
#define FORCE_FFTW_PRESERVE_INPUT 1

#include "RotationGrid.h"
#include "SphericalGrid.h"
#include "SquareGrid.h"
#include "CircularArray.h"
#include "Complex.h"

// This templated class represents the fourier coefficients of a real valued, 1D signal
// Because the input signal generating the key is assumed to be real, we know that the
// negative Fourier coefficients are just conjugates of their positive counterparts, and
// we only store the non-negative coefficients.
template< class Real=float >
class FourierKey1D
{
	int res , dim;
	Pointer( Complex< Real > ) values;
public:
	
	FourierKey1D ( void );
	~FourierKey1D( void );

	// Returns the complex dimension of the array
	int size( void ) const;
	// Returns the resolution of the signal
	int resolution( void ) const;
	// Allocates memory for the array
	int resize( int resolution , bool clr=true );

	// Clears the values of the array to 0
	void clear( void );

	// Returns a reference to the indexed array element
	Complex< Real >  operator() ( int x ) const;
	Complex< Real >& operator() ( int x );

	// Returns the square of the L2-norm of the array elements
	Real squareNorm( void ) const;

	// Reads in an array from the specified file
	int read( const char* fileName );
	int read( FILE* fp );

	// Writes out the array to the specified file
	int write( const char* fileName ) const;
	int write( FILE* fp ) const;

	// Returns the square of the L2-difference between two arrays
	static Real SquareDifference( const FourierKey1D& f1 , const FourierKey1D& f2 );

	// Returns the dot-product of two arrays
	static Real Dot( const FourierKey1D& f1 , const FourierKey1D& f2 );
};

// This templated class represents the fourier coefficients of a real valued, 2D signal
// As in the 1D case, since we assume that the original signal is real, we only store
// half the coefficients.
template< class Real=float >
class FourierKey2D
{
	int res,dim;
	Pointer( Complex< Real > ) values;
public:
	
	FourierKey2D ( void );
	~FourierKey2D( void );

	// Returns the complex dimension of the array
	int size( void ) const;
	// Returns the resolution of the signal
	int resolution( void ) const;
	// Allocates memory for the array
	int resize( int resolution , bool clr=true );

	// Clears the values of the array to 0
	void clear( void );

	// Returns a reference to the indexed array element.
	// Because we are only storing half the coefficients, we have
	// 0 <= i < resolution()
	// 0 <= j < size()
	Complex< Real >  operator() ( int i , int j ) const;
	Complex< Real >& operator() ( int i , int j );

	// Returns the square of the L2-norm of the array elements
	Real squareNorm( void ) const;

	// Reads in an array from the specified file
	int read( const char* fileName );
	int read( FILE* fp );

	// Writes out the array to the specified file
	int write( const char* fileName ) const;
	int write( FILE* fp ) const;

	// Returns the square of the L2-difference between two arrays
	static Real SquareDifference( const FourierKey2D& f1 , const FourierKey2D& f2 );

	// Returns the dot-product of two arrays
	static Real Dot( const FourierKey2D& f1 , const FourierKey2D& f2 );
};

// This templated class represents the fourier coefficients of a real valued, signal
// on the surface of the sphers.
// Since we assume that the original signal is real, we only store half the coefficients.
template<class Real=float>
class FourierKeyS2
{
	int bw;
	Pointer( Complex< Real > ) values;
public:
	FourierKeyS2 ( void );
	~FourierKeyS2( void );

	// Returns the complex dimension of the array
	int bandWidth( void ) const;
	// Returns the resolution of the signal
	int resolution( void ) const;
	// Allocates memory for the array
	int resize( int resolution , bool clr=false );

	// Clears the values of the array to 0
	void clear( void );

	// Returns a reference to the indexed array element
	// In this indexing method, "f" represents the frequency and "i"
	// represents the index of the function within the frequency:
	// 0 <= f < bandWidth()
	// 0 <= i <= f
	Complex<Real>  operator() ( int f , int i ) const;
	Complex<Real>& operator() ( int f , int i );

	// Returns the square of the L2-norm of the array elements
	Real squareNorm( void ) const;

	// Reads in an array from the specified file
	int read( const char* fileName );
	int read( FILE* fp );

	// Writes out the array to the specified file
	int write( const char* fileName ) const;
	int write( FILE* fp ) const;

	// Returns the square of the L2-difference between two arrays
	static Real SquareDifference( const FourierKeyS2& f1 , const FourierKeyS2& f2 );

	// Returns the dot-product of two arrays
	static Real Dot( const FourierKeyS2& f1 , const FourierKeyS2& f2 );

	static int Entries( int bw );
};

// This templated class represents the fourier coefficients of a real valued, signal
// on the group of 3D rotations.
// Since we assume that the original signal is real, we only store half the coefficients.
template<class Real=float>
class FourierKeySO3
{
	int bw;
	Pointer( Complex< Real > ) values;
public:
	FourierKeySO3 ( void );
	~FourierKeySO3( void );

	// Returns the complex dimension of the array
	int bandWidth( void ) const;
	// Returns the resolution of the signal
	int resolution( void ) const;
	// Allocates memory for the array
	int resize( int resolution , bool clr=true );

	// Clears the values of the array to 0
	void clear( void );

	// Returns a reference to the indexed array element
	// In this indexing method, "f" represents the frequency and the indices "i" and "j"
	// the function within the frequency:
	// 0 <= f < bandWidth()
	//  0 <= i <= f
	// -f <= j <= f
	Complex<Real>  operator() ( int f , int i , int j ) const;
	Complex<Real>& operator() ( int f , int i , int j );

	// Returns the square of the L2-norm of the array elements
	Real squareNorm( void ) const;

	// Reads in an array from the specified file
	int read( const char* fileName );
	int read( FILE* fp );

	// Writes out the array to the specified file
	int write( const char* fileName ) const;
	int write( FILE* fp ) const;

	// Returns the square of the L2-difference between two arrays
	static Real SquareDifference( const FourierKeySO3& f1 , const FourierKeySO3& f2 );

	// Returns the dot-product of two arrays
	static Real Dot( const FourierKeySO3& f1 , const FourierKeySO3& f2 );

	static int Entries( int bw );
};

// This templated class is responsible for computing the forward and inverse
// Fourier transforms of periodic functions defined either in 1D or in 2D
template<class Real=float>
class FourierTransform
{
	void _ForwardFourier1D( int res , Pointer( Real ) , Pointer( Complex< Real > ) );
	void _ForwardFourier2D( int res , Pointer( Real ) , Pointer( Complex< Real > ) );
	void _InverseFourier1D( int res , Pointer( Complex< Real > ) , Pointer( Real ) );
	void _InverseFourier2D( int res , Pointer( Complex< Real > ) , Pointer( Real ) );
public:
	// This method takes in a real valued function on the circle and computes
	// the Fourier coefficients, writing them into "key"
	int	ForwardFourier( CircularArray<Real>& g , FourierKey1D<Real>& key );

	// This method takes the Fourier coefficients of a real valued function
	// on the circle and returns the originial signal, writing it into "g"
	int InverseFourier( FourierKey1D<Real>& key , CircularArray<Real>& g );

	// This method takes in a real valued function on a 2D grid and computes
	// the Fourier coefficients, writing them into "key"
	int	ForwardFourier( SquareGrid<Real>& g , FourierKey2D<Real>& key );

	// This method takes the Fourier coefficients of a real valued function
	// on a 2D grid and returns the originial signal, writing it into "g"
	// Warning!!! This method will not preserve the input in Fourier key
	int InverseFourier( FourierKey2D<Real>& key , SquareGrid<Real>& g );

	// This static method returns the bandwidth up to which Fourier coefficients
	// will be computed for a signal of resolution "res"
	static int BandWidth( const int& res ){ return (res>>1)+1; }
};

// This templated class is responsible for computing the forward and inverse
// spherical harmonic transforms of functions defined on the sphere. It allocates
// the appropriate scratch space, based on the resolution of the grid for which
// transforms will be computed.
template<class Real=float>
class HarmonicTransform
{
	class ScratchSpace
	{
	public:
		int bw;
#if 1
		Pointer( Real ) workSpace;
		Pointer( Real ) resultSpace;
		Pointer( Real ) transposeResultSpace;
#else
		Real *workSpace,*resultSpace,*transposeResultSpace;
#endif
		Real **table , **transposeTable;
		ScratchSpace(void);
		~ScratchSpace(void);

		void resize(const int& bw);
	};
	ScratchSpace scratch;
public:
	// This method allocates the appropriate amount of scratch space, given
	// the resolution of the signals to be transformed.
	// You do not actually have to call this method, as the transforms will
	// automatically detect if the resolution of the signal doesn't match the
	// resolution of the scratch space, and will call resize if they don't.
	void resize(const int& resolution);

	// This method takes in a real valued function on a sphere and computes
	// the spherical harmonic coefficients, writing them into "key"
	int ForwardFourier(SphericalGrid<Real>& g,FourierKeyS2<Real>& key);

	// This method takes the spherical harmonic coefficients of a real valued function
	// on a sphere and returns the originial signal, writing it into "g"
	int InverseFourier(FourierKeyS2<Real>& key,SphericalGrid<Real>& g);
};

// This templated class is responsible for computing the inverse
// Wigner-D transform of functions defined on the group of 3D rotations. It
// allocates the appropriate scratch space, based on the resolution of the
// grid for which transforms will be computed.
template<class Real=float>
class WignerTransform{
	class ScratchSpace{
	public:
		int bw;
		Pointer( fftw_complex ) coeffs;
		Pointer( fftw_complex ) data;
		Pointer( fftw_complex ) workspace_cx;
		Pointer( fftw_complex ) workspace_cx2;
		Pointer( double ) workspace_re;
		fftw_plan p;
		ScratchSpace(void);
		~ScratchSpace(void);

		void resize(const int& bw);
	};
	ScratchSpace scratch;
public:
	// This method allocates the appropriate amount of scratch space, given
	// the resolution of the signals to be transformed.
	// You do not actually have to call this method, as the transforms will
	// automatically detect if the resolution of the signal doesn't match the
	// resolution of the scratch space, and will call resize if they don't.
	void resize(const int& resolution);

	// This method takes the spherical harmonic coefficients of a real valued function
	// on a sphere and returns the originial signal, writing it into "g"
	int InverseFourier(FourierKeySO3<Real>& key,RotationGrid<Real>& g);
};

#include "Fourier1D.inl"
#include "Fourier2D.inl"
#include "FourierS2.inl"
#include "FourierSO3.inl"
#endif // FOURIER_INCLUDED
