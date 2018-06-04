/*
Copyright (c) 2012, Michael Kazhdan
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
#ifndef SPARSE_MATRIX_INTERFACE_INCLUDED
#define SPARSE_MATRIX_INTERFACE_INCLUDED

#define OWN_MEMORY 1
#define USE_SIZE_T_INDEX 0
#define FORCE_TWO_BYTE_ALIGNMENT 1
#define NEW_TRANSPOSE_MULTIPLY 1
#include "Array.h"
#include <vector>


#if FORCE_TWO_BYTE_ALIGNMENT
#pragma pack(push)
#pragma pack(2)
#endif // FORCE_TWO_BYTE_ALIGNMENT
template< class T , class IndexType >
struct MatrixEntry
{
	MatrixEntry( void )       { N =-1 , Value = 0; }
	MatrixEntry( int i )      { N = i , Value = 0; }
	MatrixEntry( int n , T v ){ N = n , Value = v; }
#if USE_SIZE_T_INDEX
	size_t N;
#else // !USE_SIZE_T_INDEX
	IndexType N;
#endif // USE_SIZE_T_INDEX
	T Value;
};
#if FORCE_TWO_BYTE_ALIGNMENT
#pragma pack(pop)
#endif // FORCE_TWO_BYTE_ALIGNMENT

template< class T2 >
class ParallelSolution
{
	template< class T , class const_iterator > friend class SparseMatrixInterface;
	int _iters , _threads , _sliceDependence , _sliceCount;
#if OWN_MEMORY
	std::vector< int > __entriesPerSlice , __slicesPerThread;
#endif // OWN_MEMORY
	const int *_entriesPerSlice , *_slicesPerThread;
	std::vector< Pointer( T2 ) > _pSolution;
	std::vector< int > _startEntry;
	std::vector< std::pair< int , int > > _sliceBounds , _skirtBounds;
	void _clearSkirts( void );
	void _merge( void );
public:
	ParallelSolution( int iters , const std::vector< int >& entriesPerSlice , const std::vector< int >& slicesPerThread , int sliceDependence );
	~ParallelSolution( void );
	void SetFromArray( ConstPointer( T2 )  solution );
	void SetFromArray( const Vector< T2 >& solution ) { SetFromArray( solution+0 ); }
	void SetToArray  ( Pointer( T2 ) solution ) const;
	void SetToArray  ( Vector< T2 >& solution ) const { SetToArray( solution+0 ); }
	void clear( void );
	void synchronize( void );
	int threads( void ) const { return _threads; }
	int sliceStart( int t ) const { if( t<0 || t>=_threads ) return 0 ; return _sliceBounds[t].first;  }
	int sliceEnd  ( int t ) const { if( t<0 || t>=_threads ) return 0 ; return _sliceBounds[t].second; }
	int sliceSize( int t ) const { if( t<0 || t>=_threads ) return 0 ; return _sliceBounds[t].second - _sliceBounds[t].first; }
	int size( int t ) const { if( t<0 || t>=_threads ) return 0 ; return _startEntry[ _sliceBounds[t].second-1 ] + _entriesPerSlice[ _sliceBounds[t].second-1 ] - _startEntry[ _sliceBounds[t].first ]; }
	int size( void ) const { int sz=0 ; for( int i=0 ; i<_threads ; i++ ) sz += size( i ) ; return sz; }
	int offset( int t ) const { if( t<0 || t>_threads ) return 0 ; return _startEntry[ _sliceBounds[t].first ]; }
	void bounds( int t , int& start , int& end ) const { if( t<0 || t>=_threads ) start = end =0 ; else start = _startEntry[_skirtBounds[t].first] - _startEntry[_sliceBounds[t].first] , end = _startEntry[ _skirtBounds[t].second-1 ] + _entriesPerSlice[ _skirtBounds[t].second-1 ] - _startEntry[_sliceBounds[t].first]; }
	ConstPointer( T2 ) operator[] ( int t ) const { if( t<0 || t>=_threads ) return NullPointer< T2 >( ) ; return _pSolution[t]; }
};
template< class T2 >
class ParallelSolution2
{
	template< class T , class const_iterator > friend class SparseMatrixInterface;
	int _iters , _threads , _sliceDependence , _sliceCount , _linesPerSlice;
#if OWN_MEMORY
	std::vector< int > __entriesPerLine , __slicesPerThread;
#endif // OWN_MEMORY
	const int *_entriesPerLine , *_slicesPerThread;
	std::vector< Pointer( T2 ) > _pSolution;
	std::vector< int > _startEntry;
	std::vector< std::pair< int , int > > _sliceBounds , _skirtBounds;
	void _clearSkirts( void );
	void _merge( void );
public:
	ParallelSolution2( int iters , const std::vector< int >& entriesPerLine , int linesPerSlice , const std::vector< int >& slicesPerThread , int sliceDependence );
	~ParallelSolution2( void );
	void SetFromArray( ConstPointer( T2 )  solution );
	void SetFromArray( const Vector< T2 >& solution ){ SetFromArray( solution+0 ); }
	void SetToArray  ( Pointer( T2 ) solution ) const;
	void SetToArray  ( Vector< T2 >& solution ) const { SetToArray( solution+0 ); }
	void clear( void );
	void synchronize( void );
	int threads( void ) const { return _threads; }
	int sliceStart( int t ) const { if( t<0 || t>=_threads ) return 0 ; return _sliceBounds[t].first;  }
	int sliceEnd  ( int t ) const { if( t<0 || t>=_threads ) return 0 ; return _sliceBounds[t].second; }
	int sliceSize( int t ) const { if( t<0 || t>=_threads ) return 0 ; return _sliceBounds[t].second - _sliceBounds[t].first; }
	int size( void ) const { int sz=0 ; for( int i=0 ; i<_threads ; i++ ) sz += size( i ) ; return sz; }
	int size( int t ) const { if( t<0 || t>=_threads ) return 0 ; return _startEntry[ _sliceBounds[t].second*_linesPerSlice - 1 ] + _entriesPerLine[ _sliceBounds[t].second*_linesPerSlice - 1 ] - _startEntry[ _sliceBounds[t].first * _linesPerSlice ]; }
	void bounds( int t , int& start , int& end ) const { if( t<0 || t>=_threads ) start = end =0 ; else start = _startEntry[_skirtBounds[t].first*_linesPerSlice] - _startEntry[_sliceBounds[t].first*_linesPerSlice] , end = _startEntry[ _skirtBounds[t].second*_linesPerSlice-1 ] + _entriesPerLine[ _skirtBounds[t].second*_linesPerSlice-1 ] - _startEntry[_sliceBounds[t].first*_linesPerSlice]; }
	int offset( int t ) const { if( t<0 || t>_threads ) return 0 ; return _startEntry[ _sliceBounds[t].first * _linesPerSlice ]; }
	ConstPointer( T2 ) operator[] ( int t ) const { if( t<0 || t>=_threads ) return NullPointer< T2 >( ) ; return _pSolution[t]; }
};

enum
{
	MULTIPLY_ADD ,
	MULTIPLY_SUBTRACT ,
	MULTIPLY_CLEAR
};
template< class T , class const_iterator > class SparseMatrixInterface
{
public:
	virtual const_iterator begin( int row ) const = 0;
	virtual const_iterator end  ( int row ) const = 0;
	virtual size_t Rows   ( void )          const = 0;
	virtual size_t RowSize( size_t idx )    const = 0;

	size_t Entries( void ) const;

	double SquareNorm( void ) const;
	double SquareASymmetricNorm( void ) const;
	double SquareASymmetricNorm( int& idx1 , int& idx2 ) const;

	template< class T2 > void Multiply           (  ConstPointer( T2 ) In , Pointer( T2 ) Out ) const;
	template< class T2 > void Multiply           (       Pointer( T2 ) In , Pointer( T2 ) Out ) const { return Multiply( ( ConstPointer( T2 ) ) In , Out ); }
	template< class T2 > void Multiply           ( const Vector< T2 >& In , Vector< T2 >& Out ) const { Multiply( In+0 , Out+0 ); }
	template< class T2 > void MultiplyParallel   ( ConstPointer( T2 )  In , Pointer( T2 ) Out , int threads , int multiplyFlag=MULTIPLY_ADD ) const;
	template< class T2 > void MultiplyParallel   (      Pointer( T2 )  In , Pointer( T2 ) Out , int threads , int multiplyFlag=MULTIPLY_ADD ) const { MultiplyParallel( ( ConstPointer(T2) )( In ) , Out , threads , multiplyFlag ); }
	template< class T2 > void MultiplyParallel   ( const Vector< T2 >& In , Vector< T2 >& Out , int threads , int multiplyFlag=MULTIPLY_ADD ) const { MultiplyParallel( In+0 , Out+0 , threads , multiplyFlag ); }
#if NEW_TRANSPOSE_MULTIPLY
	template< class T2 > void MultiplyTranspose  ( ConstPointer( T2 )  In , Pointer( T2 ) Out , int outDim , int multiplyFlag , T (*TransposeFunction)( const T& )=NULL ) const;
	template< class T2 > void MultiplyTranspose  (      Pointer( T2 )  In , Pointer( T2 ) Out , int outDim , int multiplyFlag , T (*TransposeFunction)( const T& )=NULL ) const { MultiplyTranspose( ( ConstPointer( T2 ) ) In , Out , outDim , multiplyFlag , TransposeFunction ); }
	template< class T2 > void MultiplyTranspose  ( const Vector< T2 >& In , Vector< T2 >& Out , int outDim , int multiplyFlag , T (*TransposeFunction)( const T& )=NULL ) const { MultiplyTranspose( In+0 , Out+0 , outDim , multiplyFlag , TransposeFunction ); }
	template< class T2 > Vector< T2 > MultiplyTranspose( const Vector< T2 >& in , int outDim , T (*TransposeFunction)( const T& )=NULL ) const
	{
		Vector< T2 > out( outDim );
		MultiplyTranspose( in , out , outDim , MULTIPLY_CLEAR , TransposeFunction );
		return out;
	}
#else // !NEW_TRANSPOSE_MULTIPLY
	template< class T2 > void MultiplyTranspose  ( ConstPointer( T2 )  In , Pointer( T2 ) Out ) const { MultiplyTranspose( In , Out , NULL ); }
	template< class T2 > void MultiplyTranspose  (      Pointer( T2 )  In , Pointer( T2 ) Out ) const { MultiplyTranspose( ( ConstPointer( T2 ) ) In , Out ); }
	template< class T2 > void MultiplyTranspose  ( const Vector< T2 >& In , Vector< T2 >& Out ) const { MultiplyTranspose( In+0 , Out+0 ); }
	template< class T2 > void MultiplyTranspose  ( ConstPointer( T2 )  In , Pointer( T2 ) Out , T (*TransposeFunction)( const T& ) ) const;
	template< class T2 > void MultiplyTranspose  (      Pointer( T2 )  In , Pointer( T2 ) Out , T (*TransposeFunction)( const T& ) ) const { return MultiplyTranspose( ( ConstPointer( T2 ) ) In , Out , TransposeFunction ); }
	template< class T2 > void MultiplyTranspose  ( const Vector< T2 >& In , Vector< T2 >& Out , T (*TransposeFunction)( const T& ) ) const { MultiplyTranspose( In+0 , Out+0 , TransposeFunction ); }
#endif // NEW_TRANSPOSE_MULTIPLY
	template< class T2 > void BMinusMX           ( ConstPointer( T2 ) B , ConstPointer( T2 ) X , Pointer( T2 ) out ) const;
	template< class T2 > void BMinusMX           (      Pointer( T2 ) B , ConstPointer( T2 ) X , Pointer( T2 ) out ) const { return BMinusMX( ( ConstPointer( T2 ) ) B ,                        X , out ); }
	template< class T2 > void BMinusMX           ( ConstPointer( T2 ) B ,      Pointer( T2 ) X , Pointer( T2 ) out ) const { return BMinusMX(                        B , ( ConstPointer( T2 ) ) X , out ); }
	template< class T2 > void BMinusMX           (      Pointer( T2 ) B ,      Pointer( T2 ) X , Pointer( T2 ) out ) const { return BMinusMX( ( ConstPointer( T2 ) ) B , ( ConstPointer( T2 ) ) X , out ); }
	template< class T2 > void BMinusMX           ( const Vector< T2 >& B , const Vector< T2 >& X , Vector< T2 >& out ) const { BMinusMX( B+0 , X+0 , out+0 ); }

	template< class T2 , class T3 , bool unitDiagonal , bool StrippedDiagonal , bool UseSOR > void SolveGaussSeidel           ( ConstPointer( T3 ) Diagonal , ConstPointer( T2 ) b , int iters , Pointer( T2 ) Solution , T sor ) const;
	template< class T2 , class T3 , bool unitDiagonal , bool StrippedDiagonal , bool UseSOR > void SolveGaussSeidelAndResidual( ConstPointer( T3 ) Diagonal , ConstPointer( T2 ) b , int iters , Pointer( T2 ) Solution , Pointer( T2 ) Residual , T sor ) const;

	template< class T2 >
	int SolveConjugateGradient( const Vector< T2 >& b , Vector< T2 >& x , double (*SquareNorm)( T2 ) , T2 (*Dot)( T2 , T2 ) , const size_t& iters , int threads ) const;


	template< class T2 , class T3 , bool UseSOR > void SolveGaussSeidel( ConstPointer( T3 ) Diagonal , ConstPointer( T2 ) b , int iters , Pointer( T2 ) Solution , T sor , bool unitDiagonal , bool strippedDiagonal ) const
	{
		if( strippedDiagonal )
			if( unitDiagonal ) SolveGaussSeidel< T2 , T3 , true  , true  , UseSOR >( Diagonal , b , iters , Solution , sor );
			else               SolveGaussSeidel< T2 , T3 , false , true  , UseSOR >( Diagonal , b , iters , Solution , sor );
		else
			if( unitDiagonal ) SolveGaussSeidel< T2 , T3 , true  , false , UseSOR >( Diagonal , b , iters , Solution , sor );
			else               SolveGaussSeidel< T2 , T3 , false , false , UseSOR >( Diagonal , b , iters , Solution , sor );
	}
	template< class T2 , class T3 , bool UseSOR > void SolveGaussSeidelAndResidual( ConstPointer( T3 ) Diagonal , ConstPointer( T2 ) b , int iters , Pointer( T2 ) Solution , Pointer( T2 ) Residual , T sor , bool unitDiagonal , bool strippedDiagonal ) const
	{
		if( strippedDiagonal )
			if( unitDiagonal ) SolveGaussSeidelAndResidual< T2 , T3 , true  , true  , UseSOR >( Diagonal , b , iters , Solution , Residual , sor );
			else               SolveGaussSeidelAndResidual< T2 , T3 , false , true  , UseSOR >( Diagonal , b , iters , Solution , Residual , sor );
		else
			if( unitDiagonal ) SolveGaussSeidelAndResidual< T2 , T3 , true  , false , UseSOR >( Diagonal , b , iters , Solution , Residual , sor );
			else               SolveGaussSeidelAndResidual< T2 , T3 , false , false , UseSOR >( Diagonal , b , iters , Solution , Residual , sor );
	}
	template< class T2 , class T3 , bool UseSOR > void SolveGaussSeidel           ( const Vector< T3 >& Diagonal , const Vector< T2 >& b , int iters , Vector< T2 >& Solution ,                          T sor , bool unitDiagonal , bool strippedDiagonal ) const { SolveGaussSeidel           < T2 , T3 , UseSOR >( Diagonal+0 , b+0 , iters , Solution+0 ,              sor , unitDiagonal , strippedDiagonal ); }
	template< class T2 , class T3 , bool UseSOR > void SolveGaussSeidelAndResidual( const Vector< T3 >& Diagonal , const Vector< T2 >& b , int iters , Vector< T2 >& Solution , Vector< T2 >& Residual , T sor , bool unitDiagonal , bool strippedDiagonal ) const { SolveGaussSeidelAndResidual< T2 , T3 , UseSOR >( Diagonal+0 , b+0 , iters , Solution+0 , Residual+0 , sor , unitDiagonal , strippedDiagonal ); }

	template< class T2 > void MultiplyTransposeParallel( ConstPointer( T2 ) in ,                       Pointer( T2 ) out , const std::vector< int >& entriesPerLine , int linesPerSlice , const std::vector< int >& slicesPerThread , int sliceDependence ) const;
	template< class T2 > void MultiplyParallel         ( ConstPointer( T2 ) in ,                       Pointer( T2 ) out , const std::vector< int >& entriesPerLine , int linesPerSlice , const std::vector< int >& slicesPerThread ) const;
	template< class T2 > void BMinusMXParallel         ( ConstPointer( T2 ) B , ConstPointer( T2 ) X , Pointer( T2 ) out , const std::vector< int >& entriesPerLine , int linesPerSlice , const std::vector< int >& slicesPerThread ) const;

	template< class T2 > void MultiplyTransposeParallel(      Pointer( T2 ) in ,                       Pointer( T2 ) out , const std::vector< int >& entriesPerLine , int linesPerSlice , const std::vector< int >& slicesPerThread , int sliceDependence ) const { MultiplyTransposeParallel( ( ConstPointer( T2 ) ) in ,                              out , entriesPerLine , linesPerSlice , slicesPerThread , sliceDependence ); }
	template< class T2 > void MultiplyParallel         (      Pointer( T2 ) in ,                       Pointer( T2 ) out , const std::vector< int >& entriesPerLine , int linesPerSlice , const std::vector< int >& slicesPerThread )                       const { MultiplyParallel         ( ( ConstPointer( T2 ) ) in ,                              out , entriesPerLine , linesPerSlice , slicesPerThread ); }
	template< class T2 > void BMinusMXParallel         (      Pointer( T2 ) B , ConstPointer( T2 ) X , Pointer( T2 ) out , const std::vector< int >& entriesPerLine , int linesPerSlice , const std::vector< int >& slicesPerThread )                       const { BMinusMXParallel         ( ( ConstPointer( T2 ) ) B ,                        X ,    out , entriesPerLine , linesPerSlice , slicesPerThread ); }
	template< class T2 > void BMinusMXParallel         ( ConstPointer( T2 ) B ,      Pointer( T2 ) X , Pointer( T2 ) out , const std::vector< int >& entriesPerLine , int linesPerSlice , const std::vector< int >& slicesPerThread )                       const { BMinusMXParallel         (                        B , ( ConstPointer( T2 ) ) X ,    out , entriesPerLine , linesPerSlice , slicesPerThread ); }
	template< class T2 > void BMinusMXParallel         (      Pointer( T2 ) B ,      Pointer( T2 ) X , Pointer( T2 ) out , const std::vector< int >& entriesPerLine , int linesPerSlice , const std::vector< int >& slicesPerThread )                       const { BMinusMXParallel         ( ( ConstPointer( T2 ) ) B , ( ConstPointer( T2 ) ) X ,    out , entriesPerLine , linesPerSlice , slicesPerThread ); }

	template< class T2 > void MultiplyTransposeParallel( const Vector< T2 >& in ,                         Vector< T2 >& out , const std::vector< int >& entriesPerLine , int linesPerSlice , const std::vector< int >& slicesPerThread , int sliceDependence ) const { MultiplyTransposeParallel(      in+0 , out+0 , entriesPerLine , linesPerSlice , slicesPerThread , sliceDependence ); }
	template< class T2 > void MultiplyParallel         ( const Vector< T2 >& in ,                         Vector< T2 >& out , const std::vector< int >& entriesPerLine , int linesPerSlice , const std::vector< int >& slicesPerThread )                       const { MultiplyParallel         (      in+0 , out+0 , entriesPerLine , linesPerSlice , slicesPerThread ); }
	template< class T2 > void BMinusMXParallel         ( const Vector< T2 >& B  , const Vector< T2 >& X , Vector< T2 >& out , const std::vector< int >& entriesPerLine , int linesPerSlice , const std::vector< int >& slicesPerThread )                       const { BMinusMXParallel         ( B+0 , X+0 , out+0 , entriesPerLine , linesPerSlice , slicesPerThread ); }

	template< class T2 > void MultiplyTransposeParallel( ConstPointer( T2 ) in ,                       Pointer( T2 ) out , const std::vector< int >& entriesPerSlice , const std::vector< int >& slicesPerThread , int sliceDependence ) const;
	template< class T2 > void MultiplyParallel         ( ConstPointer( T2 ) in ,                       Pointer( T2 ) out , const std::vector< int >& entriesPerSlice , const std::vector< int >& slicesPerThread ) const;
	template< class T2 > void BMinusMXParallel         ( ConstPointer( T2 ) B , ConstPointer( T2 ) X , Pointer( T2 ) out , const std::vector< int >& entriesPerSlice , const std::vector< int >& slicesPerThread ) const;

	template< class T2 > void MultiplyTransposeParallel(      Pointer( T2 ) in ,                       Pointer( T2 ) out , const std::vector< int >& entriesPerSlice , const std::vector< int >& slicesPerThread , int sliceDependence ) const { MultiplyTransposeParallel( ( ConstPointer( T2 ) ) in , out , entriesPerSlice , slicesPerThread , sliceDependence ); }
	template< class T2 > void MultiplyParallel         (      Pointer( T2 ) in ,                       Pointer( T2 ) out , const std::vector< int >& entriesPerSlice , const std::vector< int >& slicesPerThread )                       const { MultiplyParallel         ( ( ConstPointer( T2 ) ) in , out , entriesPerSlice , slicesPerThread ); }
	template< class T2 > void BMinusMXParallel         (      Pointer( T2 ) B , ConstPointer( T2 ) X , Pointer( T2 ) out , const std::vector< int >& entriesPerSlice , const std::vector< int >& slicesPerThread )                       const { BMinusMXParallel         ( ( ConstPointer( T2 ) ) B ,                        X ,    out , entriesPerSlice , slicesPerThread ); }
	template< class T2 > void BMinusMXParallel         ( ConstPointer( T2 ) B ,      Pointer( T2 ) X , Pointer( T2 ) out , const std::vector< int >& entriesPerSlice , const std::vector< int >& slicesPerThread )                       const { BMinusMXParallel         (                        B , ( ConstPointer( T2 ) ) X ,    out , entriesPerSlice , slicesPerThread ); }
	template< class T2 > void BMinusMXParallel         (      Pointer( T2 ) B ,      Pointer( T2 ) X , Pointer( T2 ) out , const std::vector< int >& entriesPerSlice , const std::vector< int >& slicesPerThread )                       const { BMinusMXParallel         ( ( ConstPointer( T2 ) ) B , ( ConstPointer( T2 ) ) X ,    out , entriesPerSlice , slicesPerThread ); }

	template< class T2 > void MultiplyTransposeParallel( const Vector< T2 >& in ,                         Vector< T2 >& out , const std::vector< int >& entriesPerSlice , const std::vector< int >& slicesPerThread , int sliceDependence ) const { MultiplyTransposeParallel(      in+0 , out+0 , entriesPerSlice , slicesPerThread , sliceDependence ); }
	template< class T2 > void MultiplyParallel         ( const Vector< T2 >& in ,                         Vector< T2 >& out , const std::vector< int >& entriesPerSlice , const std::vector< int >& slicesPerThread )                       const { MultiplyParallel         (      in+0 , out+0 , entriesPerSlice , slicesPerThread ); }
	template< class T2 > void BMinusMXParallel         ( const Vector< T2 >&  B , const Vector< T2 >& X , Vector< T2 >& out , const std::vector< int >& entriesPerSlice , const std::vector< int >& slicesPerThread )                       const { BMinusMXParallel         ( B+0 , X+0 , out+0 , entriesPerSlice , slicesPerThread ); }

	template< class T2 > void MultiplyParallel( const ParallelSolution< T2 >* in , ParallelSolution< T2 >* out ) const;
	template< class T2 > void MultiplyParallel( ConstPointer( T2 ) in , ParallelSolution< T2 >* out ) const;
	template< class T2 > void MultiplyParallel( const Vector< T2 >& in , ParallelSolution< T2 >* out ) const { MultiplyParallel( in+0 , out ); }
	template< class T2 > void MultiplyTransposeParallel( const ParallelSolution< T2 >* in , ParallelSolution< T2 >* out ) const;

	template< class T2 , class T3 , bool UnitDiagonal , bool StrippedDiagonal , bool UseSOR >
	void SolveGaussSeidelParallel( ConstPointer( T3 ) Diagonal , ConstPointer( T2 ) b , int iters , ParallelSolution< T2 >* Solution , T sor ) const;

	template< class T2 , class T3 , bool UnitDiagonal , bool StrippedDiagonal , bool UseSOR >
	void SolveGaussSeidelAndResidualParallel( ConstPointer( T3 ) Diagonal , ConstPointer( T2 ) b , int iters , ParallelSolution< T2 >* Solution , Pointer( T2 ) Residual , T sor ) const;

	template< class T2 , class T3 , bool UseSOR >
	void SolveGaussSeidelParallel( ConstPointer( T3 ) Diagonal , ConstPointer( T2 ) b , int iters , ParallelSolution< T2 >* Solution , T sor ,	bool unitDiagonal , bool strippedDiagonal ) const
	{
		if( strippedDiagonal )
			if( unitDiagonal ) SolveGaussSeidelParallel< T2 , T3 , true  , true  , UseSOR >( Diagonal , b , iters , Solution , sor );
			else               SolveGaussSeidelParallel< T2 , T3 , false , true  , UseSOR >( Diagonal , b , iters , Solution , sor );
		else
			if( unitDiagonal ) SolveGaussSeidelParallel< T2 , T3 , true  , false , UseSOR >( Diagonal , b , iters , Solution , sor );
			else               SolveGaussSeidelParallel< T2 , T3 , false , false , UseSOR >( Diagonal , b , iters , Solution , sor );
	}
	template< class T2 , class T3 , bool UseSOR >
	void SolveGaussSeidelAndResidualParallel( ConstPointer( T3 ) Diagonal , ConstPointer( T2 ) b , int iters , ParallelSolution< T2 >* Solution , Pointer( T2 ) Residual ,  T sor , bool unitDiagonal , bool strippedDiagonal ) const
	{
		if( strippedDiagonal )
			if( unitDiagonal ) SolveGaussSeidelAndResidualParallel< T2 , T3 ,true  , true  , UseSOR >( Diagonal , b , iters , Solution , Residual , sor );
			else               SolveGaussSeidelAndResidualParallel< T2 , T3 ,false , true  , UseSOR >( Diagonal , b , iters , Solution , Residual , sor );
		else
			if( unitDiagonal ) SolveGaussSeidelAndResidualParallel< T2 , T3 ,true  , false , UseSOR >( Diagonal , b , iters , Solution , Residual , sor );
			else               SolveGaussSeidelAndResidualParallel< T2 , T3 ,false , false , UseSOR >( Diagonal , b , iters , Solution , Residual , sor );
	}
	template< class T2 , class T3 , bool UseSOR > void SolveGaussSeidelParallel           ( const Vector< T3 >& Diagonal , const Vector< T2 >& b , int iters , ParallelSolution< T2 >* Solution ,                          T sor , bool unitDiagonal , bool strippedDiagonal ) const { SolveGaussSeidelParallel           < T2 , T3 , UseSOR >( Diagonal+0 , b+0 , iters , Solution ,              sor , unitDiagonal , strippedDiagonal ); }
	template< class T2 , class T3 , bool UseSOR > void SolveGaussSeidelAndResidualParallel( const Vector< T3 >& Diagonal , const Vector< T2 >& b , int iters , ParallelSolution< T2 >* Solution , Vector< T2 >& Residual , T sor , bool unitDiagonal , bool strippedDiagonal ) const { SolveGaussSeidelAndResidualParallel< T2 , T3 , UseSOR >( Diagonal+0 , b+0 , iters , Solution , Residual+0 , sor , unitDiagonal , strippedDiagonal ); }

	template< class T2 > void MultiplyParallel( const ParallelSolution2< T2 >* in , ParallelSolution2< T2 >* out ) const;
	template< class T2 > void MultiplyParallel( ConstPointer( T2 ) in , ParallelSolution2< T2 >* out ) const;
	template< class T2 > void MultiplyParallel( const Vector< T2 >& in , ParallelSolution2< T2 >* out ) const { MultiplyParallel( in+0 , out ); }
	template< class T2 > void MultiplyTransposeParallel( const ParallelSolution2< T2 >* in , ParallelSolution2< T2 >* out ) const;

	template< class T2 , class T3 , bool UnitDiagonal , bool StrippedDiagonal , bool UseSOR >
	void SolveGaussSeidelParallel( ConstPointer( T3 ) Diagonal , ConstPointer( T2 ) b , int iters , ParallelSolution2< T2 >* Solution , T sor ) const;

	template< class T2 , class T3 , bool UnitDiagonal , bool StrippedDiagonal , bool UseSOR >
	void SolveGaussSeidelAndResidualParallel( ConstPointer( T3 ) Diagonal , ConstPointer( T2 ) b , int iters , ParallelSolution2< T2 >* Solution , Pointer( T2 ) Residual , T sor ) const;

	template< class T2 , class T3 , bool UnitDiagonal , bool StrippedDiagonal >
	void ResidualParallel( ConstPointer( T3 ) Diagonal , ConstPointer( T2 ) b , ParallelSolution2< T2 >* Solution , Pointer( T2 ) Residual ) const;

	template< class T2 , class T3 , bool UseSOR >
	void SolveGaussSeidelParallel( ConstPointer( T3 ) Diagonal , ConstPointer( T2 ) b , int iters , ParallelSolution2< T2 >* Solution , T sor ,	bool unitDiagonal , bool strippedDiagonal ) const
	{
		if( strippedDiagonal )
			if( unitDiagonal ) SolveGaussSeidelParallel< T2 , T3 , true  , true  , UseSOR >( Diagonal , b , iters , Solution , sor );
			else               SolveGaussSeidelParallel< T2 , T3 , false , true  , UseSOR >( Diagonal , b , iters , Solution , sor );
		else
			if( unitDiagonal ) SolveGaussSeidelParallel< T2 , T3 , true  , false , UseSOR >( Diagonal , b , iters , Solution , sor );
			else               SolveGaussSeidelParallel< T2 , T3 , false , false , UseSOR >( Diagonal , b , iters , Solution , sor );
	}
	template< class T2 , class T3 , bool UseSOR >
	void SolveGaussSeidelAndResidualParallel( ConstPointer( T3 ) Diagonal , ConstPointer( T2 ) b , int iters , ParallelSolution2< T2 >* Solution , Pointer( T2 ) Residual ,  T sor , bool unitDiagonal , bool strippedDiagonal ) const
	{
		if( strippedDiagonal )
			if( unitDiagonal ) SolveGaussSeidelAndResidualParallel< T2 , T3 , true  , true  , UseSOR >( Diagonal , b , iters , Solution , Residual , sor );
			else               SolveGaussSeidelAndResidualParallel< T2 , T3 , false , true  , UseSOR >( Diagonal , b , iters , Solution , Residual , sor );
		else
			if( unitDiagonal ) SolveGaussSeidelAndResidualParallel< T2 , T3 , true  , false , UseSOR >( Diagonal , b , iters , Solution , Residual , sor );
			else               SolveGaussSeidelAndResidualParallel< T2 , T3 , false , false , UseSOR >( Diagonal , b , iters , Solution , Residual , sor );
	}
	template< class T2 , class T3 >
	void ResidualParallel( ConstPointer( T3 ) Diagonal , ConstPointer( T2 ) b , ParallelSolution2< T2 >* Solution , Pointer( T2 ) Residual , bool unitDiagonal , bool strippedDiagonal ) const
	{
		if( strippedDiagonal )
			if( unitDiagonal ) ResidualParallel< T2 , T3 , true  , true  >( Diagonal , b , Solution , Residual );
			else               ResidualParallel< T2 , T3 , false , true  >( Diagonal , b , Solution , Residual );
		else
			if( unitDiagonal ) ResidualParallel< T2 , T3 , true  , false >( Diagonal , b , Solution , Residual );
			else               ResidualParallel< T2 , T3 , false , false >( Diagonal , b , Solution , Residual );
	}
	template< class T2 , class T3 , bool UseSOR > void SolveGaussSeidelParallel           ( const Vector< T3 >& Diagonal , const Vector< T2 >& b , int iters , ParallelSolution2< T2 >* Solution ,                          T sor ,	bool unitDiagonal , bool strippedDiagonal ) const { SolveGaussSeidelParallel           < T2 , T3 , UseSOR >( Diagonal+0 , b+0 , iters , Solution ,              sor , unitDiagonal , strippedDiagonal ); }
	template< class T2 , class T3 , bool UseSOR > void SolveGaussSeidelAndResidualParallel( const Vector< T3 >& Diagonal , const Vector< T2 >& b , int iters , ParallelSolution2< T2 >* Solution , Vector< T2 >& Residual , T sor , bool unitDiagonal , bool strippedDiagonal ) const { SolveGaussSeidelAndResidualParallel< T2 , T3 , UseSOR >( Diagonal+0 , b+0 , iters , Solution , Residual+0 , sor , unitDiagonal , strippedDiagonal ); }
	template< class T2 , class T3 , bool UseSOR > void                    ResidualParallel( const Vector< T3 >& Diagonal , const Vector< T2 >& b ,             ParallelSolution2< T2 >* Solution , Vector< T2 >& Residual ,         bool unitDiagonal , bool strippedDiagonal ) const {                    ResidualParallel< T2 , T3 , UseSOR >( Diagonal+0 , b+0 ,         Solution , Residual+0 ,       unitDiagonal , strippedDiagonal ); }

	template< class T2 , class T3 , bool UnitDiagonal , bool StrippedDiagonal , bool UseSOR >
	void SolveGaussSeidelParallel
	(
		ConstPointer( T3 ) Diagonal , ConstPointer( T2 ) b , int iters , Pointer( T2 ) Solution , T sor ,
		const std::vector< int >& entriesPerSlice , const std::vector< int >& slicesPerThread , int sliceDependence
	) const;
	template< class T2 , class T3 , bool UnitDiagonal , bool StrippedDiagonal , bool UseSOR >
	void SolveGaussSeidelAndResidualParallel
	(
		ConstPointer( T3 ) Diagonal , ConstPointer( T2 ) b , int iters , Pointer( T2 ) Solution , Pointer( T2 ) Residual , T sor ,
		const std::vector< int >& entriesPerSlice , const std::vector< int >& slicesPerThread , int sliceDependence
	) const;

	template< class T2 , class T3 , bool UseSOR >
	void SolveGaussSeidelParallel
		(
		ConstPointer( T3 ) Diagonal , ConstPointer( T2 ) b , int iters , Pointer( T2 ) Solution , T sor ,
		const std::vector< int >& entriesPerSlice , const std::vector< int >& slicesPerThread , int sliceDependence , bool unitDiagonal , bool strippedDiagonal
		) const
	{
		if( strippedDiagonal )
			if( unitDiagonal ) SolveGaussSeidelParallel< T2 , T3 , true  , true  , UseSOR >( Diagonal , b , iters , Solution , sor , entriesPerSlice , slicesPerThread , sliceDependence );
			else               SolveGaussSeidelParallel< T2 , T3 , false , true  , UseSOR >( Diagonal , b , iters , Solution , sor , entriesPerSlice , slicesPerThread , sliceDependence );
		else
			if( unitDiagonal ) SolveGaussSeidelParallel< T2 , T3 , true  , false , UseSOR >( Diagonal , b , iters , Solution , sor , entriesPerSlice , slicesPerThread , sliceDependence );
			else               SolveGaussSeidelParallel< T2 , T3 , false , false , UseSOR >( Diagonal , b , iters , Solution , sor , entriesPerSlice , slicesPerThread , sliceDependence );
	}
	template< class T2 , class T3 , bool UseSOR >
	void SolveGaussSeidelAndResidualParallel
		(
		ConstPointer( T3 ) Diagonal , ConstPointer( T2 ) b , int iters , Pointer( T2 ) Solution , Pointer( T2 ) Residual ,  T sor ,
		const std::vector< int >& entriesPerSlice , const std::vector< int >& slicesPerThread , int sliceDependence , bool unitDiagonal , bool strippedDiagonal
		) const
	{
		if( strippedDiagonal )
			if( unitDiagonal ) SolveGaussSeidelAndResidualParallel< T2 , T3 , true  , true  , UseSOR >( Diagonal , b , iters , Solution , Residual , sor , entriesPerSlice , slicesPerThread , sliceDependence );
			else               SolveGaussSeidelAndResidualParallel< T2 , T3 , false , true  , UseSOR >( Diagonal , b , iters , Solution , Residual , sor , entriesPerSlice , slicesPerThread , sliceDependence );
		else
			if( unitDiagonal ) SolveGaussSeidelAndResidualParallel< T2 , T3 , true  , false , UseSOR >( Diagonal , b , iters , Solution , Residual , sor , entriesPerSlice , slicesPerThread , sliceDependence );
			else               SolveGaussSeidelAndResidualParallel< T2 , T3 , false , false , UseSOR >( Diagonal , b , iters , Solution , Residual , sor , entriesPerSlice , slicesPerThread , sliceDependence );
	}
	template< class T2 , class T3 , bool UseSOR >
	void SolveGaussSeidelParallel
		(
		const Vector< T3 >& Diagonal , const Vector< T2 >& b , int iters , Vector< T2 >& Solution , T sor ,
		const std::vector< int >& entriesPerSlice , const std::vector< int >& slicesPerThread , int sliceDependence , bool unitDiagonal , bool strippedDiagonal
		) const { SolveGaussSeidelParallel< T2 , T3 , UseSOR >( Diagonal+0 , b+0 , iters , Solution+0 , sor , entriesPerSlice , slicesPerThread , sliceDependence , unitDiagonal , strippedDiagonal ); }
	template< class T2 , class T3 , bool UseSOR >
	void SolveGaussSeidelAndResidualParallel
		(
		const Vector< T3 >& Diagonal , const Vector< T2 >& b , int iters , Vector< T2 >& Solution , Vector< T2 >& Residual , T sor ,
		const std::vector< int >& entriesPerSlice , const std::vector< int >& slicesPerThread , int sliceDependence , bool unitDiagonal , bool strippedDiagonal
		) const { SolveGaussSeidelAndResidualParallel< T2 , T3 , UseSOR >( Diagonal+0 , b+0 , iters , Solution+0 , Residual+0 , sor , entriesPerSlice , slicesPerThread , sliceDependence , unitDiagonal , strippedDiagonal ); }
};

#include "SparseMatrixInterface.inl"

#endif // SPARSE_MATRIX_INTERFACE_INCLUDED
