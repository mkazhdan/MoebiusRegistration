/*
Copyright (c) 2006, Michael Kazhdan and Matthew Bolitho
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

#ifndef __SPARSEMATRIX_HPP
#define __SPARSEMATRIX_HPP

#include "Vector.h"
#include "SparseMatrixInterface.h"
#include "Array.h"

template <class T>
struct MatrixEntry2
{
	MatrixEntry2( void )	{ inN = outN = -1; value = 0; }
	int inN , outN;
	T value;
};
template< class T , class IndexType > class SparseMatrix : public SparseMatrixInterface< T , ConstPointer( MatrixEntry< T ,IndexType > ) >
{
	template< class T2 , class IndexType2 > friend class SparseMatrix;
	bool _contiguousMemory;
public:
	typedef SparseMatrixInterface< T , ConstPointer( MatrixEntry< T , IndexType > ) > Interface;

	int rows;
	Pointer( int ) rowSizes;
	Pointer( Pointer( MatrixEntry< T , IndexType > ) ) m_ppElements;

	SparseMatrix( );
	SparseMatrix( const SparseMatrix& M );
	template< class T2 , class IndexType2 >
	SparseMatrix( const SparseMatrix< T2 , IndexType2 >& M );
	~SparseMatrix();
	SparseMatrix< T , IndexType >& operator = (const SparseMatrix< T , IndexType >& M);
	template< class T2 , class IndexType2 >
	SparseMatrix< T , IndexType >& operator = ( const SparseMatrix< T2 , IndexType2 >& M );

	template< class T2 > Vector< T2 > operator * ( const Vector< T2 >& V ) const;

	template< class T2 , class IndexType2 >
	SparseMatrix< T , IndexType >& copy( const SparseMatrix< T2 , IndexType2 >& M , bool localIndex );

	inline ConstPointer( MatrixEntry< T , IndexType > ) begin( int row )    const { return m_ppElements[row]; }
	inline ConstPointer( MatrixEntry< T , IndexType > ) end  ( int row )    const { return m_ppElements[row]+rowSizes[row]; }
	inline size_t Rows                              ( void )       const { return rows; }
	inline size_t RowSize                           ( size_t idx ) const { return rowSizes[idx]; }

	SparseMatrix( int rows );
	void resize	( int rows );
	void SetRowSize( int row , int count );
	inline      Pointer( MatrixEntry< T , IndexType > ) operator[] ( int idx )       { return m_ppElements[idx]; }
	inline ConstPointer( MatrixEntry< T , IndexType > ) operator[] ( int idx ) const { return m_ppElements[idx]; }

	bool isContiguous( void ) const;
	void MakeContiguous( void );
	void MakeLocal( void );

	double SquareNorm(void) const;
	double ASymmetricSquareNorm( void ) const;
	double AHermitianSquareNorm( void ) const;

    /** Sets the column index of all allocated entries to -1 so that they are
     *  treated as non-existent. This is needed because SetRowSize() uses
     *  malloc instead of new and MatrixEntry's constructor is never called. */
    void invalidateEntries();

    /** Adds a scalar value to an element in the Matrix, using a new element if
     *  necessary. If no pre-allocated space for a new element exists, false is
     *  returned.
     *  WARNING: no check is done to remove entries that become zero after
     *  addition */
    bool addScalarToEntry( T s , IndexType i , IndexType j );
};
template< class T1 , class T2 , class T3 , class IndexType >
bool TransposeMultiply( const SparseMatrix< T1 , IndexType >& At , const SparseMatrix< T2 , IndexType >& B , SparseMatrix< T3 , IndexType >& out , T1 (*TransposeFunction)( const T1& )=NULL );
template< class T1 , class T2 , class T3 , class IndexType >
bool TransposeMultiply( const SparseMatrix< T1 , IndexType >& At , const SparseMatrix< T2 , IndexType >& B , SparseMatrix< T3 , IndexType >& out , int outRows , T1 (*TransposeFunction)( const T1& )=NULL );
template< class A_T , class A_const_iterator , class B_T , class B_const_iterator , class Out_T , class Out_IndexType >
bool Multiply( const SparseMatrixInterface< A_T , A_const_iterator >& A , const SparseMatrixInterface< B_T , B_const_iterator >& B , SparseMatrix< Out_T , Out_IndexType >& out , int threads = 1 );
template< class T , class In_const_iterator , class Out_IndexType >
bool Transpose( const SparseMatrixInterface< T , In_const_iterator >& At , SparseMatrix< T , Out_IndexType >& A ,               T (*TransposeFunction)( const T& )=NULL );
template< class T , class In_const_iterator , class Out_IndexType >
bool Transpose( const SparseMatrixInterface< T , In_const_iterator >& At , SparseMatrix< T , Out_IndexType >& A , int outRows , T (*TransposeFunction)( const T& )=NULL );

template< class T , class IndexType >
class SparseSymmetricMatrix : public SparseMatrix< T , IndexType >
{
public:
	template<class T2>
	Vector<T2> operator * (const Vector<T2>& V) const;
	template<class T2>
	void Multiply			( const Vector<T2>& In, Vector<T2>& Out ) const;
};
// In order to be supported:
// 1] The Matrix class has to support the Multiply(Vector<Data>,Vector<Data>) method.
// 2] The data class has to support Data*Data returning a value castable to a double

template <class Matrix>
class PivotSymmetricMatrix
{
public:
	std::vector<std::pair<Matrix,std::pair<int,int> > > rightMatrices;
	std::pair<Matrix,std::pair<int,int> > pivot;

	template <class Data>
	Vector<Data> operator * (const Vector<Data>& in) const;

	template <class Data>
	void Multiply(const Vector<Data>& in,Vector<Data>& out) const;

	void setPivot(const Matrix& M,int inDim,int outDim);
	void push(const Matrix& M,int inDim,int outDim);
	void pop(void);
};
template <class Matrix,class Data>
static int SolveConjugateGradient(const Matrix& SPD,const Vector<Data>& b,const int& iters,Vector<Data>& solution,const double eps=1e-8);
template <class Matrix,class IPS,class Real>
static int SolveConjugateGradient(const Matrix& SPD,const Vector<IPS>& b,const int& iters,Vector<IPS>& solution,const double eps=1e-8);
template <class Matrix,class IPS,class Real>
static int SolveConjugateGradient2(const Matrix& SPD,const Vector<IPS>& b,const int& iters,Vector<IPS>& solution,const double eps=1e-8);
#include "SparseMatrix.inl"

#endif /* __SPARSEMATRIX_HPP */
