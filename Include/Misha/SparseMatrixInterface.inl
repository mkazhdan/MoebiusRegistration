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

template< class T , class const_iterator > size_t SparseMatrixInterface< T , const_iterator >::Entries( void ) const
{
	size_t entries = 0;
	for( size_t i=0 ; i<Rows() ; i++ ) entries += RowSize( i );
	return entries;
}
template< class T , class const_iterator > double SparseMatrixInterface< T , const_iterator >::SquareNorm( void ) const
{
	double n=0;
	for( size_t i=0 ; i<Rows() ; i++ )
	{
		const_iterator e = end( i );
		for( const_iterator iter = begin( i ) ; iter!=e ; iter++ ) n += iter->Value * iter->Value;
	}
	return n;

}
template< class T , class const_iterator > double SparseMatrixInterface< T , const_iterator >::SquareASymmetricNorm( void ) const
{
	double n=0;
	for( size_t i=0 ; i<Rows() ; i++ )
	{
		const_iterator e = end( i );
		for( const_iterator iter1 = begin( i ) ; iter1!=e ; iter1++ )
		{
			int j = iter1->N;
			const_iterator e = end( j );
			double value = 0;
			for( const_iterator iter2 = begin( j ) ; iter2!=e ; iter2++ )
			{
				int k = iter2->N;
				if( k==i ) value += iter2->Value;
			}
			n += (iter1->Value-value) * (iter1->Value-value);
		}
	}
	return n;
}
template< class T , class const_iterator > double SparseMatrixInterface< T , const_iterator >::SquareASymmetricNorm( int& idx1 , int& idx2 ) const
{
	double n=0;
	double max=0;
	for( size_t i=0 ; i<Rows() ; i++ )
	{
		const_iterator e = end( i );
		for( const_iterator iter = begin( i ) ; iter!=e ; iter++ )
		{
			int j = iter->N;
			const_iterator e = end( j );
			double value = 0;
			for( const_iterator iter2 = begin( j ) ; iter2!=e ; iter2++ )
			{
				int k = iter2->N;
				if( k==i ) value += iter2->Value;
			}
			double temp = (iter->Value-value) * (iter->Value-value);
			n += temp;
			if( temp>=max ) idx1 = i , idx2 = j , max=temp;
		}
	}
	return n;
}
template< class T , class const_iterator >
template< class T2>
void SparseMatrixInterface< T , const_iterator >::Multiply( ConstPointer( T2 ) In , Pointer( T2 ) Out ) const
{
	ConstPointer( T2 ) in = In;
	for( size_t i=0 ; i<Rows() ; i++ )
	{
		T2 temp = Out[i];
		ConstPointer( T2 ) _in = in;
		const_iterator e = end( i );
		for( const_iterator iter = begin( i ) ; iter!=e ; iter++ ) temp += _in[ iter->N ] * iter->Value;
		Out[i] = temp;
	}
}
template< class T , class const_iterator >
template< class T2 >
void SparseMatrixInterface< T , const_iterator >::MultiplyParallel( ConstPointer( T2 ) In , Pointer( T2 ) Out , int threads , int multiplyFlag ) const
{
	ConstPointer( T2 ) in = In;
#pragma omp parallel for num_threads( threads )
	for( int i=0 ; i<(int)Rows() ; i++ )
	{
#if 1
		T2 temp;
		memset( &temp , 0 , sizeof(T2) );
		ConstPointer( T2 ) _in = in;
		const_iterator e = end( i );
		for( const_iterator iter = begin( i ) ; iter!=e ; iter++ ) temp += _in[ iter->N ] * iter->Value;
		switch( multiplyFlag )
		{
		case MULTIPLY_CLEAR:    Out[i]  = temp ; break;
		case MULTIPLY_ADD:      Out[i] += temp ; break;
		case MULTIPLY_SUBTRACT: Out[i] -= temp ; break;
		}
#else
		T2 temp;
		switch( multiplyFlag )
		{
//		case MULTIPLY_CLEAR:    temp *= 0       ; break;
		case MULTIPLY_CLEAR:    memset( &temp , 0 , sizeof(T2) ) ; break;
		case MULTIPLY_ADD:      temp  =  Out[i] ; break;
		case MULTIPLY_SUBTRACT: temp  = -Out[i] ; break;
		}
		ConstPointer( T2 ) _in = in;
		const_iterator e = end( i );
		for( const_iterator iter = begin( i ) ; iter!=e ; iter++ ) temp += _in[ iter->N ] * iter->Value;
		if( multiplyFlag==MULTIPLY_SUBTRACT ) Out[i] = -temp;
		else                                  Out[i] =  temp;
#endif
	}
}
#if NEW_TRANSPOSE_MULTIPLY
template< class T , class const_iterator >
template< class T2 >
void SparseMatrixInterface< T , const_iterator >::MultiplyTranspose( ConstPointer( T2 ) In , Pointer( T2 ) Out , int outDim , int multiplyFlag , T (*TransposeFunction)( const T& ) ) const
{
	T2* out = &Out[0];

	if     ( multiplyFlag==MULTIPLY_CLEAR    ) for( int i=0 ; i<outDim ; i++ ) out[i] *= 0;
	else if( multiplyFlag==MULTIPLY_SUBTRACT ) for( int i=0 ; i<outDim ; i++ ) out[i] = -out[i];
	for( size_t i=0 ; i<Rows() ; i++ )
	{
		T2* _out = out;
		const_iterator e = end( i );
		if( TransposeFunction ) for( const_iterator iter = begin( i ) ; iter!=e ; iter++ ) _out[ iter->N ] += In[i] * TransposeFunction( iter->Value );
		else                    for( const_iterator iter = begin( i ) ; iter!=e ; iter++ ) _out[ iter->N ] += In[i] * iter->Value;
	}
	if( multiplyFlag==MULTIPLY_SUBTRACT ) for( int i=0 ; i<outDim ; i++ ) out[i] = -out[i];
}
#else // !NEW_TRANSPOSE_MULTIPLY
template< class T , class const_iterator >
template< class T2 >
void SparseMatrixInterface< T , const_iterator >::MultiplyTranspose( ConstPointer( T2 ) In , Pointer( T2 ) Out , T (*TransposeFunction)( const T& ) ) const
{
	T2* out = &Out[0];
	for( size_t i=0 ; i<Rows() ; i++ )
	{
		T2* _out = out;
		const_iterator e = end( i );
		if( TransposeFunction ) for( const_iterator iter = begin( i ) ; iter!=e ; iter++ ) _out[ iter->N ] += In[i] * TransposeFunction( iter->Value );
		else                    for( const_iterator iter = begin( i ) ; iter!=e ; iter++ ) _out[ iter->N ] += In[i] * iter->Value;
	}
}
#endif // NEW_TRANSPOSE_MULTIPLY
template< class T , class const_iterator >
template< class T2 >
void SparseMatrixInterface< T , const_iterator >::BMinusMX( ConstPointer( T2 ) B , ConstPointer( T2 ) X , Pointer( T2 ) Out ) const
{
	ConstPointer( T2 ) x = X;
	for( size_t i=0 ; i<Rows() ; i++ )
	{
		T2 temp;
		temp *= 0;
		ConstPointer( T2 ) _x = x;
		const_iterator e = end( i );
		for( const_iterator iter = begin( i ) ; iter!=e ; iter++ ) temp += _x[ iter->N ] * iter->Value;
		Out[i] = B[i] - temp;
	}

}
template< class T , class const_iterator >
template< class T2 >
int SparseMatrixInterface< T , const_iterator >::SolveConjugateGradient( const Vector< T2 >& b , Vector< T2 >& x , double (*SquareNorm)( T2 ) , T2 (*Dot)( T2 , T2 ) , const size_t& iters , int threads ) const
{
	double eps=1e-16;
	Vector< T2 > r;
	r = b;
	MultiplyParallel( x , r , threads , MULTIPLY_SUBTRACT );
	Vector< T2 > q , d = r;
	q.resize( d.size() );
	double delta_new=0 , delta_0;
	for( int i=0 ; i<r.size() ; i++ ) delta_new += SquareNorm( r[i] );
	delta_0 = delta_new;
	if( delta_new<eps ) return 0;
	int ii;
	for( ii=0; ii<iters && delta_new>eps*delta_0 ; ii++ )
	{
		MultiplyParallel( d , q , threads , MULTIPLY_CLEAR );
        T2 dDotQ = 0 , alpha = 0;
		for( int i=0 ; i<d.size() ; i++ ) dDotQ += Dot( d[i] , q[i] );
		alpha = delta_new / dDotQ;
#pragma omp parallel for num_threads( threads )
		for( int i=0 ; i<r.size() ; i++ )
			x[i] += d[i] * alpha;
		if( !(ii%50) )
		{
			r = b;
			MultiplyParallel( x , r , threads , MULTIPLY_SUBTRACT );
		}
		else
		{
#pragma omp parallel for num_threads( threads )
			for( int i=0 ; i<r.size() ; i++ )
				r[i] = r[i] - Dot( q[i] , alpha );
		}

		double delta_old = delta_new , beta;
		delta_new = 0;
		for( int i=0 ; i<r.size() ; i++ ) delta_new += SquareNorm( r[i] );
		beta = delta_new / delta_old;
#pragma omp parallel for num_threads( threads )
		for( int i=0 ; i<d.size() ; i++ )
			d[i] = r[i] + d[i] * beta;
	}
	return ii;
}


template< class T , class const_iterator >
template< class T2 , class T3 , bool UnitDiagonal , bool StrippedDiagonal , bool UseSOR >
void SparseMatrixInterface< T , const_iterator >::SolveGaussSeidel( ConstPointer( T3 ) diagonal , ConstPointer( T2 ) b , int iters , Pointer( T2 ) solution , T sor ) const
{
	for( int i=0 ; i<iters ; i++ )
		if( UnitDiagonal )
			for( size_t j=0 ; j<Rows() ; j++ )
			{
				T2 temp;
				temp *= 0;
				Pointer( T2 ) _solution = solution;
				const_iterator e = end( j );
				for( const_iterator iter = begin( j ) ; iter!=e ; iter++ ) temp += _solution[ iter->N ] * iter->Value;
				if( UseSOR )
					if(StrippedDiagonal) solution[j]  = solution[j]*( T(1.0)-sor ) + (b[j]-temp) * sor;
					else                 solution[j] +=                              (b[j]-temp) * sor;
				else
					if(StrippedDiagonal) solution[j]  = (b[j]-temp);
					else                 solution[j] += (b[j]-temp);
			}
		else
			for( size_t j=0 ; j<Rows() ; j++ )
			{
				T2 temp;
				temp *= 0;
				Pointer( T2 ) _solution = solution;
				T dValue = T(1.) / diagonal[j];
				const_iterator e = end( j );
				for( const_iterator iter = begin( j ) ; iter!=e ; iter++ ) temp += _solution[ iter->N ] * iter->Value;
				if( UseSOR )
					if(StrippedDiagonal) solution[j]  = solution[j]*( T(1.0)-sor ) + (b[j]-temp) * dValue * sor;
					else                 solution[j] +=                              (b[j]-temp) * dValue * sor;
				else
					if(StrippedDiagonal) solution[j]  = (b[j]-temp) * dValue;
					else                 solution[j] += (b[j]-temp) * dValue;
			}
}

template< class T , class const_iterator >
template< class T2 , class T3 , bool UnitDiagonal , bool StrippedDiagonal , bool UseSOR >
void SparseMatrixInterface< T , const_iterator >::SolveGaussSeidelAndResidual( ConstPointer( T3 ) diagonal , ConstPointer( T2 ) b , int iters , Pointer( T2 ) solution , Pointer( T2 ) residual , T sor ) const
{
	for( int i=0 ; i<iters ; i++ )
		if( UnitDiagonal )
			for( size_t j=0 ; j<Rows() ; j++ )
			{
				T2 temp;
				temp *= 0;
				Pointer( T2 ) _solution = solution;
				const_iterator e = end( j );
				for( const_iterator iter = begin( j ) ; iter!=e ; iter++ ) temp += _solution[ iter->N ] * iter->Value;
				if( UseSOR )
					if( StrippedDiagonal ) solution[j]  = solution[j]*T(1.0-sor) + (b[j]-temp) * sor;
					else                   solution[j] +=                          (b[j]-temp) * sor;
				else
					if( StrippedDiagonal ) solution[j]  = (b[j]-temp);
					else                   solution[j] += (b[j]-temp);
			}
		else
			for( size_t j=0 ; j<Rows() ; j++ )
			{
				T2 temp;
				temp *= 0;
				Pointer( T2 ) _solution = solution;
				T dValue = T(1.) / diagonal[j];
				const_iterator e = end( j );
				for( const_iterator iter = begin( j ) ; iter!=e ; iter++ ) temp += _solution[ iter->N ] * iter->Value;
				if( UseSOR )
					if( StrippedDiagonal ) solution[j]  = solution[j]*T(1.0-sor) + (b[j]-temp) * dValue * sor;
					else                   solution[j] +=                          (b[j]-temp) * dValue * sor;
				else
					if( StrippedDiagonal ) solution[j]  = (b[j]-temp) * dValue;
					else                   solution[j] += (b[j]-temp) * dValue;
			}
	if( UnitDiagonal )
		for( size_t j=0 ; j<Rows() ; j++ )
		{
			T2 temp;
			temp *= 0;
			Pointer( T2 ) _solution = solution;
			const_iterator e = end( j );
			for( const_iterator iter = begin( j ) ; iter!=e ; iter++ ) temp += _solution[ iter->N ] * iter->Value;
			if( StrippedDiagonal ) residual[j] = b[j] - temp - solution[j];
			else                   residual[j] = b[j] - temp;
		}
	else
		for( size_t j=0 ; j<Rows() ; j++ )
		{
			T2 temp;
			temp *= 0;
			Pointer( T2 ) _solution = solution;
			T dValue = diagonal[j];
			const_iterator e = end( j );
			for( const_iterator iter = begin( j ) ; iter!=e ; iter++ ) temp += _solution[ iter->N ] * iter->Value;
			if( StrippedDiagonal ) residual[j] = b[j] - temp - ( solution[j] * dValue );
			else                   residual[j] = b[j] - temp;
		}
}

template< class T , class const_iterator >
template< class T2 >
void SparseMatrixInterface< T , const_iterator >::MultiplyParallel( 
	ConstPointer( T2 ) in , Pointer( T2 ) out ,
	const std::vector< int >& entriesPerSlice , const std::vector< int >& slicesPerThread ) const
{
	ConstPointer( T2 ) _in = &in[0];
	int threads = slicesPerThread.size();
	if( threads<1 ) return Multiply( in , out );
	std::vector< int > startEntry;

	startEntry.resize( entriesPerSlice.size() );
	startEntry[0] = 0;
	for( int i=1 ; i<entriesPerSlice.size() ; i++ ) startEntry[i] = startEntry[i-1] + entriesPerSlice[i-1];
#pragma omp parallel for num_threads( threads )
	for( int s=0 ; s<threads ; s++ )
	{
		int startSlice = 0 , endSlice;
		for( int i=0 ; i<s ; i++ ) startSlice += slicesPerThread[i];
		endSlice = startSlice + slicesPerThread[s];
		for( int i=startSlice ; i<endSlice ; i++ )
		{
			int start = startEntry[i];
			int stop  = startEntry[i] + entriesPerSlice[i];
			for( size_t j=start ; j<stop ; j++ )
			{
				T2 temp;
				temp *= 0;
				ConstPointer( T2 ) in = _in;
				const_iterator e = end( j );
				for( const_iterator iter = begin( j ) ; iter!=e ; iter++ ) temp += in[ iter->N ] * iter->Value;
				out[j] += temp;
			}
		}
	}
}
template< class T , class const_iterator >
template< class T2 >
void SparseMatrixInterface< T , const_iterator >::MultiplyParallel( 
	ConstPointer( T2 ) in , Pointer( T2 ) out ,
	const std::vector< int >& entriesPerLine , int linesPerSlice , const std::vector< int >& slicesPerThread ) const
{
	ConstPointer( T2 ) _in = in;
	int threads = slicesPerThread.size();
	if( threads<1 ) return Multiply( in , out );
	std::vector< int > startEntry;

	startEntry.resize( entriesPerLine.size() );
	startEntry[0] = 0;
	for( int i=1 ; i<entriesPerLine.size() ; i++ ) startEntry[i] = startEntry[i-1] + entriesPerLine[i-1];
#pragma omp parallel for num_threads( threads )
	for( int s=0 ; s<threads ; s++ )
	{
		int startSlice = 0 , endSlice;
		for( int i=0 ; i<s ; i++ ) startSlice += slicesPerThread[i];
		endSlice = startSlice + slicesPerThread[s];
		for( int i=startSlice*linesPerSlice ; i<endSlice*linesPerSlice ; i++ )
		{
			int start = startEntry[i];
			int stop  = startEntry[i] + entriesPerLine[i];
			for( size_t j=start ; j<stop ; j++ )
			{
				T2 temp;
				temp *= 0;
				ConstPointer( T2 ) in = _in;
				const_iterator e = end( j );
				for( const_iterator iter = begin( j ) ; iter!=e ; iter++ ) temp += in[ iter->N ] * iter->Value;
				out[j] += temp;
			}
		}
	}
}

template< class T , class const_iterator >
template< class T2 >
void SparseMatrixInterface< T , const_iterator >::MultiplyTransposeParallel( 
	ConstPointer( T2 ) in , Pointer( T2 ) out ,
	const std::vector< int >& entriesPerSlice , const std::vector< int >& slicesPerThread , int sliceDependence ) const
{
	T2* _out = &out[0];
	int threads = slicesPerThread.size();
	if( threads<1 ) return MultiplyTranspose( in , out );
	std::vector< int > startEntry;

	startEntry.resize( entriesPerSlice.size() );
	startEntry[0] = 0;
	for( int i=1 ; i<entriesPerSlice.size() ; i++ ) startEntry[i] = startEntry[i-1] + entriesPerSlice[i-1];

#pragma omp parallel for num_threads( threads )
	for( int s=0 ; s<threads ; s++ )
	{
		int startSlice = 0 , endSlice;
		for( int i=0 ; i<s ; i++ ) startSlice += slicesPerThread[i];
		endSlice = startSlice + slicesPerThread[s];
		startSlice += sliceDependence;
		for( int i=startSlice ; i<endSlice ; i++ )
		{
			int start = startEntry[i];
			int stop = startEntry[i] + entriesPerSlice[i];
			for( int j=start ; j<stop ; j++ )
			{
				T2* out = _out;
				const_iterator e = end( j );
				for( const_iterator iter = begin( j ) ; iter!=e ; iter++ ) out[ iter->N ] += in[ j ] * iter->Value;
			}
		}
	}
	for( int s=0 ; s<threads ; s++ )
	{
		int startSlice = 0 , endSlice;
		for( int i=0 ; i<s ; i++ ) startSlice += slicesPerThread[i];
		endSlice = startSlice + sliceDependence;
		for( int i=startSlice ; i<endSlice ; i++ )
		{
			size_t start = startEntry[i];
			size_t stop = startEntry[i] + entriesPerSlice[i];
			for( int j=start ; j<stop ; j++ )
			{
				T2* out = _out;
				const_iterator e = end( j );
				for( const_iterator iter = begin( j ) ; iter!=e ; iter++ ) out[ iter->N ] += in[ j ] * iter->Value;
			}
		}
	}
}
template< class T , class const_iterator >
template< class T2 >
void SparseMatrixInterface< T , const_iterator >::MultiplyTransposeParallel( 
	ConstPointer( T2 ) in , Pointer( T2 ) out ,
	const std::vector< int >& entriesPerLine , int linesPerSlice , const std::vector< int >& slicesPerThread , int sliceDependence ) const
{
	T2* _out = &out[0];
	int threads = slicesPerThread.size();
	if( threads<1 ) return MultiplyTranspose( in , out );
	std::vector< int > startEntry;

	startEntry.resize( entriesPerLine.size() );
	startEntry[0] = 0;
	for( int i=1 ; i<entriesPerLine.size() ; i++ ) startEntry[i] = startEntry[i-1] + entriesPerLine[i-1];

#pragma omp parallel for num_threads( threads )
	for( int s=0 ; s<threads ; s++ )
	{
		int startSlice = 0 , endSlice;
		for( int i=0 ; i<s ; i++ ) startSlice += slicesPerThread[i];
		endSlice = startSlice + slicesPerThread[s];
		startSlice += sliceDependence;
		for( int i=startSlice*linesPerSlice ; i<endSlice*linesPerSlice ; i++ )
		{
			int start = startEntry[i];
			int stop = startEntry[i] + entriesPerLine[i];
			for( int j=start ; j<stop ; j++ )
			{
				T2* out = _out;
				const_iterator e = end( j );
				for( const_iterator iter = begin( j ) ; iter!=e ; iter++ ) out[ iter->N ] += in[ j ] * iter->Value;
			}
		}
	}
	for( int s=0 ; s<threads ; s++ )
	{
		int startSlice = 0 , endSlice;
		for( int i=0 ; i<s ; i++ ) startSlice += slicesPerThread[i];
		endSlice = startSlice + sliceDependence;
		for( int i=startSlice*linesPerSlice ; i<endSlice*linesPerSlice ; i++ )
		{
			size_t start = startEntry[i];
			size_t stop = startEntry[i] + entriesPerLine[i];
			for( int j=start ; j<stop ; j++ )
			{
				T2* out = _out;
				const_iterator e = end( j );
				for( const_iterator iter = begin( j ) ; iter!=e ; iter++ ) out[ iter->N ] += in[ j ] * iter->Value;
			}
		}
	}
}

template< class T , class const_iterator >
template< class T2 >
void SparseMatrixInterface< T , const_iterator >::BMinusMXParallel( 
	ConstPointer( T2 ) B , ConstPointer( T2 ) X , Pointer( T2 ) out ,
	const std::vector< int >& entriesPerSlice , const std::vector< int >& slicesPerThread ) const
{
	ConstPointer( T2 ) x = &x[0];
	int threads = slicesPerThread.size();
	if( threads<1 ) return BMinusMX( B , X , out );
	std::vector< int > startEntry;

	startEntry.resize( entriesPerSlice.size() );
	startEntry[0] = 0;
	for( int i=1 ; i<entriesPerSlice.size() ; i++ ) startEntry[i] = startEntry[i-1] + entriesPerSlice[i-1];

#pragma omp parallel for num_threads( threads )
	for( int s=0 ; s<threads ; s++ )
	{
		int startSlice = 0 , endSlice;
		for( int i=0 ; i<s ; i++ ) startSlice += slicesPerThread[i];
		endSlice = startSlice + slicesPerThread[s];
		for( int i=startSlice ; i<endSlice ; i++ )
		{
			int start = startEntry[i];
			int stop = startEntry[i] + entriesPerSlice[i];
			for( size_t j=start ; j<stop ; j++ )
			{
				T2 temp;
				temp *= 0;
				ConstPointer( T2 ) _x = x;
				const_iterator e = end( j );
				for( const_iterator iter = begin( j ) ; iter!=e ; iter++ ) temp += _x[ iter->N ] * iter->Value;
				out[j] += B[j] - temp;
			}
		}
	}
}
// Implement with safe-zone approach from Ofir et al. (ToG '08). WARNING!!! To do this right, we would need to separate the memory spaces for the solution.
template< class T2 >
ParallelSolution< T2 >::ParallelSolution( int iters , const std::vector< int >& entriesPerSlice , const std::vector< int >& slicesPerThread , int sliceDependence )
{
	_threads = slicesPerThread.size();
	if( _threads<1 ) fprintf( stderr , "Could not allocate ParallelSolution with no threads: %d\n" , _threads ) , exit( 0 );

	_iters = iters;
#if OWN_MEMORY
	__entriesPerSlice = entriesPerSlice;
	__slicesPerThread = slicesPerThread;
	_entriesPerSlice = &__entriesPerSlice[0];
	_slicesPerThread = &__slicesPerThread[0];
#else // !OWN_MEMORY
	_entriesPerSlice = &entriesPerSlice[0];
	_slicesPerThread = &slicesPerThread[0];
#endif // OWN_MEMORY
	_sliceDependence = sliceDependence;
	_sliceCount = entriesPerSlice.size();

	_pSolution.resize( _threads );
	_startEntry.resize( _sliceCount );
	_startEntry[0] = 0;
	for( int i=1 ; i<_sliceCount ; i++ ) _startEntry[i] = _startEntry[i-1] + entriesPerSlice[i-1];

	_sliceBounds.resize( _threads );
	_skirtBounds.resize( _threads );

	for( int t=0 ; t<_threads ; t++ )
	{
		// Compute the starting and ending slices associated with the current thread.
		int startSlice = 0 , endSlice;
		_sliceBounds[t].first = 0;
		for( int i=0 ; i<t ; i++ ) _sliceBounds[t].first += slicesPerThread[i];
		_sliceBounds[t].second = _sliceBounds[t].first + slicesPerThread[t];

		_skirtBounds[t].first  = _sliceBounds[t].first  - (iters+1) * sliceDependence;
		_skirtBounds[t].second = _sliceBounds[t].second + (iters+1) * sliceDependence;
		if( _skirtBounds[t].first<0 ) _skirtBounds[t].first = 0;
		if( _skirtBounds[t].second>entriesPerSlice.size() ) _skirtBounds[t].second = entriesPerSlice.size( );

		int skirtEntryStart = _startEntry[_skirtBounds[t].first];
		int      entryStart = _startEntry[_sliceBounds[t].first];
		int      entryEnd   = _startEntry[_sliceBounds[t].second-1] + _entriesPerSlice[_sliceBounds[t].second-1];
		int skirtEntryEnd   = _startEntry[_skirtBounds[t].second-1] + _entriesPerSlice[_skirtBounds[t].second-1];
		_pSolution[t] = ( T2* ) malloc( sizeof( T2 ) * ( skirtEntryEnd - skirtEntryStart ) );
		_pSolution[t] += entryStart-skirtEntryStart;
	}
}
template< class T2 >
ParallelSolution< T2 >::~ParallelSolution( void )
{
	for( int t=0 ; t<_threads ; t++ )
	{
		int skirtEntryStart = _startEntry[_skirtBounds[t].first];
		int      entryStart = _startEntry[_sliceBounds[t].first];
		int      entryEnd   = _startEntry[_sliceBounds[t].second-1] + _entriesPerSlice[_sliceBounds[t].second-1];
		int skirtEntryEnd   = _startEntry[_skirtBounds[t].second-1] + _entriesPerSlice[_skirtBounds[t].second-1];
		_pSolution[t] -= entryStart-skirtEntryStart;
		free( _pSolution[t] );
	}
	_pSolution.clear();
}
template< class T2 >
void ParallelSolution< T2 >::SetFromArray( ConstPointer( T2 ) solution )
{
#pragma omp parallel for num_threads( _threads )
	for( int t=0 ; t<_threads ; t++ )
	{
		int skirtEntryStart = _startEntry[_skirtBounds[t].first];
		int      entryStart = _startEntry[_sliceBounds[t].first];
		int      entryEnd   = _startEntry[_sliceBounds[t].second-1] + _entriesPerSlice[_sliceBounds[t].second-1];
		int skirtEntryEnd   = _startEntry[_skirtBounds[t].second-1] + _entriesPerSlice[_skirtBounds[t].second-1];
		memcpy( _pSolution[t] - ( entryStart-skirtEntryStart ) , &solution[skirtEntryStart] , sizeof( T2 ) * ( skirtEntryEnd - skirtEntryStart ) );
	}
}
template< class T2 >
void ParallelSolution< T2 >::SetToArray( Pointer( T2 ) solution ) const
{
#pragma omp parallel for num_threads( _threads )
	for( int t=0 ; t<_threads ; t++ )
	{
		int skirtEntryStart = _startEntry[_skirtBounds[t].first];
		int      entryStart = _startEntry[_sliceBounds[t].first];
		int      entryEnd   = _startEntry[_sliceBounds[t].second-1] + _entriesPerSlice[_sliceBounds[t].second-1];
		int skirtEntryEnd   = _startEntry[_skirtBounds[t].second-1] + _entriesPerSlice[_skirtBounds[t].second-1];
		memcpy( solution + entryStart , _pSolution[t] , sizeof( T2 ) * ( entryEnd-entryStart ) );
	}
}
template< class T2 >
void ParallelSolution< T2 >::clear( void )
{
#pragma omp parallel for num_threads( _threads )
	for( int t=0 ; t<_threads ; t++ )
	{
		int skirtEntryStart = _startEntry[_skirtBounds[t].first];
		int      entryStart = _startEntry[_sliceBounds[t].first];
		int      entryEnd   = _startEntry[_sliceBounds[t].second-1] + _entriesPerSlice[_sliceBounds[t].second-1];
		int skirtEntryEnd   = _startEntry[_skirtBounds[t].second-1] + _entriesPerSlice[_skirtBounds[t].second-1];
		memset( _pSolution[t] - ( entryStart-skirtEntryStart ) , 0 , sizeof( T2 ) * ( skirtEntryEnd - skirtEntryStart ) );
	}
}
template< class T2 >
void ParallelSolution< T2 >::_clearSkirts( void )
{
#pragma omp parallel for num_threads( _threads )
	for( int t=0 ; t<_threads ; t++ )
	{
		int skirtEntryStart = _startEntry[_skirtBounds[t].first];
		int      entryStart = _startEntry[_sliceBounds[t].first];
		int      entryEnd   = _startEntry[_sliceBounds[t].second-1] + _entriesPerSlice[_sliceBounds[t].second-1];
		int skirtEntryEnd   = _startEntry[_skirtBounds[t].second-1] + _entriesPerSlice[_skirtBounds[t].second-1];
		if( entryStart - skirtEntryStart ) memset( _pSolution[t] - ( entryStart-skirtEntryStart ) , 0 , sizeof( T2 ) * ( entryStart - skirtEntryStart ) );
		if( skirtEntryEnd - entryEnd )     memset( _pSolution[t] + ( entryEnd - entryStart ) , 0 , sizeof( T2 ) * ( skirtEntryEnd - entryEnd ) );
	}
}
// WARNING!!! The synchronization  and merging code is not completely correct. If the skirts are larger than the neighbors, there will be trouble.
template< class T2 >
void ParallelSolution< T2 >::_merge( void )
{
#pragma omp parallel for num_threads( _threads-1 )
	for( int t=0 ; t<_threads-1 ; t++ )
	{
		int skirtEntryStart1 = _startEntry[_skirtBounds[t].first];
		int      entryStart1 = _startEntry[_sliceBounds[t].first];
		int      entryEnd1   = _startEntry[_sliceBounds[t].second-1] + _entriesPerSlice[_sliceBounds[t].second-1];
		int skirtEntryEnd1   = _startEntry[_skirtBounds[t].second-1] + _entriesPerSlice[_skirtBounds[t].second-1];
		int skirtEntryStart2 = _startEntry[_skirtBounds[t+1].first];
		int      entryStart2 = _startEntry[_sliceBounds[t+1].first];
		int      entryEnd2   = _startEntry[_sliceBounds[t+1].second-1] + _entriesPerSlice[_sliceBounds[t+1].second-1];
		int skirtEntryEnd2   = _startEntry[_skirtBounds[t+1].second-1] + _entriesPerSlice[_skirtBounds[t+1].second-1];

		int mergeSize;
		T2 *left , *right;

		// Merge the left side of the boundary
		mergeSize = entryStart2 - skirtEntryStart2;
		left  = _pSolution[t  ] + ( entryEnd1-entryStart1 ) - ( entryStart2-skirtEntryStart2 );
		right = _pSolution[t+1] - ( entryStart2-skirtEntryStart2 );
		for( int i=0 ; i<mergeSize ; i++ ) left[i] = right[i] = ( left[i]+right[i] );

		// Merge the right side of the boundary
		mergeSize = skirtEntryEnd1-entryEnd1;
		left  = _pSolution[t  ] + ( entryEnd1-entryStart1 );
		right = _pSolution[t+1];
		for( int i=0 ; i<mergeSize ; i++ ) left[i] = right[i] = ( left[i]+right[i] );
	}
}

template< class T2 >
void ParallelSolution< T2 >::synchronize( void )
{
#pragma omp parallel for num_threads( _threads-1 )
	for( int t=0 ; t<_threads-1 ; t++ )
	{
		int skirtEntryStart1 = _startEntry[_skirtBounds[t].first];
		int      entryStart1 = _startEntry[_sliceBounds[t].first];
		int      entryEnd1   = _startEntry[_sliceBounds[t].second-1] + _entriesPerSlice[_sliceBounds[t].second-1];
		int skirtEntryEnd1   = _startEntry[_skirtBounds[t].second-1] + _entriesPerSlice[_skirtBounds[t].second-1];
		int skirtEntryStart2 = _startEntry[_skirtBounds[t+1].first];
		int      entryStart2 = _startEntry[_sliceBounds[t+1].first];
		int      entryEnd2   = _startEntry[_sliceBounds[t+1].second-1] + _entriesPerSlice[_sliceBounds[t+1].second-1];
		int skirtEntryEnd2   = _startEntry[_skirtBounds[t+1].second-1] + _entriesPerSlice[_skirtBounds[t+1].second-1];

		memcpy( _pSolution[t+1] - ( entryStart2-skirtEntryStart2 ) , _pSolution[t] + ( entryEnd1-entryStart1 ) - ( entryStart2-skirtEntryStart2 ) , sizeof( T2 ) * ( entryStart2-skirtEntryStart2 ) );
		memcpy( _pSolution[t] + ( entryEnd1-entryStart1 ) , _pSolution[t+1] , sizeof( T2 ) * ( skirtEntryEnd1-entryEnd1 ) );
	}
}
template< class T , class const_iterator >
template< class T2 >
void SparseMatrixInterface< T , const_iterator >::MultiplyParallel( ConstPointer( T2 ) in , ParallelSolution< T2 >* out ) const
{
	const int threads = out->_threads;
	const int sliceDependence = out->_sliceDependence;
	const int sliceCount = out->_sliceCount;
	const int* startEntry = &out->_startEntry[0];
	const int* entriesPerSlice = out->_entriesPerSlice;
#pragma omp parallel for num_threads( threads )
	for( int t=0 ; t<threads ; t++ )
	{
		// Compute the starting and ending slices associated with the current thread.
		int outStartSlice = out->_sliceBounds[t].first;
		int   outEndSlice = out->_sliceBounds[t].second;
		int     outOffset = out->_startEntry[ outStartSlice ];

		T2* _out = out->_pSolution[t] - outOffset;

		for( int i=outStartSlice ; i<outEndSlice ; i++ )
		{
			int start = startEntry[i];
			int stop  = startEntry[i] + entriesPerSlice[i];
			for( size_t j=start ; j<stop ; j++ )
			{
				T2 temp;
				temp *= 0;
				ConstPointer( T2 ) _in = in;
				const_iterator e = end( j );
				for( const_iterator iter = begin( j ) ; iter!=e ; iter++ ) temp += _in[ iter->N ] * iter->Value;
				_out[j] += temp;
			}
		}
	}
}
// WARNING: In general, the next two function calls are not safe since they assume that in order to
// set out->_pSolution[t], we only need to know the values in in->_pSolution[t]!
template< class T , class const_iterator >
template< class T2 >
void SparseMatrixInterface< T , const_iterator >::MultiplyParallel( const ParallelSolution< T2 >* in , ParallelSolution< T2 >* out ) const
{
	const int threads = out->_threads;
	const int sliceDependence = out->_sliceDependence;
	const int sliceCount = out->_sliceCount;
	const int* startEntry = &out->_startEntry[0];
	const int* entriesPerSlice = out->_entriesPerSlice;
#pragma omp parallel for num_threads( threads )
	for( int t=0 ; t<threads ; t++ )
	{
		// Compute the starting and ending slices associated with the current thread.
		int  inStartSlice =  in->_sliceBounds[t].first;
		int    inEndSlice =  in->_sliceBounds[t].second;
		int outStartSlice = out->_sliceBounds[t].first;
		int   outEndSlice = out->_sliceBounds[t].second;
		int      inOffset =  in->_startEntry[  inStartSlice ];
		int     outOffset = out->_startEntry[ outStartSlice ];

		ConstPointer( T2 ) _in  =  in->_pSolution[t] -  inOffset;
		T2*       _out = out->_pSolution[t] - outOffset;

		for( int i=outStartSlice ; i<outEndSlice ; i++ )
		{
			int start = startEntry[i];
			int stop  = startEntry[i] + entriesPerSlice[i];
			for( int j=start ; j<stop ; j++ )
			{
				T2 temp;
				temp *= 0;
				ConstPointer( T2 ) in = _in;
				const_iterator e = end( j );
				for( const_iterator iter = begin( j ) ; iter!=e ; iter++ ) temp += in[ iter->N ] * iter->Value;
				_out[j] += temp;
			}
		}
	}
}
template< class T , class const_iterator >
template< class T2 >
void SparseMatrixInterface< T , const_iterator >::MultiplyTransposeParallel( const ParallelSolution< T2 >* in , ParallelSolution< T2 >* out ) const
{
	const int threads = in->_threads;
	const int sliceDependence = in->_sliceDependence;
	const int sliceCount = in->_sliceCount;
	const int* startEntry = &in->_startEntry[0];
	const int* entriesPerSlice = in->_entriesPerSlice;

	out->_clearSkirts( );
#pragma omp parallel for num_threads( threads )
	for( int t=0 ; t<threads ; t++ )
	{
		// Compute the starting and ending slices associated with the current thread.
		int  inStartSlice =  in->_sliceBounds[t].first;
		int    inEndSlice =  in->_sliceBounds[t].second;
		int outStartSlice = out->_sliceBounds[t].first;
		int   outEndSlice = out->_sliceBounds[t].second;
		int      inOffset =  in->_startEntry[  inStartSlice ];
		int     outOffset = out->_startEntry[ outStartSlice ];

		ConstPointer( T2 ) _in  =  in->_pSolution[t] -  inOffset;
		T2*       _out = out->_pSolution[t] - outOffset;

		for( int i=inStartSlice ; i<inEndSlice ; i++ )
		{
			int start = startEntry[i];
			int stop  = startEntry[i] + entriesPerSlice[i];
			for( int j=start ; j<stop ; j++ )
			{
				T2 temp = _in[ j ];
				T2* out = _out;
				const_iterator e = end( j );
				for( const_iterator iter = begin( j ) ; iter!=e ; iter++ ) out[ iter->N ] += temp * iter->Value;
			}
		}
	}
	out->_merge( );
}


template< class T , class const_iterator >
template< class T2 , class T3 , bool UnitDiagonal , bool StrippedDiagonal , bool UseSOR >
void SparseMatrixInterface< T , const_iterator >::SolveGaussSeidelParallel(	ConstPointer( T3 ) Diagonal , ConstPointer( T2 ) b , int iters , ParallelSolution< T2 >* Solution , T sor ) const
{
	const int threads = Solution->_threads;
	const int sliceDependence = Solution->_sliceDependence;
	const int sliceCount = Solution->_sliceCount;
	const int* startEntry = &Solution->_startEntry[0];
	const int* entriesPerSlice = Solution->_entriesPerSlice;
#pragma omp parallel for num_threads( threads )
	for( int t=0 ; t<threads ; t++ )
	{
		// Compute the starting and ending slices associated with the current thread.
		int startSlice = Solution->_sliceBounds[t].first;
		int   endSlice = Solution->_sliceBounds[t].second;
		int     offset = Solution->_startEntry[ startSlice ];

		Pointer( T2 ) solution = Solution->_pSolution[t] - offset;

		// When processing at index i, we will be relaxing on slices:
		// [ s , s + sliceDependence , ... , s + (iters-1)*sliceDependence ]
		for( int s=startSlice-(iters-1)*sliceDependence ; s<endSlice ; s++ )
		{

			// When processing at index j, we will be relaxing slice:
			// s + i*sliceDepdence for the (iters-i)-th time.
			for( int i=iters-1 ; i>=0 ; i-- ) 
			{
				int ss = s + i*sliceDependence;
				// Compute the starting and ending bounds of the skirt on which we are relaxing (keeping in mind that we will be doing a residual calcuation)
				int startBound = startSlice - i*sliceDependence;
				int   endBound =   endSlice + i*sliceDependence+1;
				if( startBound<0 ) startBound = 0;
				if(   endBound>sliceCount ) endBound = sliceCount;
				if( ss>=startBound && ss<endBound )
				{
					int start = startEntry[ss];
					int  stop = startEntry[ss] + entriesPerSlice[ss];
					for( int jj=start ; jj<stop ; jj++ )
					{
						T2 temp;
						temp *= 0;
						T2* _solution = solution;
						if( UnitDiagonal )
						{
							const_iterator e = end( jj );
							for( const_iterator iter = begin( jj ) ; iter!=e ; iter++ ) temp += _solution[ iter->N ] * iter->Value;
							if( UseSOR )
								if(StrippedDiagonal) solution[jj]  = solution[jj]*T(1.0-sor) + (b[jj]-temp) * sor;
								else                 solution[jj] +=                           (b[jj]-temp) * sor;
							else
								if(StrippedDiagonal) solution[jj]  = (b[jj]-temp);
								else                 solution[jj] += (b[jj]-temp);
						}
						else
						{
							T dValue = T(1.) / Diagonal[jj];
							const_iterator e = end( jj );
							for( const_iterator iter = begin( jj ) ; iter!=e ; iter++ ) temp += _solution[ iter->N ] * iter->Value;
							if( UseSOR )
								if(StrippedDiagonal) solution[jj]  = solution[jj]*T(1.0-sor) + (b[jj]-temp) * dValue * sor;
								else                 solution[jj] +=                           (b[jj]-temp) * dValue * sor;
							else
								if(StrippedDiagonal) solution[jj]  = (b[jj]-temp) * dValue;
								else                 solution[jj] += (b[jj]-temp) * dValue;
						}
					}
				}
			}
		}
	}
}
template< class T , class const_iterator >
template< class T2 , class T3 , bool UnitDiagonal , bool StrippedDiagonal , bool UseSOR >
void SparseMatrixInterface< T , const_iterator >::SolveGaussSeidelAndResidualParallel( ConstPointer( T3 ) Diagonal , ConstPointer( T2 ) b , int iters , ParallelSolution< T2 >* Solution , Pointer( T2 ) residual , T sor ) const
{
	const int threads = Solution->_threads;
	const int sliceDependence = Solution->_sliceDependence;
	const int sliceCount = Solution->_sliceCount;
	const int* startEntry = &Solution->_startEntry[0];
	const int* entriesPerSlice = Solution->_entriesPerSlice;
#pragma omp parallel for num_threads( threads )
	for( int t=0 ; t<threads ; t++ )
	{
		// Compute the starting and ending slices associated with the current thread.
		int startSlice = Solution->_sliceBounds[t].first;
		int   endSlice = Solution->_sliceBounds[t].second;
		int     offset = Solution->_startEntry[ startSlice ];

		Pointer( T2 ) solution = Solution->_pSolution[t] - offset;

		// When processing at index i, we will be relaxing on slices:
		// [ s , s + sliceDependence , ... , s + (iters-1)*sliceDependence ]

		for( int s=startSlice-iters*sliceDependence ; s<endSlice+sliceDependence ; s++ )
		{
			// When processing at index j, we will be relaxing slice:
			// s + i*sliceDepdence for the (iters-i)-th time.
			for( int i=iters-1 ; i>=0 ; i-- ) 
			{
				int ss = s + i*sliceDependence;
				// Compute the starting and ending bounds of the skirt on which we are relaxing (keeping in mind that we will be doing a residual calcuation)
				int startBound = startSlice - (i+1)*sliceDependence;
				int   endBound =   endSlice + (i+1)*sliceDependence+1;
				if( startBound<0 ) startBound = 0;
				if(   endBound>sliceCount ) endBound = sliceCount;
				if( ss>=startBound && ss<endBound )
				{
					int start = startEntry[ss];
					int  stop = startEntry[ss] + entriesPerSlice[ss];

					for( int jj=start ; jj<stop ; jj++ )
					{
						T2 temp;
						temp *= 0;
						T2* _solution = solution;
						if( UnitDiagonal )
						{
							const_iterator e = end( jj );
							for( const_iterator iter=begin( jj ) ; iter!=e ; iter++ ) temp += _solution[ iter->N ] * iter->Value;
							if( UseSOR )
								if(StrippedDiagonal) solution[jj]  = solution[jj]*T(1.0-sor) + (b[jj]-temp) * sor;
								else                 solution[jj] +=                           (b[jj]-temp) * sor;
							else
								if(StrippedDiagonal) solution[jj]  = (b[jj]-temp);
								else                 solution[jj] += (b[jj]-temp);
						}
						else
						{
							T dValue = T(1.) / Diagonal[jj];
							const_iterator e = end( jj );
							for( const_iterator iter=begin( jj ) ; iter!=e ; iter++ ) temp += _solution[ iter->N ] * iter->Value;
							if( UseSOR )
								if(StrippedDiagonal) solution[jj]  = solution[jj]*T(1.0-sor) + (b[jj]-temp) * dValue * sor;
								else                 solution[jj] +=                           (b[jj]-temp) * dValue * sor;
							else
								if(StrippedDiagonal) solution[jj]  = (b[jj]-temp) * dValue;
								else                 solution[jj] += (b[jj]-temp) * dValue;
						}
					}
				}
			}
			int ss = s-sliceDependence;
			if( ss>=startSlice && ss<endSlice )
			{
				int start = startEntry[ss];
				int  stop = startEntry[ss] + entriesPerSlice[ss];
				for( int jj=start ; jj<stop ; jj++ )
				{
					T2 temp;
					temp *= 0;
					T2* _solution = solution;
					if( UnitDiagonal )
					{
						const_iterator e = end( jj );
						for( const_iterator iter=begin( jj ) ; iter!=e ; iter++ ) temp += _solution[ iter->N ] * iter->Value;
						if( StrippedDiagonal ) residual[jj] = b[jj] - temp - solution[jj];
						else                   residual[jj] = b[jj] - temp;
					}
					else
					{
						T dValue = Diagonal[jj];
						const_iterator e = end( jj );
						for( const_iterator iter=begin( jj ) ; iter!=e ; iter++ ) temp += _solution[ iter->N ] * iter->Value;
						if( StrippedDiagonal ) residual[jj] = b[jj] - temp - ( solution[jj] * dValue );
						else                   residual[jj] = b[jj] - temp;
					}
				}
			}
		}
	}
}
template< class T2 >
ParallelSolution2< T2 >::ParallelSolution2( int iters , const std::vector< int >& entriesPerLine , int linesPerSlice , const std::vector< int >& slicesPerThread , int sliceDependence )
{
if( iters<1 ) iters = 5;
	_threads = slicesPerThread.size();
	if( _threads<1 ) fprintf( stderr , "Could not allocate ParallelSolution with no threads: %d\n" , _threads ) , exit( 0 );
	_iters = iters;
	_linesPerSlice   = linesPerSlice;
	_sliceDependence = sliceDependence;
#if OWN_MEMORY
	__entriesPerLine  = entriesPerLine;
	__slicesPerThread = slicesPerThread;
	_entriesPerLine  = &__entriesPerLine[0];
	_slicesPerThread = &__slicesPerThread[0];
#else // !OWN_MEMORY
	_entriesPerLine  = &entriesPerLine[0];
	_slicesPerThread = &slicesPerThread[0];
#endif // OWN_MEMORY
	_sliceCount = 0;
	for( int t=0 ; t<_threads ; t++ ) _sliceCount += slicesPerThread[t];
	if( _sliceCount != _linesPerSlice ) fprintf( stderr , "How come the slice count is different from the number lines per slice? %d != %d\n" , _sliceCount , _linesPerSlice );
	_pSolution.resize( _threads );
	_startEntry.resize( _sliceCount*_linesPerSlice );
	_startEntry[0] = 0;
	for( int i=1 ; i<entriesPerLine.size() ; i++ ) _startEntry[i] = _startEntry[i-1] + _entriesPerLine[i-1];
	_sliceBounds.resize( _threads );
	_skirtBounds.resize( _threads );

	for( int t=0 ; t<_threads ; t++ )
	{
		// Compute the starting and ending slices associated with the current thread.
		int startSlice = 0;
		_sliceBounds[t].first = 0;
		for( int i=0 ; i<t ; i++ ) _sliceBounds[t].first += slicesPerThread[i];
		_sliceBounds[t].second = _sliceBounds[t].first + slicesPerThread[t];

		_skirtBounds[t].first  = _sliceBounds[t].first  - (iters+1) * sliceDependence;
		_skirtBounds[t].second = _sliceBounds[t].second + (iters+1) * sliceDependence;
		if( _skirtBounds[t].first<0 ) _skirtBounds[t].first = 0;
		if( _skirtBounds[t].second>_sliceCount ) _skirtBounds[t].second = _sliceCount;

		int skirtEntryStart = _startEntry[_skirtBounds[t].first*_linesPerSlice];
		int      entryStart = _startEntry[_sliceBounds[t].first*_linesPerSlice];
		int      entryEnd   = _startEntry[_sliceBounds[t].second*_linesPerSlice-1] + _entriesPerLine[_sliceBounds[t].second*_linesPerSlice-1];
		int skirtEntryEnd   = _startEntry[_skirtBounds[t].second*_linesPerSlice-1] + _entriesPerLine[_skirtBounds[t].second*_linesPerSlice-1];
		_pSolution[t] = AllocPointer< T2 >( skirtEntryEnd - skirtEntryStart );
//		_pSolution[t] = ( T2* ) malloc( sizeof( T2 ) * ( skirtEntryEnd - skirtEntryStart ) );
		_pSolution[t] += entryStart-skirtEntryStart;
	}
}
template< class T2 >
ParallelSolution2< T2 >::~ParallelSolution2( void )
{
	for( int t=0 ; t<_threads ; t++ )
	{
		int skirtEntryStart = _startEntry[_skirtBounds[t].first*_linesPerSlice];
		int      entryStart = _startEntry[_sliceBounds[t].first*_linesPerSlice];
		int      entryEnd   = _startEntry[_sliceBounds[t].second*_linesPerSlice-1] + _entriesPerLine[_sliceBounds[t].second*_linesPerSlice-1];
		int skirtEntryEnd   = _startEntry[_skirtBounds[t].second*_linesPerSlice-1] + _entriesPerLine[_skirtBounds[t].second*_linesPerSlice-1];
		_pSolution[t] -= entryStart-skirtEntryStart;
		free( _pSolution[t] );
	}
	_pSolution.clear();
}
template< class T2 >
void ParallelSolution2< T2 >::SetFromArray( ConstPointer( T2 ) solution )
{
#pragma omp parallel for num_threads( _threads )
	for( int t=0 ; t<_threads ; t++ )
	{
		int skirtEntryStart = _startEntry[_skirtBounds[t].first*_linesPerSlice];
		int      entryStart = _startEntry[_sliceBounds[t].first*_linesPerSlice];
		int      entryEnd   = _startEntry[_sliceBounds[t].second*_linesPerSlice-1] + _entriesPerLine[_sliceBounds[t].second*_linesPerSlice-1];
		int skirtEntryEnd   = _startEntry[_skirtBounds[t].second*_linesPerSlice-1] + _entriesPerLine[_skirtBounds[t].second*_linesPerSlice-1];
		memcpy( _pSolution[t] - ( entryStart-skirtEntryStart ) , &solution[skirtEntryStart] , sizeof( T2 ) * ( skirtEntryEnd - skirtEntryStart ) );
	}
}
template< class T2 >
void ParallelSolution2< T2 >::SetToArray( Pointer( T2 ) solution ) const
{
#pragma omp parallel for num_threads( _threads )
	for( int t=0 ; t<_threads ; t++ )
	{
		int skirtEntryStart = _startEntry[_skirtBounds[t].first*_linesPerSlice];
		int      entryStart = _startEntry[_sliceBounds[t].first*_linesPerSlice];
		int      entryEnd   = _startEntry[_sliceBounds[t].second*_linesPerSlice-1] + _entriesPerLine[_sliceBounds[t].second*_linesPerSlice-1];
		int skirtEntryEnd   = _startEntry[_skirtBounds[t].second*_linesPerSlice-1] + _entriesPerLine[_skirtBounds[t].second*_linesPerSlice-1];
		memcpy( solution + entryStart , _pSolution[t] , sizeof( T2 ) * ( entryEnd-entryStart ) );
	}
}
template< class T2 >
void ParallelSolution2< T2 >::clear( void )
{
#pragma omp parallel for num_threads( _threads )
	for( int t=0 ; t<_threads ; t++ )
	{
		int skirtEntryStart = _startEntry[_skirtBounds[t].first*_linesPerSlice];
		int      entryStart = _startEntry[_sliceBounds[t].first*_linesPerSlice];
		int      entryEnd   = _startEntry[_sliceBounds[t].second*_linesPerSlice-1] + _entriesPerLine[_sliceBounds[t].second*_linesPerSlice-1];
		int skirtEntryEnd   = _startEntry[_skirtBounds[t].second*_linesPerSlice-1] + _entriesPerLine[_skirtBounds[t].second*_linesPerSlice-1];
		memset( _pSolution[t] - ( entryStart-skirtEntryStart ) , 0 , sizeof( T2 ) * ( skirtEntryEnd - skirtEntryStart ) );
	}
}
template< class T2 >
void ParallelSolution2< T2 >::_clearSkirts( void )
{
#pragma omp parallel for num_threads( _threads )
	for( int t=0 ; t<_threads ; t++ )
	{
		int skirtEntryStart = _startEntry[_skirtBounds[t].first*_linesPerSlice];
		int      entryStart = _startEntry[_sliceBounds[t].first*_linesPerSlice];
		int      entryEnd   = _startEntry[_sliceBounds[t].second*_linesPerSlice-1] + _entriesPerLine[_sliceBounds[t].second*_linesPerSlice-1];
		int skirtEntryEnd   = _startEntry[_skirtBounds[t].second*_linesPerSlice-1] + _entriesPerLine[_skirtBounds[t].second*_linesPerSlice-1];
		if( entryStart - skirtEntryStart ) memset( _pSolution[t] - ( entryStart-skirtEntryStart ) , 0 , sizeof( T2 ) * ( entryStart - skirtEntryStart ) );
		if( skirtEntryEnd - entryEnd )     memset( _pSolution[t] + ( entryEnd - entryStart ) , 0 , sizeof( T2 ) * ( skirtEntryEnd - entryEnd ) );
	}
}
// WARNING!!! The synchronization  and merging code is not completely correct. If the skirts are larger than the neighbors, there will be trouble.
template< class T2 >
void ParallelSolution2< T2 >::_merge( void )
{
#pragma omp parallel for num_threads( _threads-1 )
	for( int t=0 ; t<_threads-1 ; t++ )
	{
		int skirtEntryStart1 = _startEntry[_skirtBounds[t].first*_linesPerSlice];
		int      entryStart1 = _startEntry[_sliceBounds[t].first*_linesPerSlice];
		int      entryEnd1   = _startEntry[_sliceBounds[t].second*_linesPerSlice-1] + _entriesPerLine[_sliceBounds[t].second*_linesPerSlice-1];
		int skirtEntryEnd1   = _startEntry[_skirtBounds[t].second*_linesPerSlice-1] + _entriesPerLine[_skirtBounds[t].second*_linesPerSlice-1];
		int skirtEntryStart2 = _startEntry[_skirtBounds[t+1].first*_linesPerSlice];
		int      entryStart2 = _startEntry[_sliceBounds[t+1].first*_linesPerSlice];
		int      entryEnd2   = _startEntry[_sliceBounds[t+1].second*_linesPerSlice-1] + _entriesPerLine[_sliceBounds[t+1].second*_linesPerSlice-1];
		int skirtEntryEnd2   = _startEntry[_skirtBounds[t+1].second*_linesPerSlice-1] + _entriesPerLine[_skirtBounds[t+1].second*_linesPerSlice-1];

		int mergeSize;
//		T2* left;
//		T2* right;
		Pointer( T2 ) left;
		Pointer( T2 ) right;

		// Merge the left side of the boundary
		mergeSize = entryStart2 - skirtEntryStart2;
		left  = _pSolution[t  ] + ( entryEnd1-entryStart1 ) - ( entryStart2-skirtEntryStart2 );
		right = _pSolution[t+1] - ( entryStart2-skirtEntryStart2 );
		for( int i=0 ; i<mergeSize ; i++ ) left[i] = right[i] = ( left[i]+right[i] );

		// Merge the right side of the boundary
		mergeSize = skirtEntryEnd1-entryEnd1;
		left  = _pSolution[t  ] + ( entryEnd1-entryStart1 );
		right = _pSolution[t+1];
		for( int i=0 ; i<mergeSize ; i++ ) left[i] = right[i] = ( left[i]+right[i] );
	}
}
template< class T2 >
void ParallelSolution2< T2 >::synchronize( void )
{
#pragma omp parallel for num_threads( _threads-1 )
	for( int t=0 ; t<_threads-1 ; t++ )
	{
		int skirtEntryStart1 = _startEntry[_skirtBounds[t  ].first *_linesPerSlice  ];
		int      entryStart1 = _startEntry[_sliceBounds[t  ].first *_linesPerSlice  ];
		int      entryEnd1   = _startEntry[_sliceBounds[t  ].second*_linesPerSlice-1] + _entriesPerLine[_sliceBounds[t  ].second*_linesPerSlice-1];
		int skirtEntryEnd1   = _startEntry[_skirtBounds[t  ].second*_linesPerSlice-1] + _entriesPerLine[_skirtBounds[t  ].second*_linesPerSlice-1];
		int skirtEntryStart2 = _startEntry[_skirtBounds[t+1].first *_linesPerSlice  ];
		int      entryStart2 = _startEntry[_sliceBounds[t+1].first *_linesPerSlice  ];
		int      entryEnd2   = _startEntry[_sliceBounds[t+1].second*_linesPerSlice-1] + _entriesPerLine[_sliceBounds[t+1].second*_linesPerSlice-1];
		int skirtEntryEnd2   = _startEntry[_skirtBounds[t+1].second*_linesPerSlice-1] + _entriesPerLine[_skirtBounds[t+1].second*_linesPerSlice-1];

		memcpy( _pSolution[t+1] - ( entryStart2-skirtEntryStart2 ) , _pSolution[t] + ( entryEnd1-entryStart1 ) - ( entryStart2-skirtEntryStart2 ) , sizeof( T2 ) * ( entryStart2-skirtEntryStart2 ) );
		memcpy( _pSolution[t  ] + ( entryEnd1-entryStart1 ) , _pSolution[t+1] , sizeof( T2 ) * ( skirtEntryEnd1-entryEnd1 ) );
	}
}
template< class T , class const_iterator >
template< class T2 >
void SparseMatrixInterface< T , const_iterator >::MultiplyParallel( ConstPointer( T2 ) in , ParallelSolution2< T2 >* out ) const
{
	const int threads = out->_threads;
	const int sliceDependence = out->_sliceDependence;
	const int sliceCount = out->_sliceCount;
	const int* startEntry = &out->_startEntry[0];
	const int* entriesPerLine = out->_entriesPerLine;
	int linesPerSlice = out->_linesPerSlice;
#pragma omp parallel for num_threads( threads )
	for( int t=0 ; t<threads ; t++ )
	{
		// Compute the starting and ending slices associated with the current thread.
		int outStartLine = out->_sliceBounds[t].first  * linesPerSlice;
		int   outEndLine = out->_sliceBounds[t].second * linesPerSlice;
		int    outOffset = out->_startEntry[ outStartLine ];

//		T2* _out = out->_pSolution[t] - outOffset;
		Pointer( T2 ) _out = out->_pSolution[t] - outOffset;

		for( int i=outStartLine ; i<outEndLine ; i++ )
		{
			int start = startEntry[i];
			int stop  = startEntry[i] + entriesPerLine[i];
			for( size_t j=start ; j<stop ; j++ )
			{
				T2 temp;
				temp *= 0;
				ConstPointer( T2 ) _in = in;
				const_iterator e = end( j );
				for( const_iterator iter = begin( j ) ; iter!=e ; iter++ ) temp += _in[ iter->N ] * iter->Value;
				_out[j] += temp;
			}
		}
	}
}
// WARNING: In general, the next two function calls are not safe since they assume that in order to
// set out->_pSolution[t], we only need to know the values in in->_pSolution[t]!
template< class T , class const_iterator >
template< class T2 >
void SparseMatrixInterface< T , const_iterator >::MultiplyParallel( const ParallelSolution2< T2 >* in , ParallelSolution2< T2 >* out ) const
{
	int iters = out->_iters;
	const int threads = out->_threads;
	const int sliceDependence = out->_sliceDependence;
	const int sliceCount = out->_sliceCount;
	const int* startEntry = &out->_startEntry[0];
	const int* entriesPerLine = out->_linesPerSlice;
	int linesPerSlice = out->_linesPerSlice;
#pragma omp parallel for num_threads( threads )
	for( int t=0 ; t<threads ; t++ )
	{
		// Compute the starting and ending slices associated with the current thread.
		int  inStartSlice =  in->_sliceBounds[t].first;
		int    inEndSlice =  in->_sliceBounds[t].second;
		int outStartSlice = out->_sliceBounds[t].first;
		int   outEndSlice = out->_sliceBounds[t].second;
		int      inOffset =  in->_startEntry[  inStartSlice *  in->_linesPerSlice ];
		int     outOffset = out->_startEntry[ outStartSlice * out->_linesPerSlice ];

		ConstPointer( T2 ) _in  =  in->_pSolution[t] -  inOffset;
		T2*       _out = out->_pSolution[t] - outOffset;

		int lastStartLine  = 0;
		int lastStartSlice = outStartSlice;
		int currentSlice = lastStartSlice;
		int currentLine  = lastStartLine;

		while( 1 )
		{
			int start = startEntry[currentSlice*linesPerSlice+currentLine];
			int  stop = startEntry[currentSlice*linesPerSlice+currentLine] + entriesPerLine[currentSlice*linesPerSlice+currentLine];
			for( int jj=start ; jj<stop ; jj++ )
			{
				T2 temp;
				temp *= 0;
				ConstPointer( T2 ) in = _in;
				const_iterator e = end( jj );
				for( const_iterator iter = begin( jj ) ; iter!=e ; iter++ ) temp += in[ iter->N ] * iter->Value;
				_out[jj] += temp;
			}
			currentSlice++;
			currentLine--;
			if( currentSlice>=outEndSlice || currentLine<0 )
			{
				if( lastStartLine!=out->_linesPerSlice-1 ) lastStartLine++;
				else if( lastStartSlice!=outStartSlice-1 ) lastStartSlice++;
				else break;
				currentSlice = lastStartSlice;
				currentLine  = lastStartLine;
			}
		}
	}
}
template< class T , class const_iterator >
template< class T2 >
void SparseMatrixInterface< T , const_iterator >::MultiplyTransposeParallel( const ParallelSolution2< T2 >* in , ParallelSolution2< T2 >* out ) const
{
	const int threads = in->_threads;
	const int sliceDependence = in->_sliceDependence;
	const int sliceCount = in->_sliceCount;
	const int* startEntry = &in->_startEntry[0];
	const int* entriesPerLine = in->_entriesPerLine;
	int linesPerSlice = in->_linesPerSlice;

	out->_clearSkirts( );
#pragma omp parallel for num_threads( threads )
	for( int t=0 ; t<threads ; t++ )
	{
		// Compute the starting and ending slices associated with the current thread.
		int  inStartSlice =  in->_sliceBounds[t].first;
		int    inEndSlice =  in->_sliceBounds[t].second;
		int outStartSlice = out->_sliceBounds[t].first;
		int   outEndSlice = out->_sliceBounds[t].second;
		int      inOffset =  in->_startEntry[  inStartSlice *  in->_linesPerSlice ];
		int     outOffset = out->_startEntry[ outStartSlice * out->_linesPerSlice ];

		const Pointer( T2 ) _in  = ( const Pointer( T2 ) )( in->_pSolution[t] ) -  inOffset;
		Pointer( T2 )      _out = out->_pSolution[t] - outOffset;
//		ConstPointer( T2 ) _in  =  in->_pSolution[t] -  inOffset;
//		T2*       _out = out->_pSolution[t] - outOffset;

		int lastStartLine  = 0;
		int lastStartSlice = inStartSlice;
		int currentSlice = lastStartSlice;
		int currentLine  = lastStartLine;

		while( 1 )
		{
			int start = startEntry[currentSlice*linesPerSlice+currentLine];
			int  stop = startEntry[currentSlice*linesPerSlice+currentLine] + entriesPerLine[currentSlice*linesPerSlice+currentLine];
			for( int jj=start ; jj<stop ; jj++ )
			{
				T2 temp = _in[ jj ];
				Pointer( T2 ) out = _out;
				const_iterator e = end( jj );
				for( const_iterator iter = begin( jj ) ; iter!=e ; iter++ ) out[ iter->N ] += temp * iter->Value;
			}
			currentSlice++;
			currentLine--;
			if( currentSlice>=inEndSlice || currentLine<0 )
			{
				if( lastStartLine!=linesPerSlice-1 ) lastStartLine++;
				else if( lastStartSlice!=inEndSlice-1 ) lastStartSlice++;
				else break;
				currentSlice = lastStartSlice;
				currentLine  = lastStartLine;
			}
		}
	}
	out->_merge( );
}
template< class T , class const_iterator >
template< class T2 , class T3 , bool UnitDiagonal , bool StrippedDiagonal , bool UseSOR >
void SparseMatrixInterface< T , const_iterator >::SolveGaussSeidelParallel(	ConstPointer( T3 ) Diagonal , ConstPointer( T2 ) b , int iters , ParallelSolution2< T2 >* Solution , T sor ) const
{
	const int threads         =  Solution->_threads;
	const int sliceDependence =  Solution->_sliceDependence;
	const int sliceCount      =  Solution->_sliceCount;
	const int* startEntry     = &Solution->_startEntry[0];
	const int* entriesPerLine =  Solution->_entriesPerLine;
	const int linesPerSlice   =  Solution->_linesPerSlice;
	const int delta = (linesPerSlice+1)*sliceDependence;
#pragma omp parallel for num_threads( threads )
	for( int t=0 ; t<threads ; t++ )
	{
		// Compute the starting and ending slices associated with the current thread.
		const int startSlice = Solution->_sliceBounds[t].first;
		const int   endSlice = Solution->_sliceBounds[t].second;
#if 1
		const int     offset = startEntry[ startSlice*linesPerSlice ];
		const int startStartBound = startSlice - iters*sliceDependence;
		const int   startEndBound =   endSlice + iters*sliceDependence;
#else
//		const int     offset = Solution->_startEntry[ startSlice*linesPerSlice ];
		const int startStartBound = startSlice - (iters-1)*sliceDependence;
//		const int   startEndBound =   endSlice + (iters-1)*sliceDependence+1;
		const int   startEndBound =   endSlice + (iters-1)*sliceDependence;
#endif
		Pointer( T2 ) solution = Solution->_pSolution[t] - offset;

#if 1
		for( int s=startSlice-sliceDependence ; s<endSlice+iters*sliceDependence ; s++ )
#else
		for( int s=startSlice ; s<endSlice+(iters-1)*sliceDependence ; s++ )	// The lead slice for relaxation
#endif
		{
			const int* _startEntry     = startEntry     + s*linesPerSlice;
			const int* _entriesPerLine = entriesPerLine + s*linesPerSlice;
#if 1
			for( int l=0 ; l<linesPerSlice+iters*sliceDependence ; l++ , _startEntry++ , _entriesPerLine++ ) // The line on the lead slice for relaxation
#else
			for( int l=0 ; l<linesPerSlice+(iters-1)*sliceDependence ; l++ , _startEntry++ , _entriesPerLine++ )					// The line on the lead slice for relaxation
#endif
			{
				int ss = s;
				int ll = l;
				int startBound = startStartBound;
				int   endBound =   startEndBound;
				const int* __startEntry = _startEntry;
				const int* __entriesPerLine = _entriesPerLine;
				for( int i=0 ; i<iters ; i++ , ss-=sliceDependence , ll-=sliceDependence , startBound+=sliceDependence , endBound-=sliceDependence , __startEntry-=delta , __entriesPerLine-=delta )
				{
					if( ss>=startBound && ss<endBound && ll>=0 && ll<linesPerSlice && ss>=0 && ss<sliceCount )
					{
						int start = *__startEntry;
						int  stop = start+(*__entriesPerLine);
						for( int jj=start ; jj<stop ; jj++ )
						{
							T2 temp;
							temp *= 0;
							Pointer( T2 ) _solution = solution;
							if( UnitDiagonal )
							{
								const_iterator e = end( jj );
								for( const_iterator iter = begin( jj ) ; iter!=e ; iter++ ) temp += _solution[ iter->N ] * iter->Value;
								if( UseSOR )
									if(StrippedDiagonal) solution[jj]  = solution[jj]*(T(1.0)-sor) + (b[jj]-temp) * sor;
									else                 solution[jj] +=                             (b[jj]-temp) * sor;
								else
									if(StrippedDiagonal) solution[jj]  = (b[jj]-temp);
									else                 solution[jj] += (b[jj]-temp);
							}
							else
							{
								T dValue = T(1.) / Diagonal[jj];
								const_iterator e = end( jj );
								for( const_iterator iter = begin( jj ) ; iter!=e ; iter++ ) temp += _solution[ iter->N ] * iter->Value;
								if( UseSOR )
									if(StrippedDiagonal) solution[jj]  = solution[jj]*(T(1.0)-sor) + (b[jj]-temp) * dValue * sor;
									else                 solution[jj] +=                             (b[jj]-temp) * dValue * sor;
								else
									if(StrippedDiagonal) solution[jj]  = (b[jj]-temp) * dValue;
									else                 solution[jj] += (b[jj]-temp) * dValue;
							}
						}
					}
				}
			}
		}
	}
}

template< class T , class const_iterator >
template< class T2 , class T3 , bool UnitDiagonal , bool StrippedDiagonal , bool UseSOR >
void SparseMatrixInterface< T , const_iterator >::SolveGaussSeidelAndResidualParallel( ConstPointer( T3 ) Diagonal , ConstPointer( T2 ) b , int iters , ParallelSolution2< T2 >* Solution , Pointer( T2 ) Residual , T sor ) const
{
	const int threads         =  Solution->_threads;
	const int sliceDependence =  Solution->_sliceDependence;
	const int sliceCount      =  Solution->_sliceCount;
	const int* startEntry     = &Solution->_startEntry[0];
	const int* entriesPerLine =  Solution->_entriesPerLine;
	const int linesPerSlice   =  Solution->_linesPerSlice;
	const int delta = (linesPerSlice+1)*sliceDependence;
#pragma omp parallel for num_threads( threads )
	for( int t=0 ; t<threads ; t++ )
	{
		// Compute the starting and ending slices associated with the current thread.
		const int startSlice = Solution->_sliceBounds[t].first;
		const int   endSlice = Solution->_sliceBounds[t].second;
		const int     offset = startEntry[ startSlice*linesPerSlice ];
		const int startStartBound = startSlice - iters*sliceDependence;
//		const int   startEndBound =   endSlice + iters*sliceDependence+1;
		const int   startEndBound =   endSlice + iters*sliceDependence;
		Pointer( T2 ) solution = Solution->_pSolution[t] - offset;

		// Need to start one back so that the neighbor is that we can compute the relaxation correctly at the first slice.
		for( int s=startSlice-sliceDependence ; s<endSlice+iters*sliceDependence ; s++ )
		{
			const int* _startEntry     = startEntry     + s*linesPerSlice;
			const int* _entriesPerLine = entriesPerLine + s*linesPerSlice;
			for( int l=0 ; l<linesPerSlice+iters*sliceDependence ; l++ , _startEntry++ , _entriesPerLine++ ) // The line on the lead slice for relaxation
			{
				int ss = s;
				int ll = l;
				int startBound = startStartBound;
				int   endBound =   startEndBound;
				const int* __startEntry = _startEntry;
				const int* __entriesPerLine = _entriesPerLine;
				for( int i=0 ; i<iters ; i++ , ss-=sliceDependence , ll-=sliceDependence , startBound+=sliceDependence , endBound-=sliceDependence , __startEntry-=delta , __entriesPerLine-=delta )
				{
					if( ss>=startBound && ss<endBound && ll>=0 && ll<linesPerSlice && ss>=0 && ss<sliceCount )
					{
						int start = *__startEntry;
						int  stop = start+(*__entriesPerLine);
						for( int jj=start ; jj<stop ; jj++ )
						{
							T2 temp;
							temp *= 0;
							Pointer( T2 ) _solution = solution;
							if( UnitDiagonal )
							{
								const_iterator e = end( jj );
								for( const_iterator iter = begin( jj ) ; iter!=e ; iter++ ) temp += _solution[ iter->N ] * iter->Value;
								if( UseSOR )
									if(StrippedDiagonal) solution[jj]  = solution[jj]*T(1.0-sor) + (b[jj]-temp) * sor;
									else                 solution[jj] +=                           (b[jj]-temp) * sor;
								else
									if(StrippedDiagonal) solution[jj]  = (b[jj]-temp);
									else                 solution[jj] += (b[jj]-temp);
							}
							else
							{
								T dValue = T(1.) / Diagonal[jj];
								const_iterator e = end( jj );
								for( const_iterator iter = begin( jj ) ; iter!=e ; iter++ ) temp += _solution[ iter->N ] * iter->Value;
								if( UseSOR )
									if(StrippedDiagonal) solution[jj]  = solution[jj]*T(1.0-sor) + (b[jj]-temp) * dValue * sor;
									else                 solution[jj] +=                           (b[jj]-temp) * dValue * sor;
								else
									if(StrippedDiagonal) solution[jj]  = (b[jj]-temp) * dValue;
									else                 solution[jj] += (b[jj]-temp) * dValue;
							}
						}
					}
				}
				if( ss>=startSlice && ss<endSlice && ll>=0 && ll<linesPerSlice )
				{
					int start = *__startEntry;
					int  stop = start+(*__entriesPerLine);
					for( int jj=start ; jj<stop ; jj++ )
					{
						T2 temp;
						temp *= 0;
						Pointer( T2 ) _solution = solution;
						if( UnitDiagonal )
						{
							const_iterator e = end( jj );
							for( const_iterator iter=begin( jj ) ; iter!=e ; iter++ ) temp += _solution[ iter->N ] * iter->Value;
							if( StrippedDiagonal ) Residual[jj] = b[jj] - temp - solution[jj];
							else                   Residual[jj] = b[jj] - temp;
						}
						else
						{
							T dValue = Diagonal[jj];
							const_iterator e = end( jj );
							for( const_iterator iter=begin( jj ) ; iter!=e ; iter++ ) temp += _solution[ iter->N ] * iter->Value;
							if( StrippedDiagonal ) Residual[jj] = b[jj] - temp - ( solution[jj] * dValue );
							else                   Residual[jj] = b[jj] - temp;
						}
					}
				}
			}
		}
	}
}
template< class T , class const_iterator >
template< class T2 , class T3 , bool UnitDiagonal , bool StrippedDiagonal >
//void SparseMatrixInterface< T , const_iterator >::ResidualParallel( ConstPointer( T3 ) Diagonal , ConstPointer( T2 ) b , ParallelSolution2< T2 >* Solution , Pointer( T2 ) residual ) const
void SparseMatrixInterface< T , const_iterator >::ResidualParallel( ConstPointer( T3 ) Diagonal , ConstPointer( T2 ) b , ParallelSolution2< T2 >* Solution , Pointer( T2 ) Residual ) const
{
	const int threads         =  Solution->_threads;
	const int sliceDependence =  Solution->_sliceDependence;
	const int sliceCount      =  Solution->_sliceCount;
	const int* startEntry     = &Solution->_startEntry[0];
	const int* entriesPerLine =  Solution->_entriesPerLine;
	const int linesPerSlice   =  Solution->_linesPerSlice;
	const int delta = (linesPerSlice+1)*sliceDependence;
#pragma omp parallel for num_threads( threads )
	for( int t=0 ; t<threads ; t++ )
	{
		// Compute the starting and ending slices associated with the current thread.
		const int startSlice = Solution->_sliceBounds[t].first;
		const int   endSlice = Solution->_sliceBounds[t].second;
		const int     offset = startEntry[ startSlice*linesPerSlice ];
		Pointer( T2 ) solution = Solution->_pSolution[t] - offset;

		for( int s=startSlice ; s<endSlice ; s++ )	// The lead slice for relaxation
		{
			const int* _startEntry     = startEntry     + s*linesPerSlice;
			const int* _entriesPerLine = entriesPerLine + s*linesPerSlice;
			for( int l=0 ; l<linesPerSlice ; l++ , _startEntry++ , _entriesPerLine++ ) // The line on the lead slice for relaxation
			{
				int start = *_startEntry;
				int  stop = start+(*_entriesPerLine);
				for( int jj=start ; jj<stop ; jj++ )
				{
					T2 temp;
					temp *= 0;
					T2* _solution = solution;
					if( UnitDiagonal )
					{
						const_iterator e = end( jj );
						for( const_iterator iter=begin( jj ) ; iter!=e ; iter++ ) temp += _solution[ iter->N ] * iter->Value;
						if( StrippedDiagonal ) Residual[jj] = b[jj] - temp - solution[jj];
						else                   Residual[jj] = b[jj] - temp;
					}
					else
					{
						T dValue = Diagonal[jj];
						const_iterator e = end( jj );
						for( const_iterator iter=begin( jj ) ; iter!=e ; iter++ ) temp += _solution[ iter->N ] * iter->Value;
						if( StrippedDiagonal ) Residual[jj] = b[jj] - temp - ( solution[jj] * dValue );
						else                   Residual[jj] = b[jj] - temp;
					}
				}
			}
		}
	}
}

#if 0
// NOTE: THESE ARE BROKEN!
template< class T , class const_iterator >
template< class T2 , bool UnitDiagonal , bool StrippedDiagonal , bool UseSOR >
void SparseMatrixInterface< T , const_iterator >::SolveGaussSeidelParallel(
	const T* Diagonal , ConstPointer( T2 ) b , int iters , Pointer( T2 ) solution , T sor ,
	const std::vector< int >& entriesPerSlice , const std::vector< int >& slicesPerThread , int sliceDependence ) const
{
	ParallelSolution< T2 > pSolution( iters , entriesPerSlice , slicesPerThread , sliceDependence );
	pSolution.SetFromArray( Solution );
	SolveGaussSeidelParallel< T2 , UnitDiagonal , StrippedDiagonal , UseSOR >( Diagonal , b , iters , &pSolution , sor );
	pSolution.SetToArray( Solution );
}


template< class T , class const_iterator >
template< class T2 , bool UnitDiagonal , bool StrippedDiagonal , bool UseSOR >
void SparseMatrixInterface< T , const_iterator >::SolveGaussSeidelAndResidualParallel(
	const T* Diagonal , ConstPointer( T2 ) b , int iters , Pointer( T2 ) solution , Pointer( T2 ) residual , T sor ,
	const std::vector< int >& entriesPerSlice , const std::vector< int >& slicesPerThread , int sliceDependence ) const
{
	ParallelSolution< T2 > pSolution( iters , entriesPerSlice , slicesPerThread , sliceDependence );
	pSolution.SetFromArray( Solution );
	SolveGaussSeidelAndResidualParallel< T2 , UnitDiagonal , StrippedDiagonal , UseSOR >( Diagonal , b , iters , &pSolution , Residual , sor );
	pSolution.SetToArray( Solution );
}

template< class T , class const_iterator >
template< class T2 , class T3 , bool UnitDiagonal , bool StrippedDiagonal , bool UseSOR >
void SparseMatrixInterface< T , const_iterator >::SolveGaussSeidelParallel(
	ConstPointer( T3 ) Diagonal , ConstPointer( T2 ) b , int iters , Pointer( T2 ) solution , T sor ,
	const std::vector< int >& entriesPerSlice , const std::vector< int >& slicesPerThread , int sliceDependence ) const
{
	int threads = slicesPerThread.size();
	if( threads<1 ) return SolveGaussSeidel< T2 , T3 , UnitDiagonal , StrippedDiagonal , UseSOR >( Diagonal , b , iters , Solution , sor );
	std::vector< int > startEntry;

	startEntry.resize( entriesPerSlice.size() );
	startEntry[0] = 0;
	for( int i=1 ; i<entriesPerSlice.size() ; i++ ) startEntry[i] = startEntry[i-1] + entriesPerSlice[i-1];

#pragma omp parallel for num_threads( threads )
	for( int s=0 ; s<threads ; s++ )
	{
		// Compute the starting and ending slices associated with the current thread.
		int startSlice = 0 , endSlice;
		for( int i=0 ; i<s ; i++ ) startSlice += slicesPerThread[i];
		endSlice = startSlice + slicesPerThread[s];

		// When processing at index i, we will be relaxing on slices:
		// [ i , i + sliceDependence , ... , i + (iters-1)*sliceDependence ]
		for( int i=startSlice-iters*sliceDependence ; i<endSlice+sliceDependence ; i++ )
		{

			// When processing at index j, we will be relaxing slice:
			// i + j*sliceDepdence for the (iters-j)-th time.
			for( int j=iters-1 ; j>=0 ; j-- ) 
			{
				int ii = i + j*sliceDependence;
				// Compute the starting and ending bounds of the skirt on which we are relaxing (keeping in mind that we will be doing a residual calcuation)
				int startBound = startSlice - (j+1)*sliceDependence;
				int   endBound =   endSlice + (j+1)*sliceDependence+1;
				if( startBound<0 ) startBound = 0;
				if(   endBound>entriesPerSlice.size() ) endBound = entriesPerSlice.size();
				if( ii>=startBound && ii<endBound )
				{
					int start = startEntry[ii];
					int stop = startEntry[ii] + entriesPerSlice[ii];
					if( UnitDiagonal )
						for( size_t j=start ; j<stop ; j++ )
						{
							T2 temp;
							temp *= 0;
							T2* _Solution = Solution;
							const_iterator e = end( j );
							for( const_iterator iter = begin( j ) ; iter!=e ; iter++ ) temp += _Solution[ iter->N ] * iter->Value;
							if( UseSOR )
								if( StrippedDiagonal ) Solution[j]  = Solution[j]*(1.0-sor) + (b[j]-temp) * sor;
								else                   Solution[j] +=                         (b[j]-temp) * sor;
							else
								if( StrippedDiagonal ) Solution[j]  = (b[j]-temp);
								else                   Solution[j] += (b[j]-temp);
						}
					else
						for( size_t j=start ; j<stop ; j++ )
						{
							T2 temp;
							temp *= 0;
							T2* _Solution = Solution;
							T dValue = T(1.) / Diagonal[j];
							const_iterator e = end( j );
							for( const_iterator iter = begin( j ) ; iter!=e ; iter++ ) temp += _Solution[ iter->N ] * iter->Value;
							if( UseSOR )
								if(StrippedDiagonal) Solution[j]  = Solution[j]*(1.0-sor) + (b[j]-temp) * dValue * sor;
								else                 Solution[j] +=                         (b[j]-temp) * dValue * sor;
							else
								if(StrippedDiagonal) Solution[j]  = (b[j]-temp) * dValue;
								else                 Solution[j] += (b[j]-temp) * dValue;
						}
				}
			}
		}
	}
}
template< class T , class const_iterator >
template< class T2 , class T3 , bool UnitDiagonal , bool StrippedDiagonal , bool UseSOR >
void SparseMatrixInterface< T , const_iterator >::SolveGaussSeidelAndResidualParallel(
	ConstPointer( T3 ) Diagonal , ConstPointer( T2 ) b , int iters , Pointer( T2 ) solution , Pointer( T2 ) residual , T sor ,
	const std::vector< int >& entriesPerSlice , const std::vector< int >& slicesPerThread , int sliceDependence ) const
{
	int threads = slicesPerThread.size();
	if( threads<1 ) 
	{
		SolveGaussSeidelAndResidual< T2 , T3 , UnitDiagonal , StrippedDiagonal , UseSOR >( Diagonal , b , iters , Solution , Residual , sor );
		return;
	}
	std::vector< int > startEntry;

	startEntry.resize( entriesPerSlice.size() );
	startEntry[0] = 0;
	for( int i=1 ; i<entriesPerSlice.size() ; i++ ) startEntry[i] = startEntry[i-1] + entriesPerSlice[i-1];

#pragma omp parallel for num_threads( threads )
	for( int s=0 ; s<threads ; s++ )
	{
		// Compute the starting and ending slices associated with the current thread.
		int startSlice = 0 , endSlice;
		for( int i=0 ; i<s ; i++ ) startSlice += slicesPerThread[i];
		endSlice = startSlice + slicesPerThread[s];

		// When processing at index i, we will be relaxing on slices:
		// [ i , i + sliceDependence , ... , i + (iters-1)*sliceDependence ]
		for( int i=startSlice-iters*sliceDependence ; i<endSlice+sliceDependence ; i++ )
		{

			// When processing at index j, we will be relaxing slice:
			// i + j*sliceDepdence for the (iters-j)-th time.
			for( int j=iters-1 ; j>=0 ; j-- ) 
			{
				int ii = i + j*sliceDependence;
				// Compute the starting and ending bounds of the skirt on which we are relaxing (keeping in mind that we will be doing a residual calcuation)
				int startBound = startSlice - (j+1)*sliceDependence;
				int   endBound =   endSlice + (j+1)*sliceDependence+1;
				if( startBound<0 ) startBound = 0;
				if(   endBound>entriesPerSlice.size() ) endBound = entriesPerSlice.size();
				if( ii>=startBound && ii<endBound )
				{
					int start = startEntry[ii];
					int stop = startEntry[ii] + entriesPerSlice[ii];
					if( UnitDiagonal )
						for( size_t j=start ; j<stop ; j++ )
						{
							T2 temp;
							temp *= 0;
							T2* _Solution = Solution;
							const_iterator e = end( j );
							for( const_iterator iter=begin( j ) ; iter!=e ; iter++ ) temp += _Solution[ iter->N ] * iter->Value;
							if( UseSOR )
								if(StrippedDiagonal) Solution[j]  = Solution[j]*(1.0-sor) + (b[j]-temp) * sor;
								else                 Solution[j] +=                         (b[j]-temp) * sor;
							else
								if(StrippedDiagonal) Solution[j]  = (b[j]-temp);
								else                 Solution[j] += (b[j]-temp);
						}
					else
						for( size_t j=start ; j<stop ; j++ )
						{
							T2 temp;
							temp *= 0;
							T2* _Solution = Solution;
							T dValue = T(1.) / Diagonal[j];
							const_iterator e = end( j );
							for( const_iterator iter=begin( j ) ; iter!=e ; iter++ ) temp += _Solution[ iter->N ] * iter->Value;
							if( UseSOR )
								if(StrippedDiagonal) Solution[j]  = Solution[j]*(1.0-sor) + (b[j]-temp) * dValue * sor;
								else                 Solution[j] +=                         (b[j]-temp) * dValue * sor;
							else
								if(StrippedDiagonal) Solution[j]  = (b[j]-temp) * dValue;
								else                 Solution[j] += (b[j]-temp) * dValue;
						}
				}
			}
			int ii = i-sliceDependence;
			if( ii>=startSlice && ii<endSlice )
			{
				int start = startEntry[ii];
				int stop = startEntry[ii] + entriesPerSlice[ii];
				if( UnitDiagonal )
					for( size_t j=start ; j<stop ; j++ )
					{
						T2 temp;
						temp *= 0;
						T2* _Solution = Solution;
						const_iterator e = end(j);
						for( const_iterator iter=begin( j ) ; iter!=e ; iter++ ) temp += _Solution[ iter->N ] * iter->Value;
						if( StrippedDiagonal ) Residual[j] = b[j] - temp - Solution[j];
						else                   Residual[j] = b[j] - temp;
					}
				else
					for( size_t j=start ; j<stop ; j++ )
					{
						T2 temp;
						temp *= 0;
						T2* _Solution = Solution;
						T dValue = Diagonal[j];
						const_iterator e = end(j);
						for( const_iterator iter=begin( j ) ; iter!=e ; iter++ ) temp += _Solution[ iter->N ] * iter->Value;
						if( StrippedDiagonal ) Residual[j] = b[j] - temp - ( Solution[j] * dValue );
						else                   Residual[j] = b[j] - temp;
					}
			}
	
		}
	}
}
#endif
