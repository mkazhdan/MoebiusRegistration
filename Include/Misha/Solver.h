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
#ifndef SOLVER_H
#define SOLVER_H

#include <Eigen/Sparse>
#include <omp.h>
#include "SparseMatrix.h"
#include "Vector.h"


// Code borrowed from: https://en.wikipedia.org/wiki/Golden_section_search
template< class Real , class Functor >
std::pair< Real , Real > GoldenSectionSearch( Functor& f , Real a , Real b , Real tolerance )
{
	const static Real INVPHI   = (Real)( ( sqrt(5.) - 1 ) / 2. );
	const static Real INVPHI_2 = (Real)( ( 3 - sqrt(5.) ) / 2. );
	const static Real LOG_INVPHI = (Real)log( INVPHI );
	Real delta = b - a;
	if( delta<=tolerance ) return std::pair< Real , Real >( a , b );
	int n = (int)ceil( log(tolerance/delta) / LOG_INVPHI );
	Real c = a + INVPHI_2 * delta;
	Real d = a + INVPHI   * delta;
	Real fc = f(c) , fd = f(d);
	for( int i=0 ; i<n ; i++ )
		if( fc<fd )
		{
			b = d , d = c , fd = fc;
			delta *= INVPHI;
			c = a + INVPHI_2 * delta;
			fc = f(c);
		}
		else
		{
			a = c , c = d , fc = fd;
			delta *= INVPHI;
			d = a + INVPHI * delta;
			fd = f(d);
		}
	if( fc<fd ) return std::pair< Real , Real >( a , d );
	else        return std::pair< Real , Real >( c , b );
}

template< class Real >
class Solver
{
public:
	virtual void update( const SparseMatrix< Real , int >& M ) = 0;
	virtual void solve( const Real* b , Real* x ) = 0;
	virtual size_t dimension( void ) const = 0;
};



template< class Real >
struct EigenSolver
{
	virtual void update( const SparseMatrix< Real , int >& M ) = 0;
	virtual void solve( const Real* b , Real* x ) = 0;
	virtual size_t dimension( void ) const = 0;
};

template< class Real >
class EigenSolverCholeskyLLt : public EigenSolver< Real >
{
	typedef Eigen::SimplicialLLT< Eigen::SparseMatrix< double > > Eigen_Solver;
	typedef Eigen::VectorXd                                       Eigen_Vector;
	Eigen_Solver _solver;
	Eigen_Vector _eigenB;
	Eigen::SparseMatrix< double > _eigenM;
public:
	EigenSolverCholeskyLLt( const SparseMatrix< Real , int >& M , bool analyzeOnly=false )
	{
		_eigenM.resize( int( M.Rows() ) , int( M.Rows() ) );
		std::vector< Eigen::Triplet< double > > triplets;
		triplets.reserve( M.Entries() );
		for( int i=0 ; i<(int)M.Rows() ; i++ ) for( int j=0 ; j<(int)M.RowSize(i) ; j++ ) triplets.push_back( Eigen::Triplet< double >( i , M[i][j].N , M[i][j].Value ) );
		_eigenM.setFromTriplets( triplets.begin() , triplets.end() );
		_solver.analyzePattern( _eigenM );
		if( !analyzeOnly )
		{
			_solver.factorize( _eigenM );
			if( _solver.info()!=Eigen::Success ) fprintf( stderr , "[ERROR] EigenSolverCholeskyLLt::EigenSolverCholeskyLLt Failed to factorize matrix\n" ) , exit(0);
		}
		_eigenB.resize( M.Rows() );
	}
	void update( const SparseMatrix< Real , int >& M )
	{
#pragma omp parallel for
		for( int i=0 ; i<(int)M.Rows() ; i++ ) for( int j=0 ; j<(int)M.RowSize(i) ; j++ ) _eigenM.coeffRef( i , M[i][j].N ) = M[i][j].Value;
		_solver.factorize( _eigenM );
		switch( _solver.info() )
		{
		case Eigen::Success: break;
		case Eigen::NumericalIssue: fprintf( stderr , "[ERROR] EigenSolverCholeskyLLt::update Failed to factorize matrix (numerical issue)\n" ) , exit(0);
		case Eigen::NoConvergence:  fprintf( stderr , "[ERROR] EigenSolverCholeskyLLt::update Failed to factorize matrix (no convergence)\n" ) , exit(0);
		case Eigen::InvalidInput:   fprintf( stderr , "[ERROR] EigenSolverCholeskyLLt::update Failed to factorize matrix (invalid input)\n" ) , exit(0);
		default: fprintf( stderr , "[ERROR] EigenSolverCholeskyLLt::update Failed to factorize matrix\n" ) , exit(0);
		}
	}
	void solve( const Real* b , Real* x )
	{
#pragma omp parallel for
		for( int i=0 ; i<_eigenB.size() ; i++ ) _eigenB[i] = b[i];
		Eigen_Vector eigenX = _solver.solve( _eigenB );
#pragma omp parallel for
		for( int i=0 ; i<eigenX.size() ; i++ ) x[i] = (Real)eigenX[i];
	}
	size_t dimension( void ) const { return _eigenB.size(); }
	static void Solve( const SparseMatrix< Real , int >& M , const Real* b , Real* x ){ EigenSolverCholeskyLLt solver( M ) ; solver.solve( b , x ); }
};
template< class Real >
class EigenSolverCholeskyLDLt : public EigenSolver< Real >
{
	typedef Eigen::SimplicialLDLT< Eigen::SparseMatrix< double > > Eigen_Solver;
	typedef Eigen::VectorXd                                        Eigen_Vector;
	Eigen_Solver _solver;
	Eigen_Vector _eigenB;
public:
	EigenSolverCholeskyLDLt( const SparseMatrix< Real , int >& M , bool analyzeOnly=false )
	{
		Eigen::SparseMatrix< double > eigenM( int( M.Rows() ) , int( M.Rows() ) );
		std::vector< Eigen::Triplet<double> > triplets;
		triplets.reserve( M.Entries() );
		for( int i=0 ; i<M.Rows() ; i++ ) for( int j=0 ; j<(int)M.RowSize(i) ; j++ ) triplets.push_back( Eigen::Triplet< double >( i , M[i][j].N , M[i][j].Value ) );
		eigenM.setFromTriplets( triplets.begin() , triplets.end() );
		_solver.analyzePattern( eigenM );
		if( !analyzeOnly )
		{
			_solver.factorize( eigenM );
			if( _solver.info()!=Eigen::Success ) fprintf( stderr , "[ERROR] EigenSolverCholeskyLDLt::EigenSolverCholeskyLDLt Failed to factorize matrix\n" ) , exit(0);
		}
		_eigenB.resize( M.Rows() );
	}
	void update( const SparseMatrix< Real , int >& M )
	{
		Eigen::SparseMatrix< double > eigenM( int( M.Rows() ) , int( M.Rows() ) );
		std::vector< Eigen::Triplet<double> > triplets;
		triplets.reserve( M.Entries() );
		for( int i=0 ; i<M.Rows() ; i++ ) for( int j=0 ; j<(int)M.RowSize(i) ; j++ ) triplets.push_back( Eigen::Triplet< double >( i , M[i][j].N , M[i][j].Value ) );
		eigenM.setFromTriplets( triplets.begin() , triplets.end() );
		_solver.factorize( eigenM );
		if( _solver.info()!=Eigen::Success ) fprintf( stderr , "[ERROR] EigenSolverCholeskyLDLt::update Failed to factorize matrix\n" ) , exit(0);
	}
	void solve( const Real* b , Real* x )
	{
#pragma omp parallel for
		for( int i=0 ; i<_eigenB.size() ; i++ ) _eigenB[i] = b[i];
		Eigen_Vector eigenX = _solver.solve( _eigenB );
#pragma omp parallel for
		for( int i=0 ; i<eigenX.size() ; i++ ) x[i] = (Real)eigenX[i];
	}
	size_t dimension( void ) const { return _eigenB.size(); }
	static void Solve( const SparseMatrix< Real , int >& M , const Real* b , Real* x ){ EigenSolverCholeskyLDLt solver( M ) ; solver.solve( b , x ); }
};
#endif // SOLVER_H