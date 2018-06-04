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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <fftw3.h>
#include <math.h>

//////////////////
// FourierKey1D //
//////////////////
template<class Real> FourierKey1D<Real>::FourierKey1D(void){
	dim=res=0;
	values=NULL;
}
template<class Real> FourierKey1D<Real>::~FourierKey1D(void){
	if(values){delete[] values;}
	values=NULL;
	dim=res=0;
}
template<class Real> int FourierKey1D<Real>::read(const char* fileName){
	FILE* fp=fopen(fileName,"rb");
	if(!fp){return 0;}
	int r=read(fp);
	fclose(fp);
	return r;
}
template<class Real> int FourierKey1D<Real>::write(const char* fileName) const{
	FILE* fp=fopen(fileName,"wb");
	if(!fp){return 0;}
	int w=write(fp);
	fclose(fp);
	return w;
}
template<class Real> int FourierKey1D<Real>::read(FILE* fp){
	int resolution,r;
	r=int(fread(&resolution,sizeof(int),1,fp));
	if(!r){return 0;}
	resize(resolution);
	r=int(fread(values,sizeof(Complex<Real>),dim,fp));
	if(r==dim){return 1;}
	else{return 0;}
}
template<class Real> int FourierKey1D<Real>::write(FILE* fp) const {
	int w;
	w=int(fwrite(&res,sizeof(int),1,fp));
	if(!w){return 0;}
	w=int(fwrite(values,sizeof(Complex<Real>),dim,fp));
	if(w==dim){return 1;}
	else{return 0;}
}
template<class Real> int FourierKey1D<Real>::size(void) const{return dim;}
template<class Real> int FourierKey1D<Real>::resolution(void) const{return res;}
template<class Real> int FourierKey1D<Real>::resize( int resolution , bool clr )
{
	int d=FourierTransform<Real>::BandWidth(resolution);
	if(resolution<0){return 0;}
	else if(d!=dim){
		if(values){delete[] values;}
		values=NULL;
		dim=0;
		res=0;
		if(d)
		{
			values=new Complex<Real>[d];
			if(!values){return 0;}
			else{dim=d;}
		}
	}
	res=resolution;
	if(clr){clear();}
	return 1;
}
template<class Real> void FourierKey1D<Real>::clear(void){if(dim){memset(values,0,sizeof(Complex<Real>)*dim);}}
template<class Real> Complex<Real>& FourierKey1D<Real>::operator() ( int i )       {return values[i];}
template<class Real> Complex<Real>  FourierKey1D<Real>::operator() ( int i ) const {return values[i];}
template<class Real> Real FourierKey1D<Real>::squareNorm( void ) const{return Dot(*this,*this);}
template<class Real> Real FourierKey1D<Real>::SquareDifference(const FourierKey1D& g1,const FourierKey1D& g2){ return g1.squareNorm()+g2.squareNorm()-2*Dot(g1,g2); }
template<class Real> Real FourierKey1D<Real>::Dot( const FourierKey1D& g1 , const FourierKey1D& g2 )
{
	Real d = Real(0);
	if( g1.res!=g2.res ) fprintf( stderr , "Could not compare arrays of different sizes: %d != %d\n" , g1.dim , g2.dim ) , exit(0);
#if !FIX_SCALING
	Real n = Real(1.0/(2.0*PI));
#endif // !FIX_SCALING
	d += g1.values[0].r * g2.values[0].r;
	if( g1.res & 1 ) for( int i=1 ; i<g1.dim ; i++ ) d += ( g1.values[i] * g2.values[i].conjugate() ).r * 2;
	else
	{
		for( int i=1 ; i<g1.dim-1 ; i++ ) d += ( g1.values[i] * g2.values[i].conjugate() ).r * 2;
		d += g1.values[g1.dim-1].r * g2.values[g1.dim-1].r;
	}
#if FIX_SCALING
	return d;
#else // !FIX_SCALING
	return d*n;
#endif // FIX_SCALING
}

//////////////////////
// FourierTransform //
//////////////////////
template<> void FourierTransform< float >::_ForwardFourier1D( int res , Pointer( float ) values , Pointer( Complex< float > ) coefficients )
{
	GetAddress( coefficients );
	fftwf_plan plan = fftwf_plan_dft_r2c_1d( res , values , (fftwf_complex*)GetAddress( coefficients ) , FFTW_ESTIMATE );
	fftwf_execute( plan );
	fftwf_destroy_plan( plan );
}
template<> void FourierTransform< double >::_ForwardFourier1D( int res , Pointer( double ) values , Pointer( Complex< double > ) coefficients )
{
	fftw_plan plan = fftw_plan_dft_r2c_1d( res , values , (fftw_complex*)GetAddress( coefficients ) , FFTW_ESTIMATE );
	fftw_execute( plan );
	fftw_destroy_plan( plan );
}
template< class Real > void FourierTransform< Real >::_ForwardFourier1D( int res , Pointer( Real ) values , Pointer( Complex< Real > ) coefficients ){ fprintf( stderr , "Only float and double precision FFTs supported\n" ) , exit(0); }

template<> void FourierTransform< float >::_InverseFourier1D( int res , Pointer( Complex< float > ) coefficients , Pointer( float ) values )
{
	fftwf_plan plan = fftwf_plan_dft_c2r_1d( res , (fftwf_complex*)GetAddress( coefficients ) , values , FFTW_ESTIMATE );
	fftwf_execute( plan );
	fftwf_destroy_plan( plan );
}
template<> void FourierTransform< double >::_InverseFourier1D( int res , Pointer( Complex< double > ) coefficients , Pointer( double ) values )
{
	fftw_plan plan = fftw_plan_dft_c2r_1d( res , (fftw_complex*)GetAddress( coefficients ) , values , FFTW_ESTIMATE );
	fftw_execute( plan );
	fftw_destroy_plan( plan );
}
template< class Real > void FourierTransform< Real >::_InverseFourier1D( int res , Pointer( Complex< Real > ) coefficients , Pointer( Real ) values ){ fprintf( stderr , "Only float and double precision inverse FFTs supported\n" ) , exit(0); }

template< class Real >
int FourierTransform< Real >::ForwardFourier( CircularArray< Real >& g , FourierKey1D< Real >& key )
{
	if( key.resolution()!=g.resolution() ) key.resize( g.resolution() );

#if FORCE_FFTW_PRESERVE_INPUT
	static Real* _values = NULL;
	static int _res = 0;
	if( !_values || _res!=g.resolution() )
	{
		if( _values ) free( _values );
		_res = g.resolution();
		_values = (Real*)malloc( sizeof(Real) * _res );
	}
	memcpy( _values , &g(0) , sizeof(Real) * _res );
	_ForwardFourier1D( _res , _values , &key(0) );
#else // !FORCE_FFTW_PRESERVE_INPUT
	_ForwardFourier1D( g.resolution() , &g(0) , &key(0) );
#endif // FORCE_FFTW_PRESERVE_INPUT
#if FIX_SCALING
	Real n = Real( sqrt( 2.0 * PI ) ) / Real( g.resolution() );
#else // !FIX_SCALING
	Real n = Real( 2.0 * PI ) / Real( g.resolution() );
#endif // FIX_SCALING
	for( int i=0 ; i<key.size() ; i++ ) key(i) *= n;
	return 1;
}
template< class Real >
int FourierTransform< Real >::InverseFourier( FourierKey1D< Real >& key , CircularArray< Real >& g)
{
	if( key.resolution()!=g.resolution() ) g.resize(key.resolution());

#if FORCE_FFTW_PRESERVE_INPUT
	static Complex< Real >* _values = NULL;
	static int _res = 0;
	if( !_values || _res!=key.resolution() )
	{
		if( _values ) free( _values );
		_res = key.resolution();
		_values = ( Complex< Real >* )malloc( sizeof( Complex< Real > ) * key.size() );
	}
	memcpy( _values , &key(0) , sizeof( Complex< Real > ) * key.size() );
	_InverseFourier1D( _res , _values , &g(0) );
#else // !FORCE_FFTW_PRESERVE_INPUT
	_InverseFourier1D( g.resolution() , &key(0) , &g(0) );
#endif // FORCE_FFTW_PRESERVE_INPUT
#if FIX_SCALING
	Real n = Real(1.) / Real( sqrt( 2.0 * PI ) );
#else // !FIX_SCALING
	Real n = Real(1.) / Real( 2.0 * PI );
#endif // FIX_SCALING
	for( int i=0 ; i<g.resolution() ; i++ ) g(i) *= n;
	return 1;
}
