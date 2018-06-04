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

#include <FST_semi_memo_fftw.h>
#include <cospmls.h>
#include "fftw3.h"
#include <math.h>

//////////////////
// FourierKeyS2 //
//////////////////
template<class Real> FourierKeyS2<Real>::FourierKeyS2(void){
	bw = 0;
	values = NullPointer< Complex< Real > >();
}
template<class Real> FourierKeyS2<Real>::~FourierKeyS2(void){
	if( values ) FreePointer( values );
	values = NullPointer< Complex< Real > >();
	bw=0;
}
template<class Real> int FourierKeyS2<Real>::read(const char* fileName)
{
	FILE* fp=fopen(fileName,"rb");
	if(!fp){return 0;}
	int r=read(fp);
	fclose(fp);
	return r;
}
template<class Real> int FourierKeyS2<Real>::write(const char* fileName) const
{
	FILE* fp=fopen(fileName,"wb");
	if(!fp){return 0;}
	int w=write(fp);
	fclose(fp);
	return w;
}
template<class Real> int FourierKeyS2<Real>::read(FILE* fp){
	int b,r;
	r=int(fread(&b,sizeof(int),1,fp));
	if(!r){return 0;}
	resize(b);
	r=int(fread(values,sizeof(Complex<Real>),((bw*bw+bw)>>1),fp));
	if(r==((bw*bw+bw)>>1)){return 1;}
	else{return 0;}
}
template<class Real> int FourierKeyS2<Real>::write(FILE* fp) const {
	int w;
	w=int(fwrite(&bw,sizeof(int),1,fp));
	if(!w){return 0;}
	w=int(fwrite(values,sizeof(Complex<Real>),((bw*bw+bw)>>1),fp));
	if(w==((bw*bw+bw)>>1)){return 1;}
	else{return 0;}
}
template<class Real> int FourierKeyS2<Real>::bandWidth( void ) const{return bw;}
template<class Real> int FourierKeyS2<Real>::resolution( void ) const {return bw*2;}
template<class Real> int FourierKeyS2<Real>::resize( int resolution , bool clr )
{
	int b=resolution>>1;
	if( b<0 ) return 0;
	else if( b!=bw )
	{
		if( values ) FreePointer( values );
		values = NullPointer< Complex< Real > >();
		bw=0;
		if(b)
		{
			values = AllocPointer< Complex< Real > >( (b*b+b)>>1 );
			if( !values ) return 0;
			else bw=b;
		}
	}
	if(clr) clear();
	return 1;
}
template<class Real> void FourierKeyS2<Real>::clear(void){if(bw){memset(values,0,sizeof(Complex<Real>)*((bw*bw+bw)>>1));}}
template<class Real> Complex<Real>& FourierKeyS2<Real>::operator() ( int i , int j )       { return values[(i-j)+(j*bw)-(j*j-j)/2]; }
template<class Real> Complex<Real>  FourierKeyS2<Real>::operator() ( int i , int j ) const { return values[(i-j)+(j*bw)-(j*j-j)/2]; }
template<class Real> Real FourierKeyS2<Real>::squareNorm( void ) const { return Dot(*this,*this); }
template<class Real> Real FourierKeyS2<Real>::SquareDifference( const FourierKeyS2& g1,const FourierKeyS2& g2){ return g1.squareNorm() + g2.squareNorm() - 2*Dot(g1,g2); }
template<class Real> Real FourierKeyS2<Real>::Dot( const FourierKeyS2& g1,const FourierKeyS2& g2 )
{
	Real d = Real(0);
	int idx=0;
	if(g1.bw != g2.bw) fprintf( stderr , "Could not compare arrays of different sizes: %d != %d\n" , g1.bw , g2.bw ) , exit(0);
	for( int i=0 ; i<g1.bw ; i++ ) d+=g1.values[i].r * g2.values[i].r;
	for( int i=g1.bw ; i<((g1.bw*g1.bw+g1.bw)>>1) ; i++ ) d += ( g1.values[i]*g2.values[i].conjugate() ).r * 2;
	return d;
}
template< class Real > int FourierKeyS2< Real >::Entries( int bw ) { return (bw*bw+bw)>>1; }
/////////////////////////////////////
// HarmonicTransform::ScratchSpace //
/////////////////////////////////////
template<class Real>
HarmonicTransform<Real>::ScratchSpace::ScratchSpace(void){
	bw=0;
	workSpace = resultSpace = transposeResultSpace = NullPointer< Real >();
	table=transposeTable=NULL;
}
template<class Real>
HarmonicTransform<Real>::ScratchSpace::~ScratchSpace(void){resize(0);}
template<class Real>
void HarmonicTransform<Real>::ScratchSpace::resize(const int& b){
	if(b!=bw)
	{
		int size=b*2;
		if(workSpace)				{delete[] workSpace;}
		if(resultSpace)				{delete[] resultSpace;}
		if(transposeResultSpace)	{delete[] transposeResultSpace;}
		if(table)					{delete[] table;}
		if(transposeTable)			{delete[] transposeTable;}
		bw=0;
		workSpace = NullPointer< Real >();
		resultSpace = NullPointer< Real >() ;
		transposeResultSpace = NullPointer< Real >();
		table = NullPointer< Real* >();
		transposeTable = NullPointer< Real* >();
		if( b>0 )
		{
			bw=b;
			workSpace = AllocPointer< Real >( 4*bw*bw+36*bw );
			resultSpace = AllocPointer< Real >( Spharmonic_TableSize(bw) );
			transposeResultSpace = AllocPointer< Real >( Spharmonic_TableSize(bw) );

			table			=          Spharmonic_Pml_Table(bw,resultSpace,workSpace);
			transposeTable	=Transpose_Spharmonic_Pml_Table(table,bw,transposeResultSpace,workSpace);
		}
	}
}
///////////////////////
// HarmonicTransform //
///////////////////////
template<class Real>
void HarmonicTransform<Real>::resize(const int& resolution){scratch.resize(resolution>>1);}
template<>
int HarmonicTransform<double>::ForwardFourier(SphericalGrid<double>& g,FourierKeyS2<double>& key){
	int sz,bw;
	sz=g.resolution();
	bw=sz>>1;
	if(key.resolution()!=sz){key.resize(sz);}
	scratch.resize(bw);
	FST_semi_memo_fftw( g[0] , (fftw_complex*)&key(0,0) , sz , scratch.table , GetAddress( scratch.workSpace ) );
	return 1;
}
template<>
int HarmonicTransform<float>::ForwardFourier(SphericalGrid<float>& g,FourierKeyS2<float>& key){
	int sz,bw;
	sz=g.resolution();
	bw=sz>>1;
	if(key.resolution()!=sz){key.resize(sz);}
	scratch.resize(bw);
	FST_semi_memo_fftw(g[0],(fftwf_complex*)&key(0,0),sz,scratch.table,scratch.workSpace);
	return 1;
}
template<class Real>
int HarmonicTransform<Real>::ForwardFourier(SphericalGrid<Real>&,FourierKeyS2<Real>&){
	fprintf(stderr,"Harmonic Transform only supported for floats and doubles\n");
	return 0;
}
template<>
int HarmonicTransform<double>::InverseFourier(FourierKeyS2<double>& key,SphericalGrid<double>& g){
	if(key.resolution()!=g.resolution()){g.resize(key.resolution());}
	int bw=key.bandWidth(),sz=g.resolution();
	scratch.resize(bw);

	InvFST_semi_memo_fftw((fftw_complex*)&key(0,0),g[0],sz,scratch.transposeTable,scratch.workSpace);
	return 1;
}
template<>
int HarmonicTransform<float>::InverseFourier(FourierKeyS2<float>& key,SphericalGrid<float>& g){
	if(key.resolution()!=g.resolution()){g.resize(key.resolution());}
	int bw=key.bandWidth(),sz=g.resolution();
	scratch.resize(bw);

	InvFST_semi_memo_fftw((fftwf_complex*)&key(0,0),g[0],sz,scratch.transposeTable,scratch.workSpace);
	return 1;
}
template<class Real>
int HarmonicTransform<Real>::InverseFourier(FourierKeyS2<Real>&,SphericalGrid<Real>&){
	fprintf(stderr,"Harmonic Transform only supported for floats and doubles\n");
	return 0;
}