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
#include <soft_fftw.h>
#include <utils_so3.h>
#include <math.h>
#include "fftw3.h"

///////////////////
// FourierKeySO3 //
///////////////////
template<class Real> int FourierKeySO3<Real>::Entries( int bw ){ return (4*bw*bw*bw+3*bw*bw-bw)/6; }
template<class Real> FourierKeySO3<Real>::FourierKeySO3(void){
	bw=0;
	values = NullPointer< Complex< Real > >();
}
template<class Real> FourierKeySO3<Real>::~FourierKeySO3(void){
	if( values ) FreePointer( values );
	values = NullPointer< Complex< Real > >();
	bw=0;
}
template<class Real> int FourierKeySO3<Real>::read(const char* fileName){
	FILE* fp=fopen(fileName,"rb");
	if(!fp){return 0;}
	int r=read(fp);
	fclose(fp);
	return r;
}
template<class Real> int FourierKeySO3<Real>::write(const char* fileName) const{
	FILE* fp=fopen(fileName,"wb");
	if(!fp){return 0;}
	int w=write(fp);
	fclose(fp);
	return w;
}
template<class Real> int FourierKeySO3<Real>::read(FILE* fp){
	int b,r;
	r=int(fread(&b,sizeof(int),1,fp));
	if(!r){return 0;}
	resize(2*b);
	r=int(fread(values,sizeof(Complex<Real>),Entries(bw),fp));
	if(r==Entries(bw)){return 1;}
	else{return 0;}
}
template<class Real> int FourierKeySO3<Real>::write(FILE* fp) const {
	int w;
	w=int(fwrite(&bw,sizeof(int),1,fp));
	if(!w){return 0;}
	w=int(fwrite(values,sizeof(Complex<Real>),Entries(bw),fp));
	if(w==Entries(bw)){return 1;}
	else{return 0;}
}
template<class Real> int FourierKeySO3<Real>::bandWidth(void) const{return bw;}
template<class Real> int FourierKeySO3<Real>::resolution(void) const {return bw*2;}
template<class Real> int FourierKeySO3<Real>::resize( int resolution , bool clr )
{
	int b=resolution>>1;
	if(b<0) return 0;
	else if(b!=bw)
	{
		if( values ) FreePointer( values );
		values = NullPointer< Complex< Real > >();
		bw=0;
		if(b){
			values = AllocPointer< Complex< Real > >( Entries(b) );
			if( !values ) return 0;
			else bw=b;
		}
	}
	if(clr) clear();
	return 1;
}
template<class Real> void FourierKeySO3<Real>::clear(void){if(bw){memset(values,0,sizeof(Complex<Real>)*Entries(bw));}}
template<class Real> Complex<Real>& FourierKeySO3<Real>::operator() ( int i , int j , int k ) { return values[so3CoefLoc(j,k,i,bw)]; }
template<class Real> Complex<Real> FourierKeySO3<Real>::operator() ( int i , int j , int k ) const {
	if(j<0)
	{
		int jj = -j;
		int kk = -k;
		return values[so3CoefLoc(jj,kk,i,bw)].conjugate();
	}
	else return values[so3CoefLoc(j,k,i,bw)];
}
template<class Real> Real FourierKeySO3<Real>::squareNorm( void ) const{ return Dot(*this,*this); }
template<class Real> Real FourierKeySO3<Real>::SquareDifference( const FourierKeySO3& g1 , const FourierKeySO3& g2 ) { return g1.squareNorm() + g2.squareNorm() - 2*Dot(g1,g2); }
template<class Real> Real FourierKeySO3<Real>::Dot( const FourierKeySO3& g1 , const FourierKeySO3& g2 )
{
	Real d = Real(0);
	int idx=0;
	if(g1.bw != g2.bw) fprintf( stderr , "Could not compare arrays of different sizes: %d != %d\n" , g1.bw , g2.bw ) , exit(0);
	for( int i=0 ; i<g1.bw ; i++ ) for( int j=0 ; j<=i ; j++ )
	{
		d += g1(i,j,0).r * g2(i,j,0).r;
		for( int k=1 ; k<=i ; k++ ) d += ( g1(i,j,k) * g2(i,j,k).conjugate() + g1(i,j,-k) * g2(i,j,-k).conjugate() ).r * 2;
	}
	return d;
}
///////////////////////////////////
// WignerTransform::ScratchSpace //
///////////////////////////////////
template<class Real>
WignerTransform<Real>::ScratchSpace::ScratchSpace(void){
	bw=0;
	data = NullPointer< fftw_complex >();
	workspace_cx = NullPointer< fftw_complex >();
	workspace_cx2 = NullPointer< fftw_complex >();
	coeffs = NullPointer< fftw_complex >();
	workspace_re = NullPointer< double >();
	p=0;
}
template<class Real>
WignerTransform<Real>::ScratchSpace::~ScratchSpace( void ){ resize(0); }
template<class Real>
void WignerTransform<Real>::ScratchSpace::resize( const int& b )
{
	if( b!=bw )
	{
		int size=b*2;
		if( data )					FreePointer(data);
		if( coeffs )				FreePointer(coeffs);
		if( workspace_cx )			FreePointer(workspace_cx);
		if( workspace_cx2 )			FreePointer(workspace_cx2);
		if( workspace_re )			FreePointer( workspace_re );
		if( p )						{fftw_destroy_plan(p);}

		bw=0;
		data = NullPointer< fftw_complex >();
		workspace_cx = NullPointer< fftw_complex >();
		workspace_cx2 = NullPointer< fftw_complex >();
		coeffs = NullPointer< fftw_complex >();
		workspace_re = NullPointer< double >();
		p=0;

		if( b>0 )
		{
			bw=b;
#if 1
			data = AllocPointer< fftw_complex >( 8*bw*bw*bw );
			coeffs = AllocPointer< fftw_complex >( (4*bw*bw*bw-bw)/3 );
			workspace_cx = AllocPointer< fftw_complex >( 8*bw*bw*bw );
			workspace_cx2 = AllocPointer< fftw_complex >( 2*bw );
#else
			data=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*8*bw*bw*bw);
			coeffs=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(4*bw*bw*bw-bw)/3);
			workspace_cx=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*8*bw*bw*bw);
			workspace_cx2=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*8*bw*bw*bw);
#endif
			workspace_re = AllocPointer< double >( 12*size+size*bw );

			int na[2], inembed[2], onembed[2] ;
			int rank, howmany, istride, idist, ostride, odist ;
			howmany = size*size ;
			idist = size ;
			odist = size ;
			rank = 2 ;
			inembed[0] = size ;
			inembed[1] = size*size ;
			onembed[0] = size ;
			onembed[1] = size*size ;
			istride = 1 ;
			ostride = 1 ;
			na[0] = 1 ;
			na[1] = size ;

			p = fftw_plan_many_dft( rank, na, howmany,
				workspace_cx, inembed,
				istride, idist,
				data, onembed,
				ostride, odist,
				FFTW_FORWARD, FFTW_ESTIMATE );
		}
	}
}
/////////////////////
// WignerTransform //
/////////////////////
template<class Real>
void WignerTransform<Real>::resize( const int& resolution ){scratch.resize(resolution>>1);}
template<class Real>
int WignerTransform<Real>::InverseFourier( FourierKeySO3<Real>& key , RotationGrid<Real>& g )
{
	if( key.resolution()!=g.resolution() ) g.resize( key.resolution() );
	int bw=key.bandWidth() , sz=g.resolution();
	scratch.resize( bw );
	int m , idx=0;
	double temp;
	for( int i=0 ; i<bw ; i++ )
	{
		for( int j=0 ; j<bw ; j++ )
		{
			temp=1;
			if( i>j )
			{
				if( (i-j)%2 ){temp*=-1;}
				m=i;
			}
			else
			{
				if((j-i)%2){temp*=-1;}
				m=j;
			}
			// m = max( i , j )
			int h = bw-m;
			for( int k=0 ; k<h ; k++ )
			{
				double n=sqrt(8.0*PI*PI/(2*(m+k)+1));
				scratch.coeffs[idx][0] = key(m+k,i,j).r*n*temp;
				scratch.coeffs[idx][1] = key(m+k,i,j).i*n*temp;
				idx++;
			}
		}
		for( int j=bw-1 ; j>0 ; j-- )
		{
			if(j%2==0){temp=1;}
			else{temp=-1;}
			if( i>j )
			{
				if( (i-j)%2 ) temp*=-1;
				m=i;
			}
			else
			{
				if((j-i)%2){temp*=-1;}
				m=j;
			}
			int h=bw-m;
			for( int k=0 ; k<h ; k++ )
			{
				double n=sqrt(8.0*PI*PI/(2*(m+k)+1));
				scratch.coeffs[idx][0]=key(m+k,i,-j).r*n*temp;
				scratch.coeffs[idx][1]=key(m+k,i,-j).i*n*temp;
				idx++;
			}
		}
	}
	for( int i=bw-1 ; i>0 ; i-- )
	{
		for( int j=0 ; j<bw ; j++ )
		{
			if(i%2==0){temp=1;}
			else{temp=-1;}
			if( i>j )
			{
				if((i-j)%2){temp*=-1;}
				m=i;
			}
			else
			{
				if((j-i)%2){temp*=-1;}
				m=j;
			}
			int h=bw-m;
			for( int k=0 ; k<h ; k++ )
			{
				double n=sqrt(8.0*PI*PI/(2*(m+k)+1));
				scratch.coeffs[idx][0]= key(m+k,i,-j).r*n*temp;
				scratch.coeffs[idx][1]=-key(m+k,i,-j).i*n*temp;
				idx++;
			}
		}
		for( int j=bw-1 ; j>0 ; j-- )
		{
			if((i+j)%2==0){temp=1;}
			else{temp=-1;}
			if( i>j )
			{
				if((i-j)%2){temp*=-1;}
				m=i;
			}
			else
			{
				if((j-i)%2){temp*=-1;}
				m=j;
			}
			int h=bw-m;
			for( int k=0 ; k<h ; k++ )
			{
				double n=sqrt(8.0*PI*PI/(2*(m+k)+1));
				scratch.coeffs[idx][0]= key(m+k,i,j).r*n*temp;
				scratch.coeffs[idx][1]=-key(m+k,i,j).i*n*temp;
				idx++;
			}
		}
	}
	Inverse_SO3_Naive_fftw( bw , scratch.coeffs , scratch.data , scratch.workspace_cx , scratch.workspace_cx2 , scratch.workspace_re , &scratch.p , 1 );
	idx=0;
	for( int i=0 ; i<2*bw ; i++ ) for( int j=0 ; j<2*bw ; j++ ) for( int k=0 ; k<2*bw ; k++ ) g.values[idx++]=Real(scratch.data[2*bw*i+4*bw*bw*j+k][0]);
	return 1;
}
