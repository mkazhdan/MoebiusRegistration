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

#ifndef TIMER_INCLUDED
#define TIMER_INCLUDED
#include <sys/timeb.h>
#if !defined( WIN32 ) && !defined( _WIN64 )
#include <sys/time.h>
#endif // !WIN32 && !_WIN64
#if defined( WIN32 ) || defined( _WIN64 )
struct Timer
{
	static double Time( void )
	{
		struct _timeb t;
		_ftime( &t );
		return double(t.time)+double(t.millitm)/1000.0;
	}
	struct _timeb t;
	Timer( void ){ _ftime( &t ); }
	void reset( void ){ _ftime( &t ); }
	double elapsed( void ) const
	{
		struct _timeb _t;
		_ftime( &_t );
		return (double)(_t.time-t.time)+(double)(_t.millitm-t.millitm)/1000.;
	}
};
#else // !WIN32 && !_WIN64
struct Timer
{
	static double Time( void )
	{
		struct timeval t;
		gettimeofday(&t,NULL);
		return t.tv_sec+(double)t.tv_usec/1000000;
	}
	struct timeval t;
	Timer( void ){ gettimeofday( &t , NULL ); }
	void reset( void ){ gettimeofday( &t , NULL ); }
	double elapsed( void ) const
	{
		struct timeval _t;
		gettimeofday( &_t , NULL );
		return (double)(_t.tv_sec-t.tv_sec)+(double)(_t.tv_usec-t.tv_usec)/1000000;
	}
};
#endif // WIN32 || WIN64
#endif // TIMER_INCLUDED