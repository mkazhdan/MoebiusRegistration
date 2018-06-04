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

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>
#include <fstream>
#include <random>
#include <functional>

#include "Misha/Geometry.h"
#include "Misha/Algebra.h"
#include "Misha/Ply.h"
#include "Misha/CmdLineParser.h"
#include "Misha/Timer.h"
#include "Misha/SphericalGeometry.h"


cmdLineParameter< char* > In( "in" ) , Out( "out" );

void Usage( const char* ex )
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input mesh>\n" , In.name );
	printf( "\t[--%s <output spherical mesh>]\n" , Out.name );
}

cmdLineReadable* params[] = { &In , &Out , NULL };

int main( int argc , char* argv[] )
{
	cmdLineParse( argc-1 , argv+1 , params );
	if( !In.set )
	{
		Usage( argv[0] );
		return EXIT_FAILURE;
	}

	SphericalGeometry::Mesh< float > mesh;
	std::vector< Point3D< float > > vertices , colors;
	mesh.read( In.value , vertices , colors );

	std::vector< PlyColorVertex< float > > _vertices( vertices.size() );
	for( int i=0 ; i<_vertices.size() ; i++ ) _vertices[i].point =  mesh.vertices[i] , _vertices[i].color = colors[i];
	if( Out.set ) PlyWriteTriangles( Out.value , _vertices , mesh.triangles , PlyColorVertex< float >::WriteProperties , PlyColorVertex< float >::WriteComponents , PLY_BINARY_NATIVE );
	return EXIT_SUCCESS;
}
