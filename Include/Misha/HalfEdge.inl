/*
Copyright (c) 2012, Michael Kazhdan and Julian Panetta
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
#include <unordered_map>
#include "HalfedgeDictionary.hh"

template< class VertexData , class EdgeData , class FacetData >
const typename HalfEdgeMesh< VertexData , EdgeData , FacetData >::E*
HalfEdgeMesh< VertexData , EdgeData , FacetData >::_previous
(
 const E* e
) const
{
    const E* _e = e;
    size_t idx = e - &halfedges[0];
    while( _e->_nextEdgeIndex!=idx ) _e = &halfedges[_e->_nextEdgeIndex];
    return _e;
}
template< class VertexData , class EdgeData , class FacetData >
typename HalfEdgeMesh< VertexData , EdgeData , FacetData >::E*
HalfEdgeMesh< VertexData , EdgeData , FacetData >::_previous
( 
 E* e
)
{
    E* _e = e;
    int idx = e - &halfedges[0];
    while( _e->_nextEdgeIndex!=idx ) _e = &halfedges[_e->_nextEdgeIndex];
    return _e;
}

template< class VertexData , class EdgeData , class FacetData >
inline const typename HalfEdgeMesh< VertexData , EdgeData , FacetData >::E*
HalfEdgeMesh< VertexData , EdgeData , FacetData >::_next
(
 const E* e
) const
{
    return &halfedges[e->_nextEdgeIndex];
}
template< class VertexData , class EdgeData , class FacetData >
inline typename HalfEdgeMesh< VertexData , EdgeData , FacetData >::E*
HalfEdgeMesh< VertexData , EdgeData , FacetData >::_next
( 
 E* e
)
{
    return &halfedges[e->_nextEdgeIndex];
}
template< class VertexData , class EdgeData , class FacetData >
inline const typename HalfEdgeMesh< VertexData , EdgeData , FacetData >::E*
HalfEdgeMesh< VertexData , EdgeData , FacetData >::_opposite
(
 const E* e
) const
{
    return &halfedges[e->_oppositeEdgeIndex];
}
template< class VertexData , class EdgeData , class FacetData >
inline typename HalfEdgeMesh< VertexData , EdgeData , FacetData >::E*
HalfEdgeMesh< VertexData , EdgeData , FacetData >::_opposite
( 
 E* e
)
{
    return &halfedges[e->_oppositeEdgeIndex];
}
template< class VertexData , class EdgeData , class FacetData >
inline const typename HalfEdgeMesh< VertexData , EdgeData , FacetData >::F*
HalfEdgeMesh< VertexData , EdgeData , FacetData >::_facet
(
 const E* e
) const
{
    return &facets[e->_facetIndex];
}
template< class VertexData , class EdgeData , class FacetData >
inline typename HalfEdgeMesh< VertexData , EdgeData , FacetData >::F*
HalfEdgeMesh< VertexData , EdgeData , FacetData >::_facet
( 
 E* e
)
{
    return &facets[e->_facetIndex];
}

template< class VertexData , class EdgeData , class FacetData >
inline const typename HalfEdgeMesh< VertexData , EdgeData , FacetData >::E*
HalfEdgeMesh< VertexData , EdgeData , FacetData >::_halfedge
(
 const F* f
) const
{
    return &halfedges[f->_edgeIndex];
}
template< class VertexData , class EdgeData , class FacetData >
inline typename HalfEdgeMesh< VertexData , EdgeData , FacetData >::E*
HalfEdgeMesh< VertexData , EdgeData , FacetData >::_halfedge
( 
 F* f
)
{
    return &halfedges[f->_edgeIndex];
}

template< class VertexData , class EdgeData , class FacetData >
inline const typename HalfEdgeMesh< VertexData , EdgeData , FacetData >::E*
HalfEdgeMesh< VertexData , EdgeData , FacetData >::_out_halfedge
(
 const V* v
) const
{
    return &halfedges[v->_outEdgeIndex];
}
template< class VertexData , class EdgeData , class FacetData >
inline typename HalfEdgeMesh< VertexData , EdgeData , FacetData >::E*
HalfEdgeMesh< VertexData , EdgeData , FacetData >::_out_halfedge
( 
 V* v
)
{
    return &halfedges[v->_outEdgeIndex];
}

template< class VertexData , class EdgeData , class FacetData >
inline const typename HalfEdgeMesh< VertexData , EdgeData , FacetData >::V*
HalfEdgeMesh< VertexData , EdgeData , FacetData >::_end_vertex
(
 const E* e
) const
{
    return &vertices[e->_endVertexIndex];
}
template< class VertexData , class EdgeData , class FacetData >
inline typename HalfEdgeMesh< VertexData , EdgeData , FacetData >::V*
HalfEdgeMesh< VertexData , EdgeData , FacetData >::_end_vertex
( 
 E* e
)
{
    return &vertices[e->_endVertexIndex];
}

template< class  VertexData , class  EdgeData , class  FacetData , class _VertexData , class _EdgeData , class _FacetData >
void CopyTopology( const HalfEdgeMesh< VertexData , EdgeData , FacetData >& inMesh , HalfEdgeMesh< _VertexData , _EdgeData , _FacetData >& outMesh )
{
#if TEST_HE_MESH_VALIDITY
    for( int i=0 ; i<inMesh.halfedge_size() ; i++ ) if( !inMesh.is_valid( inMesh.halfedge(i) ) ) fprintf( stderr , "Invalid in halfedge %d\n" , i ) , exit( 0 );
    for( int i=0 ; i<  inMesh.vertex_size() ; i++ ) if( !inMesh.is_valid(   inMesh.vertex(i) ) ) fprintf( stderr , "Invalid in vertex %d\n"   , i ) , exit( 0 );
    for( int i=0 ; i<   inMesh.facet_size() ; i++ ) if( !inMesh.is_valid(    inMesh.facet(i) ) ) fprintf( stderr , "Invalid in facet  %d\n"   , i ) , exit( 0 );
#endif // TEST_HE_MESH_VALIDITY

    outMesh.vertices.resize ( inMesh.vertex_size() );
    outMesh.halfedges.resize( inMesh.halfedge_size() );
    outMesh.facets.resize   ( inMesh.facet_size() );
    for( int i=0 ; i< outMesh.vertices.size() ; i++ ) outMesh.vertices[i]._outEdgeIndex = inMesh.vertices[i]._outEdgeIndex;
    for( int i=0 ; i<   outMesh.facets.size() ; i++ )   outMesh.facets[i]._edgeIndex    = inMesh.facets[i]._edgeIndex;
    for( int i=0 ; i<outMesh.halfedges.size() ; i++ )
    {
        const typename HalfEdgeMesh<  VertexData ,  EdgeData ,  FacetData >::E* in  = &inMesh.halfedges[i];
              typename HalfEdgeMesh< _VertexData , _EdgeData , _FacetData >::E* out = &outMesh.halfedges[i];
        out->_nextEdgeIndex = in->_nextEdgeIndex;
        out->_oppositeEdgeIndex = in->_oppositeEdgeIndex;
        out->_facetIndex = in->_facetIndex;
        out->_endVertexIndex = in->_endVertexIndex;
    }

#if TEST_HE_MESH_VALIDITY
    for( int i=0 ; i<outMesh.halfedges.size() ; i++ ) if( !outMesh.is_valid( &outMesh.halfedges[i] ) ) fprintf( stderr , "Invalid out halfedge %d\n" , i ) , exit( 0 );
    for( int i=0 ; i< outMesh.vertices.size() ; i++ ) if( !outMesh.is_valid(  &outMesh.vertices[i] ) ) fprintf( stderr , "Invalid out vertex %d\n"   , i ) , exit( 0 );
    for( int i=0 ; i<   outMesh.facets.size() ; i++ ) if( !outMesh.is_valid(    &outMesh.facets[i] ) ) fprintf( stderr , "Invalid out facet  %d\n"   , i ) , exit( 0 );
#endif // TEST_HE_MESH_VALIDITY
}

template< class VertexData , class EdgeData , class FacetData >
template< class _VertexData , class _EdgeData , class _FacetData >
HalfEdgeMesh< VertexData , EdgeData , FacetData >&
HalfEdgeMesh< VertexData , EdgeData , FacetData >::operator = ( const HalfEdgeMesh< _VertexData , _EdgeData , _FacetData >& mesh )
{
    CopyTopology( mesh , *this );
    return *this;
}

template< class VertexData , class EdgeData , class FacetData >
template< class Polygon >
bool HalfEdgeMesh< VertexData , EdgeData , FacetData >::_SetHalfEdgeData(
        const std::vector< Polygon >& polygons, std::vector< V >& vertices, 
        std::vector< E >& halfedges, std::vector< F >& facets,
        unsigned int (*PSize)(const Polygon &) ,
		void (*startFunction)( void ) , void (*endFunction )( void )
		)
{
	if( startFunction ) startFunction();

	unsigned int maxVIndex = 0;
    for (unsigned int i = 0; i < polygons.size(); i++)  {
        for (unsigned int v = 0; v < PSize(polygons[i]); v++)   {
            if (polygons[i][v] > maxVIndex) {
                maxVIndex = polygons[i][v];
            }
        }
    }
    vertices.resize( maxVIndex+1 );
    facets.resize( polygons.size() );

    int eCount = 0;
    for( size_t i=0 ; i<polygons.size() ; i++ ) eCount += PSize( polygons[i] );
    halfedges.resize( eCount );
    HalfedgeDictionary heDict((int)vertices.size(), (int)halfedges.size());

    eCount = 0;
    // First Pass (add all the half-edges except those on the boundary)
    for( size_t i=0 ; i<polygons.size() ; i++ )
    {
        int startEdgeIndex = eCount;
        int pSize = PSize( polygons[i] );
        if (pSize < 3) fprintf(stderr, "[Warning] Could not add polygon with %d vertices: %d\n", pSize, (int) i);
        else
        {
            int thisVertexIndex = polygons[i][0];
            for( int v=0 ; v<pSize ; v++ ) 
            {
                E& edge = halfedges[eCount];
                edge._facetIndex = (int)i;
                int v1 = (v+1)%pSize;
                int nextVertexIndex = polygons[i][v1];
                edge._endVertexIndex = nextVertexIndex;
                edge._nextEdgeIndex = startEdgeIndex + v1;
                heDict.insert(thisVertexIndex, nextVertexIndex, eCount);
                vertices[ thisVertexIndex ]._outEdgeIndex = eCount;

                // Force the face to point to the last edge
                facets[i]._edgeIndex = eCount;
                thisVertexIndex = nextVertexIndex;
                eCount++;
            }
        }
    }

    // Second Pass. For 1. filling in information of opposite edges and 
    //                  2. linking boundary-vertices to the "opposite edges" of boundary-edges
    //                  [Note] When we say an edge is a boundary-edge, we mean its facePointer is NULL (negative).
    for( size_t i=0 ; i<polygons.size() ; i++ )
    {
        int pSize = PSize( polygons[i] );
        for( int v=0 ; v<pSize ; v++ )
        {
            int thisVertexIndex = polygons[i][ v         ];
            int nextVertexIndex = polygons[i][(v+1)%pSize];
            size_t thisEdgeIndex = facets[i]._edgeIndex + v - pSize + 1;
            // Try to find the opposite edge (next->this)
            unsigned int heIdx = heDict.halfedge(nextVertexIndex, thisVertexIndex);
            if (heIdx != HalfedgeDictionary::INVALID_HALFEDGE)  {
                halfedges[thisEdgeIndex]._oppositeEdgeIndex = heIdx;
            }
            else {
                // Create a new boundary edge opposite thisEdge
                E edge;
                edge._facetIndex                            = -1;
                edge._nextEdgeIndex                         = -1; // to be set later
                edge._endVertexIndex                        = thisVertexIndex;
                edge._oppositeEdgeIndex                     = (int)thisEdgeIndex;
                halfedges[thisEdgeIndex]._oppositeEdgeIndex = (int)halfedges.size();
                halfedges.push_back(edge);
            }
        }
    }

    // Link each boundary edge to its next edge
    // (Boundary edges are in halfEdges[eCount:end])
    for( size_t i=eCount ; i<halfedges.size() ; i++ )
    {
        size_t currentEdgeIndex = halfedges[i]._oppositeEdgeIndex;
        // Circulate CCW around endVertex until we find the other outEdge on the
        // boundary
        // (CCW circulation around a vertex is the operation:
        //  currentEdge->prev->opposite)
        while (halfedges[currentEdgeIndex]._facetIndex != -1)   {
            // Walk CCW around polygon to find the previous half edge
            int startEdgeIndex = (int)currentEdgeIndex;
            while (halfedges[ currentEdgeIndex ]._nextEdgeIndex != startEdgeIndex)
                currentEdgeIndex = halfedges[ currentEdgeIndex ]._nextEdgeIndex;
            // Advance to the previous edge's opposite
            currentEdgeIndex = halfedges[ currentEdgeIndex ]._oppositeEdgeIndex;
        }
        halfedges[i]._nextEdgeIndex = (int)currentEdgeIndex;
    }

	if( endFunction ) endFunction();
    return true;
}

template< class VertexData , class EdgeData , class FacetData >
void
HalfEdgeMesh< VertexData , EdgeData , FacetData >::_swap_vertices
(
 V* v1 ,
 V* v2
)
{
    if( v1==v2 ) return;

    // Adjust the pointers from the halfedges into the vertices
    {
        E *e , *s;
        e = s = _out_halfedge( v1 );
        do
        {
            e = _opposite( e );
            e->_endVertexIndex = v2 - &vertices[0];
            e = _next( e );
        }
        while( e!=s );

        e = s = _out_halfedge( v2 );
        do
        {
            e = _opposite( e );
            e->_endVertexIndex = v1 - &vertices[0];
            e = _next( e );
        }
        while( e!=s );
    }

    // Copy the vertex data
    {
		if( swap_vertices_function )
		{
			swap_vertices_function( v1 , v2 );
			std::swap( v1->_outEdgeIndex , v2->_outEdgeIndex );
		}
		else
		{
			V temp = *v1;
			*v1  = *v2;
			*v2 = temp;
		}
    }
}
template< class VertexData , class EdgeData , class FacetData >
void HalfEdgeMesh< VertexData , EdgeData , FacetData >::_remove_vertex ( int v)
{
    if( v<0 || v>=vertices.size() ) return;
    if( v==vertices.size()-1 )
    {
        vertices.pop_back();
        return;
    }

    V* _v1 = &vertices[vertices.size()-1];
    V* _v2 = &vertices[v];

    // Adjust the tip pointers for halfedges incident on the vertex moved to
    // index v.
    {
        E *e , *s;
        e = s = _out_halfedge( _v1 );
        do
        {
            e = _opposite( e );
            e->_endVertexIndex = _v2 - &vertices[0];
            e = _next( e );
        }
        while( e!=s );
    }

    // Copy the vertex data
    *_v2 = *_v1;
    vertices.pop_back();
}

template< class VertexData , class EdgeData , class FacetData >
void
HalfEdgeMesh< VertexData , EdgeData , FacetData >::_swap_facets
(
 F* f1 ,
 F* f2 
)
{
    if( f1==f2 || !f1 || !f2 ) return;

    E *e , *next;
    // Adjust the pointers from the halfedges into the facets
    {
        next = e = _halfedge( f1 );
        do
        {
            next->_facetIndex = f2 - &facets[0];
            next = _next( next );
        }
        while( next!=e );
        next = e = _halfedge( f2 );
        do
        {
            next->_facetIndex = f1 - &facets[0];
            next = _next( next );
        }
        while( next!=e );
    }

    // Copy the facet data
    {
		if( swap_facets_function )
		{
			swap_facets_function( f1 , f2 );
			std::swap( f1->_edgeIndex , f2->_edgeIndex );
		}
		else
		{
			F temp = *f1;
			*f1 = *f2;
			*f2 = temp;
		}
    }
}
template< class VertexData , class EdgeData , class FacetData >
void
HalfEdgeMesh< VertexData , EdgeData , FacetData >::_remove_facet
(
 int f 
)
{
    if( f<0 || f>=facets.size() ) return;
    if( f==facets.size()-1 )
    {
        facets.pop_back();
        return;
    }
    F* _f1 = &facets[facets.size()-1];
    F* _f2 = &facets[f];

    E *e , *next;
    // Adjust the pointers from the halfedges into the facets
    {
        next = e = _halfedge( _f1 );
        do {
            next->_facetIndex = _f2 - &facets[0];
            next = _next( next );
        } while( next!=e );
    }

    // Copy the facet data
    *_f2 = *_f1;
    facets.pop_back();
}
template< class VertexData , class EdgeData , class FacetData >
void
HalfEdgeMesh< VertexData , EdgeData , FacetData >::_swap_halfedges
(
 E* e1 ,
 E* e2
)
{
    if( e1==e2 ) return;

    E *o1 = _opposite( e1 ) , *o2 = _opposite( e2 );

    // Adjust the pointers from the facets into the halfedges
    {
        E* next;
        int idx1 = e1 - &halfedges[0] , idx2 = e2 - &halfedges[0];
        next = e1;
        do
        {
            if( next->_facetIndex!=-1 && facets[next->_facetIndex]._edgeIndex==idx1 ) facets[next->_facetIndex]._edgeIndex = idx2;
            next = _next( next );
        }
        while( next!=e1 );
        next = e2;
        do
        {
            if( next->_facetIndex!=-1 && facets[next->_facetIndex]._edgeIndex==idx2 ) facets[next->_facetIndex]._edgeIndex = idx1;
            next = _next( next );
        }
        while( next!=e2 );
    }

    // Adjust the pointers from the vertices into the halfedges
    {
        int idx1 = e1 - &halfedges[0] , idx2 = e2 - &halfedges[0];
        if( vertices[o1->_endVertexIndex]._outEdgeIndex==idx1 ) vertices[o1->_endVertexIndex]._outEdgeIndex = idx2;
        if( vertices[o2->_endVertexIndex]._outEdgeIndex==idx2 ) vertices[o2->_endVertexIndex]._outEdgeIndex = idx1;
    }

    // Adjust the pointers from the halfedges into the halfedges
    {
        E * p1 = _previous( e1 ) , *p2 = _previous( e2 );

        int idx1 = e1 - &halfedges[0] , idx2 = e2 - &halfedges[0];
        o1->_oppositeEdgeIndex = idx2 , o2->_oppositeEdgeIndex = idx1;
        p1->_nextEdgeIndex     = idx2 , p2->_nextEdgeIndex     = idx1;
    }

    // Copy the halfedge data
    {
		if( swap_halfedges_function )
		{
			swap_halfedges_function( e1 , e2 );
			std::swap( e1->_endVertexIndex    , e2->_endVertexIndex    );
			std::swap( e1->_facetIndex        , e2->_facetIndex        );
			std::swap( e1->_nextEdgeIndex     , e2->_nextEdgeIndex     );
			std::swap( e1->_oppositeEdgeIndex , e2->_oppositeEdgeIndex );
		}
		else
		{
			E temp = *e1;
			*e1 = *e2;
			*e2 = temp;
		}
    }
}
template< class VertexData , class EdgeData , class FacetData >
void
HalfEdgeMesh< VertexData , EdgeData , FacetData >::_remove_halfedge( int e )
{
    if( e<0 || e>=halfedges.size() ) return;
    if( e==halfedges.size()-1 )
    {
        halfedges.pop_back();
        return;
    }
    E* _e1 = &halfedges[halfedges.size()-1];
    E* _e2 = &halfedges[e];
    E* o1 = _opposite( _e1 );
    // Adjust the pointers from the facets into the halfedges
    {
        E* next;
        next = _e1;
        int idx1 = _e1 - &halfedges[0] , idx2 = _e2 - &halfedges[0];
        do
        {
            if( next->_facetIndex!=-1 && facets[next->_facetIndex]._edgeIndex==idx1 ) facets[next->_facetIndex]._edgeIndex = idx2;
            next = _next( next );
        }
        while( next!=_e1 );
    }

    // Adjust the pointers from the vertices into the halfedges
    {
        if( vertices[o1->_endVertexIndex]._outEdgeIndex==(_e1-&halfedges[0]) ) vertices[o1->_endVertexIndex]._outEdgeIndex = _e2-&halfedges[0];
    }

    // Adjust the pointers from the halfedges into the halfedges
    {
        E * p1 = _previous( _e1 );

        o1->_oppositeEdgeIndex = _e2 - &halfedges[0];
        p1->_nextEdgeIndex     = _e2 - &halfedges[0];
    }

    // Copy the halfedge data
    *_e2 = *_e1;

    halfedges.pop_back();
}

template< class VertexData , class EdgeData , class FacetData >
void
HalfEdgeMesh< VertexData , EdgeData , FacetData >::print_facet_vertices( Facet_const_handle f , bool new_line ) const
{
    if( !f )
    {
        printf( "NULL" );
        if( new_line ) printf( "\n" );
        return;
    }

    E *e = _halfedge( f ) , *next = _halfedge( f );
    printf( "[%d] {%d %d}" , index( f ) , index( f.opposite().end_vertex() ) , index( f.halfedge().end_vertex() ) );
    do {
        printf( " %d" , next->_endVertexIndex );
        next = _next( next );
    } while( next!=e );
    if( new_line ) printf( "\n" );
}
template< class VertexData , class EdgeData , class FacetData >
void
HalfEdgeMesh< VertexData , EdgeData , FacetData >::print_facet_halfedges( Facet_const_handle f , bool new_line ) const
{
    if( !f )
    {
        printf( "NULL" );
        if( new_line ) printf( "\n" );
        return;
    }
    E *e = _halfedge( f ) , *next = _halfedge( f );
    printf( "[%d] {%d %d}" , index( f ) , index( f.halfedge() ) , index( f.halfedge().opposite() ) );
    do {
        printf( " %d" , next-&halfedges[0] );
        next = _next( next );
    }
    while( next!=e );
    if( new_line ) printf( "\n" );
}
template< class VertexData , class EdgeData , class FacetData >
void
HalfEdgeMesh< VertexData , EdgeData , FacetData >::print_halfedge( Halfedge_const_handle e , bool new_line ) const
{
    printf( "{%d} %d %d (%d)" , index( e ) , index( e.opposite().end_vertex() ) , index( e.end_vertex() ) );
    if( new_line ) printf( "\n" );
}
template< class VertexData , class EdgeData , class FacetData >
typename HalfEdgeMesh< VertexData , EdgeData , FacetData >::Facet_handle
HalfEdgeMesh< VertexData , EdgeData , FacetData >::merge_facets
(
 Halfedge_handle e,
 void (*merge)( const FacetData* , FacetData* )
)
{
    if( !e ) return Facet_handle();
    E *e1 = &halfedges[e._idx] , *e2 = _opposite( e1 );
    F *f1 = _facet( e1 ) , *f2 = _facet( e2 );
    if( !f1 || !f2 ) return Facet_handle();

    if( f1==f2 )
    {
        int idx1 = e1-&halfedges[0] , idx2 = e2-&halfedges[0];
        if( e1->_nextEdgeIndex==idx2 && e2->_nextEdgeIndex==idx1 )
        {
            // Although we can remove this edge, it's existence indicates that we have a facet with two diconnected halfedge loops (almost)
            fprintf( stderr , "[Warning] Isolated edge in winged-edge mesh\n" );
            return Facet_handle();
        }
        else if( e1->_nextEdgeIndex==idx2 || e2->_nextEdgeIndex==idx1 )
        {
            E *e;
            V *v;
            _swap_halfedges( e1 , &halfedges[ halfedges.size()-2 ] );
            if( (e2-&halfedges[0])==halfedges.size()-2 ) e2 = e1;
            e1 = &halfedges[ halfedges.size()-2 ];
            _swap_halfedges( e2 , &halfedges[ halfedges.size()-1 ] );
            e2 = &halfedges[ halfedges.size()-1 ];
            int idx1 = e1-&halfedges[0] , idx2 = e2-&halfedges[0];
            if( e1->_nextEdgeIndex==idx2 )
            {
                v = &vertices[e1->_endVertexIndex];                                     // The vertex that will be removed
                vertices[e2->_endVertexIndex]._outEdgeIndex = e2->_nextEdgeIndex;       // Make sure the vertex doesn't point to the removed edge
                f1->_edgeIndex = e2->_nextEdgeIndex;                                    // Make sure the facet doesn't point to the removed edge
                _previous( e1 )->_nextEdgeIndex = e2->_nextEdgeIndex;                   // Adjust the edge pointers
            }
            else
            {
                v = &vertices[e2->_endVertexIndex];                                     // The vertex that will be removed
                vertices[e1->_endVertexIndex]._outEdgeIndex = e1->_nextEdgeIndex;       // Make sure the vertex doesn't point to the removed edge
                f1->_edgeIndex = e1->_nextEdgeIndex;                                    // Make sure the facet doesn't point to the removed edge
                _previous( e2 )->_nextEdgeIndex = e1->_nextEdgeIndex;               // Adjust the edge pointers
            }
            _swap_vertices( v , & vertices[ vertices.size()-1 ] ) , v  = &vertices[ vertices.size()-1 ];
            vertices.pop_back();
            halfedges.pop_back();
            halfedges.pop_back();

            return Facet_handle( f1-&facets[0] , this );
        }
        else return Facet_handle();
    }
    // Swap f1 and f2 into the back of the facet array
    {
        _swap_facets( f1 , &facets[ facets.size()-2 ] );
        if( (f2-&facets[0])==facets.size()-2 ) f2 = f1;
        f1 = &facets[ facets.size()-2 ];
        _swap_facets( f2 , &facets[ facets.size()-1 ] );
        f2 = &facets[ facets.size()-1 ];
    }
    // Swap e1 and e2 into the back of the halfedge array
    {
        _swap_halfedges( e1 , &halfedges[ halfedges.size()-2 ] );
        if( (e2-&halfedges[0])==halfedges.size()-2 ) e2 = e1;
        e1 = &halfedges[ halfedges.size()-2 ];
        _swap_halfedges( e2 , &halfedges[ halfedges.size()-1 ] );
        e2 = &halfedges[ halfedges.size()-1 ];
    }

    // Make sure that no vertices point to the edges that will be removed
    {
        int idx1 = e1-&halfedges[0] , idx2=e2-&halfedges[0];
        if( vertices[e2->_endVertexIndex]._outEdgeIndex==idx1 ) vertices[e2->_endVertexIndex]._outEdgeIndex = e2->_nextEdgeIndex;
        if( vertices[e1->_endVertexIndex]._outEdgeIndex==idx2 ) vertices[e1->_endVertexIndex]._outEdgeIndex = e1->_nextEdgeIndex;
    }

    // Adjust the edge connectivity of the new face and choose a representative edge for the new face
    {
        E *p1 = _previous( e1 ) , *p2 = _previous( e2 );
        p1->_nextEdgeIndex = e2->_nextEdgeIndex;
        p2->_nextEdgeIndex = e1->_nextEdgeIndex;
        e = Halfedge_handle( p1-&halfedges[0] , this );
    }

    // Copy the new face into f1
    {
		f1->_edgeIndex = e._idx;
		if( merge ) merge( f2 , f1 );
    }

    // Set the face pointers of the edges along the new face
    {
        E* _e   = &halfedges[e._idx];
        E* next = &halfedges[e._idx];
        do
        {
            next->_facetIndex = f1-&facets[0];
            next = _next( next );
        }
        while( next!=_e );
    }

    // Remove facet f2
    {
        facets.pop_back();
    }

    // Remove edges e1 and e2
    {
        halfedges.pop_back();
        halfedges.pop_back();
    }
    return Facet_handle( f1-&facets[0] , this );
}
template< class VertexData , class EdgeData , class FacetData >
bool
HalfEdgeMesh< VertexData , EdgeData , FacetData >::merge_facets_full
(
 bool (*toMerge)( Halfedge_const_handle ) ,
 void (*mergeFacetsFunction)( const FacetData* , const FacetData* , FacetData* )
)
{
    bool merged = false;
    std::vector< bool > validEdges ( halfedges.size() , true );
    std::vector< bool > validFacets(  2*facets.size() , true );
    std::vector< bool > validVerts (  vertices.size() , true );
    for( int i=0 ; i<halfedges.size() ; i++ )
    {
        E* e = &halfedges[i];
        E *e1 = e , *e2 = _opposite( e1 );
        if( validEdges[i] && toMerge( e ) )
        {
            F *f1 = _facet( e1 ) , *f2 = _facet( e2 );
            if( !f1 || !f2 ) continue;
            if( f1==f2 )
            {
                int idx1 = e1-&halfedges[0] , idx2 = e2-&halfedges[0];
                if( e1->_nextEdgeIndex==idx2 && e2->_nextEdgeIndex==idx1 )
                {
                    // Although we can remove this edge, it's existence indicates that we have a facet with two disconnected halfedge loops (almost)
                    fprintf( stderr , "[Warning] Isolated edge in winged-edge mesh\n" );
                    continue;
                }
                else if( e1->_nextEdgeIndex==idx2 || e2->_nextEdgeIndex==idx1 )
                {
                    V *v;
                    if( e1->_nextEdgeIndex==idx2 )
                    {
                        v = _end_vertex( e1 );                                                  // The vertex that will be removed
                        vertices[e2->_endVertexIndex]._outEdgeIndex = e2->_nextEdgeIndex;       // Make sure the vertex doesn't point to the removed edge
                        f1->_edgeIndex = e2->_nextEdgeIndex;                                    // Make sure the facet doesn't point to the removed edge
                        _previous( e1 )->_nextEdgeIndex = e2->_nextEdgeIndex;                   // Adjust the edge pointers
                    }
                    else
                    {
                        v = _end_vertex( e1 );                                                  // The vertex that will be removed
                        vertices[e1->_endVertexIndex]._outEdgeIndex = e1->_nextEdgeIndex;       // Make sure the vertex doesn't point to the removed edge
                        f1->_edgeIndex = e1->_nextEdgeIndex;                                    // Make sure the facet doesn't point to the removed edge
                        _previous( e2 )->_nextEdgeIndex = e1->_nextEdgeIndex;                   // Adjust the edge pointers
                    }
                    validVerts[ index(v) ] = validEdges[ index(e1) ] = validEdges[ index(e2) ] = false;
                    merged = true;
                }
            }
            else
            {
                // Make sure that no vertices point to the edges that will be removed
                {
                    int idx1 = e1-&halfedges[0] , idx2 = e2-&halfedges[0];
                    if( vertices[e2->_endVertexIndex]._outEdgeIndex==idx1 ) vertices[e2->_endVertexIndex]._outEdgeIndex = e2->_nextEdgeIndex;
                    if( vertices[e1->_endVertexIndex]._outEdgeIndex==idx2 ) vertices[e1->_endVertexIndex]._outEdgeIndex = e1->_nextEdgeIndex;
                }

                // Adjust the edge connectivity of the new face and choose a representative edge for the new face
                {
                    E *p1 = _previous( e1 ) , *p2 = _previous( e2 );
                    p1->_nextEdgeIndex = e2->_nextEdgeIndex;
                    p2->_nextEdgeIndex = e1->_nextEdgeIndex;
                    e = p1;
                }

                // Copy the new face into f1
                {
                    F newFace;
                    newFace._edgeIndex = e-&halfedges[0];
                    if( mergeFacetsFunction ) mergeFacetsFunction( f1 , f2 , &newFace );
                    *f1 = newFace;
                }

                // Set the face pointers of the edges along the new face
                {
                    E* next = e;
                    do
                    {
                        next->_facetIndex = f1-&facets[0];
                        next = _next( next );
                    }
                    while( next!=e );
                }

                validFacets[ index(f2) ] = validEdges[ index(e1) ] = validEdges[ index(e2) ] = false;
                merged = true;
            }
        }
    }
    for( int i=halfedges.size()-1 ; i>=0 ; i-- ) if(  !validEdges[i] ) _remove_halfedge( i );
    for( int i=   facets.size()-1 ; i>=0 ; i-- ) if( !validFacets[i] ) _remove_facet( i );
    for( int i= vertices.size()-1 ; i>=0 ; i-- ) if(  !validVerts[i] ) _remove_vertex( i );
    return merged;
}

template< class VertexData , class EdgeData , class FacetData >
bool
HalfEdgeMesh< VertexData , EdgeData , FacetData >::remove_facet
(
 Facet_handle f
)
{
    if( !f ) return false;
    // Swap fto the back of the facet array
    {
        _swap_facets( f , &facets[ facets.size()-1 ] );
        f._idx = facets.size()-1;
    }
    // Set all edges pointing to the facet to NULL

    {
        E* next = _halfedge( f );
        E* start = next;
        do
        {
            next->_facetIndex = -1;
            next = _next( next );
        }
        while( next!=start );
    }
    // Remove facet f
    {
        facets.pop_back();
    }

    return true;
}

template< class VertexData , class EdgeData , class FacetData >
typename HalfEdgeMesh< VertexData , EdgeData , FacetData >::Halfedge_handle
HalfEdgeMesh< VertexData , EdgeData , FacetData >::merge_halfedges
(
 Halfedge_handle e_handle ,
 void (*merge)( const EdgeData* , const EdgeData* , EdgeData* )
)
{
    E* e = &halfedges[e_handle._idx];
    E* o = _opposite( e );
    if( !e_handle ) return Halfedge_handle();                                                                               // Make sure the halfedge is valid
    if( e->_nextEdgeIndex==(e-&halfedges[0]) || halfedges[e->_nextEdgeIndex]._nextEdgeIndex==(e-&halfedges[0]) || halfedges[halfedges[e->_nextEdgeIndex]._nextEdgeIndex]._nextEdgeIndex==(e-&halfedges[0]) ) return Halfedge_handle();      // Make sure facets have at least four halfedges
    if( o->_nextEdgeIndex==(o-&halfedges[0]) || halfedges[o->_nextEdgeIndex]._nextEdgeIndex==(o-&halfedges[0]) || halfedges[halfedges[o->_nextEdgeIndex]._nextEdgeIndex]._nextEdgeIndex==(o-&halfedges[0]) ) return Halfedge_handle();      // Make sure facets have at least four halfedges
    if( o->_endVertexIndex==halfedges[e->_nextEdgeIndex]._endVertexIndex ) return Halfedge_handle();                        // Make sure you don't end up with an edge has the same vertex as both endpoints
    if( e->_nextEdgeIndex==e->_oppositeEdgeIndex ) return Halfedge_handle();                                        // Make sure we're not merging the halfedge with its opposite
    if( halfedges[e->_oppositeEdgeIndex]._nextEdgeIndex==e->_oppositeEdgeIndex || halfedges[halfedges[e->_oppositeEdgeIndex]._nextEdgeIndex]._nextEdgeIndex==e->_oppositeEdgeIndex ) return Halfedge_handle();
    if( (e-&halfedges[0])!=halfedges[halfedges[halfedges[e->_nextEdgeIndex]._oppositeEdgeIndex]._nextEdgeIndex]._oppositeEdgeIndex ) return Halfedge_handle();  // Make sure the the in-between vertex has valence two
    E *e2 = _next( e );
    E *o2 = _opposite( e2 );
    if( e==o2 ) return Halfedge_handle();
    if( e2->_nextEdgeIndex==(o2-&halfedges[0]) ) return Halfedge_handle();

    // Swap e2 and o2 into the back of the halfedge array
    {
        _swap_halfedges( e2 , &halfedges[ halfedges.size()-2 ] );
        if( (o2-&halfedges[0])==halfedges.size()-2 ) o2 = e2;
        e2 = &halfedges[ halfedges.size()-2 ];

        _swap_halfedges( o2 , &halfedges[ halfedges.size()-1 ] );
        if( (e2-&halfedges[0])==halfedges.size()-1 ) e2 = o2;
        o2 = &halfedges[ halfedges.size()-1 ];
    }

    E* o1 = _next( o2 );
    E* e1 = _opposite( o1 );
    E* n = _next( e2 );
    E* p = _previous( o2 );
    V* v = _end_vertex( e1 );

    // Swap v into the back of the vertex array
    {
        _swap_vertices( v , &vertices[ vertices.size()-1 ] );
        v = &vertices[ vertices.size()-1 ];
    }

    // Adjust any facets pointing to e2 or o2
    {
        if( e2->_facetIndex!=-1 && facets[e2->_facetIndex]._edgeIndex==(e2-&halfedges[0]) ) facets[e2->_facetIndex]._edgeIndex = e1-&halfedges[0];
        if( o2->_facetIndex!=-1 && facets[o2->_facetIndex]._edgeIndex==(o2-&halfedges[0]) ) facets[o2->_facetIndex]._edgeIndex = o1-&halfedges[0];
    }
    // Adjust any vertex pointing to e2 or o2
    {
        if( vertices[e2->_endVertexIndex]._outEdgeIndex==(o2-&halfedges[0]) ) vertices[e2->_endVertexIndex]._outEdgeIndex = o1-&halfedges[0];
    }

    // Adjust the vertex pointers
    {
        e1->_endVertexIndex = e2->_endVertexIndex;
    }

    // Adjust the edge connectivity of the merged edges
    {
        e1->_nextEdgeIndex = n-&halfedges[0];
        p->_nextEdgeIndex = o1-&halfedges[0];
    }
    // Copy the new edges into e1 and o1
    {
        E newHalfedge;
        newHalfedge = *e1;
        if( merge ) merge( e1 , e2 , &newHalfedge );
        *e1 = newHalfedge;

        newHalfedge = *o1;
        if( merge ) merge( o1 , o2 , &newHalfedge );
        *o1 = newHalfedge;
    }
    // Remove vertex v
    {
        vertices.pop_back();
    }

    // Remove edges e2 and o2
    {
        halfedges.pop_back();
        halfedges.pop_back();
    }
    return Halfedge_handle( e1-&halfedges[0] , this );
}

template< class VertexData , class EdgeData , class FacetData >
void
HalfEdgeMesh< VertexData , EdgeData , FacetData >::sort_facets
(
 int (*facet_compare_function)( const void* , const void* )
)
{
    qsort( &facets[0] , facets.size() , sizeof( F ) , facet_compare_function );
    for( int i=0 ; i<facets.size() ; i++ )
    {
        E *e , *next;
        next = e = _halfedge( &facets[i] );
        do
        {
            next->_facetIndex = i;
            next = _next( next );
        }
        while( next!=e );
    }
}
template< class VertexData , class EdgeData , class FacetData >
void
HalfEdgeMesh< VertexData , EdgeData , FacetData >::sort_vertices
(
 int (*vertex_compare_function)( const void* , const void* )
)
{
    qsort( &vertices[0] , vertices.size() , sizeof( V ) , vertex_compare_function );
    for( int i=0 ; i<vertices.size() ; i++ ) halfedges[halfedges[vertices[i]._outEdgeIndex]._oppositeEdgeIndex]._endVertexIndex = i;
}
template< class VertexData , class EdgeData , class FacetData >
bool
HalfEdgeMesh< VertexData , EdgeData , FacetData >::is_valid
(
 bool aggressive ,
 bool testEdges
)
const
{
    printf( "Testing mesh validity\n" );
    bool v , e , f , fullE;
    v = e = f = fullE = true;
    for( int i=0 ; i< vertices.size() ; i++ ) if( !is_valid(  &vertices[i] , aggressive ) ) fprintf( stderr , "Invalid vertex   %d\n" , i ) , v = false;
    for( int i=0 ; i<halfedges.size() ; i++ ) if( !is_valid( &halfedges[i] , aggressive ) ) fprintf( stderr , "Invalid halfedge %d\n" , i ) , e = false;
    for( int i=0 ; i<   facets.size() ; i++ ) if( !is_valid(    &facets[i] , aggressive ) ) fprintf( stderr , "Invalid facet    %d\n" , i ) , f = false;
    if( testEdges )
    {
        std::unordered_map< long long , int > eTable;
        for( int i=0 ; i<halfedges.size() ; i++ )
        {
            const E* e = &halfedges[i];
            std::pair< const V* , const V* > verts = e->vertices();
            long long key = HalfEdgeKey( index( verts.first ) , index( verts.second ) );
            if( eTable.find(key)!=eTable.end() )
            {
                fprintf( stderr , "redundant edge failure\n" );
                eTable[key]++;
                fullE = false;
            }
            else eTable[key] = 1;
        }
    }
    return v & e & f & fullE;
}

template< class VertexData , class EdgeData , class FacetData >
bool
HalfEdgeMesh< VertexData , EdgeData , FacetData >::is_valid
(
 Vertex_const_handle v ,
 bool aggressive
)
const
{
    if( index(v)<0 || index(v)>=vertices.size() ){ fprintf( stderr , "vertex index out of bounds\n" ) ; return false; }
    if( v._outEdgeIndex<0 || v._outEdgeIndex>=halfedges.size() ){ fprintf( stderr , "vertex edge index out of bounds\n" ) ; return false; }
    const E* e = _out_halfedge( v );
    const E* s = _out_halfedge( v );
    do
    {
        if( halfedges[e->_oppositeEdgeIndex]._endVertexIndex!=v._idx ){ fprintf( stderr , "vertex to edge to vertex failure\n" ) ; return false; }
        e = &halfedges[halfedges[e->_oppositeEdgeIndex]._nextEdgeIndex];
    }
    while( e!=s );
    return true;
}
template< class VertexData , class EdgeData , class FacetData >
bool
HalfEdgeMesh< VertexData , EdgeData , FacetData >::is_valid
(
 Facet_const_handle f ,
 bool aggressive
)
const
{
    if( index(f)<0 || index(f)>=facets.size() ){ fprintf( stderr , "facet index out of bounds\n" ) ; return false; }
    if( f._idx<0 || f._idx>=halfedges.size() ){ fprintf( stderr , "facet edge index out of bounds\n" ) ; return false; }
    const E* e = _halfedge( f );
    const E* s = _halfedge( f );
    int count = 0;
    do
    {
        if( e._facetIndex!=(f-&facets[0]) ){ fprintf( stderr , "facet to edge to facet failure\n" ) ; return false; }
        e = _next( e );
        count++;
    }
    while( e!=s );
    if( aggressive && count<3 ){ fprintf( stderr , "facet with less than three edges\n" ) ; return false; }
    return true;
}
template< class VertexData , class EdgeData , class FacetData >
bool
HalfEdgeMesh< VertexData , EdgeData , FacetData >::is_valid
(
 Halfedge_const_handle e ,
 bool aggressive
)
const
{
    if( index(e)<0 || index(e)>=halfedges.size() ){ fprintf( stderr , "halfedge index out of bounds\n" ) ; return false; }
    if( aggressive )
    {
        if( e->_endVertexIndex==halfedges[e->_oppositeEdgeIndex]._endVertexIndex ){ fprintf( stderr , "halfedge endpoints failure\n" ) ; return false; }
//      if( e->nextEdgePointer==e->oppositeEdgePointer ){ fprintf( stderr , "halfedge to next to opposite failure\n" ) ; return false; }
    }
    if( e._facetIndex!=-1 && (e._facetIndex<0 || e._facetIndex=facets.size() ) ){ fprintf( stderr , "halfedge opposite index out of bounds\n" ) ; return false; }
    if( e._endVertexIndex<0 || e._endVertexIndex>= vertices.size() ){ fprintf( stderr , "halfedge to end-vertex index out of bounds\n" ) ; return false; }
    if( e._oppositeEdgeIndex<0 || e._oppositeEdgeIndex>=halfedges.size() ){ fprintf( stderr , "halfedge to opposite index out of bounds\n" ) ; return false; }
    if( e._nextEdgeIndex<0 || e._nextEdgeIndex>=halfedges.size() ){ fprintf( stderr , "halfedge to next index out of bounds\n" ) ; return false; }

    if( &halfedges[halfedges[e._oppositeEdgeIndex]._oppositeEdgeIndex]!=e ){ fprintf( stderr , "halfedge to opposite to opposite failure\n" ) ; return false; }
    if( _previous( e )->_nextEdgeIndex!=(e-&halfedges[0]) ){ fprintf( stderr , "halfedge to previous to next failure\n" ) ; return false; }
    if( _previous( &halfedges[e._nextEdgeIndex] )!=e ){ fprintf( stderr , "halfedge to next to previous failure\n" ) ; return false; }
    return true;
}
template< class VertexData , class EdgeData , class FacetData >
int HalfEdgeMesh< VertexData , EdgeData , FacetData >::vertex_valence
(
 Vertex_const_handle v
)
const
{
    int count = 0;
    Halfedge_around_vertex_const_circulator s( v ) , e( v );
    do
    {
        count++;
        e++;
    }
    while( s!=e );
    return count;
}
template< class VertexData , class EdgeData , class FacetData >
int HalfEdgeMesh< VertexData , EdgeData , FacetData >::facet_size
(
 Facet_const_handle f
)
const
{
    int count = 0;
    Halfedge_around_facet_const_circulator s( f ) , e( f );
    do
    {
        count++;
        e++;
    }
    while( s!=e );
    return count;
}
