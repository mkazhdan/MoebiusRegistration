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

///////////////////////////////////////////////////////////////////////////////
//  HalfEdge.h
////////////////////////////////////////////////////////////////////////////////
//
//  This is the second version of our half-edge datastructure.
//  We completely reimplemented this. The main features are:
//
//  1. No hash table anymore. Hopefully this will gain us some speedup.
//
//  2. Every half-edge has an opposiate edge now. If on the boundary, that edge's
//	   triangle pointer will be NULL (negative).
//
//  3. Given any boundary-edge, we can now walk along the whole boundary easily 
//     by querying the next edge (nextEdgePointer) repeatedly. Eventually, we 
//     we will be back to the original edge.
//
//  4. For the sake of flexibility, each edge/face/vertex now has a property 
//     container (EdgeData/FaceData/VertexData).
//
//  5. Since each edge points to its next edge, walking around a vertex by going:
//	        next-edge, opposite-edge , next-edge, opposite-edge, next, etc.
//     goes in the opposite direction of the face vertices
//     Also, circulators around vertices on the boundary WILL travel around
//     every vertex in the one-ring.
//
//  
// DZ: The names of member functions compatible with CGAL are marked
// Best not to change these name, and if any changes re made to these functions,
// best to check semantic compatibility with CGAL definitions
// Note: this is not a complete CGAL halfedge interface;
// an attempt was made to minimize changes needed to convert CGAL code
////////////////////////////////////////////////////////////////////////////////
#ifndef HALF_EDGE_INCLUDED
#define HALF_EDGE_INCLUDED

#define TEST_HE_MESH_VALIDITY 0

#include "Geometry.h"
#include <vector>

////////////////////////////////////////////////////////////////////////////////
//  [ HalfEdgeMesh ]
/*! A container for storing the half-edge representation of a mesh
//  @tparam     VertexData      Data and methods on the vertices
//  @tparam     EdgeData        Data and methods on the edges
//  @tparam     FacetData       Data and methods on the faces
*///////////////////////////////////////////////////////////////////////////////
template< class VertexData , class EdgeData , class FacetData >
class HalfEdgeMesh
{
	template<class inVertexData, class inEdgeData, class inFacetData,
             class outVertexData, class outEdgeData, class outFacetData>
    friend void CopyTopology(const HalfEdgeMesh< inVertexData, inEdgeData,
            inFacetData >& inMesh, HalfEdgeMesh< outVertexData, outEdgeData,
            outFacetData >& outMesh);

	// Forward Declarations
protected:
	class V;
	class E;
	class F;
public:
	class Halfedge_around_vertex_circulator;
	class Halfedge_around_vertex_const_circulator;
	class Halfedge_handle;
	class Halfedge_const_handle;
	class Vertex_handle;
	class Vertex_const_handle;
	class Facet_handle;
	class Facet_const_handle;
protected:

	////////////////////////////////////////////////////////////////////////////////
	//  [ _const_handle ]
	/*! A container for storing a handle to mesh elements
	//  @tparam     ConstHandle     The type of handle this will be
	*///////////////////////////////////////////////////////////////////////////////
	template< class ConstHandle >
	class _const_handle
	{
		friend class HalfEdgeMesh;
		int _idx;
		const HalfEdgeMesh* _mesh;
		_const_handle( int idx , const HalfEdgeMesh* mesh ){ _idx=idx , _mesh=mesh; }
	public:
		_const_handle( void ) { _idx=-1 , _mesh=NULL; }
		operator bool( void ) const { return _mesh!=NULL && _idx!=-1; }
		bool operator == ( const ConstHandle& h ) const { return (h._idx==_idx && h._mesh==_mesh) || (h._mesh==NULL && _mesh==NULL); }
		bool operator != ( const ConstHandle& h ) const { return (h._idx!=_idx || h._mesh!=_mesh) && (h._mesh!=NULL || _mesh!=NULL); }
		_const_handle& operator++ ( void ) { if( _mesh && _idx!=-1 ) _idx++ ; return *this; }
		_const_handle& operator-- ( void ) { if( _mesh && _idx!=-1 ) _idx-- ; return *this; }
		_const_handle  operator++ ( int  ) { ConstHandle old = *this; ++(*this); return old; }
		_const_handle  operator-- ( int  ) { ConstHandle old = *this; --(*this); return old; }
        int index() const { return _idx; }
	};
	////////////////////////////////////////////////////////////////////////////////
	//  [ _handle ]
	/*! A container for storing a handle to mesh elements
	//  @tparam     Handle          The type of handle this will be
	//  @tparam     ConstHandle     The constant version of the handle
	*///////////////////////////////////////////////////////////////////////////////
	template< class Handle , class ConstHandle >
	class _handle
	{
		friend class HalfEdgeMesh;
		int _idx;
		HalfEdgeMesh* _mesh;
		_handle( int idx , HalfEdgeMesh* mesh ){ _idx=idx , _mesh=mesh; }
	public:
		_handle( void ) { _idx=-1 , _mesh=NULL; }
		operator bool( void ) const { return _mesh!=NULL && _idx!=-1; }
		operator ConstHandle( void ) const { return ConstHandle( _idx , _mesh ); }
		bool operator == ( const Handle& h ) const { return (h._idx==_idx && h._mesh==_mesh) || (h._mesh==NULL && _mesh==NULL); }
		bool operator != ( const Handle& h ) const { return (h._idx!=_idx || h._mesh!=_mesh) && (h._mesh!=NULL || _mesh!=NULL); }
		_handle& operator++ ( void ) { if( _mesh && _idx!=-1 ) _idx++ ; return *this; }
		_handle& operator-- ( void ) { if( _mesh && _idx!=-1 ) _idx-- ; return *this; }
		_handle  operator++ ( int  ) { Handle old = *this; ++(*this); return old; }
		_handle  operator-- ( int  ) { Handle old = *this; --(*this); return old; }
        int index() const { return _idx; }
	};
public:
	////////////////////////////////////////////////////////////////////////////////
	//  [ Halfedge_const_handle ]
	//  [ Halfedge_handle       ]
	/*! Handles for supporting halfedges
	*///////////////////////////////////////////////////////////////////////////////
	class Halfedge_handle : public _handle<Halfedge_handle, Halfedge_const_handle>
	{
#ifdef WIN32
		// Got errors:
		// "HalfEdgeMesh<VertexData,EdgeData,FacetData>::_handle<Handle,ConstHandle>::_mesh' is not accessible from  'HalfEdgeMesh<VertexData,EdgeData,FacetData>::_handle<Handle,ConstHandle>"
		// "HalfEdgeMesh<VertexData,EdgeData,FacetData>::_handle<Handle,ConstHandle>::_idx' is not accessible from  'HalfEdgeMesh<VertexData,EdgeData,FacetData>::_handle<Handle,ConstHandle>"
#else // !WIN32
        using _handle<Halfedge_handle, Halfedge_const_handle>::_mesh;
        using _handle<Halfedge_handle, Halfedge_const_handle>::_idx;
#endif // WIN32
		friend class HalfEdgeMesh;
        friend class Vertex_handle;

		Halfedge_handle(int idx, HalfEdgeMesh *mesh)
            : _handle<Halfedge_handle, Halfedge_const_handle>(idx, mesh)
        { }
	public:
        typedef E value_type;

        // JP: Superclass is automatically called for constructors without an
        // argument
		// Halfedge_handle( void ) : _handle( ){ }
		Halfedge_handle(void) { }
        // JP: changed to return E instead of EdgeData
	    E* operator->( void ) const { return &_mesh->halfedges[_idx]; }
		E& operator* ()       const { return  _mesh->halfedges[_idx]; }

		Halfedge_handle halfedge  ( void ) const { return *this; }
		Halfedge_handle next      ( void ) const { return Halfedge_handle( _mesh->halfedges[_idx]._nextEdgeIndex , _mesh ); }
		Halfedge_handle prev      ( void ) const { return next().next(); }
		Halfedge_handle opposite  ( void ) const { return Halfedge_handle( _mesh->halfedges[_idx]._oppositeEdgeIndex , _mesh ); }
		Vertex_handle vertex      ( void ) const { return Vertex_handle  ( _mesh->halfedges[_idx]._endVertexIndex , _mesh ); }
		Vertex_handle end_vertex  ( void ) const { return Vertex_handle  ( _mesh->halfedges[_idx]._endVertexIndex , _mesh ); }
		Vertex_handle start_vertex( void ) const { return Vertex_handle  ( _mesh->halfedges[ _mesh->halfedges[_idx]._oppositeEdgeIndex]._endVertexIndex , _mesh ); }
		Facet_handle facet        ( void ) const { return Facet_handle   ( _mesh->halfedges[_idx]._facetIndex , _mesh ); }
        using _handle<Halfedge_handle, Halfedge_const_handle>::index;

        ////////////////////////////////////////////////////////////////////////
        /*! Changes this halfedge's opposite halfedge
        //  @param[in]  opposite halfedge to which this halfedge should point
        *///////////////////////////////////////////////////////////////////////
        void set_opposite(Halfedge_handle opp)
        {
            assert(_mesh == opp._mesh);
            _mesh->halfedges[_idx]._oppositeEdgeIndex = opp._idx;
        }

		bool is_border( void ) const { return _mesh->halfedges[_idx]._facetIndex==-1; }
		bool is_border_edge( void ) const { return _mesh->halfedges[_idx]._facetIndex==-1 || _mesh->halfedges[_mesh->halfedges[_idx]._oppositeEdgeIndex]._facetIndex==-1; }
        // Arbitrary priority assignment
        bool hasPriority(void) const { return _idx < opposite()._idx; }

		std::pair< Vertex_handle , Vertex_handle > vertices( void ) const
		{
			return std::pair< Vertex_handle , Vertex_handle >(
				Vertex_handle( _mesh->halfedges[ _mesh->halfedges[_idx]._oppositeEdgeIndex ]._endVertexIndex , _mesh ) ,
				Vertex_handle( _mesh->halfedges[_idx]._endVertexIndex , _mesh ) );
		}
		inline long long key( void ) const { return ( ( (long long) _mesh->halfedges[ _mesh->halfedges[_idx]._oppositeEdgeIndex ]._endVertexIndex )<<32 ) | ( (long long) _mesh->halfedges[_idx]._endVertexIndex ); }
	};
	class Halfedge_const_handle : public _const_handle< Halfedge_const_handle >
	{
#ifdef WIN32
#else // !WIN32
        using _const_handle<Halfedge_const_handle>::_mesh;
        using _const_handle<Halfedge_const_handle>::_idx;
#endif // WIN32
		friend class HalfEdgeMesh;

		Halfedge_const_handle(int idx, const HalfEdgeMesh *mesh)
            : _const_handle<Halfedge_const_handle>(idx, mesh)
        { }
	public:
        typedef E value_type;

        // JP: Superclass is automatically called for constructors without an
        // argument
		// Halfedge_const_handle( void ) : _const_handle( ){ }
		Halfedge_const_handle(void) { }
        // JP: changed to return E instead of EdgeData
		const E* operator->( void ) const { return &_mesh->halfedges[_idx]; }
		const E& operator* ()       const { return  _mesh->halfedges[_idx]; }

		Halfedge_const_handle halfedge  ( void ) const { return *this; }
		Halfedge_const_handle next      ( void ) const { return Halfedge_const_handle( _mesh->halfedges[_idx]._nextEdgeIndex , _mesh ); }
		Halfedge_const_handle prev      ( void ) const { return next().next(); }
		Halfedge_const_handle opposite  ( void ) const { return Halfedge_const_handle( _mesh->halfedges[_idx]._oppositeEdgeIndex , _mesh ); }
		Vertex_const_handle vertex      ( void ) const { return Vertex_const_handle( _mesh->halfedges[_idx]._endVertexIndex , _mesh ); }
		Vertex_const_handle end_vertex  ( void ) const { return Vertex_const_handle( _mesh->halfedges[_idx]._endVertexIndex , _mesh ); }
		Vertex_const_handle start_vertex( void ) const { return Vertex_const_handle( _mesh->halfedges[ _mesh->halfedges[_idx]._oppositeEdgeIndex]._endVertexIndex , _mesh ); }
		Facet_const_handle facet        ( void ) const { return Facet_const_handle( _mesh->halfedges[_idx]._facetIndex , _mesh ); }
        using _const_handle<Halfedge_const_handle>::index;

		bool is_border( void ) const { return _mesh->halfedges[_idx]._facetIndex==-1; }
		bool is_border_edge( void ) const { return _mesh->halfedges[_idx]._facetIndex==-1 || _mesh->halfedges[_mesh->halfedges[_idx]._oppositeEdgeIndex]._facetIndex==-1; }
        // Arbitrary priority assignment
        bool hasPriority(void) const { return _idx < opposite()._idx; }

		std::pair< Vertex_const_handle , Vertex_const_handle > vertices( void ) const
		{ 
			return std::pair< Vertex_const_handle , Vertex_const_handle >(
				Vertex_const_handle( _mesh->halfedges[ _mesh->halfedges[_idx]._oppositeEdgeIndex ]._endVertexIndex , _mesh ) ,
				Vertex_const_handle( _mesh->halfedges[_idx]._endVertexIndex , _mesh ) );
		}
		inline long long key( void ) const { return ( ( (long long) _mesh->halfedges[ _mesh->halfedges[_idx]._oppositeEdgeIndex ]._endVertexIndex )<<32 ) | ( (long long) _mesh->halfedges[_idx]._endVertexIndex ); }
	};

	////////////////////////////////////////////////////////////////////////////////
	//  [ Vertex_const_handle ]
	//  [ Vertex_handle       ]
	/*! Handles for supporting vertices
	*///////////////////////////////////////////////////////////////////////////////
	class Vertex_handle : public _handle< Vertex_handle , Vertex_const_handle >
	{
#ifdef WIN32
#else // !WIN32
        using _handle<Vertex_handle, Vertex_const_handle>::_mesh;
        using _handle<Vertex_handle, Vertex_const_handle>::_idx;
#endif // WIN32
		friend class HalfEdgeMesh;

		Vertex_handle(int idx, HalfEdgeMesh *mesh)
            : _handle<Vertex_handle, Vertex_const_handle>(idx, mesh)
        { }
	public:
        typedef V value_type;

        // JP: Superclass is automatically called for constructors without an
        // argument
		// Vertex_handle( void ) : _handle( ){ }
		Vertex_handle(void) { }
        // JP: changed to return V instead of VertexData
		V* operator->( void ) const { return &_mesh->vertices[_idx]; }
		V& operator* ()       const { return  _mesh->vertices[_idx]; }

        /*! Gets the halfedge incident on this vertex */
		Halfedge_handle halfedge( void ){ return Halfedge_handle( _mesh->halfedges[ _mesh->vertices[_idx]._outEdgeIndex ]._oppositeEdgeIndex , _mesh ); }

		Halfedge_handle out_halfedge( void ) const { return Halfedge_handle( _mesh->vertices[_idx]._outEdgeIndex , _mesh ); }
		Halfedge_handle  in_halfedge( void ) const { return Halfedge_handle( _mesh->halfedges[ _mesh->vertices[_idx]._outEdgeIndex ]._oppositeEdgeIndex , _mesh ); }
        using _handle<Vertex_handle, Vertex_const_handle>::index;

        ////////////////////////////////////////////////////////////////////////
        /*! Updates the incident edge information
        //  NOTE: changes not only out_halfedge but also in_halfedge
        //  (changed to out_halfedge.opposite())
        //  @param[in]  outHalfedge new incident halfedge
        *///////////////////////////////////////////////////////////////////////
        void set_out_halfedge(Halfedge_handle outHalfedge)
        {
            assert(_mesh == outHalfedge._mesh);
            assert(*this == outHalfedge.opposite().vertex());
            _mesh->vertices[_idx]._outEdgeIndex = outHalfedge._idx;
        }

		/*! Gets a circulator of halfedges around this vertex (clockwise). */
		Halfedge_around_vertex_circulator vertex_begin( void ) const { return Halfedge_around_vertex_circulator( *this ); }
	};
	class Vertex_const_handle : public _const_handle< Vertex_const_handle >
	{
#ifdef WIN32
#else // !WIN32
        using _const_handle<Vertex_const_handle>::_mesh;
        using _const_handle<Vertex_const_handle>::_idx;
#endif // WIN32
		friend class HalfEdgeMesh;

		Vertex_const_handle(int idx, const HalfEdgeMesh *mesh)
            : _const_handle<Vertex_const_handle>(idx, mesh)
        { }
	public:
        typedef V value_type;

        // JP: Superclass is automatically called for constructors without an
        // argument
		// Vertex_const_handle( void ) : _const_handle( ){ }
		Vertex_const_handle(void) { }
        // JP: changed to return V instead of VertexData
		const V* operator->( void ) const { return &_mesh->vertices[_idx]; }
		const V& operator* ()       const { return  _mesh->vertices[_idx]; }

        /*! Gets the halfedge incident on this vertex */
		Halfedge_const_handle halfedge( void ) const { return Halfedge_const_handle( _mesh->halfedges[ _mesh->vertices[_idx]._outEdgeIndex ]._oppositeEdgeIndex , _mesh ); }

		// What's going on here? Why is the default half-edge returned by a vertex the incoming one?
		Halfedge_const_handle out_halfedge( void ) const { return Halfedge_const_handle( _mesh->vertices[_idx]._outEdgeIndex , _mesh ); }
		Halfedge_const_handle  in_halfedge( void ) const { return Halfedge_const_handle( _mesh->halfedges[ _mesh->vertices[_idx]._outEdgeIndex ]._oppositeEdgeIndex , _mesh ); }
        using _const_handle<Vertex_const_handle>::index;

		Halfedge_around_vertex_const_circulator vertex_begin( void ) const { return Halfedge_around_vertex_const_circulator( *this ); }

		operator bool( void ) const {
            return _mesh!=NULL && ((unsigned int) _idx) < _mesh->vertices.size();
        }
	};

	////////////////////////////////////////////////////////////////////////////////
	//  [ Facet_const_handle ]
	//  [ Facet_handle       ]
	/*! Handles for supporting facets
	*///////////////////////////////////////////////////////////////////////////////
	class Facet_handle : public _handle< Facet_handle , Facet_const_handle >
	{
#ifdef WIN32
#else // !WIN32
        using _handle<Facet_handle, Facet_const_handle>::_mesh;
        using _handle<Facet_handle, Facet_const_handle>::_idx;
#endif // WIN32
		friend class HalfEdgeMesh;

		Facet_handle(int idx, HalfEdgeMesh *mesh)
            : _handle<Facet_handle, Facet_const_handle>(idx, mesh)
        { }
	public:
        typedef F value_type;

        // JP: Superclass is automatically called for constructors without an
        // argument
		// Facet_handle( void ) : _handle( ){ }
		Facet_handle(void) { }
        // JP: Changed to return F instead of FacetData
		F* operator->( void ) const { return &_mesh->facets[_idx]; }
		F& operator* ()       const { return  _mesh->facets[_idx]; }
		Halfedge_handle halfedge( void ) const { return Halfedge_handle( _mesh->facets[_idx]._edgeIndex , _mesh ); }
        using _handle<Facet_handle, Facet_const_handle>::index;
	};
	class Facet_const_handle : public _const_handle< Facet_const_handle >
	{
#ifdef WIN32
#else // !WIN32
        using _const_handle<Facet_const_handle>::_mesh;
        using _const_handle<Facet_const_handle>::_idx;
#endif // WIN32
		friend class HalfEdgeMesh;

		Facet_const_handle(int idx, const HalfEdgeMesh *mesh)
            : _const_handle<Facet_const_handle>(idx, mesh)
        { }
	public:
        typedef F value_type;

        // JP: Superclass is automatically called for constructors without an
        // argument
		// Facet_const_handle( void ) : _const_handle( ){ }
		Facet_const_handle(void) { }
        // JP: Changed to return F instead of FacetData
		const F* operator->( void ) const { return &_mesh->facets[_idx]; }
		const F& operator* ()       const { return  _mesh->facets[_idx]; }
		Halfedge_const_handle halfedge( void ) const { return Halfedge_const_handle( _mesh->facets[_idx]._edgeIndex , _mesh ); }
        using _const_handle<Facet_const_handle>::index;
	};

protected:
	////////////////////////////////////////////////////////////////////////////////
	//  [ const_around_facet_circulator ]
	//  [       around_facet_circulator ]
	//  Circulators for traversing the around a facet.
	//  In principal, the circulator only supports forward access (++).
	//  However, _only_ for triangle meshes we support backward access (--).
    //  Note: halfedge is "inward" (its tip is the vertex around which we
    //        circulate)
	////////////////////////////////////////////////////////////////////////////////
	template< class HandleType , HandleType (Halfedge_handle::*MemberFunction)( void ) const >
	struct around_facet_circulator : public HandleType
	{
	protected:
		Halfedge_handle e;
	public:
		around_facet_circulator( void ){ e = Halfedge_handle(); }
		around_facet_circulator( Facet_handle f ){ e = f.halfedge() ; *( (HandleType*)this ) = ((&e)->*MemberFunction)( ); }
		around_facet_circulator( Halfedge_handle e ){ this->e = e  ; *( (HandleType*)this ) = ((&e)->*MemberFunction)( ) ; }
		around_facet_circulator& operator++ ( void ){ e = e.next() ; *( (HandleType*)this ) = ((&e)->*MemberFunction)( ) ; return *this; }
		around_facet_circulator& operator-- ( void ){ e = e.prev() ; *( (HandleType*)this ) = ((&e)->*MemberFunction)( ) ; return *this; }
		around_facet_circulator  operator++ ( int  ) { around_facet_circulator old = *this; ++(*this); return old; }
		around_facet_circulator  operator-- ( int  ) { around_facet_circulator old = *this; --(*this); return old; }
	};

	template< class HandleType , HandleType (Halfedge_const_handle::*MemberFunction)( void ) const >
	struct const_around_facet_circulator : public HandleType
	{
	protected:
		Halfedge_const_handle e;
	public:
		const_around_facet_circulator( void ){ e = Halfedge_const_handle(); }
		const_around_facet_circulator( Facet_const_handle f ){ e = f.halfedge() ; *( (HandleType*)this ) = ((&e)->*MemberFunction)( ); }
		const_around_facet_circulator( Halfedge_const_handle e ){ this->e = e ; *( (HandleType*)this ) = ((&e)->*MemberFunction)( ); }
		const_around_facet_circulator& operator++ ( void ) { e = e.next() ; *( (HandleType*)this ) = ((&e)->*MemberFunction)() ; return *this; }
		const_around_facet_circulator& operator-- ( void ) { e = e.prev() ; *( (HandleType*)this ) = ((&e)->*MemberFunction)() ; return *this; }
		const_around_facet_circulator  operator++ ( int  ) { const_around_facet_circulator old = *this ; ++(*this) ; return old; }
		const_around_facet_circulator  operator-- ( int  ) { const_around_facet_circulator old = *this ; --(*this) ; return old; }
	};

	////////////////////////////////////////////////////////////////////////////////
	//  [ const_around_vertex_circulator ]
	//  [       around_vertex_circulator ]
	//  Circulators for traversing the around a vertex in CW order.
	//  In principal, the circulator only supports forward access (++).
	//  However, _only_ for triangle meshes we support backward access (--).
    //  Note: halfedge is "inward" (its tip is the vertex around which we
    //        circulate)
	////////////////////////////////////////////////////////////////////////////////
	// Basic circulators to go around a vertex
	template< class HandleType , HandleType (Halfedge_handle::*MemberFunction)( void ) const >
	struct around_vertex_circulator : public HandleType
	{
	protected:
		Halfedge_handle e;
	public:
		around_vertex_circulator( void ){ e = Halfedge_handle(); }
		around_vertex_circulator( Vertex_handle v ){ e = v.halfedge() ; *( (HandleType*)this ) = ((&e)->*MemberFunction)( ) ; }
		around_vertex_circulator( Halfedge_handle e ){ this->e = e  ; *( (HandleType*)this ) = ((&e)->*MemberFunction)( ) ; }
		around_vertex_circulator& operator++ ( void ) { e = e.next().opposite() ; *( (HandleType*)this ) = ((&e)->*MemberFunction)() ; return *this; }
		around_vertex_circulator& operator-- ( void ) { e = e.opposite().prev() ; *( (HandleType*)this ) = ((&e)->*MemberFunction)() ; return *this; }
		around_vertex_circulator  operator++ ( int ) { around_vertex_circulator old = *this; ++(*this); return old; }
		around_vertex_circulator  operator-- ( int ) { around_vertex_circulator old = *this; --(*this); return old; }
	};

	template< class HandleType , HandleType (Halfedge_const_handle::*MemberFunction)( void ) const >
	struct const_around_vertex_circulator : public HandleType
	{
	protected:
		Halfedge_const_handle e;
	public:
		const_around_vertex_circulator( void ){ e = Halfedge_const_handle(); }
		const_around_vertex_circulator( Vertex_const_handle v ){ e = v.halfedge() ; *( (HandleType*)this ) = ((&e)->*MemberFunction)( ) ; }
		const_around_vertex_circulator( Halfedge_const_handle e ){ this->e = e  ; *( (HandleType*)this ) = ((&e)->*MemberFunction)( ) ; }
	    const_around_vertex_circulator& operator++ ( void ) { e = e.next().opposite() ; *( (HandleType*)this ) = ((&e)->*MemberFunction)() ; return *this; }
	    const_around_vertex_circulator& operator-- ( void ) { e = e.opposite().prev() ; *( (HandleType*)this ) = ((&e)->*MemberFunction)() ; return *this; }
		const_around_vertex_circulator  operator++ ( int ) { const_around_vertex_circulator old = *this; ++(*this); return old; }
		const_around_vertex_circulator  operator-- ( int ) { const_around_vertex_circulator old = *this; --(*this); return old; }
	};

public:
	////////////////////////////////////////////////////////////////////////////////
	//  [ Vertex_const_around_facet_circulator ]
	//  [       Vertex_around_facet_circulator ]
	//  Circulators for traversing the vertices in a facet in CCW order.
	//  In principal, the circulator only supports forward access (++).
	//  However, _only_ for triangle meshes we support backward access (--).
	////////////////////////////////////////////////////////////////////////////////
	struct Vertex_around_facet_circulator : public around_facet_circulator< Vertex_handle , &Halfedge_handle::vertex >
	{
		friend class HalfEdgeMesh;
	public:
		Vertex_around_facet_circulator( void )
            : around_facet_circulator<Vertex_handle, &Halfedge_handle::vertex>( ) { }
		Vertex_around_facet_circulator(Facet_handle f)
            : around_facet_circulator<Vertex_handle, &Halfedge_handle::vertex>(f) { }
		Vertex_around_facet_circulator(Halfedge_handle e)
            : around_facet_circulator<Vertex_handle, &Halfedge_handle::vertex>(e) { }
	};
	struct Vertex_around_facet_const_circulator : public const_around_facet_circulator< Vertex_const_handle , &Halfedge_const_handle::vertex >
	{
		friend class HalfEdgeMesh;
	public:
		Vertex_around_facet_const_circulator(void)
            : const_around_facet_circulator<Vertex_const_handle, &Halfedge_const_handle::vertex>( )
        { }
		Vertex_around_facet_const_circulator(Facet_const_handle f)
            : const_around_facet_circulator<Vertex_const_handle, &Halfedge_const_handle::vertex>(f)
        { }
		Vertex_around_facet_const_circulator( Halfedge_const_handle e )
            : const_around_facet_circulator<Vertex_const_handle, &Halfedge_const_handle::vertex>( e )
        { }
	};

	////////////////////////////////////////////////////////////////////////////////
	//  [ Halfedge_around_const_facet_circulator ]
	//  [       Halfedge_around_facet_circulator ]
	//  Circulators for traversing the half-edges in a facet in CCW order.
	//  In principal, the circulator only supports forward access (++).
	//  However, _only_ for triangle meshes we support backward access (--).
	////////////////////////////////////////////////////////////////////////////////
	struct Halfedge_around_facet_circulator : public around_facet_circulator< Halfedge_handle , &Halfedge_handle::halfedge >
	{
		friend class HalfEdgeMesh;
	public:
		Halfedge_around_facet_circulator( void )
            : around_facet_circulator<Halfedge_handle, &Halfedge_handle::halfedge>( )
        { }
		Halfedge_around_facet_circulator( Facet_handle f )
            : around_facet_circulator<Halfedge_handle, &Halfedge_handle::halfedge>(f)
        { }
		Halfedge_around_facet_circulator( Halfedge_handle e )
            : around_facet_circulator<Halfedge_handle, &Halfedge_handle::halfedge>(e)
        { }
	};
	struct Halfedge_around_facet_const_circulator : public const_around_facet_circulator< Halfedge_const_handle , &Halfedge_const_handle::halfedge >
	{
		friend class HalfEdgeMesh;
	public:
		Halfedge_around_facet_const_circulator( void )
            : const_around_facet_circulator<Halfedge_const_handle, &Halfedge_const_handle::halfedge>( )
        { }
		Halfedge_around_facet_const_circulator( Facet_const_handle f )
            : const_around_facet_circulator<Halfedge_const_handle, &Halfedge_const_handle::halfedge>(f)
        { }
		Halfedge_around_facet_const_circulator( Halfedge_const_handle e )
            : const_around_facet_circulator<Halfedge_const_handle, &Halfedge_const_handle::halfedge>(e)
        { }
	};

	////////////////////////////////////////////////////////////////////////////////
	//  [ const_Halfedge_around_vertex_circulator ]
	//  [       Halfedge_around_vertex_circulator ]
	//  Circulators for traversing the half-edges coming around a vertex in CW order.
	////////////////////////////////////////////////////////////////////////////////
	struct Halfedge_around_vertex_const_circulator : public const_around_vertex_circulator< Halfedge_const_handle , &Halfedge_const_handle::halfedge >
	{
	protected:
		using const_around_vertex_circulator<Halfedge_const_handle, &Halfedge_const_handle::halfedge>::e;
		friend class HalfEdgeMesh;
	public:
		Halfedge_around_vertex_const_circulator( void )
            : const_around_vertex_circulator<Halfedge_const_handle, &Halfedge_const_handle::halfedge>( )
        { }
		Halfedge_around_vertex_const_circulator( Vertex_const_handle v )
            : const_around_vertex_circulator<Halfedge_const_handle, &Halfedge_const_handle::halfedge>(v)
        { }
		Halfedge_around_vertex_const_circulator( Halfedge_const_handle e )
            : const_around_vertex_circulator<Halfedge_const_handle, &Halfedge_const_handle::halfedge>(e)
        { }
	};
	struct Halfedge_around_vertex_circulator : public around_vertex_circulator< Halfedge_handle , &Halfedge_handle::halfedge >
	{
	  //	protected:
	public: // DZ making this public as in CGAL
		using around_vertex_circulator<Halfedge_handle, &Halfedge_handle::halfedge>::e;
		friend class HalfEdgeMesh;
	public:
		Halfedge_around_vertex_circulator( void )
            : around_vertex_circulator<Halfedge_handle , &Halfedge_handle::halfedge>( )
        { }
		Halfedge_around_vertex_circulator( Vertex_handle v )
            : around_vertex_circulator<Halfedge_handle , &Halfedge_handle::halfedge>(v)
        { }
		Halfedge_around_vertex_circulator( Halfedge_handle e )
            : around_vertex_circulator<Halfedge_handle , &Halfedge_handle::halfedge>(e)
        { }

        /** Allow typecast to Halfedge_around_vertex_const_circulator */
        operator Halfedge_around_vertex_const_circulator() const
        {
            return Halfedge_around_vertex_const_circulator(e);
        }
	};

	////////////////////////////////////////////////////////////////////////////////
	//  [       Vertex_around_vertex_const_circulator ]
	//  [       Vertex_around_vertex_circulator ]
	//  Circulators for traversing the vertices in the one-ring of a vertex in CW order.
	////////////////////////////////////////////////////////////////////////////////
	struct Vertex_around_vertex_const_circulator : public const_around_vertex_circulator< Vertex_const_handle , &Halfedge_const_handle::start_vertex > //CGAL name
	{
	protected:
		using const_around_vertex_circulator<Vertex_const_handle, &Halfedge_const_handle::start_vertex>::e;
		friend class HalfEdgeMesh;
	public:
		Vertex_around_vertex_const_circulator( void )
            : const_around_vertex_circulator<Vertex_const_handle, &Halfedge_const_handle::start_vertex>( )
        { }
		Vertex_around_vertex_const_circulator( Vertex_const_handle v )
            : const_around_vertex_circulator<Vertex_const_handle, &Halfedge_const_handle::start_vertex>(v)
        { }
		Vertex_around_vertex_const_circulator( Halfedge_const_handle e )
            : const_around_vertex_circulator<Vertex_const_handle, &Halfedge_const_handle::start_vertex>(e)
        { }
	};
	struct Vertex_around_vertex_circulator : public around_vertex_circulator< Vertex_handle , &Halfedge_handle::start_vertex >  //CGAL name
	{
	protected:
    using around_vertex_circulator<Vertex_handle, &Halfedge_handle::start_vertex>::e;
		friend class HalfEdgeMesh;
	public:
		Vertex_around_vertex_circulator( void )
            : around_vertex_circulator<Vertex_handle, &Halfedge_handle::start_vertex>( )
        { }
		Vertex_around_vertex_circulator( Vertex_handle v )
            : around_vertex_circulator<Vertex_handle, &Halfedge_handle::start_vertex>(v)
        { }
		Vertex_around_vertex_circulator( Halfedge_handle e )
            : around_vertex_circulator<Vertex_handle, &Halfedge_handle::start_vertex>(e)
        { }
        /** Allow typecast to Vertex_around_vertex_const_circulator */
        operator Vertex_around_vertex_const_circulator() const
        {
            return Vertex_around_vertex_const_circulator(e);
        }
	};

	////////////////////////////////////////////////////////////////////////////////
	//  [ E ]
	//  Edge in the half-edge mesh
	//  Note that the implementation of previous is only correct for triangle meshes
	////////////////////////////////////////////////////////////////////////////////
protected:
	class E: public EdgeData
	{
    public:
        typedef typename HalfEdgeMesh::Facet_const_handle    Facet_const_handle;
        typedef typename HalfEdgeMesh::Facet_handle          Facet_handle;
        typedef typename HalfEdgeMesh::Halfedge_const_handle Halfedge_const_handle;
        typedef typename HalfEdgeMesh::Halfedge_handle       Halfedge_handle;
        typedef typename HalfEdgeMesh::Vertex_const_handle   Vertex_const_handle;
        typedef typename HalfEdgeMesh::Vertex_handle         Vertex_handle;
    private:
		template<class inVertexData,  class inEdgeData,  class inFacetData, class outVertexData, class outEdgeData, class outFacetData>
		friend void CopyTopology( const HalfEdgeMesh< inVertexData , inEdgeData , inFacetData >& inMesh , HalfEdgeMesh< outVertexData , outEdgeData , outFacetData >& outMesh );
		friend class HalfEdgeMesh;
	protected:
		int _nextEdgeIndex;			// Pointer to the next edge on the face.
		int _oppositeEdgeIndex;		// Pointer to the opposite edge.
		int _facetIndex;			// Pointer to the facet tied to this half-edge or NULL if this is a boundary edge.
		int _endVertexIndex;		// Pointer to the vertex at the end of the edge.
	};

	////////////////////////////////////////////////////////////////////////////////
	//  [ F ]
	//  Face in the half-edge mesh
	////////////////////////////////////////////////////////////////////////////////
	// Assume that the face has a polygon boundary that is a single connected component.
	class F: public FacetData
	{
    public:
        typedef typename HalfEdgeMesh::Halfedge_const_handle Halfedge_const_handle;
    private:
		template<class inVertexData,  class inEdgeData,  class inFacetData, class outVertexData, class outEdgeData, class outFacetData>
		friend void CopyTopology( const HalfEdgeMesh< inVertexData , inEdgeData , inFacetData >& inMesh , HalfEdgeMesh< outVertexData , outEdgeData , outFacetData >& outMesh );
		friend class HalfEdgeMesh;
	protected:
		int _edgeIndex;		// Index of an edge on the face

	};

	////////////////////////////////////////////////////////////////////////////////
	//  [ V ]
	//  Vertex in the half-edge mesh
	////////////////////////////////////////////////////////////////////////////////
	class V: public VertexData
	{
    public:
        typedef typename HalfEdgeMesh::Halfedge_const_handle Halfedge_const_handle;
    private:
		template<class inVertexData,  class inEdgeData,  class inFacetData, class outVertexData, class outEdgeData, class outFacetData>
		friend void CopyTopology( const HalfEdgeMesh< inVertexData , inEdgeData , inFacetData >& inMesh , HalfEdgeMesh< outVertexData , outEdgeData , outFacetData >& outMesh );
		friend class HalfEdgeMesh;
	protected:
		int _outEdgeIndex;		// Index of a half-edge exiting the vertex.
	};
public:
	// This method copies over the topology from one half-edge mesh to an other.
	template< class _VertexData , class _EdgeData , class _FacetData >
	HalfEdgeMesh& operator = ( const HalfEdgeMesh< _VertexData , _EdgeData , _FacetData >& mesh );

	// These two methods take in a mesh and construct the half-edge topology for it.
	// These methods assume:
	// 1] That the mesh is topologically manifold (with boundary).
	// 2] That the triangles of the mesh are oriented in CCW order.
	// These methods guarantee:
	// 1] The vertices will appear in the same order that they are given.
	// 2] The facets will appear in the same order that they are given.
	// 3] The edge pointed to by a face will be the last edge in the original face.
	// Note that since the vertex list is not given, the number of vertices allocated
	// is the largest vertex index (plus one).
	
	bool SetHalfEdgeData( const std::vector<      TriangleIndex >& triangles , void (*start)(void)=NULL , void (*end)( void )=NULL ) { return _SetHalfEdgeData( triangles , vertices , halfedges , facets , TriangleSize , start , end ); }
	bool SetHalfEdgeData( const std::vector< std::vector< int > >& polygons  , void (*start)(void)=NULL , void (*end)( void )=NULL ) { return _SetHalfEdgeData( polygons  , vertices , halfedges , facets , PolygonSize  , start , end ); }
	template< class P >
	bool SetHalfEdgeData( const std::vector< P >& polygons , unsigned int (*PSize)( const P& ) , void (*start)(void)=NULL , void (*end)( void )=NULL ) { return _SetHalfEdgeData( polygons , vertices , halfedges , facets , PSize , start , end ); }

	////////////////////////////////////////////////////////////////////////////
	// Methods for accessing the vertex, half-edge, and facet data by index
	////////////////////////////////////////////////////////////////////////////
	size_t   vertex_size( void ) const { return  vertices.size(); }
	size_t halfedge_size( void ) const { return halfedges.size(); }
	size_t    facet_size( void ) const { return    facets.size(); }
	Vertex_const_handle     vertex( size_t v ) const { return   Vertex_const_handle( v , this ); }
	Vertex_handle           vertex( size_t v )       { return         Vertex_handle( v , this ); }
	Halfedge_const_handle halfedge( size_t e ) const { return Halfedge_const_handle( e , this ); }
	Halfedge_handle       halfedge( size_t e )       { return       Halfedge_handle( e , this ); }
	Facet_const_handle       facet( size_t f ) const { return    Facet_const_handle( (int)f , this ); }
	Facet_handle             facet( size_t f )       { return          Facet_handle( f , this ); }
	// Methods for getting the index of a vertex, half-edge, and facet handle
	size_t index(   Vertex_const_handle v ) const { return v._idx; }
	size_t index( Halfedge_const_handle e ) const { return e._idx; }
	size_t index(    Facet_const_handle f ) const { return f._idx; }

	// Methods for asserting whether a vertex handle is in the mesh
	Vertex_const_handle vertex( Vertex_const_handle v ) const { assert(index(v) < vertices.size()); return v; }
	Vertex_handle       vertex( Vertex_handle       v )       { assert(index(v) < vertices.size()); return v; }

	////////////////////////////////////////////////////////////////////////////////
	//  [       Vertex_const_iterator ]
	//  [       Vertex_iterator ]
	//  Iterators for traversing the vertices in the mesh
	////////////////////////////////////////////////////////////////////////////////
	typedef Vertex_const_handle       Vertex_const_iterator;
	typedef Vertex_handle             Vertex_iterator;

	////////////////////////////////////////////////////////////////////////////////
	//  [       Halfedge_const_iterator ]
	//  [       Halfedge_iterator ]
	//  Iterators for traversing the half-edges in the mesh
	////////////////////////////////////////////////////////////////////////////////
	typedef Halfedge_const_handle       Halfedge_const_iterator;
	typedef Halfedge_handle             Halfedge_iterator;

	////////////////////////////////////////////////////////////////////////////////
	//  [       Facet_const_iterator ]
	//  [       Facet_iterator ]
	//  Iterators for traversing the facets in the mesh
	////////////////////////////////////////////////////////////////////////////////
	typedef Facet_const_handle       Facet_const_iterator;
	typedef Facet_handle             Facet_iterator;

	// Methods for getting the begining and ending of the iterators
        // CGAL names
	Vertex_const_iterator    vertices_begin ( void ) const { return  vertices.begin(); }  
		Vertex_iterator      vertices_begin ( void )       { return  vertices.begin(); }  
	Halfedge_const_iterator halfedges_begin ( void ) const { return halfedges.begin(); }   
		Halfedge_iterator   halfedges_begin ( void )       { return halfedges.begin(); }  
	Facet_const_iterator       facets_begin ( void ) const { return    facets.begin(); } 
		Facet_iterator         facets_begin ( void )       { return    facets.begin(); } 
	Vertex_const_iterator      vertices_end ( void ) const { return  vertices.end  (); }
		Vertex_iterator        vertices_end ( void )       { return  vertices.end  (); }
	Halfedge_const_iterator   halfedges_end ( void ) const { return halfedges.end  (); } 
		Halfedge_iterator     halfedges_end ( void )       { return halfedges.end  (); } 
	Facet_const_iterator         facets_end ( void ) const { return    facets.end  (); } 
		Facet_iterator           facets_end ( void )       { return    facets.end  (); } 

	// Methods for getting the begining the circulators
	// these are NOT in CGAL; CGAL convention is to retrieve begin
	// circulators from the corresponding element
	Vertex_around_facet_const_circulator vertex_around_facet_begin( int fIndex ) const { return Vertex_around_facet_const_circulator( Facet_const_handle( fIndex , this ) ); }
	Vertex_around_facet_circulator       vertex_around_facet_begin( int fIndex )       { return Vertex_around_facet_circulator      (       Facet_handle( fIndex , this ) ); }
	Vertex_around_facet_const_circulator vertex_around_facet_begin( Facet_const_handle f ) const { return Vertex_around_facet_const_circulator( f ); }
	Vertex_around_facet_circulator       vertex_around_facet_begin( Facet_handle       f )       { return Vertex_around_facet_circulator      ( f ); }

	Halfedge_around_facet_const_circulator halfedge_around_facet_begin( int fIndex ) const { return Halfedge_around_facet_const_circulator( Facet_const_handle( fIndex , this ) ); }
	Halfedge_around_facet_circulator       halfedge_around_facet_begin( int fIndex )       { return Halfedge_around_facet_circulator      (       Facet_handle( fIndex , this ) ); }
	Halfedge_around_facet_const_circulator halfedge_around_facet_begin( Facet_const_handle f ) const { return Halfedge_around_facet_const_circulator( f ); }
	Halfedge_around_facet_circulator       halfedge_around_facet_begin( Facet_handle       f )       { return Halfedge_around_facet_circulator      ( f ); }

	Halfedge_around_vertex_const_circulator halfedge_around_vertex_begin( int vIndex ) const { return Halfedge_around_vertex_const_circulator( Vertex_const_handle( vIndex , this ) ); }
	Halfedge_around_vertex_circulator       halfedge_around_vertex_begin( int vIndex )       { return Halfedge_around_vertex_circulator      (       Vertex_handle( vIndex , this ) ); }
	Halfedge_around_vertex_const_circulator halfedge_around_vertex_begin( Vertex_const_handle v ) const { return Halfedge_around_vertex_const_circulator( v ); }
	Halfedge_around_vertex_circulator       halfedge_around_vertex_begin( Vertex_handle       v )       { return Halfedge_around_vertex_circulator      ( v ); }

	Vertex_around_vertex_const_circulator vertex_around_vertex_begin( int vIndex ) const { return Vertex_around_vertex_const_circulator( Vertex_const_handle( vIndex , this ) ); }
	Vertex_around_vertex_circulator       vertex_around_vertex_begin( int vIndex )       { return Vertex_around_vertex_circulator      (       Vertex_handle( vIndex , this ) ); }
	Vertex_around_vertex_const_circulator vertex_around_vertex_begin( Vertex_const_handle v ) const { return Vertex_around_vertex_const_circulator( v->halfedge() ); }
	Vertex_around_vertex_circulator       vertex_around_vertex_begin( Vertex_handle       v )       { return Vertex_around_vertex_circulator      ( v->halfedge() ); }

	static inline Halfedge_handle       Opposite     ( Halfedge_handle       e ) { return e.opposite(); }
	static inline Halfedge_const_handle Opposite     ( Halfedge_const_handle e ) { return e.opposite(); }
	static inline    Facet_handle       Facet        ( Halfedge_handle       e ) { return e.facet(); }
	static inline    Facet_const_handle Facet        ( Halfedge_const_handle e ) { return e.facet(); }
	static inline    Facet_handle       OppositeFacet( Halfedge_handle       e ) { return e.opposite().facet(); }
	static inline    Facet_const_handle OppositeFacet( Halfedge_const_handle e ) { return e.opposite().facet(); }

	// This method will try to remove a facet.
	bool remove_facet( Facet_handle f );
	bool remove_facet( int fIndex ) { return remove_facet( &facets[fIndex] ); }

	// This method will try to merge the two facets on either side of the edge into a single facet and will return a pointer to the new face.
	// The method will fail (return NULL) if:
	// 1] The two facets on either side of the edge are the same and neither of the two vertices has valence one,
	// 2] One of the facets is NULL.
	Facet_handle merge_facets( Halfedge_handle e , void (*merge)( const FacetData* , FacetData* )=NULL );
	Facet_handle merge_facets(        int eIndex , void (*merge)( const FacetData* , FacetData* )=NULL ) { return merge_facets( halfedge( eIndex ) , merge ); }

	bool merge_facets_full( bool (*toMerge)( Halfedge_const_handle ) , void (*mergeFacets)( const FacetData* , const FacetData* , FacetData* )=NULL );

	// This method will try to merge an edge with its next edge.
	// The method will fail (return NULL) if:
	// 1] The vertex between the two does not have valence 2
	// 2] The next and opposite halfedges are the same (i.e. the vertex at the end of the halfedge has valence one)
	// 3] The facets tied to the halfedge and its opposite have less than four halfedges (i.e. removing would create a facet with two or fewer halfedges)
	// 4] The merged edge would have the same endpoints
	Halfedge_handle merge_halfedges( Halfedge_handle e , void (*merge)( const EdgeData* , const EdgeData* , EdgeData* )=NULL );
	Halfedge_handle merge_halfedges(        int eIndex , void (*merge)( const EdgeData* , const EdgeData* , EdgeData* )=NULL ) { return merge_halfedges( halfedge( eIndex ) , merge ); }

	void sort_facets  ( int (*facet_compare_function )( const void* , const void* ) );
	void sort_vertices( int (*vertex_compare_function)( const void* , const void* ) );

	void print_facet_vertices ( Facet_const_handle f , bool new_line=true ) const;
	void print_facet_halfedges( Facet_const_handle f , bool new_line=true ) const;
	void print_halfedge( Halfedge_const_handle e , bool new_line=true ) const;
	bool is_valid( Vertex_const_handle   v , bool aggressive=false ) const;
	bool is_valid( Halfedge_const_handle e , bool aggressive=false ) const;
	bool is_valid( Facet_const_handle    f , bool aggressive=false ) const;
	// Note that failure with testEdges=true does not indicate a bad mesh. Only that there are more than two halfedges with the same start and end vertices.
	bool is_valid( bool aggressive=false , bool testEdges=false ) const;

	int vertex_valence( Vertex_const_handle v ) const;
	int facet_size( Facet_const_handle v ) const;

	static void (*swap_vertices_function )( VertexData* , VertexData* );
	static void (*swap_halfedges_function)( EdgeData*   , EdgeData*   );
	static void (*swap_facets_function   )( FacetData*  , FacetData*  );
protected:

	std::vector< V > vertices;	// The vertices in the mesh
	std::vector< E > halfedges;	// The half-edges in the mesh
	std::vector< F > facets;	// The facets in the mesh
	const E*            _previous( const E* ) const;
	      E*            _previous(       E* );
	inline const E*         _next( const E* ) const;
	inline       E*         _next(       E* );
	inline const E*     _opposite( const E* ) const;
	inline       E*     _opposite(       E* );
	inline const F*        _facet( const E* ) const;
	inline       F*        _facet(       E* );
	inline const E*     _halfedge( const F* ) const;
	inline       E*     _halfedge(       F* );
	inline const E* _out_halfedge( const V* ) const;
	inline       E* _out_halfedge(       V* );
	inline const V*   _end_vertex( const E* ) const;
	inline       V*   _end_vertex(       E* );

	void _swap_vertices ( V* v1 , V* v2 );
	void _swap_facets   ( F* f1 , F* f2 );
	void _swap_halfedges( E* e1 , E* e2 );
	void _remove_vertex  ( int v );
	void _remove_facet   ( int f );
	void _remove_halfedge( int e );

	static unsigned int TriangleSize( const TriangleIndex& tri ) { return 3; }
	static unsigned int PolygonSize ( const std::vector< int >& poly ) { return poly.size(); }
	// Transforms a mesh into a half-edge data structure
	template< class Polygon >
	static bool _SetHalfEdgeData
		(
		const std::vector< Polygon >& polygons ,
		std::vector< V >& vertices , 
		std::vector< E >& halfedges , 
		std::vector< F >& facets ,
		unsigned int (*)( const Polygon& ) , 
		void (*startFunction)( void ) , void (*endFunction)( void )
		);
};

template< class  VertexData , class  EdgeData , class  FacetData , class _VertexData , class _EdgeData , class _FacetData >
void CopyTopology(const HalfEdgeMesh< VertexData,  EdgeData,  FacetData > &inMesh,
                        HalfEdgeMesh<_VertexData, _EdgeData, _FacetData > &outMesh);

template< class VertexData , class EdgeData , class FacetData >
void (*HalfEdgeMesh< VertexData , EdgeData , FacetData >::swap_facets_function   )( FacetData*  , FacetData*  ) = NULL;
template< class VertexData , class EdgeData , class FacetData >
void (*HalfEdgeMesh< VertexData , EdgeData , FacetData >::swap_halfedges_function)( EdgeData*   , EdgeData*   ) = NULL;
template< class VertexData , class EdgeData , class FacetData >
void (*HalfEdgeMesh< VertexData , EdgeData , FacetData >::swap_vertices_function )( VertexData* , VertexData* ) = NULL;
#include "HalfEdge.inl"
#endif // HALF_EDGE_INCLUDED
