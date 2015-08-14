cimport dilap.core.vector as dpv
cimport dilap.core.pointset as dps

cdef class triangulation:

    cdef public dps.pointset points
    cdef public list triangles
    cdef public int tricnt
    cdef public dict eg_tri_lookup
    cdef public list ghosts
    cdef public int ghostcnt
    cdef public dict eg_ghost_lookup
    cdef bint locally_delaunay(self,int u,int v)
    cdef tuple flip_edge(self,int u,int v)
    cdef void add_triangle(self,int u,int v,int w)
    cdef void delete_triangle(self,int u,int v,int w)
    cdef void add_ghost(self,int u,int v)
    cdef void delete_ghost(self,int u,int v)
    cdef void dig_cavity(self,int u,int v,int w)
    cdef void insert_vertex(self,int u,int v,int w,int x)
    cdef void insert_ghost_vertex(self,int u,int v,int w,str x)
    cdef void ghost_border(self,list edges)


