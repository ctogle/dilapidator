cimport dilap.core.vector as dpv
cimport dilap.core.quaternion as dpq
cimport dilap.core.tools as dpr
cimport dilap.core.pointset as dps

cdef class triangulation:

    cdef public dpv.vector p0
    cdef public dpv.vector pn
    cdef public dps.pointset points
    cdef public list triangles
    cdef public int tricnt
    cdef public dict eg_tri_lookup
    cdef public list ghosts
    cdef public int ghostcnt
    cdef public dict eg_ghost_lookup

    cpdef skinny_triangle(self,int u,int v,int w,float h)
    cdef int point_on_boundary(self,int u)
    cdef bint locally_delaunay(self,int u,int v)
    cdef int adjacent(self,int u,int v)

    cdef tuple flip_edge(self,int u,int v)
    cdef void add_triangle(self,int u,int v,int w)
    cdef void delete_triangle(self,int u,int v,int w)
    cdef void add_ghost(self,int u,int v)
    cdef void delete_ghost(self,int u,int v)
    cdef void dig_cavity(self,int u,int v,int w)
    cdef void insert_vertex(self,int u,int v,int w,int x)
    cdef void insert_ghost_vertex(self,int u,int v,int w,int x)

cdef float PI

cdef tuple triangulate_c(tuple ebnd,tuple ibnds,float hmin)





