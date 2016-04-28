cimport dilap.geometry.tools as gtl
import dilap.geometry.tools as gtl
cimport dilap.geometry.pointset as dps

cimport dilap.geometry.vec3 as dpv
from dilap.geometry.vec3 cimport vec3
from dilap.geometry.quat cimport quat

stuff = 'hi'



cdef class triangulation:

    cdef public vec3 p0
    cdef public vec3 pn
    cdef public dps.pointset points
    cdef public list triangles
    cdef public int tricnt
    cdef public dict eg_tri_lookup
    cdef public list ghosts
    cdef public int ghostcnt
    cdef public dict eg_ghost_lookup

    cdef int point_on_boundary(self,int u)
    cdef bint locally_delaunay(self,int u,int v)
    cdef int adjacent(self,int u,int v)
    cdef tuple flip_edge(self,int u,int v)
    cdef void insert_vertex(self,int u,int v,int w,int x)
    cdef void insert_ghost_vertex(self,int u,int v,int w,int x)
    cdef void add_triangle(self,int u,int v,int w)
    cdef void delete_triangle(self,int u,int v,int w)
    cdef void add_ghost(self,int u,int v)
    cdef void delete_ghost(self,int u,int v)
    cdef void dig_cavity(self,int u,int v,int w)






