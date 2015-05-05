#imports
# cython: profile=True
#cimport cython
cimport dp_vector as dpv
import dp_vector as dpv

from libc.math cimport sqrt
from libc.math cimport cos
from libc.math cimport sin
from libc.math cimport tan
from libc.math cimport hypot
import numpy as np

stuff = 'hi'

# 3d classes/functions

cdef float epsilon = 0.0000001

cdef class ray:

    def __cinit__(self,dpv.vector origin,dpv.vector direction):
        self.origin = origin 
        self.direction = direction
        self.cast = dpv.xhat.copy().flip()

    def __str__(self):
        strr = 'ray:' + str((self.origin,self.direction))
        return strr

    cpdef ray copy(self):
        cdef ray new = ray(self.origin.copy(),self.direction.copy())
        return new

    # tests if this ray intersects a triangle (assuming culling)
    cpdef bint intersect_tri(self,dpv.vector v0,dpv.vector v1,dpv.vector v2):
        cdef dpv.vector e1 = v1 - v0
        cdef dpv.vector e2 = v2 - v0
        cdef dpv.vector T  = self.origin - v0
        cdef dpv.vector P  = dpv.cross(self.direction,e2)
        cdef dpv.vector Q  = dpv.cross(T,e1)
        cdef float denom = dpv.dot(P,e1)
        cdef float t = dpv.dot(Q,e2)
        cdef float u = dpv.dot(P,T)
        cdef float v = dpv.dot(Q,self.direction)

        self.cast = dpv.xhat.copy().flip()
        if denom < epsilon:return 0
        if u < 0.0 or u > denom:return 0
        if v < 0.0 or u + v > denom:return 0
        self.cast.x = t/denom
        self.cast.y = u/denom
        self.cast.z = v/denom
        return 1

# return the indices of triangles in triangles which are hit by ray r
cdef list intersect_filter_c(ray r,list triangles):
    cdef list hits = []
    cdef int tcnt = len(triangles)
    for tdx in range(tcnt):
        tri = triangles[tdx]
        hit = r.intersect_tri(*tri)
        if hit:hits.append(tdx)
    return hits

cpdef list intersect_filter(ray r,list triangles):
    return intersect_filter_c(r,triangles)

# return tuples of index,cast vector for triangles which are hit by ray r 
cdef list intersect_hits_c(ray r,list triangles):
    cdef list hits = []
    cdef int tcnt = len(triangles)
    for tdx in range(tcnt):
        tri = triangles[tdx]
        hit = r.intersect_tri(*tri)
        if hit:hits.append((tdx,r.cast.copy()))
    return hits

cpdef list intersect_hits(ray r,list triangles):
    return intersect_hits_c(r,triangles)


