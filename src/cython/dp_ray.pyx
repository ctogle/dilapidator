#imports
# cython: profile=True
#cimport cython
cimport dp_vector as dpv
import dp_vector as dpv
cimport dp_bbox as dbb
import dp_bbox as dbb

from libc.math cimport sqrt
from libc.math cimport cos
from libc.math cimport sin
from libc.math cimport tan
from libc.math cimport hypot
cimport numpy as np
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

# return tuple of index,cast vector for the closest triangle to r.origin
cdef tuple intersect_hits_closest_c(ray r,list triangles):
    cdef int tcnt = len(triangles)
    cdef float closest = 100000000000.0
    cdef int hitdex = -1
    cdef dpv.vector hitcast = dpv.nxhat.copy()
    for tdx in range(tcnt):
        tri = triangles[tdx]
        hit = r.intersect_tri(*tri)
        if hit:
            if r.cast.x < closest:
                closest = r.cast.x
                hitdex = tdx
                hitcast = r.cast.copy()
    return (hitdex,hitcast)

cpdef tuple intersect_hits_closest(ray r,list triangles):
    return intersect_hits_closest_c(r,triangles)

# currently hard coded to assume direction == dpv.nzhat...
cdef tuple ray_grid_c(dpv.vector direction,dbb.bbox bb,list faces,float dx):
    cdef double [:] xax = np.arange(bb.x.x,bb.x.y,dx)
    cdef double [:] yax = np.arange(bb.y.x,bb.y.y,dx)
    cdef float z = bb.z.y + 1
    cdef list pts = [(x,y) for x in xax for y in yax]
    cdef list raygrid = [ray(dpv.vector(x,y,z),dpv.nzhat) for x,y in pts]
    cdef list hitfaces = []
    cdef list hitcasts = []
    for zray in raygrid:
        hf,hc = intersect_hits_closest(zray,faces)
        if not hf == -1:
            hitfaces.append(hf)
            hitcasts.append(hc)
    return hitfaces,hitcasts

cpdef tuple ray_grid(dpv.vector direction,dbb.bbox bb,list faces,float dx):
    return ray_grid_c(direction,bb,faces,dx)







