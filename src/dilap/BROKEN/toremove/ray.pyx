#imports
# cython: profile=True
#cimport cython
cimport dilap.core.tools as dpr
import dilap.core.tools as dpr
cimport dilap.core.vector as dpv
import dilap.core.vector as dpv
#import dp_vector as dpv
cimport dilap.core.bbox as dbb
#import dp_bbox as dbb

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

    cdef public dpv.vector origin
    cdef public dpv.vector direction
    cdef public dpv.vector cast

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
        cdef dpv.vector e1 = dpv.v1_v2(v0,v1)
        cdef dpv.vector e2 = dpv.v1_v2(v0,v2)
        cdef dpv.vector T  = self.origin - v0
        cdef dpv.vector P  = dpv.cross(self.direction,e2)
        cdef dpv.vector Q  = dpv.cross(T,e1)
        cdef float denom = dpr.near(dpv.dot(P,e1),0)
        cdef float t = dpr.near(dpv.dot(Q,e2),0)
        cdef float u = dpr.near(dpv.dot(P,T),0)
        cdef float v = dpr.near(dpv.dot(Q,self.direction),0)
        self.cast = dpv.xhat.copy().flip()
        #if denom < epsilon:return 0
        if abs(denom) < epsilon:return 0
        t /= denom
        u /= denom
        v /= denom
        if t < 0.0:return 0
        if u < 0.0 or u > 1:return 0
        if v < 0.0 or u + v > 1:return 0
        self.cast.x = t
        self.cast.y = u
        self.cast.z = v
        return 1

    # tests if this ray intersects a plane 
    # given point in plane r0 and plane normal n
    #
    # currently just checks if parametric line across ray intersects...
    cpdef bint intersect_plane(self,dpv.vector r0,dpv.vector n):
        cdef dpv.vector p0 = self.origin.copy()
        cdef dpv.vector p1 = self.origin.copy().translate(self.direction)
        cdef float denom = n.dot(p1 - p0)
        cdef float numer = n.dot(r0 - p0)
        cdef float t
        if denom == 0:
            print('ray is parallel to the plane')
            # check if p0 or p1 is in the plane
            return 0
        else:t = numer/denom
        self.cast.x = self.origin.x + t*self.direction.x
        self.cast.y = self.origin.y + t*self.direction.y
        self.cast.z = self.origin.z + t*self.direction.z
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







