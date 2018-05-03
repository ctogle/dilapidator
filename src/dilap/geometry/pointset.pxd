#cimport dilap.core.tools as dpr

from dilap.geometry.vec3 cimport vec3
from dilap.geometry.quat cimport quat

stuff = 'hi'










__doc__ = '''dilapidator\'s implementation of a pointset'''
# dilapidators implementation of a pointset (for quats,vec3s,vec2s,etc.)
cdef class pointset:

    cdef public list ps
    cdef public int pcnt

    cdef list gpscp_c(self,rng)
    cdef list gps_c(self,rng)
    cdef list gpsset_c(self)
    cdef int ap_c(self, np)
    cdef list aps_c(self, list nps)
    cdef int np_c(self, np, e)
    cdef list nps_c(self, list nps, float e)
    cdef int fp_c(self, p, e)
    cdef list fps_c(self, ps, float e)
    cdef bint disjoint_c(self,pointset o)
    cdef pointset trn_c(self,vec3 v)
    cdef pointset rot_c(self,quat q)
    cdef pointset scl_c(self,vec3 o)
    cdef pointset uscl_c(self,float f)

    cpdef list gpscp(self,rng)
    cpdef list gps(self,rng)
    cpdef list gpsset(self)
    cpdef int ap(self, np)
    cpdef list aps(self, list nps)
    cpdef int np(self, np, e)
    cpdef list nps(self, list nps, float e)
    cpdef int fp(self, p, e)
    cpdef list fps(self, ps, float e)
    cpdef bint disjoint(self,pointset o)
    cpdef pointset trn(self,vec3 v)
    cpdef pointset rot(self,quat q)
    cpdef pointset scl(self,vec3 o)
    cpdef pointset uscl(self,float f)
        

 







