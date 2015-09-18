cimport dilap.core.tools as dpr

from dilap.geometry.vec3 cimport vec3

stuff = 'hi'










# dilapidators implementation of a quaternion in R3
cdef class quat:

    cdef public float w
    cdef public float x    
    cdef public float y
    cdef public float z

    cdef quat av_c(self,float a,vec3 v)
    cdef quat uu_c(self,vec3 x,vec3 y)
    cdef quat cp_c(self)
    cdef bint isnear_c(self,quat o)
    cdef float mag2_c(self)
    cdef float mag_c(self)
    cdef quat nrm_c(self)
    cdef quat flp_c(self)
    cdef quat scl_c(self,float s)
    cdef quat cnj_c(self)
    cdef quat rot_c(self,quat o)

    cpdef quat av(self,float a,vec3 v)
    cpdef quat uu(self,vec3 x,vec3 y)
    cpdef quat cp(self)
    cpdef bint isnear(self,quat o)
    cpdef float mag2(self)
    cpdef float mag(self)
    cpdef quat nrm(self)
    cpdef quat flp(self)
    cpdef quat scl(self,float s)
    cpdef quat cnj(self)
    cpdef quat rot(self,quat o)










