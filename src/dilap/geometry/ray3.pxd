cimport dilap.core.tools as dpr

from dilap.geometry.vec3 cimport vec3
from dilap.geometry.quat cimport quat

stuff = 'hi'










# dilapidators implementation of a ray in R3
cdef class ray3:

    cdef public vec3 o
    cdef public vec3 d
    cdef public vec3 c

    cdef ray3 cp_c(self)
    cdef bint isnear_c(self,ray3 o)
    cdef bint hittri_c(self,vec3 t1,vec3 t2,vec3 t3)
    cdef list hittris_c(self,list tris)
    cdef int hittri_close_c(self,list tris)
    cdef bint hitpln_c(self,vec3 r,vec3 n)

    cpdef ray3 cp(self)
    cpdef bint isnear(self,ray3 o)
    cpdef bint hittri(self,vec3 t1,vec3 t2,vec3 t3)
    cpdef list hittris(self,list tris)
    cpdef int hittri_close(self,list tris)
    cpdef bint hitpln(self,vec3 r,vec3 n)









