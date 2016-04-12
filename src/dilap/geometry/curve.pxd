cimport dilap.geometry.tools as dpr

from dilap.geometry.vec3 cimport vec3
from dilap.geometry.quat cimport quat

stuff = 'hi'










# dilapidators implementation of a curve in R3
cdef class curve:

    cdef public int meth
    cdef public int segs
    cdef public vec3 tail
    cdef public vec3 tip
    cdef public list pts
    cdef public list nms
    cdef public list tns

    cdef curve cp_c(self)
    cdef curve cpxy_c(self)
    cdef curve cln_c(self)
    cdef void calcone_c(self,vec3 p1,vec3 p2)
    cdef curve calc_c(self)

    cpdef curve cp(self)
    cpdef curve cpxy(self)
    cpdef curve cln(self)
    cpdef calcone(self,vec3 p1,vec3 p2)
    cpdef curve calc(self)










