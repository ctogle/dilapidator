cimport dilap.core.tools as dpr

from dilap.geometry.vec3 cimport vec3
from dilap.geometry.quat cimport quat

stuff = 'hi'










# dilapidator\'s implementation of a transform object in R3'''
cdef class tform:

    cdef public vec3 pos    
    cdef public quat rot
    cdef public vec3 scl

    cdef tform cp_c(self)
    cdef tform true_c(self,parent)

    cpdef tform cp(self)
    cpdef tform true(self,parent)









