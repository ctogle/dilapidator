cimport dp_vector as dpv

stuff = 'hi'

cdef class ray:
    cdef public dpv.vector origin
    cdef public dpv.vector direction
    cdef public dpv.vector cast
    cpdef ray copy(self)
    cpdef bint intersect_tri(self,dpv.vector v0,dpv.vector v1,dpv.vector v2)
    cpdef bint intersect_plane(self,dpv.vector r0,dpv.vector n)


