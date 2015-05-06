cimport dp_vector as dpv

stuff = 'hi'

cdef class bbox:
    cdef public dpv.vector2d x    
    cdef public dpv.vector2d y
    cdef public dpv.vector2d z


