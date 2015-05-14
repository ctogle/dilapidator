cimport dp_vector as dpv

stuff = 'hi'

cdef class bbox:
    cdef public dpv.vector2d x    
    cdef public dpv.vector2d y
    cdef public dpv.vector2d z
    cdef public list corners
    cdef public list edgenorms
    cdef public int edgecount
    cdef public float radius
    cdef public dpv.vector center
    cpdef bbox _consume_x(self,dpv.vector2d proj)
    cpdef bbox _consume_y(self,dpv.vector2d proj)
    cpdef bbox _consume_z(self,dpv.vector2d proj)
    cpdef bbox _consume(self,bbox other)


