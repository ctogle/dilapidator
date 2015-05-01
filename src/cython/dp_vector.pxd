

stuff = 'hi'

cdef class vector2d:
    cdef public float x    
    cdef public float y
    cpdef vector2d copy(self)
    cpdef list to_list(self)
    cpdef scale(self, vector2d sv)
    cpdef translate(self, vector2d sv)

cdef vector2d midpoint2d_c(vector2d v1, vector2d v2)

cdef class vector:
    cdef public float x    
    cdef public float y
    cdef public float z
    cpdef vector copy(self)
    cpdef scale(self, vector sv)
    cpdef vector translate(self, vector sv)
    cpdef vector normalize(self)
    cpdef float magnitude(self)
    cpdef list to_list(self)

cdef int find_closest_xy_c(vector one, list bunch, int bcnt, float close_enough)
cdef float distance_xy_c(vector v1, vector v2)
cdef float distance_c(vector v1, vector v2)
cdef float dot_c(vector v1, vector v2)
cdef void translate_coords_c(list coords, vector t)
cdef float angle_from_xaxis_xy_c(vector v)
cdef float angle_from_xaxis_c(vector v)

cdef vector cross_c(vector v1, vector v2)
cdef vector v1_v2_c(vector v1, vector v2)
cdef vector vzip_c(vector v1, vector v2)
cdef vector midpoint_c(vector v1, vector v2)
cdef vector com(list coords)




