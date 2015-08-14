cimport dilap.core.vector as dpv

stuff = 'hi'

cdef class quaternion:
    cdef public float w
    cdef public float x    
    cdef public float y
    cdef public float z
    cpdef quaternion copy(self)
    cpdef list to_list(self)
    cpdef tuple to_tuple(self)
    cpdef float magnitude(self)
    cpdef quaternion normalize(self)
    cpdef quaternion flip(self)
    cpdef quaternion conjugate(self)
    cpdef quaternion multiply(self, quaternion q)
    cpdef quaternion translate(self, quaternion sv)
    cpdef quaternion translate_w(self, float tv)
    cpdef quaternion translate_x(self, float tv)
    cpdef quaternion translate_y(self, float tv)
    cpdef quaternion translate_z(self, float tv)
    cpdef quaternion scale(self, quaternion sv)
    cpdef quaternion scale_w(self, float sw)
    cpdef quaternion scale_x(self, float sx)
    cpdef quaternion scale_y(self, float sy)
    cpdef quaternion scale_z(self, float sz)
    cpdef quaternion scale_u(self, float s)
    cpdef quaternion rotate(self, quaternion q)
    cpdef dpv.vector rotate_vector(self, dpv.vector v)
    cpdef quaternion rotate_x(self, float zang)
    cpdef quaternion rotate_y(self, float zang)
    cpdef quaternion rotate_z(self, float zang)

    #cpdef vector reciprocate(self)
    #cpdef vector normalize(self)

#cdef float distance_xy_c(vector v1, vector v2)
#cdef float distance_c(vector v1, vector v2)
#cdef void translate_coords_x_c(list coords, float tv)
#cdef void translate_coords_y_c(list coords, float tv)
#cdef void translate_coords_z_c(list coords, float tv)
#cdef void translate_coords_c(list coords, vector t)
#cdef void scale_coords_c(list coords, vector t)
#cdef void rotate_x_coords_c(list coords, float ang)
#cdef void rotate_y_coords_c(list coords, float ang)
#cdef void rotate_z_coords_c(list coords, float ang)
#cdef float angle_from_xaxis_xy_c(vector v)
#cdef float angle_from_xaxis_c(vector v)
#cdef float angle_between_xy_c(vector v1, vector v2)
#cdef float angle_between_c(vector v1, vector v2)

cdef quaternion zero_c()
# THIS IS INCOMLETE



