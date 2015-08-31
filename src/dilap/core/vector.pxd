cimport dilap.core.quaternion as dpq
#cimport dp_bbox as dbb

#import dp_bbox as dbb

stuff = 'hi'

cdef class vector2d:
    cdef public float x    
    cdef public float y
    cpdef bint near(self,vector2d other)
    cpdef vector2d copy(self)
    cpdef list to_list(self)
    cpdef tuple to_tuple(self)
    cpdef vector2d flip(self)
    cpdef vector2d reciprocate(self)
    cpdef vector2d translate(self, vector2d tv)
    cpdef vector2d translate_x(self, float tx)
    cpdef vector2d translate_y(self, float ty)
    cpdef vector2d scale(self, vector2d sv)
    cpdef vector2d scale_x(self, float sx)
    cpdef vector2d scale_y(self, float sy)

cdef vector2d midpoint2d_c(vector2d v1, vector2d v2)

cdef class vector:
    cdef public float x    
    cdef public float y
    cdef public float z
    cpdef bint neighborhood(self,vector other,float epsilon)
    cpdef bint near_xy(self,vector other)
    cpdef bint near(self,vector other)
    cpdef int nearest(self,list others)
    cpdef vector2d xy2d(self)
    cpdef vector2d xz2d(self)
    cpdef vector2d yz2d(self)
    cpdef vector xy(self)
    cpdef vector xz(self)
    cpdef vector yz(self)
    cpdef vector copy(self)
    cpdef vector rotate(self, dpq.quaternion q)
    cpdef vector rotate_x(self, float zang)
    cpdef vector rotate_y(self, float zang)
    cpdef vector rotate_z(self, float zang)
    cpdef vector translate_x(self, float tv)
    cpdef vector translate_y(self, float tv)
    cpdef vector translate_z(self, float tv)
    cpdef vector translate(self, vector sv)
    cpdef vector scale(self, vector sv)
    cpdef vector scale_u(self, float s)
    cpdef vector project_plane(self, vector , vector n)
    cpdef vector cross(self, vector v)
    cpdef float dot(self, vector v)
    cpdef vector flip(self)
    cpdef vector reciprocate(self)
    cpdef vector normalize(self)
    cpdef float magnitude2(self)
    cpdef float magnitude(self)
    cpdef list to_list(self)
    cpdef tuple to_tuple(self)
    cpdef vector linterpolate(self,vector other,float delta)
    cpdef vector in_basis(self,vector b1,vector b2,vector b3)

cdef vector xhat
cdef vector yhat
cdef vector zhat
cdef vector nxhat
cdef vector nyhat
cdef vector nzhat

cdef int find_closest_xy_c(vector one, list bunch, int bcnt, float close_enough)
cdef float distance_xy_c(vector v1, vector v2)
cdef float distance_c(vector v1, vector v2)
cdef float distance_to_edge_c(vector pt, vector e1, vector e2, vector nm)
cdef float dot_c(vector v1, vector v2)
cdef void translate_coords_x_c(list coords, float tv)
cdef void translate_coords_y_c(list coords, float tv)
cdef void translate_coords_z_c(list coords, float tv)
cdef list translate_coords_c(list coords, vector t)
cdef void scale_coords_c(list coords, vector t)
cdef void rotate_x_coords_c(list coords, float ang)
cdef void rotate_y_coords_c(list coords, float ang)
cdef void rotate_z_coords_c(list coords, float ang)
cdef list edge_normals_xy_c(list verts)
cdef float distance_to_border_xy_c(vector pt,list border)
#cdef float angle_from_xaxis_xy_c(vector v)
#cdef float angle_from_xaxis_c(vector v)
#cdef float angle_between_xy_c(vector v1, vector v2)
#cdef float angle_between_c(vector v1, vector v2)
#cdef float signed_angle_between_c(vector v1, vector v2, vector n)
#cdef float signed_angle_between_xy_c(vector v1, vector v2)

cdef vector2d flip2d_c(vector2d f)

cdef vector zero_c()
cdef vector one_c()
cdef vector flip_c(vector f)
cdef vector normalize_c(vector v)
cdef vector cross_c(vector v1, vector v2)
cdef vector v1_v2_xy_c(vector v1, vector v2)
cdef vector v1_v2_c(vector v1, vector v2)
cdef vector vzip_c(vector v1, vector v2)
cdef vector midpoint_c(vector v1, vector v2)
cdef vector com(list coords)
cdef vector com_weighted(list coords,list weights)
cdef vector barymetric_to_world_c(float u,float v,vector v0,vector v1,vector v2)
cdef vector2d project_coords_c(list coords, vector axis)

cdef list line_normals_c(list verts)



