cimport dilap.geometry.vec3 as dpv
from dilap.geometry.vec3 cimport vec3
from dilap.geometry.quat cimport quat
from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat

cdef float epsilon_c
cdef float epsilonsq_c

cdef float maxfloat_c

stuff = 'hi'



cdef bint isnear_c(float a,float b)
cdef float near_c(float a,float b)
cdef float rad_c(float deg)
cdef float deg_c(float rad)
cdef float clamp_c(float v,float f,float c)
cdef float wrap_c(float v,float f,float c)
cdef bint inrng_c(float a,float b,float c)
cdef float adist_c(float a1,float a2)
cdef tuple circumscribe_tri_c(vec3 p1,vec3 p2,vec3 p3)
#cdef float ang_c(vec3 v1,vec3 v2)
#cdef float ang_xy_c(vec3 v1,vec3 v2)
#cdef float sang_xy_c(vec3 v1,vec3 v2)
#cdef float xang_xy_c(vec3 v)
#cdef tuple bary_xy_c(vec3 pt,vec3 a,vec3 b,vec3 c)

cdef bint onseg_xy_c(vec3 p,vec3 s1,vec3 s2)
cdef bint inseg_xy_c(vec3 p,vec3 s1,vec3 s2)
cdef float orient2d_c(vec3 a,vec3 b,vec3 c)
cdef float incircle_c(vec3 a,vec3 b,vec3 c,vec3 d)
cdef bint intri_xy_c(vec3 pt,vec3 a,vec3 b,vec3 c)
cdef bint onconcave_xy_c(vec3 pt,tuple py)
cdef bint inconcave_xy_c(vec3 pt,tuple poly)
cdef bint segs_isect_perp_c(vec3 s11,vec3 s12,vec3 s21,vec3 s22)
cdef bint segs_isect_int_c(vec3 s11,vec3 s12,vec3 s21,vec3 s22)
cdef bint polyinpoly_c(tuple p1,tuple p2)
cdef bint poly_isect_c(tuple p1,tuple p2)
#cdef list pline_c(vec3 s,vec3 e,int n)
#cdef list pring_c(float r,int n)
cdef tuple rot_poly_c(tuple polygon,quat q)
cdef quat q_to_xy_c(vec3 v)
cdef vec3 nrm_c(vec3 c1,vec3 c2,vec3 c3)
cdef vec3 poly_nrm_c(tuple poly)
cdef int winding_c(vec3 pt,tuple py)



#cdef bint cyclic_permutation_c(seq1,seq2)
#cdef float distance_to_line_c(dpv.vector pt,dpv.vector e1,dpv.vector e2,dpv.vector nm)
#cdef float angle_from_xaxis_xy_c(dpv.vector v)
#cdef float signed_angle_between_c(dpv.vector v1,dpv.vector v2,dpv.vector n)
#cdef float signed_angle_between_xy_c(dpv.vector v1,dpv.vector v2)
#cdef vec3 tangent_c(vec3 c1,vec3 c2,vec3 c3)
#cdef tuple translate_polygon_c(tuple polygon,vec3 tv)
#cdef tuple rotate_x_polygon_c(tuple polygon,float a)
#cdef list square_c(float l,float w,p = ?,a = ?)
#cdef bint inside_c(vec3 pt,list corners)
#cdef bint inside_circle_c(dpv.vector pt,dpv.vector c,float r)
#cdef tuple circumscribe_tri_c(dpv.vector p1,dpv.vector p2,dpv.vector p3)
#cdef bint segments_intersect_noncolinear_c(dpv.vector s11,dpv.vector s12,dpv.vector s21,dpv.vector s22)
#cdef bint insegment_xy_c(dpv.vector p,dpv.vector s1,dpv.vector s2)
#cdef tuple barycentric_xy_c(dpv.vector pt,dpv.vector a,dpv.vector b,dpv.vector c)
#cdef dpv.vector2d barycentric_c(dpv.vector pt,dpv.vector a,dpv.vector b,dpv.vector c)
#cdef bint intriangle_xy_c(dpv.vector pt,dpv.vector a,dpv.vector b,dpv.vector c)
#cdef bint intriangle_c(dpv.vector pt,dpv.vector a,dpv.vector b,dpv.vector c)
#cdef bint inconvex_c(dpv.vector pt,tuple poly)
#cdef bint inconcave_xy_c(dpv.vector pt,tuple poly)
#cdef bint concaves_contains_c(tuple p1,tuple p2)
#cdef dpv.vector find_y_apex_c(list pts)
#cdef dpv.vector find_x_apex_c(list pts)
#cdef dpv.vector sweep_search_c(list pts,dpv.vector center,tangent = ?)
#cdef list pts_to_convex_xy_c(list pts)
#cdef list inflate_c(list convex,float radius)
#cdef list offset_faces_c(list faces,int offset)
#cdef float orient2d_c(dpv.vector a,dpv.vector b,dpv.vector c)
#cdef float orient3d_c(dpv.vector a,dpv.vector b,dpv.vector c,dpv.vector d)
#cdef float incircle_c(dpv.vector a,dpv.vector b,dpv.vector c,dpv.vector d)
#cdef float insphere_c(dpv.vector a,dpv.vector b,dpv.vector c,dpv.vector d,dpv.vector e)








