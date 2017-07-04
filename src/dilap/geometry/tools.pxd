from dilap.geometry.vec3 cimport vec3
from dilap.geometry.quat cimport quat
from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat

cdef float epsilon_c
cdef float epsilonsq_c

cdef float maxfloat_c

stuff = 'hi'



cdef bint isnear_c(float a,float b,float c)
cdef float near_c(float a,float b,float c)
cdef float rad_c(float deg)
cdef float deg_c(float rad)
cdef float clamp_c(float v,float f,float c)
cdef float wrap_c(float v,float f,float c)
cdef bint inrng_c(float a,float b,float c)
cdef float adist_c(float a1,float a2)
cdef tuple circumscribe_tri_c(vec3 p1,vec3 p2,vec3 p3)
cdef float orient2d_c(vec3 a,vec3 b,vec3 c)
cdef float incircle_c(vec3 a,vec3 b,vec3 c,vec3 d)
cdef tuple rot_poly_c(tuple polygon,quat q)
cdef quat q_to_xy_c(vec3 v)
cdef vec3 nrm_c(vec3 c1,vec3 c2,vec3 c3)








