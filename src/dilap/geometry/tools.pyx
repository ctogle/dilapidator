#imports
# cython: profile=True
#cimport cython

from dilap.geometry.vec3 cimport vec3
from dilap.geometry.quat cimport quat
from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat

import sys,math,numpy
import matplotlib.pyplot as plt

PI = numpy.pi
PI2 = PI/2.0
PI4 = PI/4.0
twoPI = PI*2.0
threePI = PI*3.0
threePI2 = PI*3.0/2.0
threePI4 = PI*3.0/4.0

epsilon   = 0.0001
epsilonsq = epsilon*epsilon
cdef float epsilon_c   = 0.0001
cdef float epsilonsq_c = epsilon*epsilon

maxfloat = sys.float_info.max
cdef float maxfloat_c = maxfloat

stuff = 'hi'

__doc__ = '''General purpose tool functions...'''

###############################################################################
### c space
###############################################################################

# if a is within c of b, return True
# else return False
cdef bint isnear_c(float a,float b):
    cdef d = a-b
    d *= d
    if d < epsilon_c:return 1
    else:return 0

# if a is within c of b, return b
# else return a
cdef float near_c(float a,float b):
    cdef d = a-b
    d *= d
    if d < epsilon_c:return b
    else:return a

# convert an angle from degrees to radians
cdef float rad_c(float deg):return PI*deg/180.0

# convert an angle from radians to degrees
cdef float deg_c(float rad):return 180.0*rad/PI

# keep the value val bounded by f and c by flooring
cdef float clamp_c(float v,float f,float c):
    if v < f:return f
    elif v > c:return c
    else:return v

# keep the value val bounded by f and c by wrapping around
cdef float wrap_c(float v,float f,float c):
    period = c - f
    while v < f:v += period
    while v > c:v -= period
    else:return v                                  

# is a on the interior of (a,b) given an error of d
cdef bint inrng_c(float a,float b,float c):
    cdef float r = near_c(near_c(a,b),c)
    #cdef bint inr = r >= b and r <= c
    cdef bint inr = r > b and r < c
    return inr

# return the shortest angular distance between two angles
cdef float adist_c(float a1,float a2):
    cdef float da = wrap_c(a1-a2,0.0,twoPI)
    return da if da < PI else twoPI - da

# given 3 points and a plane, determine the center and radius
# of the circumcenter found in the plane plane by projecting p1,p2,p3
# in the plane
cdef tuple circumscribe_tri_c(vec3 p1,vec3 p2,vec3 p3):
    cdef vec3 cp1 = p1.cpxy()
    cdef vec3 cp2 = p2.cpxy()
    cdef vec3 cp3 = p3.cpxy()
    cdef vec3 e1 = cp3.tovxy_c(cp1)
    cdef vec3 e2 = cp3.tovxy_c(cp2)
    cdef float th = e1.angxy_c(e2)
    cdef float cr = cp1.dxy_c(cp2)/(2*numpy.sin(th))
    cdef vec3 cp = e2.cp_c().uscl_c(e1.mag2_c())-e1.cp_c().uscl_c(e2.mag2_c())
    cdef vec3 fp = p3+cp.crs_c(e1.crs_c(e2)).uscl_c(1.0/(2.0*(e1.crs_c(e2).mag2())))
    return fp,cr

# return the signed area of the triangle created 
# by the vectors a-c,b-c
# return 0 if a,b,c are colinear
# this assumes theyre in the xy plane!!!
cdef float orient2d_c(vec3 a,vec3 b,vec3 c):
    cdef float m11 = a.x-c.x
    cdef float m12 = a.y-c.y
    cdef float m21 = b.x-c.x
    cdef float m22 = b.y-c.y
    cdef float det = m11*m22-m12*m21
    return near_c(det,0)

# determine if d is inside the circumcircle of the triangle a,b,c
cdef float incircle_c(vec3 a,vec3 b,vec3 c,vec3 d):
    cdef float m11 = a.x-d.x
    cdef float m12 = a.y-d.y
    cdef float m13 = m11*m11 + m12*m12
    cdef float m21 = b.x-d.x
    cdef float m22 = b.y-d.y
    cdef float m23 = m21*m21 + m22*m22
    cdef float m31 = c.x-d.x
    cdef float m32 = c.y-d.y
    cdef float m33 = m31*m31 + m32*m32
    cdef float det1 = m11*(m22*m33-m23*m32)
    cdef float det2 = m12*(m21*m33-m23*m31)
    cdef float det3 = m13*(m21*m32-m22*m31)
    cdef float inc = near_c(det1 - det2 + det3,0)
    return inc

# rotate a polygon: (extbnd,(holes...)) by a quaternion q
cdef tuple rot_poly_c(tuple polygon,quat q):
    cdef tuple ebnd
    cdef tuple ibnds
    cdef tuple ibnd
    ebnd,ibnds = polygon
    cdef int elen = len(ebnd)
    cdef int ex
    cdef int islen = len(ibnds)
    cdef int ibx
    cdef int ilen 
    cdef int ix
    for ex in range(elen):
        ebnd[ex].rot(q)
    for ibx in range(islen):
        ibnd = ibnds[ibx]
        ilen = len(ibnd)
        for ix in range(ilen):
            ibnd[ix].rot(q)
    return polygon

# given a vector, return a quaternion that rotates it to coincide with zhat
cdef quat q_to_xy_c(vec3 v):
    cdef quat prot
    if v.isnear(vec3(0,0,-1)):prot = quat(1,0,0,0).av_c(PI,vec3(1,0,0))
    elif not v.isnear(vec3(0,0,1)):prot = quat(1,0,0,0).uu_c(v,vec3(0,0,1))
    else:prot = quat(1,0,0,0)
    return prot

# return a vector normal to the plane containing c1,c2,c3
# returns 0 if c1,c2,c3 are colinear
cdef vec3 nrm_c(vec3 c1,vec3 c2,vec3 c3):
    cdef vec3 c1c2 = c1.tov_c(c2).nrm_c()
    cdef vec3 c2c3 = c2.tov_c(c3).nrm_c()
    cdef vec3 cn = c1c2.crs_c(c2c3).nrm_c()
    return cn

# return a vector tanget to the plane containing c1,c2,c3
cdef vec3 tng_c(vec3 c1,vec3 c2,vec3 c3):
    cdef vec3 tn = c1.tov_c(c2).nrm_c()
    return tn

###############################################################################
### python space
###############################################################################

# if a is within c of b, return True
# else return False
cpdef bint isnear(float a,float b):
    '''determine if a is within a neighborhood c of b'''
    return isnear_c(a,b)

# if a is within c of b, return b
# else return a
cpdef float near(float a,float b):
    '''effectively round a to b if within a neighborhood c'''
    return near_c(a,b)

# convert an angle from degrees to radians
cpdef float rad(float deg):
    '''convert an angle from degrees to radians'''
    return rad_c(deg)

# convert an angle from radians to degrees
cpdef float deg(float rad):
    '''convert an angle from radians to degrees'''
    return deg_c(rad)

# keep the value val bounded by f and c by flooring
cpdef float clamp(float v,float f,float c):
    '''clamp a float between two floats using closed boundaries'''
    return clamp_c(v,f,c)

# keep the value val bounded by f and c by wrapping around
cpdef float wrap(float v,float f,float c):
    '''clamp a float between two floats using periodic boundaries'''
    return wrap_c(v,f,c)

# is a on the interior of (a,b) given an error of d
cpdef bint inrng(float a,float b,float c):
    '''determine if a value is on an open interval'''
    return inrng_c(a,b,c)

# return the shortest angular distance between two angles
cpdef float adist(float a1,float a2):
    '''find the angular distance between two angles'''
    return adist_c(a1,a2)

# given 3 points and a plane, determine the center and radius
# of the circumcenter found in the plane plane by projecting p1,p2,p3
# in the plane
cpdef tuple circumscribe_tri(vec3 p1,vec3 p2,vec3 p3):
    '''return the circumcenter and circumradius of a triangle'''
    return circumscribe_tri_c(p1,p2,p3)

# return the signed area of the triangle created 
# by the vectors a-c,b-c
# return 0 if a,b,c are colinear
# this assumes theyre in the xy plane!!!
cpdef float orient2d(vec3 a,vec3 b,vec3 c):
    '''determine whether a,b,c form a positively oriented triangle'''
    return orient2d_c(a,b,c)

# determine if d is inside the circumcircle of the triangle a,b,c
cpdef float incircle(vec3 a,vec3 b,vec3 c,vec3 d):
    '''determine if d is inside the circumcircle of the triangle a,b,c'''
    return incircle_c(a,b,c,d)

# rotate a polygon: (extbnd,(holes...)) by a quaternion q
cpdef tuple rot_poly(tuple polygon,quat q):
    '''rotate a concave polygon and its holes by a quaternion'''
    return rot_poly_c(polygon,q)

# given a vector, return a quaternion that rotates it to coincide with zhat
cpdef quat q_to_xy(vec3 v):
    '''provide a quaternion which rotates a vector onto zhat'''
    return q_to_xy_c(v)

# return a vector normal to the plane containing c1,c2,c3
# returns 0 if c1,c2,c3 are colinear
cpdef vec3 nrm(vec3 c1,vec3 c2,vec3 c3):
    '''find a vector normal to the plane containing c1,c2,c3'''
    return nrm_c(c1,c2,c3)

# return a vector tanget to the plane containing c1,c2,c3
cpdef vec3 tng(vec3 c1,vec3 c2,vec3 c3):
    '''find a vector tangent to the plane containing c1,c2,c3'''
    return tng_c(c1,c2,c3)

def lexicographic(unops):
    ufn = unops[:]
    ops = []
    while ufn:
        xmin = 1e10
        ymin = 1e10
        ux = None
        for x in range(len(ufn)):
            u = ufn[x]
            if u.x < xmin:
                ux = x
                xmin = u.x
            elif u.x == xmin and u.y < ymin:
                ux = x
                xmin = u.x
                ymin = u.y
        ops.append(ufn.pop(ux))
    return ops
    







