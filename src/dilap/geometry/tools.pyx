#imports
# cython: profile=True
#cimport cython

cimport dilap.geometry.vec3 as dpv
from dilap.geometry.vec3 cimport vec3
from dilap.geometry.quat cimport quat
from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat

#import dilap.core.base as db
#cimport dilap.core.ray as dr
#import dilap.core.ray as dr

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
    #if abs(a-b) < epsilon:return 1
    cdef d = a-b
    d *= d
    if d < epsilonsq_c:return 1
    else:return 0

# if a is within c of b, return b
# else return a
cdef float near_c(float a,float b):
    #if abs(a-b) < c:return b
    cdef d = a-b
    d *= d
    if d < epsilonsq_c:return b
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

# compute the center of mass for a set of vectors
cdef vec3 com_c(ps):
    cdef int pcnt = len(ps)
    cdef int px
    cdef float pcntf = float(pcnt)
    cdef vec3 n = vec3(0,0,0)
    for px in range(pcnt):n.trn_c(ps[px])
    return n.uscl_c(1.0/pcntf)

# find the angle between two vectors 
cdef float ang_c(vec3 v1,vec3 v2):
    cdef vec3 n1 = v1.cp().nrm()
    cdef vec3 n2 = v2.cp().nrm()
    cdef float n12dot = n1.x*n2.x + n1.y*n2.y + n1.z*n2.z
    cdef float ang = 0.0
    if   isnear_c(n12dot, 1.0):pass
    elif isnear_c(n12dot,-1.0):ang = PI
    else:ang = numpy.arccos(n12dot)
    return ang                    

# find the signed angle between two vectors, given a vector pointing "up"
cdef float sang_c(vec3 v1,vec3 v2,vec3 n):
    cdef vec3 n1 = v1.cp_c().nrm_c()
    cdef vec3 n2 = v2.cp_c().nrm_c()
    cdef vec3 vn = n1.crs_c(n2)
    cdef float vdot = vn.x*n.x + vn.y*n.y + vn.z*n.z
    cdef float n12dot = n1.x*n2.x + n1.y*n2.y + n1.z*n2.z
    cdef float ang = 0.0
    if   isnear_c(n12dot, 1.0):pass
    elif isnear_c(n12dot,-1.0):ang = PI
    else:ang = numpy.arccos(n12dot)
    if vdot < 0.0:ang *= -1.0
    return ang                    

# find the angle between the xy projections of two vectors 
cdef float ang_xy_c(vec3 v1,vec3 v2):
    cdef vec3 n1 = v1.xy().normalize()
    cdef vec3 n2 = v2.xy().normalize()
    cdef float n12dot = n1.x*n2.x + n1.y*n2.y
    cdef float ang = 0.0
    if   isnear_c(n12dot,-1):pass
    elif isnear_c(n12dot, 1):ang = PI
    else:ang = numpy.arccos(n12dot)
    return ang                    

# find the signed angle between the xy projections of two vectors
cdef float sang_xy_c(vec3 v1,vec3 v2):
    cdef vec3 n1 = v1.xy().nrm()
    cdef vec3 n2 = v2.xy().nrm()
    cdef float n12dot = n1.x*n2.x + n1.y*n2.y
    cdef float ang = 0.0
    cdef float vn = n1.x*n2.y-n1.y*n2.x
    if   isnear_c(n12dot,-1):pass
    elif isnear_c(n12dot, 1):ang = PI
    else:ang = numpy.arccos(n12dot)
    if vn < 0.0:ang *= -1.0
    return ang                    

# find the positive angle between the xy projection of a vector and the x axis
cdef float xang_xy_c(vec3 v):
    cdef vec3 nv = v.xy_c().nrm_c()
    cdef float ang = 0.0
    if   isnear_c(nv.x, 1):pass
    elif isnear_c(nv.x,-1):ang = PI
    else:ang = numpy.arccos(nv.x)
    if nv.y < 0.0:ang = twoPI - ang
    return ang

# calculate the barycentric coordinates of the point pt for the triangle abc
# assume all points are in the xy plane 
cdef tuple bary_xy_c(vec3 pt,vec3 a,vec3 b,vec3 c): 
    cdef float v0x =  c.x-a.x
    cdef float v0y =  c.y-a.y
    cdef float v1x =  b.x-a.x
    cdef float v1y =  b.y-a.y
    cdef float v2x = pt.x-a.x
    cdef float v2y = pt.y-a.y
    cdef float dot00 = v0x*v0x + v0y*v0y
    cdef float dot01 = v0x*v1x + v0y*v1y
    cdef float dot02 = v0x*v2x + v0y*v2y
    cdef float dot11 = v1x*v1x + v1y*v1y
    cdef float dot12 = v1x*v2x + v1y*v2y
    cdef float denom = (dot00 * dot11 - dot01 * dot01)
    cdef float invdenom = 1.0 / denom
    cdef float u = (dot11 * dot02 - dot01 * dot12) * invdenom
    cdef float v = (dot00 * dot12 - dot01 * dot02) * invdenom
    return u,v

# determine if a point lies on the interior of a line segment or its endpoints
cdef bint onseg_xy_c(vec3 p,vec3 s1,vec3 s2):
    if not orient2d_c(p,s1,s2) == 0:return 0
    if p.isnear(s1) or p.isnear(s2):return 1
    else:return inseg_xy(p,s1,s2)

# determine if a point lies on the interior of a line segment
cdef bint inseg_xy_c(vec3 p,vec3 s1,vec3 s2):
    if not orient2d_c(p,s1,s2) == 0:return 0
    if p.isnear(s1) or p.isnear(s2):return 0
    #if (s1.x != s2.x): # S is not  vertical
    if not isnear_c(s1.x,s2.x): # S is not  vertical
        if (s1.x <= p.x and p.x <= s2.x):return 1
        if (s1.x >= p.x and p.x >= s2.x):return 1
    else: # S is vertical, so test y  coordinate
        if (s1.y <= p.y and p.y <= s2.y):return 1
        if (s1.y >= p.y and p.y >= s2.y):return 1
    return 0

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

# determine if the point pt is inside the triangle abc
# assume all points are in the xy plane 
cdef bint intri_xy_c(vec3 pt,vec3 a,vec3 b,vec3 c):
    cdef float u,v
    u,v = bary_xy_c(pt,a,b,c)
    if u > 0 or abs(u) < epsilon:
        if v > 0 or abs(v) < epsilon:
            if 1-u-v > 0 or abs(1-u-v) < epsilon:
                return 1
    return 0

# determine if the point pt is inside the concave polygon poly or on its boundary
cdef bint onconcave_xy_c(vec3 pt,tuple py):
    cdef int wn = 0
    cdef int px
    cdef int pcnt = len(py)
    cdef float read
    for px in range(pcnt):
        read = orient2d_c(pt,py[px-1],py[px])
        if read == 0:
            if onseg_xy_c(pt,py[px-1],py[px]):return 1
        if py[px-1].y <= pt.y:
            if py[px].y > pt.y:
                if read > 0:wn += 1
        else:
            if py[px].y <= pt.y:
                if read < 0:wn -= 1
    return not wn == 0

# determine if the point pt is inside the concave polygon poly
cdef bint inconcave_xy_c(vec3 pt,tuple poly):
    return not winding_c(pt,poly) == 0

# given line segment s1, line segment s2
# does s1 overlap the interior or s2?
# a segment is a tuple of two points
cdef bint segs_isect_perp_c(vec3 s11,vec3 s12,vec3 s21,vec3 s22):
    p,q = s11,s21
    r = s11.tov_c(s12)
    s = s21.tov_c(s22)
    qmp = q-p
    rcs = r.crs_c(s)
    rcsmag = rcs.mag_c()
    qmpcr = qmp.crs_c(r)
    qmpcrmag = qmpcr.mag_c()
    rmag2 = r.mag2_c()
    if isnear_c(rcsmag,0) and isnear_c(qmpcrmag,0):return 0
    elif isnear_c(rcsmag,0) and not isnear_c(qmpcrmag,0):return 0
    elif not isnear_c(rcs.z,0):
        u = near_c(near_c(       qmpcr.z/rcs.z,0),1)
        t = near_c(near_c(qmp.crs_c(s).z/rcs.z,0),1)
        if (u == 0 or u == 1) and (t == 0 or t == 1):
            #if include_endpoints:return 1
            #else:return 0
            return 0
        if not inrng_c(u,0,1) or not inrng_c(t,0,1):return 0
        else:return 1

# given line segment s1, line segment s2
# does s1 overlap the interior or s2?
# a segment is a tuple of two points
#cdef bint segments_intersect_interior_c(dpv.vector s11,dpv.vector s12,dpv.vector s21,dpv.vector s22):
cdef bint segs_isect_int_c(vec3 s11,vec3 s12,vec3 s21,vec3 s22):
    p,q = s11,s21
    r = s11.tov_c(s12)
    s = s21.tov_c(s22)
    qmp = q-p
    rcs = r.crs_c(s)
    rcsmag = rcs.mag_c()
    qmpcr = qmp.crs_c(r)
    qmpcrmag = qmpcr.mag_c()
    rmag2 = r.mag2_c()
    if isnear_c(rcsmag,0) and isnear_c(qmpcrmag,0):
        t0 = near_c(near_c(qmp.dot_c(r)/rmag2,0),1)
        t1 = near_c(near_c(t0+s.dot_c(r)/rmag2,0),1)
        if near_c(s.dot_c(r),0) < 0:t0,t1 = t1,t0
        if inrng_c(t0,0,1) or inrng_c(t1,0,1):return 1
        elif inrng_c(0,t0,t1) and inrng_c(1,t0,t1):return 1
        elif (t0 == 0 and t1 == 1) or (t0 == 1 and t1 == 0):return 1
        else:return 0
    elif isnear_c(rcsmag,0) and not isnear_c(qmpcrmag,0):return 0
    elif not isnear_c(rcs.z,0):
        u = near_c(near_c(       qmpcr.z/rcs.z,0),1)
        t = near_c(near_c(qmp.crs_c(s).z/rcs.z,0),1)
        if (u == 0 or u == 1) and (t == 0 or t == 1):
            #if include_endpoints:return 1
            #else:return 0
            return 0
        if not inrng_c(u,0,1) or not inrng_c(t,0,1):return 0
        else:return 1

# given concave polygon p1, concave polygon p2
# does p1 contain the entire interior of p2?
# a polygon is a tuple of points
# NOTE: this assumes com(p2) is inside p2!
cdef bint polyinpoly_c(tuple p1,tuple p2):
    cdef int p2cnt = len(p2)
    cdef int px
    cdef vec3 i2 = com_c(p2)
    if not inconcave_xy_c(i2,p1):return 0
    for px in range(p2cnt):
        if not onconcave_xy_c(p2[px],p1):
            return 0
    return 1

# given concave polygon p1, concave polygon p2
# determine if facets of p1 intersect facets of p2
# a polygon is a tuple of points
cdef bint poly_isect_c(tuple p1,tuple p2):
    cdef bint isegsectfound = 0
    cdef int p1cnt = len(p1)
    cdef int p2cnt = len(p2)
    cdef int px
    cdef int py
    for px in range(p1cnt):
        if isegsectfound:break
        for py in range(p2cnt):
            #if segments_intersect_noncolinear_c(p1[px-1],p1[px],p2[py-1],p2[py]):
            if segs_isect_int_c(p1[px-1],p1[px],p2[py-1],p2[py]):
                isegsectfound = 1
                break
    return isegsectfound

# given start point s, end point e, and n segments, 
# return a colinear set of points equally spaced between s and e
cdef list pline_c(vec3 s,vec3 e,int n):
    cdef list line = [s.cp_c()]
    cdef vec3 tn = s.tov_c(e)
    cdef float l = tn.mag_c()
    cdef int x
    cdef vec3 nxt
    tn.nrm_c()
    tn.uscl_c(float(l)/(n+1))
    for x in range(n):
        print('wtfff',line)
        #nxt = line[-1].cp_c().trn_c(tn)
        nxt = line[x].cp().trn(tn)
        line.append(nxt)
    line.append(e.cp_c())
    return line

# return a ring of points of radius r with n corners
cdef list pring_c(float r,int n):
    cdef vec3 st = vec3(0,0,0).xtrn_c(r)
    cdef vec3 nv
    cdef float alpha = PI*(2.0/n)
    cdef list points = []
    cdef int x
    for x in range(n):
        nv = st.cp_c().zrot_c(x*alpha)
        points.append(nv)
    return points

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

# return a vector normal to the polygon poly
cdef vec3 poly_nrm_c(tuple poly):
    cdef vec3 zero = vec3(0,0,0)
    cdef vec3 c1c2 
    cdef vec3 c2c3
    cdef vec3 pn = zero
    cdef int x = 0
    while pn.isnear_c(zero):
        c1,c2,c3 = poly[x-2],poly[x-1],poly[x]
        c1c2 = c1.tov_c(c2)
        c2c3 = c2.tov_c(c3)
        cang = ang_c(c1c2,c2c3)
        if cang > 0.1 and cang < PI-0.1:
            pn = c1c2.crs_c(c2c3).nrm_c()
        x += 1
    return pn.nrm_c()

# determine the winding number of a point wrt a polygon
cdef int winding_c(vec3 pt,tuple py):
    cdef int wn = 0
    cdef int px
    cdef int pcnt = len(py)
    cdef float read
    for px in range(pcnt):
        read = orient2d_c(pt,py[px-1],py[px])
        if read == 0:return 0
        if py[px-1].y <= pt.y:
            if py[px].y > pt.y:
                if read > 0:wn += 1
        else:
            if py[px].y <= pt.y:
                if read < 0:wn -= 1
    return wn

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

# compute the center of mass for a set of vectors
cpdef vec3 com(ps):
    '''compute the center of mass for a set of vectors'''
    return com_c(ps)

# find the angle between two vectors 
cpdef float ang(vec3 v1,vec3 v2):
    '''find the angle between two vectors'''
    return ang_c(v1,v2)

# find the signed angle between two vectors, given a vector pointing "up"
cpdef float sang(vec3 v1,vec3 v2,vec3 n):
    '''find the signed angle between two vectors, given a vector pointing "up"'''
    return sang_c(v1,v2,n)

# find the angle between the xy projections of two vectors 
cpdef float ang_xy(vec3 v1,vec3 v2):
    '''determine the angle between two vectors projected into the xy plane'''
    return ang_xy_c(v1,v2)

# find the signed angle between the xy projections of two vectors
cpdef float sang_xy(vec3 v1,vec3 v2):
    '''find the signed angle between the xy projections of two vectors'''
    return sang_xy_c(v1,v2)

# find the positive angle between the xy projection of a vector and the x axis
cpdef float xang_xy(vec3 v):
    '''find the angle between the xy projection of a vector and the x axis'''
    return xang_xy_c(v)

# calculate the barycentric coordinates of the point pt for the triangle abc
# assume all points are in the xy plane
cpdef tuple bary_xy(vec3 pt,vec3 a,vec3 b,vec3 c): 
    '''calculate a barycentric representation of a point relative to a triangle in the xy plane'''
    return bary_xy_c(pt,a,b,c)

# determine if a point lies on the interior of a line segment or its endpoints
cpdef bint onseg_xy(vec3 p,vec3 s1,vec3 s2):
    '''determine if a point lies on the interior of a line segment or its endpoints'''
    return onseg_xy_c(p,s1,s2)

# determine if a point p lies on the interior of a line segment s1,s2
cpdef bint inseg_xy(vec3 p,vec3 s1,vec3 s2):
    '''determine if a point lies on the interior of a line segment'''
    return inseg_xy_c(p,s1,s2)

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

# determine if the point pt is inside the triangle abc
# assume all points are in the xy plane
cpdef bint intri_xy(vec3 pt,vec3 a,vec3 b,vec3 c):
    '''determine if a point lies within a triangle in the xy plane'''
    return intri_xy_c(pt,a,b,c)

# determine if the point pt is inside the concave polygon poly or on its boundary
cpdef bint onconcave_xy(vec3 pt,tuple poly):
    '''determine if a point is inside a concave polygon or on its boundary'''
    return onconcave_xy_c(pt,poly)

# determine if the point pt is inside the concave polygon poly
cpdef bint inconcave_xy(vec3 pt,tuple poly):
    '''determine if a point is inside a concave polygon'''
    return inconcave_xy_c(pt,poly)

# given line segment s1, line segment s2
# does s1 overlap the interior or s2?
# a segment is a tuple of two points
cpdef bint segs_isect_perp(vec3 s11,vec3 s12,vec3 s21,vec3 s22):
    '''determine if two line segments intersect or not'''
    return segs_isect_perp_c(s11,s12,s21,s22)

# given line segment s1, line segment s2
# does s1 overlap the interior or s2?
# a segment is a tuple of two points
#cdef bint segments_intersect_interior_c(dpv.vector s11,dpv.vector s12,dpv.vector s21,dpv.vector s22):
cpdef bint segs_isect_int(vec3 s11,vec3 s12,vec3 s21,vec3 s22):
    '''determine if two line segments intersect or not'''
    return segs_isect_int_c(s11,s12,s21,s22)

# given concave polygon p1, concave polygon p2
# does p1 overlap the interior p2?
# a polygon is a tuple of points
# NOTE: this assumes com(p2) is inside p2!
cpdef bint polyinpoly(tuple p1,tuple p2):
    '''determine if one concave polygon overlaps the interior of another'''
    return polyinpoly_c(p1,p2)

# given concave polygon p1, concave polygon p2
# determine if facets of p1 intersect facets of p2
# a polygon is a tuple of points
cpdef bint poly_isect(tuple p1,tuple p2):
    '''determine if the bounds of one concave polygon intersect those of another'''
    return poly_isect_c(p1,p2)

# return a ring of points of radius r with n corners
cpdef list pline(vec3 s,vec3 e,int n):
    '''return a line of n points starting with s, ending with e'''
    return pline_c(s,e,n)

# return a ring of points of radius r with n corners
cpdef list pring(float r,int n):
    '''return a ring of points of radius r with n corners'''
    return pring_c(r,n)

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

# return a vector normal to the polygon poly
cpdef vec3 poly_nrm(tuple poly):
    '''find a vector normal to a polygon'''
    return poly_nrm_c(poly)

# determine the winding number of a point wrt a polygon
cpdef int winding(vec3 pt,tuple py):
    '''calculate the winding number of a point wrt a polygon'''
    return winding_c(pt,py)









