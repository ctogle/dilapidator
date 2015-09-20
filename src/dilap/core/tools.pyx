#imports
# cython: profile=True
#cimport cython

#import dilap.core.base as db
cimport dilap.core.vector as dpv
cimport dilap.core.quaternion as dpq
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

__doc__ = '''General purpose tool functions...'''

###############################################################################
### c space
###############################################################################

###############################################################################
### 1d inputs/problems
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
cdef float clamp_periodic_c(float v,float f,float c):
    period = c - f
    while v < f:v += period
    while v > c:v -= period
    else:return v                                  

# is a on the interior of (a,b) given an error of d
cdef bint inrange_c(float a,float b,float c):
    cdef float r = near_c(near_c(a,b),c)
    #cdef bint inr = r >= b and r <= c
    cdef bint inr = r > b and r < c
    return inr

# return the shortest angular distance between two angles
cdef float adist_c(float a1,float a2):
    cdef float da = clamp_periodic(a1-a2,0.0,twoPI)
    return da if da < PI else twoPI - da

###############################################################################
###############################################################################

###############################################################################
### sequence inputs/problems
###############################################################################

# return the index of the smallest value in values
cdef int locate_smallest_c(list values):
    cdef int vcnt = len(values)
    cdef int vx
    cdef float v
    cdef float sv = values[0]
    cdef int si = 0
    for vx in range(1,vcnt):
        v = values[vx]
        if v < sv:
            si = vx
            sv = v
    return si

# return the index of the largest value in values
cdef int locate_largest_c(list values):
    cdef int vcnt = len(values)
    cdef int vx
    cdef float v
    cdef float sv = values[0]
    cdef int si = 0
    for vx in range(1,vcnt):
        v = values[vx]
        if v > sv:
            si = vx
            sv = v
    return si

# return the ordered indices for a list of values (ascending)
cdef list order_ascending_c(list values):
    cdef list xrng = [x for x in range(len(values))]
    cdef list od = list(list(zip(*sorted(zip(values,xrng))))[1])
    #od.reverse()
    return od

# is seq1 a cyclic permutation of seq2?
cdef bint cyclic_permutation_c(seq1,seq2):
    cdef int s1cnt = len(seq1)
    cdef int s2cnt = len(seq2)
    cdef int pmdx
    if not s1cnt == s2cnt:return 0
    for pmdx in range(s1cnt):
        perm = seq1[pmdx:] + seq1[:pmdx]
        if perm == seq2:return 1
    return 0

###############################################################################
###############################################################################

###############################################################################
### 2d problems
###############################################################################

# find the angle between the xy projections of two vectors 
cdef float angle_between_xy_c(dpv.vector v1,dpv.vector v2):
    cdef dpv.vector n1 = v1.xy().normalize()
    cdef dpv.vector n2 = v2.xy().normalize()
    cdef float n12dot = n1.x*n2.x + n1.y*n2.y
    cdef float ang = 0.0
    if   isnear_c(n12dot,-1):pass
    elif isnear_c(n12dot, 1):ang = PI
    else:ang = numpy.arccos(n12dot)
    return ang                    

# find the signed angle between the xy projections of two vectors
cdef float signed_angle_between_xy_c(dpv.vector v1,dpv.vector v2):
    cdef dpv.vector n1 = v1.xy().normalize()
    cdef dpv.vector n2 = v2.xy().normalize()
    cdef float n12dot = n1.x*n2.x + n1.y*n2.y
    cdef float ang = 0.0
    cdef float vn = n1.x*n2.y-n1.y*n2.x
    if   isnear_c(n12dot,-1):pass
    elif isnear_c(n12dot, 1):ang = PI
    else:ang = numpy.arccos(n12dot)
    if vn < 0.0:ang *= -1.0
    return ang                    

# find the positive angle between the xy projection of a vector and the x axis
cdef float angle_from_xaxis_xy_c(dpv.vector v):
    cdef dpv.vector nv = v.xy().normalize()
    cdef float ang = 0.0
    if   isnear_c(nv.x, 1):pass
    elif isnear_c(nv.x,-1):ang = PI
    else:ang = numpy.arccos(nv.x)
    if nv.y < 0.0:ang = twoPI - ang
    return ang

# determine if a point lies on the interior of a line segment
cdef bint insegment_xy_c(dpv.vector p,dpv.vector s1,dpv.vector s2):
    if not orient2d_c(p,s1,s2) == 0:return 0
    if p.near(s1) or p.near(s2):return 0
    #if (s1.x != s2.x): # S is not  vertical
    if not isnear_c(s1.x,s2.x): # S is not  vertical
        if (s1.x <= p.x and p.x <= s2.x):return 1
        if (s1.x >= p.x and p.x >= s2.x):return 1
    else: # S is vertical, so test y  coordinate
        if (s1.y <= p.y and p.y <= s2.y):return 1
        if (s1.y >= p.y and p.y >= s2.y):return 1
    return 0

# determine if a point lies on the interior of a line segment or its endpoints
cdef bint onsegment_xy_c(dpv.vector p,dpv.vector s1,dpv.vector s2):
    if not orient2d_c(p,s1,s2) == 0:return 0
    if p.near(s1) or p.near(s2):return 1
    else:return insegment_xy(p,s1,s2)

# return the signed area of the triangle created 
# by the vectors a-c,b-c
# return 0 if a,b,c are colinear
# this assumes theyre in the xy plane!!!
cdef float orient2d_c(dpv.vector a,dpv.vector b,dpv.vector c):
    cdef float m11 = a.x-c.x
    cdef float m12 = a.y-c.y
    cdef float m21 = b.x-c.x
    cdef float m22 = b.y-c.y
    cdef float det = m11*m22-m12*m21
    return near_c(det,0)

###############################################################################
###############################################################################

###############################################################################
### 3d problems
###############################################################################

# find the angle between two vectors 
cdef float angle_between_c(dpv.vector v1,dpv.vector v2):
    cdef dpv.vector n1 = v1.copy().normalize()
    cdef dpv.vector n2 = v2.copy().normalize()
    cdef float n12dot = n1.x*n2.x + n1.y*n2.y + n1.z*n2.z
    cdef float ang = 0.0
    if   isnear_c(n12dot, 1.0):pass
    elif isnear_c(n12dot,-1.0):ang = PI
    else:ang = numpy.arccos(n12dot)
    return ang                    

# find the signed angle between two vectors, given a vector pointing "up"
cdef float signed_angle_between_c(dpv.vector v1,dpv.vector v2,dpv.vector n):
    cdef dpv.vector n1 = v1.copy().normalize()
    cdef dpv.vector n2 = v2.copy().normalize()
    cdef dpv.vector vn = dpv.cross_c(n1,n2)
    cdef float vdot = vn.x*n.x + vn.y*n.y + vn.z*n.z
    cdef float n12dot = n1.x*n2.x + n1.y*n2.y + n1.z*n2.z
    cdef float ang = 0.0
    if   isnear_c(n12dot, 1.0):pass
    elif isnear_c(n12dot,-1.0):ang = PI
    else:ang = numpy.arccos(n12dot)
    if vdot < 0.0:ang *= -1.0
    return ang                    

# compute the distance from pt to the line intersecting e1,e2 along nm
### CONSIDER BROKEN!! ONLY GOOD FOR PERP CASE
cdef float distance_to_line_c(dpv.vector pt,dpv.vector e1,
                              dpv.vector e2,dpv.vector nm):
    cdef float eproj11 = dpv.dot_c(e1,nm)
    cdef float eproj12 = dpv.dot_c(e2,nm)
    cdef float pproj   = dpv.dot_c(pt,nm)
    cdef dpv.vector2d eproj = dpv.vector2d(
        min(eproj11,eproj12),max(eproj11,eproj12))
    return abs(eproj.x - pproj)

# given a pt and a boundary, 
# return the distance to the nearest edge in the boundary
cdef float distance_to_border_c(dpv.vector pt,list border):
    cdef list edgenorms = normals_c(border,dpv.zhat)
    cdef list dists = []
    cdef dpv.vector e1,e2,norm
    cdef int ecnt = len(border)
    cdef int edx
    for edx in range(len(border)):
        e1,e2 = border[edx-1],border[edx]
        norm = edgenorms[edx-1]
        dists.append(distance_to_line(pt,e1,e2,norm))
    cdef float distance = min(dists)
    return distance

# return tangent vectors between points in a list
cdef list tangents_c(list verts):
    cdef dpv.vector v1,v2,tang
    cdef float dx,dy,dv
    cdef list tangs = []
    cdef int vcnt = len(verts)
    cdef int vdx
    for vdx in range(vcnt):
        v1,v2 = verts[vdx-1],verts[vdx]
        tang = dpv.v1_v2_c(v1,v2).normalize()
        tangs.append(tang)
    tangs.append(tangs.pop(0))
    return tangs

# return normal vectors between points in a list, given an "up" direction
cdef list normals_c(list verts,dpv.vector nm):
    cdef dpv.vector norm
    cdef list norms = []
    cdef list tangs = tangents_c(verts)
    cdef int vcnt = len(tangs)
    cdef int vdx
    for vdx in range(vcnt):
        v1,v2 = verts[vdx-1],verts[vdx]
        norm = tangs[vdx].cross(nm)
        norms.append(norm)
    return norms

###############################################################################
###############################################################################

###############################################################################
### in-place transformations on collections of points
###############################################################################

# revolve pt around the edge segment e1,e2 by ang
cdef dpv.vector revolve_about_line_c(dpv.vector pt,dpv.vector e1,dpv.vector e2,float ang):
    cdef dpv.vector etn = dpv.v1_v2_c(e1,e2)
    cdef dpq.quaternion q = dpq.q_from_av_c(ang,etn)
    pt.translate(e1.flip()).rotate(q).translate(e1.flip())
    return pt

# rotate a list of vectors by a quaternion q
cdef list rotate_coords_c(list ps,dpq.quaternion q):
    cdef int pcnt = len(ps)
    cdef int px
    cdef dpv.vector p
    for px in range(pcnt):
        p = <dpv.vector>ps[px]
        p.rotate(q)
    return ps

### UNTESTED
### UNTESTED
### UNTESTED
# rotate a list of segments by a quaternion q
cdef list rotate_segments_c(list segments,dpq.quaternion q):
    cdef int slen = len(segments)
    cdef int ex
    for ex in range(slen):
        s1,s2 = segments[ex]
        s1.rotate(q)
        #s2.rotate(q)
    return segments

# copy a polygon: (extbnd,(holes...)) 
cdef tuple copy_polygon_c(tuple polygon):
    cdef list eb = [x.copy() for x in polygon[0]]
    cdef list ibnds = []
    for ibnd in polygon[1]:
        ibnds.append(tuple([x.copy() for x in ibnd]))
    return (tuple(eb),tuple(ibnds))

# translate a polygon: (extbnd,(holes...)) by vector tv
cdef tuple translate_polygon_c(tuple polygon,dpv.vector tv):
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
        ebnd[ex].translate(tv)
    for ibx in range(islen):
        ibnd = ibnds[ibx]
        ilen = len(ibnd)
        for ix in range(ilen):
            ibnd[ix].translate(tv)
    return polygon

# rotate a polygon: (extbnd,(holes...)) by a quaternion q
cdef tuple rotate_polygon_c(tuple polygon,dpq.quaternion q):
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
        ebnd[ex].rotate(q)
    for ibx in range(islen):
        ibnd = ibnds[ibx]
        ilen = len(ibnd)
        for ix in range(ilen):
            ibnd[ix].rotate(q)
    return polygon

# rotate a polygon: (extbnd,(holes...)) by float a around xhat
cdef tuple rotate_x_polygon_c(tuple polygon,float a):
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
        ebnd[ex].rotate_x(a)
    for ibx in range(islen):
        ibnd = ibnds[ibx]
        ilen = len(ibnd)
        for ix in range(ilen):
            ibnd[ix].rotate_x(a)
    return polygon

# translate a polygon: (extbnd,(holes...)) by float a around zhat
cdef tuple rotate_z_polygon_c(tuple polygon,float a):
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
        ebnd[ex].rotate_z(a)
    for ibx in range(islen):
        ibnd = ibnds[ibx]
        ilen = len(ibnd)
        for ix in range(ilen):
            ibnd[ix].rotate_z(a)
    return polygon
### UNTESTED
### UNTESTED
### UNTESTED

###############################################################################
###############################################################################


###############################################################################
### DISORGANIZED!!
###############################################################################

### UNTESTED
### UNTESTED
### UNTESTED

# given a vector, return a quaternion that rotates it to coincide with zhat
cdef dpq.quaternion q_to_xy_c(dpv.vector v):
    cdef dpq.quaternion prot
    if v.near(dpv.nzhat):prot = dpq.q_from_av_c(PI,dpv.xhat)
    elif not v.near(dpv.zhat):prot = dpq.q_from_uu_c(v,dpv.zhat)
    else:prot = dpq.zero_c()
    return prot

# return a vector normal to the plane containing c1,c2,c3
# returns 0 if c1,c2,c3 are colinear
cdef dpv.vector normal_c(dpv.vector c1,dpv.vector c2,dpv.vector c3):
    cdef dpv.vector c1c2 = dpv.v1_v2_c(c1,c2).normalize()
    cdef dpv.vector c2c3 = dpv.v1_v2_c(c2,c3).normalize()
    cdef dpv.vector cn = c1c2.cross(c2c3).normalize()
    return cn

# return a vector normal to the polygon poly
cdef dpv.vector polygon_normal_c(tuple poly):
    cdef dpv.vector zero = dpv.zero_c()
    cdef dpv.vector c1c2 
    cdef dpv.vector c2c3
    cdef dpv.vector pn = zero
    cdef int x = 0
    while pn.near(zero):
        c1,c2,c3 = poly[x-2],poly[x-1],poly[x]
        c1c2 = dpv.v1_v2_c(c1,c2)
        c2c3 = dpv.v1_v2_c(c2,c3)
        cang = angle_between_c(c1c2,c2c3)
        if cang > 0.1 and cang < PI-0.1:
            pn = c1c2.cross(c2c3).normalize()
        x += 1
    return pn.normalize()

# return a vector tanget to the plane containing c1,c2,c3
cdef dpv.vector tangent_c(dpv.vector c1,dpv.vector c2,dpv.vector c3):
    cdef dpv.vector tn = dpv.v1_v2_c(c1,c2).normalize()
    return tn

# given start point s, end point e, and n segments, 
# return a colinear set of points equally spaced between s and e
cdef list point_line_c(dpv.vector s,dpv.vector e,int n):
    cdef list line = [s.copy()]
    cdef dpv.vector tn = dpv.v1_v2_c(s,e)
    cdef float l = tn.magnitude()
    cdef int x
    cdef dpv.vector new
    tn.normalize()
    tn.scale_u(float(l)/n)
    for x in range(n):
        new = line[-1].copy().translate(tn)
        line.append(new)
    return line

# return a ring of points of radius r with n corners
cdef list point_ring_c(float r,int n):
    cdef dpv.vector st = dpv.zero_c().translate_x(r)
    cdef dpv.vector nv
    cdef float alpha = PI*(2.0/n)
    cdef list points = []
    cdef int x
    for x in range(n):
        nv = st.copy().rotate_z(x*alpha)
        points.append(nv)
    return points

# return a square of length,width l,w at position p and zrot a
cdef list square_c(float l,float w,p = None,a = None):
    cdef float l2 = l/2.0
    cdef float w2 = w/2.0
    cdef list cs = [
        dpv.vector(-l2,-w2,0),dpv.vector( l2,-w2,0), 
        dpv.vector( l2, w2,0),dpv.vector(-l2, w2,0)]
    if not a is None:dpv.rotate_z_coords_c(cs,a)
    if not p is None:dpv.translate_coords_c(cs,p)
    return cs

# consider xy project of polygon corners and xy projection of pt
# return 0 is pt is outside of polygon in this projection
cdef bint inside_c(dpv.vector pt,list corners):
    cdef float x = pt.x
    cdef float y = pt.y
    cdef int n = len(corners)
    cdef bint ins = 0
    cdef dpv.vector p1 = corners[0]
    cdef dpv.vector p2
    cdef float p1x
    cdef float p1y
    cdef float p2x
    cdef float p2y
    cdef int i
    p1x,p1y = p1.x,p1.y
    for i in range(n+1):
        p2 = corners[i % n]
        p2x,p2y = p2.x,p2.y
        if y > min(p1y,p2y):
            if y <= max(p1y,p2y):
                if x <= max(p1x,p2x):
                    if p1y != p2y:xints = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xints:ins = 1 - ins
        p1x,p1y = p2x,p2y
    return ins

# return 1 if pt is inside a sphere of radius r centered at c
# otherwise return 0
cdef bint inside_circle_c(dpv.vector pt,dpv.vector c,float r):
    cdef bint ins = not dpv.distance_c(pt,c) > r
    return ins

# given 3 points and a plane, determine the center and radius
# of the circumcenter found in the plane plane by projecting p1,p2,p3
# in the plane
cdef tuple circumscribe_tri_c(dpv.vector p1,dpv.vector p2,dpv.vector p3):
    cdef dpv.vector cp1 = p1.xy()
    cdef dpv.vector cp2 = p2.xy()
    cdef dpv.vector cp3 = p3.xy()
    cdef dpv.vector e1 = dpv.v1_v2_xy_c(cp3,cp1)
    cdef dpv.vector e2 = dpv.v1_v2_xy_c(cp3,cp2)
    cdef float th = angle_between_xy_c(e1,e2)
    cdef float cr = dpv.distance_xy_c(cp1,cp2)/(2*numpy.sin(th))
    cdef dpv.vector cp = e2.copy().scale_u(
        e1.magnitude2())-e1.copy().scale_u(e2.magnitude2())
    cdef dpv.vector fp = p3+cp.cross(e1.cross(e2)).scale_u(
                    1.0/(2.0*(e1.cross(e2).magnitude2())))
    return fp,cr

# given line segment s1, line segment s2
# does s1 overlap the interior or s2?
# a segment is a tuple of two points
cdef bint segments_intersect_noncolinear_c(dpv.vector s11,dpv.vector s12,dpv.vector s21,dpv.vector s22):
    p,q = s11,s21
    r = dpv.v1_v2_c(s11,s12)
    s = dpv.v1_v2_c(s21,s22)
    qmp = q-p
    rcs = r.cross(s)
    rcsmag = rcs.magnitude()
    qmpcr = qmp.cross(r)
    qmpcrmag = qmpcr.magnitude()
    rmag2 = r.magnitude2()
    if isnear_c(rcsmag,0) and isnear_c(qmpcrmag,0):return 0
    elif isnear_c(rcsmag,0) and not isnear_c(qmpcrmag,0):return 0
    elif not isnear_c(rcs.z,0):
        u = near_c(near_c(       qmpcr.z/rcs.z,0),1)
        t = near_c(near_c(qmp.cross(s).z/rcs.z,0),1)
        if (u == 0 or u == 1) and (t == 0 or t == 1):
            #if include_endpoints:return 1
            #else:return 0
            return 0
        if not inrange_c(u,0,1) or not inrange_c(t,0,1):return 0
        else:return 1

# given line segment s1, line segment s2
# does s1 overlap the interior or s2?
# a segment is a tuple of two points
cdef bint segments_intersect_interior_c(dpv.vector s11,dpv.vector s12,dpv.vector s21,dpv.vector s22):
    p,q = s11,s21
    r = dpv.v1_v2_c(s11,s12)
    s = dpv.v1_v2_c(s21,s22)
    qmp = q-p
    rcs = r.cross(s)
    rcsmag = rcs.magnitude()
    qmpcr = qmp.cross(r)
    qmpcrmag = qmpcr.magnitude()
    rmag2 = r.magnitude2()
    if isnear_c(rcsmag,0) and isnear_c(qmpcrmag,0):
        t0 = near_c(near_c(dpv.dot_c(qmp,r)/rmag2,0),1)
        t1 = near_c(near_c(t0 + dpv.dot_c(s,r)/rmag2,0),1)
        if near_c(dpv.dot_c(s,r),0) < 0:t0,t1 = t1,t0
        if inrange_c(t0,0,1) or inrange_c(t1,0,1):return 1
        elif inrange_c(0,t0,t1) and inrange_c(1,t0,t1):return 1
        elif (t0 == 0 and t1 == 1) or (t0 == 1 and t1 == 0):return 1
        else:return 0
    elif isnear_c(rcsmag,0) and not isnear_c(qmpcrmag,0):return 0
    elif not isnear_c(rcs.z,0):
        u = near_c(near_c(       qmpcr.z/rcs.z,0),1)
        t = near_c(near_c(qmp.cross(s).z/rcs.z,0),1)
        if (u == 0 or u == 1) and (t == 0 or t == 1):
            #if include_endpoints:return 1
            #else:return 0
            return 0
        if not inrange_c(u,0,1) or not inrange_c(t,0,1):return 0
        else:return 1

# calculate the barycentric coordinates of the point pt for the triangle abc
# assume all points are in the xy plane 
cdef tuple barycentric_xy_c(dpv.vector pt,dpv.vector a,dpv.vector b,dpv.vector c): 
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
    #if denom == 0:
    #    print('holla',a.__str__(),b.__str__(),c.__str__(),pt.__str__())
    #    raise ValueError
    cdef float invdenom = 1.0 / denom
    cdef float u = (dot11 * dot02 - dot01 * dot12) * invdenom
    cdef float v = (dot00 * dot12 - dot01 * dot02) * invdenom
    return u,v

# determine if the point pt is inside the triangle abc
# assume all points are in the xy plane 
cdef bint intriangle_xy_c(dpv.vector pt,dpv.vector a,dpv.vector b,dpv.vector c):
    cdef float u,v
    u,v = barycentric_xy_c(pt,a,b,c)
    if u > 0 or abs(u) < epsilon:
        if v > 0 or abs(v) < epsilon:
            if 1-u-v > 0 or abs(1-u-v) < epsilon:
                return 1
    return 0

# calculate the barycentric coordinates of the point pt for the triangle abc
# NOTE: DOES THE POINT NEED TO BE PROJECTED INTO THE PLANE OF THE TRIANGLE????
cdef dpv.vector2d barycentric_c(dpv.vector pt,dpv.vector a,dpv.vector b,dpv.vector c): 
    cdef dpv.vector v0 = dpv.v1_v2_c(a,c)
    cdef dpv.vector v1 = dpv.v1_v2_c(a,b)
    cdef dpv.vector v2 = dpv.v1_v2_c(a,pt)
    cdef float dot00 = v0.x*v0.x + v0.y*v0.y + v0.z*v0.z  
    cdef float dot01 = v0.x*v1.x + v0.y*v1.y + v0.z*v1.z  
    cdef float dot02 = v0.x*v2.x + v0.y*v2.y + v0.z*v2.z  
    cdef float dot11 = v1.x*v1.x + v1.y*v1.y + v1.z*v1.z  
    cdef float dot12 = v1.x*v2.x + v1.y*v2.y + v1.z*v2.z  
    cdef float invdenom = 1.0 / (dot00 * dot11 - dot01 * dot01)
    cdef float u = (dot11 * dot02 - dot01 * dot12) * invdenom
    cdef float v = (dot00 * dot12 - dot01 * dot02) * invdenom
    cdef dpv.vector2d bary = dpv.vector2d(u,v)
    return bary

# determine if the point pt is inside the triangle abc
cdef bint intriangle_c(dpv.vector pt,dpv.vector a,dpv.vector b,dpv.vector c):
    cdef dpv.vector2d bary = barycentric_c(pt,a,b,c)
    cdef bint i1 = near_c(bary.x,0) >= 0
    cdef bint i2 = near_c(bary.y,0) >= 0
    cdef bint i3 = near_c(1-bary.x-bary.y,0) >= 0
    cdef bint ins = i1 and i2 and i3
    return ins

##### DO WORK
# determine if the point pt is inside the convex polygon poly
cdef bint inconvex_c(dpv.vector pt,tuple poly):
    raise NotImplemented
    # why the hell doenst this work...
    cdef int ecnt = len(poly)
    cdef int edx
    for edx in range(ecnt):
        p0,p1 = poly[edx-1],poly[edx]
        sarea = orient2d_c(p1,p0,pt)
        if near_c(sarea,0) > 0:return 1
    return 0
#####

# determine if the point pt is inside the concave polygon poly or on its boundary
cdef bint onconcave_xy_c(dpv.vector pt,tuple py):
    cdef int wn = 0
    cdef int px
    cdef int pcnt = len(py)
    cdef float read
    for px in range(pcnt):
        read = orient2d_c(pt,py[px-1],py[px])
        if read == 0:
            if onsegment_xy_c(pt,py[px-1],py[px]):return 1
        if py[px-1].y <= pt.y:
            if py[px].y > pt.y:
                if read > 0:wn += 1
        else:
            if py[px].y <= pt.y:
                if read < 0:wn -= 1
    return not wn == 0

# determine if the point pt is inside the concave polygon poly
cdef bint inconcave_xy_c(dpv.vector pt,tuple poly):
    return not winding_c(pt,poly) == 0

# determine the winding number of a point wrt a polygon
cdef int winding_c(dpv.vector pt,tuple py):
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

# given concave polygon p1, concave polygon p2
# does p1 contain the entire interior of p2?
# a polygon is a tuple of points
# NOTE: this assumes com(p2) is inside p2!
cdef bint concaves_contains_c(tuple p1,tuple p2):
    cdef int p2cnt = len(p2)
    cdef int px
    cdef dpv.vector i2 = dpv.com(list(p2))
    if not inconcave_xy_c(i2,p1):return 0
    for px in range(p2cnt):
        if not onconcave_xy_c(p2[px],p1):
            return 0
    return 1

# given concave polygon p1, concave polygon p2
# determine if facets of p1 intersect facets of p2
# a polygon is a tuple of points
cdef bint concaves_intersect_c(tuple p1,tuple p2):
    cdef bint isegsectfound = 0
    cdef int p1cnt = len(p1)
    cdef int p2cnt = len(p2)
    cdef int px
    cdef int py
    for px in range(p1cnt):
        if isegsectfound:break
        for py in range(p2cnt):
            #if segments_intersect_noncolinear_c(p1[px-1],p1[px],p2[py-1],p2[py]):
            if segments_intersect_interior_c(p1[px-1],p1[px],p2[py-1],p2[py]):
                isegsectfound = 1
                break
    return isegsectfound

# return the point with the highest y value
cdef dpv.vector find_y_apex_c(list pts):
    cdef dpv.vector hi = pts[0]
    cdef int pcnt = len(pts)
    cdef int px
    for px in range(1,pcnt):
        pt = pts[px]
        if pt.y > hi.y:hi = pt
    return hi

# return the point with the highest x value
cdef dpv.vector find_x_apex_c(list pts):
    cdef dpv.vector hi = pts[0]
    cdef int pcnt = len(pts)
    cdef int px
    for px in range(1,pcnt):
        pt = pts[px]
        if pt.x > hi.x:hi = pt
    return hi

# perform a sweep search of a set a points with respect to another point
cdef dpv.vector sweep_search_c(list pts,dpv.vector center,tangent = None):
    cdef dpv.vector offset = center.copy().flip()
    cdef float tangent_rot
    cdef float tpang
    dpv.translate_coords_c(pts,offset)
    if not tangent is None:
        tangent_rot = angle_from_xaxis_xy_c(tangent)
        dpv.rotate_z_coords_c(pts,-tangent_rot)
    cdef dpv.vector which = center
    cdef float pang = twoPI
    cdef int pcnt = len(pts)
    cdef int adx
    cdef dpv.vector pt
    for adx in range(pcnt):
        pt = pts[adx]
        #if pt is center:continue
        if pt.near(center):continue
        tpang = angle_from_xaxis_xy_c(pt)
        if tpang < pang:
            pang = tpang
            which = pt
    if not tangent is None:dpv.rotate_z_coords_c(pts,tangent_rot)
    dpv.translate_coords_c(pts,offset.flip())
    return which

# find a convex hull in the xy-plane for a set of points
# NOTE:pts should not be a colinear set!
cdef list pts_to_convex_xy_c(list pts):
    cdef dpv.vector new = find_x_apex_c(pts)
    cdef list shape = []
    tang = None
    while not new in shape:
        shape.append(new)
        if len(shape) > 1:
            tang = dpv.v1_v2_c(shape[-2],shape[-1]).rotate_z(-0.001)
        new = sweep_search_c(pts,new,tang)
    return shape

# move the points of a convex (or star) polygon along local normals
cdef list inflate_c(list convex,float radius):
    cdef list enorms = normals_c(convex,dpv.zhat)
    cdef int ccnt = len(convex)
    cdef int cdx
    cdef dpv.vector lead
    cdef dpv.vector rear
    cdef dpv.vector norm
    for cdx in range(ccnt):
        lead = enorms[cdx]
        rear = enorms[cdx-1]
        norm = dpv.midpoint_c(lead,rear).normalize()
        convex[cdx].translate(norm.scale_u(radius))
    return convex

# add index offset to a list of faces
cdef list offset_faces_c(list faces,int offset):
    cdef int fcnt = len(faces)
    cdef int fdx
    cdef list fa
    for fdx in range(fcnt):
        fa = faces[fdx]
        #tfcnt = len(fa)
        #for tfdx in range(tfcnt):
        #    fa[tfdx] += offset
        fa[0] += offset
        fa[1] += offset
        fa[2] += offset
    return faces

# return the signed volume of the parallelpiped created
# by the vectors a-d,b-d,c-d
# return 0 if a,b,c,d are coplanar
cdef float orient3d_c(dpv.vector a,dpv.vector b,dpv.vector c,dpv.vector d):
    cdef float m11 = a.x-d.x
    cdef float m12 = a.y-d.y
    cdef float m13 = a.z-d.z
    cdef float m21 = b.x-d.x
    cdef float m22 = b.y-d.y
    cdef float m23 = b.z-d.z
    cdef float m31 = c.x-d.x
    cdef float m32 = c.y-d.y
    cdef float m33 = c.z-d.z
    cdef float det = m11*(m22*m33-m23*m32)-m12*(m21*m33-m23*m31)+m13*(m21*m32-m22*m31)
    return near_c(det,0)

# determine if d is inside the circumcircle of the triangle a,b,c
cdef float incircle_c(dpv.vector a,dpv.vector b,dpv.vector c,dpv.vector d):
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
    cdef float incirc = near_c(det1 - det2 + det3,0)
    return incirc

### THIS NEEDS MORE WORK
### THIS NEEDS MORE WORK
### THIS NEEDS MORE WORK
# let a,b,c,d be such that orient3d(a,b,c,d) is nonnegative
# return > 0 if e is inside sphere passing through a,b,c,d
# return < 0 if e is outside sphere passing through a,b,c,d
# return 0 if e is on the surface of the sphere passing through a,b,c,d
# return 0 if all five points are coplanar
cdef float insphere_c(dpv.vector a,dpv.vector b,dpv.vector c,dpv.vector d,dpv.vector e):
    m11 = a.x-e.x;m12 = a.y-e.y;m13 = a.z-e.z
    m14 = m11*m11 + m12*m12 + m13*m13
    m21 = b.x-e.x;m22 = b.y-e.y;m23 = b.z-e.z
    m24 = m21*m21 + m22*m22 + m23*m23
    m31 = c.x-e.x;m32 = c.y-e.y;m33 = c.z-e.z
    m34 = m31*m31 + m32*m32 + m33*m33
    m41 = d.x-e.x;m42 = d.y-e.y;m43 = d.z-e.z
    m44 = m41*m41 + m42*m42 + m43*m43
    det1 = m11*(m22*(m33*m44-m34*m43)-m23*(m32*m44-m34*m42)+m24*(m32*m43-m33*m42))
    det2 = m12*(m21*(m33*m44-m34*m43)-m23*(m31*m44-m34*m41)+m24*(m31*m43-m33*m41))
    det3 = m13*(m21*(m32*m44-m34*m42)-m22*(m31*m44-m34*m41)+m24*(m31*m42-m32*m41))
    det4 = m14*(m21*(m32*m43-m33*m42)-m22*(m31*m43-m33*m41)+m23*(m31*m42-m32*m41))
    insphr = near_c(det1 - det2 + det3 - det4,0)
    return insphr

###############################################################################
### python space
###############################################################################

# if a is within c of b, return b
# else return a
cpdef float near(float a,float b):
    '''effectively round a to b if within a neighborhood c'''
    return near_c(a,b)

# if a is within c of b, return True
# else return False
cpdef bint isnear(float a,float b):
    '''determine if a is within a neighborhood c of b'''
    return isnear_c(a,b)

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
cpdef float clamp_periodic(float v,float f,float c):
    '''clamp a float between two floats using periodic boundaries'''
    return clamp_periodic_c(v,f,c)

# is a on the interior of (a,b) given an error of d
cpdef bint inrange(float a,float b,float c):
    '''determine if a value is on an open interval'''
    return inrange_c(a,b,c)

# return the shortest angular distance between two angles
cpdef float adist(float a1,float a2):
    '''find the angular distance between two angles'''
    return adist_c(a1,a2)

# return the index of the smallest value in values
cpdef int locate_smallest(list values):
    '''locate the smallest value in a list of values'''
    return locate_smallest_c(values)

# return the index of the largest value in values
cpdef int locate_largest(list values):
    '''locate the largest value in a list of values'''
    return locate_largest_c(values)

# return the ordered indices for a list of values (ascending)
cpdef list order_ascending(list values):
    '''determine the ascending ordering of a list of values'''
    return order_ascending_c(values)

# is seq1 a cyclic permutation of seq2?
cpdef bint cyclic_permutation(seq1,seq2):
    '''determine if one sequence is a cyclic permutation of another'''
    return cyclic_permutation_c(seq1,seq2)

# find the angle between the xy projections of two vectors 
cpdef float angle_between_xy(dpv.vector v1,dpv.vector v2):
    '''determine the angle between two vectors projected into the xy plane'''
    return angle_between_xy_c(v1,v2)

# find the signed angle between the xy projections of two vectors
cpdef float signed_angle_between_xy(dpv.vector v1,dpv.vector v2):
    '''find the signed angle between the xy projections of two vectors'''
    return signed_angle_between_xy_c(v1,v2)

# find the positive angle between the xy projection of a vector and the x axis
cpdef float angle_from_xaxis_xy(dpv.vector v):
    '''find the angle between the xy projection of a vector and the x axis'''
    return angle_from_xaxis_xy_c(v)

# determine if a point lies on the interior of a line segment or its endpoints
cpdef bint onsegment_xy(dpv.vector p,dpv.vector s1,dpv.vector s2):
    '''determine if a point lies on the interior of a line segment or its endpoints'''
    return onsegment_xy_c(p,s1,s2)

# determine if a point p lies on the interior of a line segment s1,s2
cpdef bint insegment_xy(dpv.vector p,dpv.vector s1,dpv.vector s2):
    '''determine if a point lies on the interior of a line segment'''
    return insegment_xy_c(p,s1,s2)

# return the signed area of the triangle created 
# by the vectors a-c,b-c
# return 0 if a,b,c are colinear
# this assumes theyre in the xy plane!!!
cpdef float orient2d(dpv.vector a,dpv.vector b,dpv.vector c):
    '''determine whether a,b,c form a positively oriented triangle'''
    return orient2d_c(a,b,c)

# find the angle between two vectors 
cpdef float angle_between(dpv.vector v1,dpv.vector v2):
    '''find the angle between two vectors'''
    return angle_between_c(v1,v2)

# find the signed angle between two vectors, given a vector pointing "up"
cpdef float signed_angle_between(dpv.vector v1,dpv.vector v2,dpv.vector n):
    '''find the signed angle between two vectors, given a vector pointing "up"'''
    return signed_angle_between_c(v1,v2,n)

# compute the distance from pt to the line intersecting e1,e2 along nm
cpdef float distance_to_line(dpv.vector pt,dpv.vector e1,dpv.vector e2,dpv.vector nm):
    '''compute the distance from a point to an edge segment along a unit vector'''
    return distance_to_line_c(pt,e1,e2,nm)

# given a pt and a boundary, 
# return the distance to the nearest edge in the boundary
cpdef float distance_to_border(dpv.vector pt,list border):
    ''' find the distance from a point to a convex boundary containing it'''
    return distance_to_border_c(pt,border)

# return tangent vectors between points in a list
cpdef list tangents(list verts):
    '''get tangent vectors between points in a list'''
    return tangents_c(verts)

# return normal vectors between points in a list, given an "up" direction
cpdef list normals(list verts,dpv.vector nm):
    '''get normal vectors between points in a list, given an "up" direction'''
    return normals_c(verts,nm)

# revolve pt around the edge segment e1,e2 by ang
cpdef dpv.vector revolve_about_line(dpv.vector pt,dpv.vector e1,dpv.vector e2,float ang):
    '''compute the distance from a point to an edge segment along a unit vector'''
    return revolve_about_line_c(pt,e1,e2,ang)

# rotate a list of vectors by a quaternion q
cpdef list rotate_coords(list coords,dpq.quaternion q):
    '''rotate a list of vectors by a quaternion'''
    return rotate_coords_c(coords,q)

# given a vector, return a quaternion that rotates it to coincide with zhat
cpdef dpq.quaternion q_to_xy(dpv.vector v):
    '''provide a quaternion which rotates a vector onto zhat'''
    return q_to_xy_c(v)

# return a vector normal to the plane containing c1,c2,c3
# returns 0 if c1,c2,c3 are colinear
cpdef dpv.vector normal(dpv.vector c1,dpv.vector c2,dpv.vector c3):
    '''find a vector normal to the plane containing c1,c2,c3'''
    return normal_c(c1,c2,c3)

# return a vector normal to the polygon poly
cpdef dpv.vector polygon_normal(tuple poly):
    '''find a vector normal to a polygon'''
    return polygon_normal_c(poly)

# return a vector tanget to the plane containing c1,c2,c3
cpdef dpv.vector tangent(dpv.vector c1,dpv.vector c2,dpv.vector c3):
    '''find a vector tangent to the plane containing c1,c2,c3'''
    return tangent_c(c1,c2,c3)

# return a ring of points of radius r with n corners
cpdef list point_line(dpv.vector s,dpv.vector e,int n):
    '''return a line of n points starting with s, ending with e'''
    return point_line_c(s,e,n)

# return a ring of points of radius r with n corners
cpdef list point_ring(float r,int n):
    '''return a ring of points of radius r with n corners'''
    return point_ring_c(r,n)

# return a square of length,width l,w at position p and zrot a
cpdef list square(float l,float w,p = None,a = None):
    '''return a rectangle centered and oriented about a point in the xy plane'''
    return square_c(l,w,p,a)

# consider xy project of polygon corners and xy projection of pt
# return 0 is pt is outside of polygon in this projection
cpdef bint inside(dpv.vector pt,list corners):
    '''determine whether a polygon contains a point'''
    return inside_c(pt,corners)

# return 1 if pt is inside a sphere of radius r centered at c
# otherwise return 0
cpdef bint inside_circle(dpv.vector pt,dpv.vector c,float r):
    '''determine whether a point is contained by a circle'''
    return inside_circle_c(pt,c,r)

# given 3 points and a plane, determine the center and radius
# of the circumcenter found in the plane plane by projecting p1,p2,p3
# in the plane
cpdef tuple circumscribe_tri(dpv.vector p1,dpv.vector p2,dpv.vector p3):
    '''return the circumcenter and circumradius of a triangle'''
    return circumscribe_tri_c(p1,p2,p3)
                                            
# given line segment s1, line segment s2
# does s1 overlap the interior or s2?
# a segment is a tuple of two points
cpdef bint segments_intersect_noncolinear(dpv.vector s11,dpv.vector s12,dpv.vector s21,dpv.vector s22):
    '''determine if two line segments intersect or not'''
    return segments_intersect_noncolinear_c(s11,s12,s21,s22)

# calculate the barycentric coordinates of the point pt for the triangle abc
# assume all points are in the xy plane
cpdef dpv.vector2d barycentric_xy(dpv.vector pt,dpv.vector a,dpv.vector b,dpv.vector c): 
    '''calculate a barycentric representation of a point relative to a triangle in the xy plane'''
    return barycentric_xy_c(pt,a,b,c)

# determine if the point pt is inside the triangle abc
# assume all points are in the xy plane
cpdef bint intriangle_xy(dpv.vector pt,dpv.vector a,dpv.vector b,dpv.vector c):
    '''determine if a point lies within a triangle in the xy plane'''
    return intriangle_xy_c(pt,a,b,c)

# calculate the barycentric coordinates of the point pt for the triangle abc
cpdef dpv.vector2d barycentric(dpv.vector pt,dpv.vector a,dpv.vector b,dpv.vector c): 
    '''calculate a barycentric representation of a point relative to a triangle'''
    return barycentric_c(pt,a,b,c)

# determine if the point pt is inside the triangle abc
cpdef bint intriangle(dpv.vector pt,dpv.vector a,dpv.vector b,dpv.vector c):
    '''determine if a point lies within a triangle'''
    return intriangle_c(pt,a,b,c)
  
# determine if the point pt is inside the convex polygon poly
cpdef bint inconvex(dpv.vector pt,tuple poly):
    '''determine if a pt lies within a convex polygon'''
    return inconvex_c(pt,poly)

# determine if the point pt is inside the concave polygon poly or on its boundary
cpdef bint onconcave_xy(dpv.vector pt,tuple poly):
    '''determine if a point is inside a concave polygon or on its boundary'''
    return onconcave_xy_c(pt,poly)

# determine if the point pt is inside the concave polygon poly
cpdef bint inconcave_xy(dpv.vector pt,tuple poly):
    '''determine if a point is inside a concave polygon'''
    return inconcave_xy_c(pt,poly)

# determine the winding number of a point wrt a polygon
cpdef int winding(dpv.vector pt,tuple py):
    '''calculate the winding number of a point wrt a polygon'''
    return winding_c(pt,py)

# given concave polygon p1, concave polygon p2
# does p1 overlap the interior p2?
# a polygon is a tuple of points
# NOTE: this assumes com(p2) is inside p2!
cpdef bint concaves_contains(tuple p1,tuple p2):
    '''determine if one polygon overlaps the interior of another'''
    return concaves_contains_c(p1,p2)

# given concave polygon p1, concave polygon p2
# determine if facets of p1 intersect facets of p2
# a polygon is a tuple of points
cpdef bint concaves_intersect(tuple p1,tuple p2):
    '''determine if the bounds of one polygon intersect those of another'''
    return concaves_intersect_c(p1,p2)

# return the point with the highest y value
cpdef dpv.vector find_y_apex(list pts):
    '''find the point which maximizes the y coordinate within a set of points'''
    return find_y_apex_c(pts)

# return the point with the highest x value
cpdef dpv.vector find_x_apex(list pts):
    '''find the point which maximizes the x coordinate within a set of points'''
    return find_x_apex_c(pts)

# perform a sweep search of a set a points with respect to another point
cpdef dpv.vector sweep_search(list pts,dpv.vector center,tangent = None):
    '''perform a sweep search on a set of points with respect to a point'''
    return sweep_search_c(pts,center,tangent)

# find a convex hull in the xy-plane for a set of points
# NOTE:pts should not be a colinear set!
cpdef list pts_to_convex_xy(list pts):
    '''generate a convex hull in the xy-plane for a set of points'''
    return pts_to_convex_xy_c(pts)

# move the points of a convex (or star) polygon along local normals
cpdef list inflate(list convex,float radius):
    '''inflate a convex polygon along local normals by a radius'''
    return inflate_c(convex,radius)

# add index offset to a list of faces
cpdef list offset_faces(list faces,int offset):
    '''apply an index offset to a list of triangles'''
    return offset_faces_c(faces,offset)

# rotate a list of segments by a quaternion q
cpdef list rotate_segments(list segments,dpq.quaternion q):
    '''rotate a group of segments by a quaternion'''
    return rotate_segments_c(segments,q)

# copy a polygon: (extbnd,(holes...)) 
cpdef tuple copy_polygon(tuple polygon):
    '''create a copy of a concave polygon with holes'''
    return copy_polygon_c(polygon)

# translate a polygon: (extbnd,(holes...)) by vector tv
cpdef tuple translate_polygon(tuple polygon,dpv.vector tv):
    '''translate a concave polygon and its holes by a vector'''
    return translate_polygon_c(polygon,tv)

# rotate a polygon: (extbnd,(holes...)) by a quaternion q
cpdef tuple rotate_polygon(tuple polygon,dpq.quaternion q):
    '''rotate a concave polygon and its holes by a quaternion'''
    return rotate_polygon_c(polygon,q)

# rotate a polygon: (extbnd,(holes...)) by float a about zhat
cpdef tuple rotate_z_polygon(tuple polygon,float a):
    '''rotate a concave polygon and its holes by a float about zhat'''
    return rotate_z_polygon_c(polygon,a)

# rotate a polygon: (extbnd,(holes...)) by float a about xhat
cpdef tuple rotate_x_polygon(tuple polygon,float a):
    '''rotate a concave polygon and its holes by a float about xhat'''
    return rotate_x_polygon_c(polygon,a)

# return the signed volume of the parallelpiped created
# by the vectors a-d,b-d,c-d
# return 0 if a,b,c,d are coplanar
cpdef float orient3d(dpv.vector a,dpv.vector b,dpv.vector c,dpv.vector d):
    '''determine whether a,b,c,d form a positively oriented tetrahedron'''
    return orient3d_c(a,b,c,d)

# determine if d is inside the circumcircle of the triangle a,b,c
cpdef float incircle(dpv.vector a,dpv.vector b,dpv.vector c,dpv.vector d):
    '''determine if d is inside the circumcircle of the triangle a,b,c'''
    return incircle_c(a,b,c,d)

# let a,b,c,d be such that orient3d(a,b,c,d) is nonnegative
# return > 0 if e is inside sphere passing through a,b,c,d
# return < 0 if e is outside sphere passing through a,b,c,d
# return 0 if e is on the surface of the sphere passing through a,b,c,d
# return 0 if all five points are coplanar
cpdef float insphere(dpv.vector a,dpv.vector b,dpv.vector c,dpv.vector d,dpv.vector e):
    '''determine if e is inside the circumsphere of the tetrahedron a,b,c,d'''
    return insphere_c(a,b,c,d,e)









