#imports
# cython: profile=True
#cimport cython

from libc.math cimport sqrt

import numpy as np
 



stuff = 'hi'

# 2d classes/functions

cdef class vector2d:

    def __cinit__(self,float x,float y):
        self.x = x
        self.y = y

    def __str__(self):
        strr = 'vector2d:' + str((self.x,self.y))
        return strr

    cpdef vector2d copy(self):
        cdef vector2d new = vector2d(self.x,self.y)
        return new

    cpdef list to_list(self):
        cdef list new = [self.x,self.y,self.z]
        return new

    cpdef translate(self, vector2d tv):
        self.x += tv.x
        self.y += tv.y

    cpdef scale(self, vector2d sv):
        self.x *= sv.x
        self.y *= sv.y

cdef vector2d midpoint2d_c(vector2d v1, vector2d v2):
    cdef float x = (v1.x + v2.x)/2.0
    cdef float y = (v1.y + v2.y)/2.0
    cdef vector2d new = vector2d(x,y)
    return new

cpdef vector2d midpoint2d(vector2d v1, vector2d v2):
    cdef vector2d pt = midpoint2d_c(v1, v2)
    return pt

# 3d classes/functions

cdef class vector:

    def __cinit__(self,float x,float y,float z):
        self.x = x
        self.y = y
        self.z = z

    def __str__(self):
        strr = 'vector:' + str((self.x,self.y,self.z))
        return strr

    cpdef vector copy(self):
        cdef vector new = vector(self.x,self.y,self.z)
        return new

    cpdef list to_list(self):
        cdef list new = [self.x,self.y,self.z]
        return new

    cpdef float magnitude(self):
        cdef float xx = self.x**2
        cdef float yy = self.y**2
        cdef float zz = self.z**2
        cdef float ss = sqrt(xx + yy + zz)
        return ss

    cpdef vector normalize(self):
        cdef float mag = self.magnitude()
        if not mag == 0.0:
            self.x /= mag
            self.y /= mag
            self.z /= mag
        return self

    cpdef vector translate(self, vector tv):
        self.x += tv.x
        self.y += tv.y
        self.z += tv.z
        return self

    cpdef scale(self, vector sv):
        self.x *= sv.x
        self.y *= sv.y
        self.z *= sv.z

cdef float dot_c(vector v1, vector v2):
    return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z

cpdef float dot(vector v1, vector v2):
    return dot_c(v1,v2)

cdef vector cross_c(vector v1, vector v2):
    cdef float cx = v1.y*v2.z-v1.z*v2.y
    cdef float cy = v1.z*v2.x-v1.x*v2.z
    cdef float cz = v1.x*v2.y-v1.y*v2.x
    cdef vector res = vector(cx,cy,cz)
    return res

cpdef vector cross(vector v1, vector v2):
    return cross_c(v1,v2)

cdef vector v1_v2_c(vector v1, vector v2):
    cdef float dx = v2.x - v1.x
    cdef float dy = v2.y - v1.y
    cdef float dz = v2.z - v1.z
    cdef vector new = vector(dx,dy,dz)
    return new

cpdef vector v1_v2(vector v1, vector v2):
    cdef vector pt = v1_v2_c(v1, v2)
    return pt

cdef vector vzip_c(vector v1, vector v2):
    cdef float x = v1.x*v2.x
    cdef float y = v1.y*v2.y
    cdef float z = v1.z*v2.z
    return vector(x,y,z)

cpdef vector vzip(vector v1, vector v2):
    cdef vector pt = vzip_c(v1, v2)
    return pt

cdef vector midpoint_c(vector v1, vector v2):
    cdef float x = (v1.x + v2.x)/2.0
    cdef float y = (v1.y + v2.y)/2.0
    cdef float z = (v1.z + v2.z)/2.0
    cdef vector new = vector(x,y,z)
    return new

cpdef vector midpoint(vector v1, vector v2):
    cdef vector pt = midpoint_c(v1, v2)
    return pt

cdef void translate_coords_c(list coords, vector t):
    cdef int ccnt = len(coords)
    cdef int cdx
    cdef vector coo
    for cdx in range(ccnt):
        coo = <vector>coords[cdx]
        coo.translate(t)

cpdef translate_coords(list coords, vector t):
    translate_coords_c(coords,t)

cdef vector com(list coords):
    cdef int ccnt = len(coords)
    cdef float ccntf = float(ccnt)
    cdef int cdx = 0
    cdef float x = 0.0
    cdef float y = 0.0
    cdef float z = 0.0
    cdef vector coo
    cdef vector new
    for cdx in range(ccnt):
        coo = <vector>coords[cdx]
        x += coo.x
        y += coo.y
        z += coo.z
    new = vector(x/ccntf,y/ccntf,z/ccntf)
    return new

cpdef vector center_of_mass(list coords):
    return com(coords)

cdef float distance_xy_c(vector v1, vector v2):
    cdef float dx = v2.x - v1.x
    cdef float dy = v2.y - v1.y
    cdef float ds = sqrt(dx**2 + dy**2)
    return ds

cpdef float distance_xy(vector v1, vector v2):
    return distance_xy_c(v1, v2)

cdef float distance_c(vector v1, vector v2):
    cdef vector v1v2 = v1_v2(v1,v2)
    cdef float mag = v1v2.magnitude()
    return mag

cpdef float distance(vector v1, vector v2):
    return distance_c(v1, v2)

cdef int find_closest_xy_c(vector one,list bunch,int bcnt,float close_enough):
    cdef float nearest = 100000000.0
    cdef float ds = nearest
    cdef int bdx
    cdef int ndx = 0
    cdef vector which
    for bdx in range(bcnt):
        which = <vector>bunch[bdx]
        ds = distance_xy_c(one,which)
        if ds < nearest:
            nearest = ds
            ndx = bdx
            if ds <= close_enough:
                return ndx
    return ndx

cpdef int find_closest_xy(vector one,list bunch,int bcnt,float close_enough):
    return find_closest_xy_c(one,bunch,bcnt,close_enough)

cdef float angle_from_xaxis_xy_c(vector v):
    cdef float oldz = v.z
    cdef float angle = angle_from_xaxis_c(v)
    v.z = oldz
    return angle

cpdef float angle_from_xaxis_xy(vector v):
    return angle_from_xaxis_xy_c(v)

cdef xhat = vector(1,0,0)
cdef yhat = vector(0,1,0)
cdef zhat = vector(0,0,1)

cdef float angle_from_xaxis_c(vector v):
    v.normalize()
    cdef float xproj = dot_c(v,xhat)
    cdef float yproj = dot_c(v,yhat)
    cdef float ang = np.arccos(xproj)
    if yproj < 0.0: ang = 2.0*np.pi - ang
    return ang

cpdef float angle_from_xaxis(vector v):
    return angle_from_xaxis_c(v)

cpdef bint test_afxa():
    cdef vector v1 = vector(10,0,0)
    cdef vector v2 = vector(10,10,0)
    cdef vector v3 = vector(10,10,10)
    cdef vector v4 = vector(-10,10,0)
    cdef vector v5 = vector(-10,0,0)
    cdef vector v6 = vector(-10,-10,0)
    cdef vector v7 = vector(0,-10,0)
    cdef vector v8 = vector(10,-10,0)
    t1 = str(angle_from_xaxis_c(v1))
    t2 = str(angle_from_xaxis_c(v2))
    t3 = str(angle_from_xaxis_c(v3))
    t4 = str(angle_from_xaxis_c(v4))
    t5 = str(angle_from_xaxis_c(v5))
    t6 = str(angle_from_xaxis_c(v6))
    t7 = str(angle_from_xaxis_c(v7))
    t8 = str(angle_from_xaxis_c(v8))
    print 'ANGLETEST'
    print '\t'.join([t1,t2,t3,t4,t5,t6,t7,t8])

cpdef bint inside(vector pt, list corners):
    poly = [(c.x,c.y) for c in corners]
    x,y = pt.x,pt.y
    n = len(poly)
    inside = False

    p1x,p1y = poly[0]
    for i in range(n+1):
        p2x,p2y = poly[i % n]
        if y > min(p1y,p2y):
            if y <= max(p1y,p2y):
                if x <= max(p1x,p2x):
                    if p1y != p2y:
                        xints = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xints:
                        inside = not inside
        p1x,p1y = p2x,p2y

    return inside



