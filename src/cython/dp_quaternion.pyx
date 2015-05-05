#imports
# cython: profile=True
#cimport cython
cimport dp_vector as dpv
import dp_vector as dpv

from libc.math cimport sqrt
from libc.math cimport cos
from libc.math cimport sin
from libc.math cimport tan
from libc.math cimport hypot
import numpy as np

stuff = 'hi'

# 3d classes/functions

cdef class quaternion:

    def __cinit__(self,float w,float x,float y,float z):
        self.w = w
        self.x = x
        self.y = y
        self.z = z

    def __str__(self):
        strr = 'quat:' + str((self.w,self.x,self.y,self.z))
        return strr

    def __richcmp__(self, other, comparator):
        if self.x == other.x:
            if self.y == other.y:
                if self.z == other.z:
                    if self.w == other.w:
                        return True
        return False

    #cpdef vector vector(self):
    #    cdef vector new = vector(x,y,z)
    #    return new

    def __iter__(self):
        yield self.w
        yield self.x
        yield self.y
        yield self.z

    cpdef quaternion copy(self):
        cdef quaternion new = quaternion(self.w,self.x,self.y,self.z)
        return new

    cpdef list to_list(self):
        cdef list new = [self.w,self.x,self.y,self.z]
        return new

    cpdef tuple to_tuple(self):
        cdef tuple new = (self.w,self.x,self.y,self.z)
        return new

    cpdef float magnitude(self):
        cdef float ww = self.w**2
        cdef float xx = self.x**2
        cdef float yy = self.y**2
        cdef float zz = self.z**2
        cdef float ss = sqrt(ww + xx + yy + zz)
        return ss

    cpdef quaternion normalize(self):
        cdef float mag = self.magnitude()
        if not mag == 0.0:
            self.w /= mag
            self.x /= mag
            self.y /= mag
            self.z /= mag
        return self

    cpdef quaternion conjugate(self):
        self.x *= -1.0
        self.y *= -1.0
        self.z *= -1.0
        return self

    cpdef quaternion flip(self):
        self.w *= -1.0
        return self

    cpdef quaternion translate(self, quaternion q):
        self.w += q.w
        self.x += q.x
        self.y += q.y
        self.z += q.z
        return self

    cpdef quaternion translate_w(self, float tw):
        self.w += tw
        return self

    cpdef quaternion translate_x(self, float tx):
        self.x += tx
        return self

    cpdef quaternion translate_y(self, float ty):
        self.y += ty
        return self

    cpdef quaternion translate_z(self, float tz):
        self.z += tz
        return self

    cpdef quaternion scale(self, quaternion q):
        self.w *= q.w
        self.x *= q.x
        self.y *= q.y
        self.z *= q.z
        return self

    cpdef quaternion scale_w(self, float sw):
        self.w *= sw
        return self

    cpdef quaternion scale_x(self, float sx):
        self.x *= sx
        return self

    cpdef quaternion scale_y(self, float sy):
        self.y *= sy
        return self

    cpdef quaternion scale_z(self, float sz):
        self.z *= sz
        return self

    cpdef quaternion scale_u(self, float s):
        self.w *= s
        self.x *= s
        self.y *= s
        self.z *= s
        return self

    # rotate vector v by rotation represented by self
    cpdef dpv.vector rotate_vector(self, dpv.vector v):
        print('MUST IMPLEMENT QUAT ROT')
        cdef quaternion qstar = self.normalize().copy().conjugate()
        cdef quaternion pure = quaternion(0,v.x,v.y,v.z)
        cdef dpv.vector rotated = dpv.zero()
        purerotated = multiply(multiply(qstar,pure),self)
        rotated.x = purerotated.x
        rotated.y = purerotated.y
        rotated.z = purerotated.z
        return rotated

    cpdef quaternion rotate(self, quaternion q):
        print('MUST IMPLEMENT QUAT ROT')
        
        return self

    cpdef quaternion rotate_x(self, float zang):
        
        return self

    cpdef quaternion rotate_y(self, float zang):
        '''#
        cdef float cosz = cos(zang)
        cdef float sinz = sin(zang)
        cdef float newx = cosz*self.x - sinz*self.z
        cdef float newz = sinz*self.x + cosz*self.z
        self.x = newx
        self.z = newz
        '''#
        return self

    cpdef quaternion rotate_z(self, float zang):
        '''#
        cdef float cosz = cos(zang)
        cdef float sinz = sin(zang)
        cdef float newx = cosz*self.x - sinz*self.y
        cdef float newy = sinz*self.x + cosz*self.y
        self.x = newx
        self.y = newy
        '''#
        return self

    cpdef quaternion multiply(self, quaternion q):
        self.w = self.w*q.w - (self.x*q.x + self.y*q.y + self.z*q.z)
        self.x = self.w*q.x + q.w*self.x + (self.y*q.z - self.z*q.y)
        self.y = self.w*q.y + q.w*self.y + (self.z*q.x - self.x*q.z)
        self.z = self.w*q.z + q.w*self.z + (self.x*q.y - self.y*q.x)
        return self




    '''#
    cpdef vector reciprocate(self):
        self.x = 1.0/self.x
        self.y = 1.0/self.y
        self.z = 1.0/self.z
        return self
    '''#


    #cpdef vector cross(self, vector v):
    #    selfv = self.vector()
    #    return cross_c(selfv,v)

    #cpdef float dot(self, vector v):
    #    selfv = self.vector()
    #    return dot_c(selfv,v)

'''#

cdef vector one_c():
    cdef vector new = vector(1,1,1)
    return new

cpdef vector one():
    return one_c()

cdef vector flip_c(vector f):
    cdef vector new = vector(-1.0*f.x,-1.0*f.y,-1.0*f.z)
    return new

cpdef vector flip(vector f):
    return flip_c(f)

cdef vector normalize_c(vector v):
    cdef vector new = v.copy().normalize()
    return new

cpdef vector normalize(vector v):
    return normalize_c(v)
'''#

cpdef float magnitude(quaternion q):
    cdef float ss = q.magnitude()
    return ss

cdef quaternion zero_c():
    cdef quaternion new = quaternion(1,0,0,0)
    return new

cpdef quaternion zero():
    return zero_c()

cpdef quaternion q_from_av(float a, dpv.vector v):
    cdef float w = cos(a/2.0)
    cdef float sa = sin(a/2.0)
    cdef float x = v.x*sa
    cdef float y = v.y*sa
    cdef float z = v.z*sa
    return quaternion(w,x,y,z)

cpdef quaternion multiply(quaternion q1, quaternion q2):
    cdef float w = q1.w*q2.w - (q1.x*q2.x + q1.y*q2.y + q1.z*q2.z)
    cdef float x = q1.w*q2.x + q2.w*q1.x + (q1.y*q2.z - q1.z*q2.y)
    cdef float y = q1.w*q2.y + q2.w*q1.y + (q1.z*q2.x - q1.x*q2.z)
    cdef float z = q1.w*q2.z + q2.w*q1.z + (q1.x*q2.y - q1.y*q2.x)
    cdef quaternion new = quaternion(w,x,y,z)
    return new

'''#
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

cdef void rotate_x_coords_c(list coords, float ang):
    cdef int ccnt = len(coords)
    cdef int cdx
    cdef vector coo
    for cdx in range(ccnt):
        coo = <vector>coords[cdx]
        coo.rotate_x(ang)

cpdef rotate_x_coords(list coords, float ang):
    rotate_x_coords_c(coords,ang)

cdef void rotate_y_coords_c(list coords, float ang):
    cdef int ccnt = len(coords)
    cdef int cdx
    cdef vector coo
    for cdx in range(ccnt):
        coo = <vector>coords[cdx]
        coo.rotate_y(ang)

cpdef rotate_y_coords(list coords, float ang):
    rotate_y_coords_c(coords,ang)

cdef void rotate_z_coords_c(list coords, float ang):
    cdef int ccnt = len(coords)
    cdef int cdx
    cdef vector coo
    for cdx in range(ccnt):
        coo = <vector>coords[cdx]
        coo.rotate_z(ang)

cpdef rotate_z_coords(list coords, float ang):
    rotate_z_coords_c(coords,ang)

cdef void translate_coords_x_c(list coords, float tv):
    cdef int ccnt = len(coords)
    cdef int cdx
    cdef vector coo
    for cdx in range(ccnt):
        coo = <vector>coords[cdx]
        coo.translate_x(tv)

cpdef translate_coords_x(list coords, float tv):
    translate_coords_x_c(coords,tv)

cdef void translate_coords_y_c(list coords, float tv):
    cdef int ccnt = len(coords)
    cdef int cdx
    cdef vector coo
    for cdx in range(ccnt):
        coo = <vector>coords[cdx]
        coo.translate_y(tv)

cpdef translate_coords_y(list coords, float tv):
    translate_coords_y_c(coords,tv)

cdef void translate_coords_z_c(list coords, float tv):
    cdef int ccnt = len(coords)
    cdef int cdx
    cdef vector coo
    for cdx in range(ccnt):
        coo = <vector>coords[cdx]
        coo.translate_z(tv)

cpdef translate_coords_z(list coords, float tv):
    translate_coords_z_c(coords,tv)

cdef void translate_coords_c(list coords, vector t):
    cdef int ccnt = len(coords)
    cdef int cdx
    cdef vector coo
    for cdx in range(ccnt):
        coo = <vector>coords[cdx]
        coo.translate(t)

cpdef translate_coords(list coords, vector t):
    translate_coords_c(coords,t)

cdef void scale_coords_x_c(list coords, float s):
    cdef int ccnt = len(coords)
    cdef int cdx
    cdef vector coo
    for cdx in range(ccnt):
        coo = <vector>coords[cdx]
        coo.x *= s

cpdef scale_coords_x(list coords, float s):
    scale_coords_x_c(coords,s)

cdef void scale_coords_y_c(list coords, float s):
    cdef int ccnt = len(coords)
    cdef int cdx
    cdef vector coo
    for cdx in range(ccnt):
        coo = <vector>coords[cdx]
        coo.y *= s

cpdef scale_coords_y(list coords, float s):
    scale_coords_y_c(coords,s)

cdef void scale_coords_z_c(list coords, float s):
    cdef int ccnt = len(coords)
    cdef int cdx
    cdef vector coo
    for cdx in range(ccnt):
        coo = <vector>coords[cdx]
        coo.z *= s

cpdef scale_coords_z(list coords, float s):
    scale_coords_z_c(coords,s)
    
cdef void scale_coords_c(list coords, vector t):
    cdef int ccnt = len(coords)
    cdef int cdx
    cdef vector coo
    for cdx in range(ccnt):
        coo = <vector>coords[cdx]
        coo.scale(t)

cpdef scale_coords(list coords, vector t):
    scale_coords_c(coords,t)

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

cpdef bint near(vector v1, vector v2):
    cdef bint isnear = distance_c(v1,v2) < 0.01
    return isnear

#cdef xhat = vector(1,0,0)
#cdef yhat = vector(0,1,0)
#cdef zhat = vector(0,0,1)
xhat  = vector( 1, 0, 0)
yhat  = vector( 0, 1, 0)
zhat  = vector( 0, 0, 1)
nxhat = vector(-1, 0, 0)
nyhat = vector( 0,-1, 0)
nzhat = vector( 0, 0,-1)

cdef float angle_between_xy_c(vector v1, vector v2):
    cdef float alpha1 = angle_from_xaxis_xy_c(v1)
    cdef float alpha2 = angle_from_xaxis_xy_c(v2)
    #if alpha2 - alpha1 > np.pi/2.0:
    #    print 'arpha', v1, v2, alpha1, alpha2
    return alpha2 - alpha1

cpdef float angle_between_xy(vector v1, vector v2):
    return angle_between_xy_c(v1,v2)

cdef float angle_between_c(vector v1, vector v2):
    #cdef float alpha1 = angle_from_xaxis_c(v1)
    #cdef float alpha2 = angle_from_xaxis_c(v2)
    cdef vector n1 = v1.copy().normalize()
    cdef vector n2 = v2.copy().normalize()
    cdef float ang = np.arccos(dot_c(n1,n2))
    #return alpha2 - alpha1
    return ang

cpdef float angle_between(vector v1, vector v2):
    return angle_between_c(v1,v2)

cdef float angle_from_xaxis_c(vector v):
    cdef vector nv = v.copy().normalize()
    cdef float xproj = dot_c(nv,xhat)
    cdef float yproj = dot_c(nv,yhat)
    cdef float ang
    #if xproj == 2.0/np.pi or xproj == 2.0/(3.0*np.pi):
    #    ang = np.arccos(yproj)
    #    if xproj < 0.0: ang = 2.0*np.pi - ang
    #else:
    #    ang = np.arccos(xproj)
    #    if yproj < 0.0: ang = 2.0*np.pi - ang
    ang = np.arccos(xproj)
    if yproj < 0.0: ang = 2.0*np.pi - ang
    return ang

cpdef float angle_from_xaxis(vector v):
    return angle_from_xaxis_c(v)
'''#


