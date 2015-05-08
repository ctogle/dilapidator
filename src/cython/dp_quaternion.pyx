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
    cdef dpv.vector d = v.copy().normalize()
    cdef float x = d.x*sa
    cdef float y = d.y*sa
    cdef float z = d.z*sa
    return quaternion(w,x,y,z)

cpdef quaternion multiply(quaternion q1, quaternion q2):
    cdef float w = q1.w*q2.w - (q1.x*q2.x + q1.y*q2.y + q1.z*q2.z)
    cdef float x = q1.w*q2.x + q2.w*q1.x + (q1.y*q2.z - q1.z*q2.y)
    cdef float y = q1.w*q2.y + q2.w*q1.y + (q1.z*q2.x - q1.x*q2.z)
    cdef float z = q1.w*q2.z + q2.w*q1.z + (q1.x*q2.y - q1.y*q2.x)
    cdef quaternion new = quaternion(w,x,y,z)
    return new

# this needs more work
def quat_test():

    i = quaternion(0,1,0,0)
    j = quaternion(0,0,1,0)
    k = quaternion(0,0,0,1)
    print('qtest ij = :',i.copy().multiply(j).__str__())

    x = dpv.xhat.copy()
    a = np.pi/2
    v = dpv.yhat.copy()
    q = q_from_av(a,v)

    #z = q.rotate_vector(x)
    z = x.copy().rotate(q)
    print('quat_test:')
    print(x.__str__(),'rotated by',a.__str__())
    print('about axis',v.__str__(),'yields',z.__str__())


