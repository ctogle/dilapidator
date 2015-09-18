# cython: profile=True
#cimport cython

cimport dilap.core.tools as dpr

from dilap.geometry.vec3 cimport vec3

from libc.math cimport sqrt
from libc.math cimport sqrt
from libc.math cimport cos
from libc.math cimport sin
from libc.math cimport tan
from libc.math cimport hypot
import numpy as np

stuff = 'hi'










# dilapidators implementation of a quaternion in R3
cdef class quat:

    ###########################################################################
    ### python object overrides ###############################################
    ###########################################################################

    def __str__(self):return 'quat:'+str(tuple(self))
    def __iter__(self):yield self.w;yield self.x;yield self.y;yield self.z

    # should return a new quat representing rot by self and then o
    #def __add__(self,o):return vec3(self.x+o.x,self.y+o.y,self.z+o.z)

    # should return a new quat representing rot by self and then o-inv
    #def __sub__(self,o):return vec3(self.x-o.x,self.y-o.y,self.z-o.z)

    # should return a new quat representing this one rotated by o
    #def __mul__(self,o):return vec3(self.x*o.x,self.y*o.y,self.z*o.z)
    #cpdef quaternion multiply(self, quaternion q):
    #    self.w = self.w*q.w - (self.x*q.x + self.y*q.y + self.z*q.z)
    #    self.x = self.w*q.x + q.w*self.x + (self.y*q.z - self.z*q.y)
    #    self.y = self.w*q.y + q.w*self.y + (self.z*q.x - self.x*q.z)
    #    self.z = self.w*q.z + q.w*self.z + (self.x*q.y - self.y*q.x)
    #    return self

    def __is_equal(self,o):return self.isnear(o)
    def __richcmp__(x,y,op):
        if op == 2:return x.__is_equal(y)
        else:assert False

    ###########################################################################

    ###########################################################################
    ### c methods #############################################################
    ###########################################################################

    def __cinit__(self,float w,float x,float y,float z):
        self.w = w
        self.x = x
        self.y = y
        self.z = z

    # given an angle a and axis v, modify to represent 
    # a rotation by a around v and return self
    # NOTE: negative a still maps to positive self.w
    cdef quat av_c(self,float a,vec3 v):
        cdef float a2 = a/2.0
        cdef float sa = sin(a2)
        cdef float vm = v.mag()
        self.w = cos(a2)
        self.x = v.x*sa/vm
        self.y = v.y*sa/vm
        self.z = v.z*sa/vm
        return self

    # modify to represent a rotation 
    # between two vectors and return self
    cdef quat uu_c(self,vec3 x,vec3 y):
        cdef float a = x.ang_c(y)
        cdef vec3 v = x.crs_c(y)
        return self.av_c(a,v)

    # return an independent copy of this quaternion
    cdef quat cp_c(self):
        cdef quat n = quat(self.w,self.x,self.y,self.z)
        return n

    # is quat o within a very small neighborhood of self
    cdef bint isnear_c(self,quat o):
        cdef float dw = (self.w-o.w)
        if dw*dw > dpr.epsilonsq_c:return 0
        cdef float dx = (self.x-o.x)
        if dx*dx > dpr.epsilonsq_c:return 0
        cdef float dy = (self.y-o.y)
        if dy*dy > dpr.epsilonsq_c:return 0
        cdef float dz = (self.z-o.z)
        if dz*dz > dpr.epsilonsq_c:return 0
        return 1

    # return the squared magintude of self
    cdef float mag2_c(self):
        cdef float w2 = self.w*self.w
        cdef float x2 = self.x*self.x
        cdef float y2 = self.y*self.y
        cdef float z2 = self.z*self.z
        cdef float m2 = w2 + x2 + y2 + z2
        return m2

    # return the magintude of self
    cdef float mag_c(self):
        return sqrt(self.mag2_c())

    # normalize and return self
    cdef quat nrm_c(self):
        cdef float m = self.mag_c()
        if m == 0.0:return self
        else:return self.scl_c(1.0/m)

    # flip the direction of and return self
    cdef quat flp_c(self):
        self.w *= -1.0
        return self

    # multiply each component by a scalar of and return self
    cdef quat scl_c(self,float s):
        self.w *= s
        self.x *= s
        self.y *= s
        self.z *= s
        return self

    # conjugate and return self
    cdef quat cnj_c(self):
        self.x *= -1.0
        self.y *= -1.0
        self.z *= -1.0
        return self

    # given quat o, rotate self so that self represents
    # a rotation by self and then q (q * self)
    cdef quat rot_c(self,quat o):
        print('MUST IMPLEMENT QUAT ROT')
        raise NotImplemented
        '''#
        cdef rotw = q.w*self.w - q.x*self.x - q.y*self.y - q.z*self.z
        cdef rotx = q.w*self.x + q.x*self.w + q.y*self.z - q.z*self.y
        cdef roty = q.w*self.y - q.x*self.z + q.y*self.w + q.z*self.x
        cdef rotz = q.w*self.z + q.x*self.y - q.y*self.x + q.z*self.w
        if self.magnitude() < 0.1:rotw,rotx,roty,rotz = q.to_tuple()
        self.w = rotw
        self.x = rotx
        self.y = roty
        self.z = rotz
        '''#
        return self

    ###########################################################################

    ###########################################################################
    ### python wrappers for c methods #########################################
    ###########################################################################

    # given an angle a and axis v, modify to represent 
    # a rotation by a around v and return self
    cpdef quat av(self,float a,vec3 v):
        '''modify to represent a rotation about a vector by an angle'''
        return self.av_c(a,v)

    # modify to represent a rotation 
    # between two vectors and return self
    cpdef quat uu(self,vec3 x,vec3 y):
        '''modify to represent a rotation from one vector to another'''
        return self.uu_c(x,y)

    # return an independent copy of this quaternion
    cpdef quat cp(self):
        '''create an independent copy of this quaternion'''
        return self.cp_c()

    # is quat o within a very small neighborhood of self
    cpdef bint isnear(self,quat o):
        '''determine if a point is numerically close to another'''
        return self.isnear_c(o)
    
    # return the squared magintude of self
    cpdef float mag2(self):
        '''compute the squared magnitude of this quaternion'''
        return self.mag2_c()

    # return the magintude of self
    cpdef float mag(self):
        '''compute the magnitude of this quaternion'''
        return self.mag_c()

    # normalize and return self
    cpdef quat nrm(self):
        '''normalize this quaternion'''
        return self.nrm_c()

    # flip the direction of and return self
    cpdef quat flp(self):
        '''flip the direction of rotation represented by this quaternion'''
        return self.flp_c()

    # multiply each component by a scalar of and return self
    cpdef quat scl(self,float s):
        '''multiply components of this point by a scalar'''
        return self.scl_c(s)

    # conjugate and return self
    cpdef quat cnj(self):
        '''conjugate this quaternion'''
        return self.cnj_c()

    # given quat o, rotate self so that self represents
    # a rotation by self and then q (q * self)
    cpdef quat rot(self,quat o):
        '''rotate this quaternion by another quaternion'''
        return self.rot_c(o)

    ###########################################################################










