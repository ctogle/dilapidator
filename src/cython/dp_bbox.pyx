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

import matplotlib.pyplot as plt

stuff = 'hi'

# 3d classes/functions

cdef class bbox:

    def __str__(self):
        xstr = (self.x.x,self.x.y)
        ystr = (self.y.x,self.y.y)
        zstr = (self.z.x,self.z.y)
        strr = 'bbox:' + str((xstr,ystr,zstr))
        return strr

    def __cinit__(self,dpv.vector2d x,dpv.vector2d y,dpv.vector2d z):
        self.x = x
        self.y = y
        self.z = z

    # modify self.x to encompass proj
    cpdef bbox _consume_x(self,dpv.vector2d proj):
        if self.x.x > proj.x:
            self.x.x = proj.x
        if self.x.y < proj.y:
            self.x.y = proj.y
        return self

    # modify self.y to encompass proj
    cpdef bbox _consume_y(self,dpv.vector2d proj):
        if self.y.x > proj.x:
            self.y.x = proj.x
        if self.y.y < proj.y:
            self.y.y = proj.y
        return self

    # modify self.z to encompass proj
    cpdef bbox _consume_z(self,dpv.vector2d proj):
        if self.z.x > proj.x:
            self.z.x = proj.x
        if self.z.y < proj.y:
            self.z.y = proj.y
        return self

    # modify self to encompass other
    cpdef bbox _consume(self,bbox other):
        self._consume_x(other.x)
        self._consume_y(other.y)
        self._consume_z(other.z)
        return self

cdef bbox zero_c():
    cdef dpv.vector2d z = dpv.zero2d()
    cdef bbox new = bbox(z.copy(),z.copy(),z.copy())
    return new

cpdef bbox zero():
    return zero_c()


