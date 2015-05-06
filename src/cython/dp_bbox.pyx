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

    def __cinit__(self,dpv.vector2d x,dpv.vector2d y,dpv.vector2d z):
        self.x = x
        self.y = y
        self.z = z

    def __str__(self):
        xstr = (self.x.x,self.x.y)
        ystr = (self.y.x,self.y.y)
        zstr = (self.z.x,self.z.y)
        strr = 'bbox:' + str((xstr,ystr,zstr))
        return strr


