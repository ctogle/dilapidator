#imports
# cython: profile=True
#cimport cython

from dilap.geometry.vec3 cimport vec3
from dilap.geometry.quat cimport quat

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

__doc__ = '''General purpose tool functions...'''

###############################################################################
### c space
###############################################################################

# compute the center of mass for a set of vectors
cdef vec3 com_c(ps):
    cdef int pcnt = len(ps)
    cdef int px
    cdef float pcntf = float(pcnt)
    cdef vec3 n = vec3(0,0,0)
    for px in range(pcnt):n.trn_c(ps[px])
    return n.uscl_c(1.0/pcntf)

###############################################################################
### python space
###############################################################################

# compute the center of mass for a set of vectors
cpdef vec3 com(ps):
    '''compute the center of mass for a set of vectors'''
    return com_c(ps)










