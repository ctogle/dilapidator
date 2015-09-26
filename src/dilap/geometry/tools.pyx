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

stuff = 'hi'

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

###############################################################################
### python space
###############################################################################

# compute the center of mass for a set of vectors
cpdef vec3 com(ps):
    '''compute the center of mass for a set of vectors'''
    return com_c(ps)

# return a ring of points of radius r with n corners
cpdef list pline(vec3 s,vec3 e,int n):
    '''return a line of n points starting with s, ending with e'''
    return pline_c(s,e,n)









