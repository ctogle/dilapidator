# cython: profile=True
#cimport cython

cimport dilap.geometry.tools as dpr

from dilap.geometry.vec3 cimport vec3
from dilap.geometry.quat cimport quat

from libc.math cimport sqrt
from libc.math cimport cos
from libc.math cimport sin
from libc.math cimport tan
from libc.math cimport hypot

cimport numpy
import numpy

stuff = 'hi'



__doc__ = '''dilapidator\'s implementation of a ray in R3'''
# dilapidators implementation of a ray in R3
cdef class ray3:

    ###########################################################################
    ### python object overrides ###############################################
    ###########################################################################

    def __str__(self):return 'ray:'+self.o.__str__()+','+self.d.__str__()
    #def __iter__(self):yield self.x;yield self.y;yield self.z
    #def __add__(self,o):return vec3(self.x+o.x,self.y+o.y,self.z+o.z)
    #def __sub__(self,o):return vec3(self.x-o.x,self.y-o.y,self.z-o.z)
    #def __mul__(self,o):return vec3(self.x*o.x,self.y*o.y,self.z*o.z)
    def __is_equal(self,o):return self.isnear(o)
    def __richcmp__(x,y,op):
        if op == 2:return x.__is_equal(y)
        else:assert False

    ###########################################################################

    ###########################################################################
    ### c methods #############################################################
    ###########################################################################

    def __cinit__(self,vec3 o,vec3 d):
        self.o = o 
        self.d = d.nrm_c()
        self.c = vec3(-1,0,0)

    # return an independent copy of this ray
    cdef ray3 cp_c(self):
        cdef ray3 n = ray3(self.o.cp_c(),self.d.cp_c())
        return n

    # is ray3 o within a very small neighborhood of self (both self.o and self.d)
    cdef bint isnear_c(self,ray3 o):
        return self.o.isnear_c(o.o) and self.d.isnear_c(o.d)

    # perform a ray triangle intersection test:
    # if intesecting, set self.c to contain the barycentric coords of the hit 
    # and the distance from the hit to self.o and return 1, otherwise return 0
    cdef bint hittri_c(self,vec3 t1,vec3 t2,vec3 t3):
        cdef vec3 e1 = t1.tov_c(t2)
        cdef vec3 e2 = t1.tov_c(t3)
        cdef vec3 T  = t1.tov_c(self.o)
        cdef vec3 P  = self.d.crs_c(e2)
        cdef vec3 Q  = T.crs_c(e1)
        cdef float d = dpr.near_c(P.dot_c(e1),0,dpr.epsilon_c)
        cdef float t = dpr.near_c(Q.dot_c(e2),0,dpr.epsilon_c)
        cdef float u = dpr.near_c(P.dot_c(T),0,dpr.epsilon_c)
        cdef float v = dpr.near_c(Q.dot_c(self.d),0,dpr.epsilon_c)
        if d*d < dpr.epsilonsq_c:return 0
        else:d = 1/d
        t *= d;u *= d;v *= d
        if (t < 0) or (u < 0 or u > 1) or (v < 0 or u + v > 1):
            self.c.x = -1;self.c.y = 0;self.c.z = 0
            return 0
        else:
            self.c.x = t;self.c.y = u;self.c.z = v
            return 1

    # given triangles, return the indices 
    # of the subset which are hit by ray r
    cdef list hittris_c(self,list tris):
        cdef int tcnt = len(tris)
        cdef int tx
        cdef vec3 t1,t2,t3
        cdef list hxs = []
        for tx in range(tcnt):
            t1,t2,t3 = tris[tx]
            if self.hittri_c(t1,t2,t3):
                hxs.append(tx)
        return hxs

    # return index for the closest triangle 
    # to self.o which this ray intersects
    cdef int hittri_close_c(self,list tris):
        cdef int tcnt = len(tris)
        cdef int tx,hx
        cdef float close = dpr.maxfloat_c
        cdef vec3 t1,t2,t3,hcst
        for tx in range(tcnt):
            t1,t2,t3 = tris[tx]
            if self.hittri_c(t1,t2,t3):
                if self.c.x < close:
                    close = self.c.x
                    hx = tx
                    hcst = self.c.cp_c()
        return hx

    # perform a ray plane intersection test:
    # if intesecting, set self.c to contain the R3 location 
    # of the hit and return 1, otherwise return 0
    cdef bint hitpln_c(self,vec3 r,vec3 n):
        cdef vec3 p0 = self.o.cp_c()
        cdef vec3 p1 = self.o.cp_c().trn_c(self.d)
        cdef float den = n.dot_c(p0.tov_c(p1))
        cdef float num = n.dot_c(p0.tov_c(r))
        if dpr.isnear_c(den,0,dpr.epsilon_c):return 0
        else:t = num/den
        if t < 0:
            self.c.x = -1;self.c.y = 0;self.c.z = 0
            return 0
        else:
            self.c.x = self.o.x + t*self.d.x
            self.c.y = self.o.y + t*self.d.y
            self.c.z = self.o.z + t*self.d.z
            return 1

    ###########################################################################

    ###########################################################################
    ### python wrappers for c methods #########################################
    ###########################################################################

    # return an independent copy of this ray
    cpdef ray3 cp(self):
        '''create an independent copy of this ray'''
        return self.cp_c()

    # is ray3 o within a very small neighborhood of self (both self.o and self.d)
    cpdef bint isnear(self,ray3 o):
        '''determine if a ray is near in origin and direction to another ray'''
        return self.isnear_c(o)

    # perform a ray triangle intersection test:
    # if intesecting, set self.c to contain the barycentric coords of the hit 
    # and the distance from the hit to self.o and return 1, otherwise return 0
    cpdef bint hittri(self,vec3 t1,vec3 t2,vec3 t3):
        '''perform a hit test on a triangle, given its corners'''
        return self.hittri_c(t1,t2,t3)

    # given triangles, return the indices 
    # of the subset which are hit by ray r
    cpdef list hittris(self,list tris):
        '''perform a hit test on a set of triangles'''
        return self.hittris_c(tris)

    # return index for the closest triangle 
    # to self.o which this ray intersects
    cpdef int hittri_close(self,list tris):
        '''perform a hit test on a set of triangles to find the nearest hit'''
        return self.hittri_close_c(tris)

    # perform a ray plane intersection test:
    # if intesecting, set self.c to contain the R3 location 
    # of the hit and return 1, otherwise return 0
    cpdef bint hitpln(self,vec3 r,vec3 n):
        '''perform a hit test on a plane, given its normal and point in it'''
        return self.hitpln_c(r,n)

    ###########################################################################










