# cython: profile=True
#cimport cython

cimport dilap.core.tools as dpr

from dilap.geometry.quat cimport quat

from libc.math cimport sqrt
from libc.math cimport cos
from libc.math cimport sin
from libc.math cimport tan
from libc.math cimport hypot

cimport numpy
import numpy
 
stuff = 'hi'










__doc__ = '''dilapidator\'s implementation of a vector in R3'''
# dilapidators implementation of a vector in R3
cdef class vec3:

    ###########################################################################
    ### python object overrides ###############################################
    ###########################################################################

    def __str__(self):return 'vec3:'+str(tuple(self))
    def __iter__(self):yield self.x;yield self.y;yield self.z
    def __add__(self,o):return vec3(self.x+o.x,self.y+o.y,self.z+o.z)
    def __sub__(self,o):return vec3(self.x-o.x,self.y-o.y,self.z-o.z)
    def __mul__(self,o):return self.cp().mul(o)
    def __is_equal(self,o):return self.isnear(o)
    # could use <,> for lexicographic ordering?
    def __richcmp__(x,y,op):
        if op == 2:return x.__is_equal(y)
        else:assert False

    ###########################################################################

    ###########################################################################
    ### c methods #############################################################
    ###########################################################################

    def __cinit__(self,float x,float y,float z):
        self.x = x
        self.y = y
        self.z = z

    # return an independent copy of this point
    cdef vec3 cp_c(self):
        cdef vec3 n = vec3(self.x,self.y,self.z)
        return n

    # return an independent copy of the xy projection of this point
    cdef vec3 cpxy_c(self):
        cdef vec3 n = vec3(self.x,self.y,0.0)
        return n

    # return the R3 euclidean distance between self and vec3 o
    cdef float d_c(self,vec3 o):
        cdef float dx = self.x - o.x
        cdef float dy = self.y - o.y
        cdef float dz = self.z - o.z
        cdef float ds = sqrt(dx*dx + dy*dy + dz*dz)
        return ds

    # return the R3 euclidean distance between the xy projections of self and vec3 o
    cdef float dxy_c(self,vec3 o):
        cdef float dx = self.x - o.x
        cdef float dy = self.y - o.y
        cdef float ds = sqrt(dx*dx + dy*dy)
        return ds

    # return the angle between self and vec3 o
    cdef float ang_c(self,vec3 o):
        cdef float sm = self.mag_c()
        cdef float om = o.mag_c()
        cdef float sod = (self.x*o.x + self.y*o.y + self.z*o.z)/(sm*om)
        cdef float a
        if   dpr.isnear_c(sod, 1.0):a = 0.0
        elif dpr.isnear_c(sod,-1.0):a = dpr.PI
        else:a = numpy.arccos(sod)
        return a

    # return the angle between the xy projections of self and vec3 o
    cdef float angxy_c(self,vec3 o):
        return self.cpxy_c().ang_c(o.cpxy_c())

    # return the dot product of self and vec3 o
    cdef float dot_c(self,vec3 o):
        return self.x*o.x + self.y*o.y + self.z*o.z

    # return the cross product of self and vec3 o
    cdef vec3 crs_c(self,vec3 o):
        cdef float cx = self.y*o.z-self.z*o.y
        cdef float cy = self.z*o.x-self.x*o.z
        cdef float cz = self.x*o.y-self.y*o.x
        cdef vec3 n = vec3(cx,cy,cz)
        return n

    # project into a plane and return self
    cdef vec3 prj_c(self,vec3 r,vec3 n):
        cdef float d = (self.x-r.x)*n.x+(self.y-r.y)*n.y+(self.z-r.z)*n.z
        cdef vec3 pj = n.cp_c().scl_c(-d)
        return self.trn_c(pj)

    # 1-1 multiplication by vec3 o
    cdef vec3 mul_c(self,vec3 o):
        #cdef vec3 n = vec3(self.x*o.x,self.y*o.y,self.z*o.z)
        #return n
        self.x *= o.x;self.y *= o.y;self.z *= o.z
        return self

    # is vec3 o within an open ball of raidus e centered at self
    cdef bint inneighborhood_c(self,vec3 o,float e):
        cdef float d = self.d_c(o)
        if d < e:return 1
        else:return 0

    # is vec3 o within a very small neighborhood of self
    cdef bint isnear_c(self,vec3 o):
        cdef float dx = (self.x-o.x)
        if dx*dx > dpr.epsilonsq_c:return 0
        cdef float dy = (self.y-o.y)
        if dy*dy > dpr.epsilonsq_c:return 0
        cdef float dz = (self.z-o.z)
        if dz*dz > dpr.epsilonsq_c:return 0
        return 1

    # is vec3 o within a very small neighborhood of self in the xy plane
    cdef bint isnearxy_c(self,vec3 o):
        cdef float dx = (self.x-o.x)
        if dx*dx > dpr.epsilonsq_c:return 0
        cdef float dy = (self.y-o.y)
        if dy*dy > dpr.epsilonsq_c:return 0
        return 1

    # return the squared magintude of self
    cdef float mag2_c(self):
        cdef float x2 = self.x*self.x
        cdef float y2 = self.y*self.y
        cdef float z2 = self.z*self.z
        cdef float m2 = x2 + y2 + z2
        return m2

    # return the magintude of self
    cdef float mag_c(self):
        return sqrt(self.mag2_c())

    # normalize and return self
    cdef vec3 nrm_c(self):
        cdef float m = self.mag_c()
        if m == 0.0:return self
        else:return self.scl_c(1.0/m)

    # translate self by vec3 o
    cdef vec3 trn_c(self,vec3 o):
        self.x += o.x
        self.y += o.y
        self.z += o.z
        return self

    # translate self in the x direction by d
    cdef vec3 xtrn_c(self,float d):
        self.x += d
        return self

    # translate self in the y direction by d
    cdef vec3 ytrn_c(self,float d):
        self.y += d
        return self

    # translate self in the z direction by d
    cdef vec3 ztrn_c(self,float d):
        self.z += d
        return self

    # multiply each component by a scalar of and return self
    cdef vec3 scl_c(self,float s):
        self.x *= s
        self.y *= s
        self.z *= s
        return self

    # scale self.x by scalar s
    cdef vec3 xscl_c(self,float s):
        self.x *= s
        return self

    # scale self.y by scalar s
    cdef vec3 yscl_c(self,float s):
        self.y *= s
        return self

    # scale self.z by scalar s
    cdef vec3 zscl_c(self,float s):
        self.z *= s
        return self

    # rotate by a quaternion q and return self
    cdef vec3 rot_c(self,quat q):
        if dpr.isnear_c(q.w,0):return self
        cdef float row1x = q.w**2 + q.x**2 - q.y**2 - q.z**2
        cdef float row1y = 2*(q.x*q.y - q.w*q.z)
        cdef float row1z = 2*(q.x*q.z + q.w*q.y)
        cdef float row2x = 2*(q.x*q.y + q.w*q.z)
        cdef float row2y = q.w**2 - q.x**2 + q.y**2 - q.z**2
        cdef float row2z = 2*(q.y*q.z - q.w*q.x)
        cdef float row3x = 2*(q.x*q.z - q.w*q.y)
        cdef float row3y = 2*(q.y*q.z + q.w*q.x)
        cdef float row3z = q.w**2 - q.x**2 - q.y**2 + q.z**2
        cdef float nx = self.x*row1x + self.y*row1y + self.z*row1z
        cdef float ny = self.x*row2x + self.y*row2y + self.z*row2z
        cdef float nz = self.x*row3x + self.y*row3y + self.z*row3z
        self.x = nx;self.y = ny;self.z = nz
        return self

    # rotate around the x axis by an angle a and return self
    # NOTE: this is slower than using rot with the appropriate quaternion
    cdef vec3 xrot_c(self,float a):
        cdef float ca = cos(a)
        cdef float sa = sin(a)
        cdef float ny = ca*self.y - sa*self.z
        cdef float nz = sa*self.y + ca*self.z
        self.y = ny;self.z = nz
        return self

    # rotate around the y axis by an angle a and return self
    # NOTE: this is slower than using rot with the appropriate quaternion
    cdef vec3 yrot_c(self,float a):
        cdef float ca = cos(a)
        cdef float sa = sin(a)
        cdef float nx =  ca*self.x + sa*self.z
        cdef float nz = -sa*self.x + ca*self.z
        self.x = nx;self.z = nz
        return self

    # rotate around the z axis by an angle a and return self
    # NOTE: this is slower than using rot with the appropriate quaternion
    cdef vec3 zrot_c(self,float a):
        cdef float ca = cos(a)
        cdef float sa = sin(a)
        cdef float nx = ca*self.x - sa*self.y
        cdef float ny = sa*self.x + ca*self.y
        self.x = nx;self.y = ny
        return self

    # flip the direction of and return self
    cdef vec3 flp_c(self):
        return self.scl_c(-1.0)

    # return a vector point from self to vec3 o
    cdef vec3 tov_c(self,vec3 o):
        cdef float dx = o.x - self.x
        cdef float dy = o.y - self.y
        cdef float dz = o.z - self.z
        cdef vec3 n = vec3(dx,dy,dz)
        return n

    # return a vector at the midpoint between self and vec3 o
    cdef vec3 mid_c(self,vec3 o):
        return self.lerp_c(o,0.5)

    # linearly interpolate between self and vec3 o proportionally to ds
    cdef vec3 lerp_c(self,vec3 o,float ds):
        cdef float dx = self.x + (o.x - self.x)*ds
        cdef float dy = self.y + (o.y - self.y)*ds
        cdef float dz = self.z + (o.z - self.z)*ds
        cdef vec3 n = vec3(dx,dy,dz)
        return n

    ###########################################################################

    ###########################################################################
    ### python wrappers for c methods #########################################
    ###########################################################################

    # return an independent copy of this point
    cpdef vec3 cp(self):
        '''create an independent copy of this point'''
        return self.cp_c()

    # return an independent copy of the xy projection of this point
    cpdef vec3 cpxy(self):
        '''create an independent copy of the xy projection of this point'''
        return self.cpxy_c()

    # return the R3 euclidean distance between self and vec3 o
    cpdef float d(self,vec3 o):
        '''determine the R3 euclidean distance between this point and another'''
        return self.d_c(o)

    # return the R3 euclidean distance between the xy projections of self and vec3 o
    cpdef float dxy(self,vec3 o):
        '''determine the R3 euclidean distance between the xy projections of this point and another'''
        return self.dxy_c(o)

    # return the angle between self and vec3 o
    cpdef float ang(self,vec3 o):
        '''determine the angle between this point and another'''
        return self.ang_c(o)

    # return the angle between the xy projections of self and vec3 o
    cpdef float angxy(self,vec3 o):
        '''determine the angle between the xy projections of this point and another'''
        return self.angxy_c(o)

    # return the dot product of self and vec3 o
    cpdef float dot(self,vec3 o):
        '''compute the dot product of this vector and another'''
        return self.dot_c(o)

    # return the cross product of self and vec3 o
    cpdef vec3 crs(self,vec3 o):
        '''generate the cross product of this vector and another'''
        return self.crs_c(o)

    # project into a plane and return self
    cpdef vec3 prj(self,vec3 r,vec3 n):
        '''project this point into a plane'''
        return self.prj_c(r,n)

    # 1-1 multiplication by vec3 o
    cpdef vec3 mul(self,vec3 o):
        '''1-1 multiplication by another vector'''
        return self.mul_c(o)

    # is vec3 o within an open ball of raidus e centered at self
    cpdef bint inneighborhood(self,vec3 o,float e):
        '''determine if a point lies in an open ball centered at this point'''
        return self.inneighborhood_c(o,e)

    # is vec3 o within a very small neighborhood of self
    cpdef bint isnear(self,vec3 o):
        '''determine if a point is numerically close to another'''
        return self.isnear_c(o)

    # is vec3 o within a very small neighborhood of self in the xy plane
    cpdef bint isnearxy(self,vec3 o):
        '''determine if a point is numerically close to another in the xy plane'''
        return self.isnearxy_c(o)

    # return the squared magintude of self
    cpdef float mag2(self):
        '''compute the squared magnitude of this point'''
        return self.mag2_c()

    # return the magintude of self
    cpdef float mag(self):
        '''compute the magnitude of this point'''
        return self.mag_c()

    # normalize and return self
    cpdef vec3 nrm(self):
        '''normalize this vector'''
        return self.nrm_c()

    # translate self by vec3 o
    cpdef vec3 trn(self,vec3 o):
        '''translate this point by a vector'''
        return self.trn_c(o)

    # translate self in the x direction by d
    cpdef vec3 xtrn(self,float d):
        '''translate this point in the x direction by a distance'''
        return self.xtrn_c(d)

    # translate self in the y direction by d
    cpdef vec3 ytrn(self,float d):
        '''translate this point in the y direction by a distance'''
        return self.ytrn_c(d)

    # translate self in the z direction by d
    cpdef vec3 ztrn(self,float d):
        '''translate this point in the z direction by a distance'''
        return self.ztrn_c(d)

    # multiply each component by a scalar of and return self
    cpdef vec3 scl(self,float s):
        '''multiply components of this point by a scalar'''
        return self.scl_c(s)

    # scale self.x by scalar s
    cpdef vec3 xscl(self,float s):
        '''multiply x component of this point by a scalar'''
        return self.xscl_c(s)

    # scale self.y by scalar s
    cpdef vec3 yscl(self,float s):
        '''multiply y component of this point by a scalar'''
        return self.yscl_c(s)

    # scale self.z by scalar s
    cpdef vec3 zscl(self,float s):
        '''multiply z component of this point by a scalar'''
        return self.zscl_c(s)

    # rotate by a quaternion q and return self
    cpdef vec3 rot(self,quat q):
        '''rotate by a quaternion and return self'''
        return self.rot_c(q)

    # rotate around the x axis by an angle a and return self
    cpdef vec3 xrot(self,float a):
        '''rotate around the x axis by an angle and return self'''
        return self.xrot_c(a)

    # rotate around the y axis by an angle a and return self
    # NOTE: this is slower than using rot with the appropriate quaternion
    cpdef vec3 yrot(self,float a):
        '''rotate around the y axis by an angle and return self'''
        return self.yrot_c(a)

    # rotate around the z axis by an angle a and return self
    # NOTE: this is slower than using rot with the appropriate quaternion
    cpdef vec3 zrot(self,float a):
        '''rotate around the z axis by an angle and return self'''
        return self.zrot_c(a)

    # flip the direction of and return self
    cpdef vec3 flp(self):
        '''multiply components of this point by -1.0'''
        return self.flp_c()

    # return a vector point from self to vec3 o
    cpdef vec3 tov(self,vec3 o):
        '''return a vector which points from this point to another'''
        return self.tov_c(o)

    # return a vector at the midpoint between self and vec3 o
    cpdef vec3 mid(self,vec3 o):
        '''return the midpoint between this point and another'''
        return self.mid_c(o)

    # linearly interpolate between self and vec3 o proportionally to ds
    cpdef vec3 lerp(self,vec3 o,float ds):
        '''create a new point linearly interpolated between this point and another'''
        return self.lerp_c(o,ds)

    ###########################################################################

###########################################################################

# functions to quickly generate R3 basis vectors and their flips
cdef vec3 x_c():return vec3(1,0,0)
cdef vec3 y_c():return vec3(0,1,0)
cdef vec3 z_c():return vec3(0,0,1)
cdef vec3 nx_c():return vec3(-1,0,0)
cdef vec3 ny_c():return vec3(0,-1,0)
cdef vec3 nz_c():return vec3(0,0,-1)

# get a copy of xhat
cpdef vec3 x():
    '''create an x axis basis vector'''
    return x_c()

# get a copy of yhat
cpdef vec3 y():
    '''create an y axis basis vector'''
    return y_c()

# get a copy of zhat
cpdef vec3 z():
    '''create an z axis basis vector'''
    return z_c()

# get a copy of negative xhat
cpdef vec3 nx():
    '''create an negative x axis basis vector'''
    return nx_c()

# get a copy of negative yhat
cpdef vec3 ny():
    '''create an negative y axis basis vector'''
    return ny_c()

# get a copy of negative zhat
cpdef vec3 nz():
    '''create an negative z axis basis vector'''
    return nz_c()

###########################################################################

    








