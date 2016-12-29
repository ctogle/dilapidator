# cython: profile=True
#cimport cython

#cimport dilap.core.tools as dpr

import dilap.geometry.tools as gtl
cimport dilap.geometry.tools as gtl

from dilap.geometry.quat cimport quat

from libc.math cimport sqrt
from libc.math cimport cos
from libc.math cimport sin
from libc.math cimport tan
from libc.math cimport hypot

cimport numpy
import numpy
 
import copyreg

stuff = 'hi'



__doc__ = '''dilapidator\'s implementation of a vector in R3'''



# dilapidators implementation of a vector in R3
cdef class vec3:

    ###########################################################################
    ### python object overrides ###############################################
    ###########################################################################

    def __str__(self):return 'vec3:'+str(tuple(self))
    def __repr__(self):return 'vec3:'+str(tuple(self))
    def __iter__(self):yield self.x;yield self.y;yield self.z
    def __add__(self,o):return vec3(self.x+o.x,self.y+o.y,self.z+o.z)
    def __sub__(self,o):return vec3(self.x-o.x,self.y-o.y,self.z-o.z)
    def __mul__(self,o):return self.cp().scl(o)
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

    # return an independent reciprocated copy of this point
    cdef vec3 cpr_c(self):
        cdef vec3 n = vec3(1.0/self.x,1.0/self.y,1.0/self.z)
        return n

    # return an independent flipped copy of this point
    cdef vec3 cpf_c(self):
        cdef vec3 n = vec3(-self.x,-self.y,-self.z)
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
        if sm*om == 0:
            print('ang_c asserting fake zero')
            return 0.0
        cdef float sod = (self.x*o.x + self.y*o.y + self.z*o.z)/(sm*om)
        cdef float a
        if   gtl.isnear_c(sod, 1.0):a = 0.0
        elif gtl.isnear_c(sod,-1.0):a = numpy.pi
        else:a = numpy.arccos(sod)
        return a

    # return the signed angle between self and vec3 o given z-up vec3 n
    cdef float sang_c(self,vec3 o,vec3 n):
        cdef vec3 n1 = self.cp_c().nrm_c()
        cdef vec3 n2 = o.cp_c().nrm_c()
        cdef vec3 vn = n1.crs_c(n2)
        cdef float vd = vn.x*n.x + vn.y*n.y + vn.z*n.z
        cdef float sod = n1.x*n2.x + n1.y*n2.y + n1.z*n2.z
        cdef float a = 0.0
        if   gtl.isnear_c(sod, 1.0):pass
        elif gtl.isnear_c(sod,-1.0):a = numpy.pi
        else:a = numpy.arccos(sod)
        if vd < 0.0:a *= -1.0
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
        cdef vec3 pj = n.cp_c().uscl_c(-d)
        return self.trn_c(pj)

    # project a set of vectors onto this one, returning the min/max 
    cdef tuple prjps_c(self,ps):
        cdef pcnt = len(ps)
        cdef int px
        cdef float pj = ps[0].dot(self)
        cdef float pmin = pj
        cdef float pmax = pj
        for px in range(1,pcnt):
            pj = ps[px].dot(self)
            if pj < pmin:pmin = pj
            if pj > pmax:pmax = pj
        return pmin,pmax

    # return the u,v barycentric coordinates of self given 3 corners of a triangle
    # everything is assumed to be in the xy plane
    cdef tuple baryxy_c(self,vec3 a,vec3 b,vec3 c): 
        cdef float v0x,v0y,v1x,v1y,v2x,v2y
        cdef float dot00,dot01,dot02,dot11,dot12,u,v,invdenom
        v0x,v0y = c.x-a.x,c.y-a.y
        v1x,v1y = b.x-a.x,b.y-a.y
        v2x,v2y = self.x-a.x,self.y-a.y
        dot00 = v0x*v0x + v0y*v0y
        dot01 = v0x*v1x + v0y*v1y
        dot02 = v0x*v2x + v0y*v2y
        dot11 = v1x*v1x + v1y*v1y
        dot12 = v1x*v2x + v1y*v2y
        denom = (dot00 * dot11 - dot01 * dot01)
        if denom == 0:
            print('colinear triangle?',a,b,c,self,denom)
            #import dilap.core.plotting as dtl
            #import matplotlib.pyplot as plt
            #ax = dtl.plot_axes_xy(100)
            #ax = dtl.plot_polygon_xy((a,b,c),ax,lw = 3,col = 'b')
            #ax = dtl.plot_point_xy(self,ax,mk = 's')
            #plt.show()
            u,v = self.baryxy_c(c,a,b)
            print('not colinear afterall',u,v)
        else:
            invdenom = 1.0 / denom
            u = (dot11 * dot02 - dot01 * dot12) * invdenom
            v = (dot00 * dot12 - dot01 * dot02) * invdenom
        return u,v

    # is vec3 o within an open ball of raidus e centered at self
    cdef bint inneighborhood_c(self,vec3 o,float e):
        cdef float d = self.d_c(o)
        #if gtl.near(d,e) < e:return 1
        if d < e:return 1
        else:return 0

    # is vec3 o within a very small neighborhood of self
    cdef bint isnear_c(self,vec3 o):
        cdef float dx = (self.x-o.x)
        if dx*dx > gtl.epsilonsq_c:return 0
        cdef float dy = (self.y-o.y)
        if dy*dy > gtl.epsilonsq_c:return 0
        cdef float dz = (self.z-o.z)
        if dz*dz > gtl.epsilonsq_c:return 0
        return 1

    # is vec3 o within a very small neighborhood of self in the xy plane
    cdef bint isnearxy_c(self,vec3 o):
        cdef float dx = (self.x-o.x)
        if dx*dx > gtl.epsilonsq_c:return 0
        cdef float dy = (self.y-o.y)
        if dy*dy > gtl.epsilonsq_c:return 0
        return 1

    # is self within the edges between any two adjacent points in ps
    cdef bint inbxy_c(self,ps):
        cdef int wn = 0
        cdef int px
        cdef int pcnt = len(ps)
        cdef float read
        for px in range(pcnt):
            if self.onsxy_c(ps[px-1],ps[px],1):return 0
            read = gtl.orient2d_c(self,ps[px-1],ps[px])
            if ps[px-1].y <= self.y:
                if ps[px].y > self.y:
                    if read > 0:wn += 1
            else:
                if ps[px].y <= self.y:
                    if read < 0:wn -= 1
        #if wn > 0:return 1
        if not wn == 0:return 1
        else:return 0

    # determine if the point pt is inside the triangle abc
    # assume all points are in the xy plane 
    cdef bint intrixy_c(self,vec3 a,vec3 b,vec3 c):
        cdef float u,v
        if gtl.orient2d_c(a,b,c) == 0:
            return 0
            '''#
            ab,bc,ca = a.d(b),b.d(c),c.d(a)
            if ab > bc and ab > ca:
                if self.d_c(a.mid_c(b)) < ab/2.0:return 1
                else:return 0
            elif bc > ca and bc > ab:
                if self.d_c(b.mid_c(c)) < bc/2.0:return 1
                else:return 0
            elif ca > ab and ca > bc:
                if self.d_c(c.mid_c(a)) < ca/2.0:return 1
                else:return 0
            '''#
        else:
            u,v = self.baryxy_c(a,b,c)
            if u > 0 or abs(u) < gtl.epsilon_c:
                if v > 0 or abs(v) < gtl.epsilon_c:
                    if 1-u-v > 0 or abs(1-u-v) < gtl.epsilon_c:
                        return 1
            return 0

    # is self on an edge between any two points
    cdef bint onsxy_c(self,vec3 s1,vec3 s2,bint ie = 0):
        if not gtl.orient2d_c(self,s1,s2) == 0:return 0
        if self.isnear(s1) or self.isnear(s2):
            if ie:return 1
            else:return 0
        #if (s1.x != s2.x): # S is not  vertical
        if not gtl.isnear_c(s1.x,s2.x): # S is not  vertical
            if (s1.x <= self.x and self.x <= s2.x):return 1
            if (s1.x >= self.x and self.x >= s2.x):return 1
        else: # S is vertical, so test y  coordinate
            if (s1.y <= self.y and self.y <= s2.y):return 1
            if (s1.y >= self.y and self.y >= s2.y):return 1
        return 0

    # is self on an edge between any two adjacent points in ps
    cdef bint onbxy_c(self,ps):
        cdef int px
        cdef int pcnt = len(ps)
        cdef vec3 p1,p2
        for px in range(pcnt):
            p1,p2 = ps[px-1],ps[px]
            if p1.isnear_c(self):return 1
            if self.onsxy_c(p1,p2,1):return 1
            #if gtl.inseg_xy_c(self,p1,p2):return 1
        return 0

    # is self on the boundary or any holes of a polygon
    cdef bint onpxy_c(self,py):
        eb,ibs = py
        if self.onbxy(eb):return 1
        for ib in ibs:
            if self.onbxy(ib):return 1
        return 0

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
        else:return self.uscl_c(1.0/m)

    # translate self by vec3 o
    cdef vec3 trn_c(self,vec3 o):
        self.x += o.x
        self.y += o.y
        self.z += o.z
        return self

    # translate a set of vec3s by self
    cdef vec3 trnos_c(self,os):
        cdef int ocnt = len(os)
        cdef int ox
        cdef vec3 o
        for ox in range(ocnt):
            o = os[ox]
            o.trn(self)
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

    # 1-1 multiplication by vec3 o
    cdef vec3 scl_c(self,vec3 o):
        self.x *= o.x;self.y *= o.y;self.z *= o.z
        return self

    # multiply each component by a scalar of and return self
    cdef vec3 uscl_c(self,float s):
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

    # rotate a set of points around self
    cdef vec3 fulc_c(self,quat q,pts):
        cdef int pcnt = len(pts)
        cdef int px
        cdef vec3 pt
        for px in range(pcnt):
            pt = pts[px]
            pt.trn_c(self.flp_c()).rot_c(q).trn_c(self.flp_c())
        return self

    # return a clockwise ordered version of a set of vectors
    # fst : a vector whose angle is considered to be minimized
    cdef list cwoxy_c(self,vec3 fst,ps):
        cdef vec3 ftn = self.tov_c(fst)
        cdef vec3 o,otn
        cdef quat q
        cdef list os = []
        cdef int pcnt = len(ps)
        cdef int px
        for px in range(pcnt):
            o = ps[px]
            otn = self.tov_c(o)
            q = quat(0,0,0,0).uu(ftn,otn)
            print('cwoxy',o,q)
            os.append(o)
        return os

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
        return self.uscl_c(-1.0)

    # return a vector point from self to vec3 o
    cdef vec3 tov_c(self,vec3 o):
        cdef float dx = o.x - self.x
        cdef float dy = o.y - self.y
        cdef float dz = o.z - self.z
        cdef vec3 n = vec3(dx,dy,dz)
        return n

    # return a vector point from self to vec3 o in the xy plane
    cdef vec3 tovxy_c(self,vec3 o):
        cdef float dx = o.x - self.x
        cdef float dy = o.y - self.y
        cdef vec3 n = vec3(dx,dy,0)
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

    # generate a polyline between self and vec3 o with n points between ends
    # self and o are not modified nor contained in the result
    cdef list pline_c(self,vec3 o,int n):
        #cdef float s = self.d_c(o)
        cdef float t
        cdef int x
        cdef list line = []
        for x in range(n):
            t = (x+1.0)/(n+1.0)
            line.append(self.lerp_c(o,t))
            #print('t',t,self,line[-1],o)
        return line

    # return the n points appearing on a spline between self and o
    #   st is the tangent vector as self
    #   ot is the tangent vector at o (the end point)
    cdef list spline_c(self,vec3 o,vec3 st,vec3 ot,int n):
        n += 1
        cdef vec3 p2 = self.cp().trn(st)
        cdef vec3 p3 = o.cp().trn(ot)
        cdef list cox = [self.x,p2.x,p3.x,o.x]
        cdef list coy = [self.y,p2.y,p3.y,o.y]
        cdef list coz = [self.z,p2.z,p3.z,o.z]
        cdef list tim = [0.0,1.0,2.0,3.0]
        cdef float alpha = 1.0/2.0

        self.ptime_c([self,p2,p3,o],tim,alpha)

        cox = self.catmull_rom_c(cox,tim,n)[1:-1]
        coy = self.catmull_rom_c(coy,tim,n)[1:-1] 
        coz = self.catmull_rom_c(coz,tim,n)[1:-1] 

        cdef list spline = [vec3(*i) for i in zip(cox,coy,coz)]
        return spline

    cdef list ptime_c(self,list ps,list time,float alpha):
        cdef float total = 0.0
        cdef int idx
        cdef list v1v2
        for idx in range(1,4):
            total += ps[idx-1].tov(ps[idx]).mag2()**(alpha)
            time[idx] = total

    # FIXXXXX
    cdef list catmull_rom_c(self,list P,list T,int tcnt):
        cdef int j
        cdef int t
        cdef float tt
        cdef float p
        cdef list spl = P[:1]
        for j in range(1, len(P)-2):  # skip the ends
            for t in range(tcnt):  # t: 0 .1 .2 .. .9
                tt = float(t)/tcnt
                tt = T[1] + tt*(T[2]-T[1])
                p = self.spline__FIXME_c(
                    [P[j-1], P[j], P[j+1], P[j+2]],
                    [T[j-1], T[j], T[j+1], T[j+2]],tt)
                spl.append(p)
        spl.extend(P[-2:])
        return spl

    # FIXXXXX
    cdef float spline__FIXME_c(self,list p,list time,float t):
        L01 = p[0] * (time[1] - t) / (time[1] - time[0]) + p[1] * (t - time[0]) / (time[1] - time[0])
        L12 = p[1] * (time[2] - t) / (time[2] - time[1]) + p[2] * (t - time[1]) / (time[2] - time[1])
        L23 = p[2] * (time[3] - t) / (time[3] - time[2]) + p[3] * (t - time[2]) / (time[3] - time[2])
        L012 = L01 * (time[2] - t) / (time[2] - time[0]) + L12 * (t - time[0]) / (time[2] - time[0])
        L123 = L12 * (time[3] - t) / (time[3] - time[1]) + L23 * (t - time[1]) / (time[3] - time[1])
        C12 = L012 * (time[2] - t) / (time[2] - time[1]) + L123 * (t - time[1]) / (time[2] - time[1])
        return C12





    # return a ring of points of radius r with n corners
    cdef list pring_c(self,float r,int n):
        cdef float alpha = numpy.pi*(2.0/n)
        cdef float sr = r/cos(alpha/2.0)
        cdef vec3 st = vec3(0,0,0).xtrn_c(sr)
        cdef vec3 nv
        cdef list points = []
        cdef int x
        for x in range(n):
            nv = st.cp_c().zrot_c(x*alpha-alpha/2.0)
            points.append(nv.trn_c(self))
        return points

    # return a rectangle of dims l by w, centered at self
    cdef list sq_c(self,float l,float w):
        cdef float hl = l/2.0
        cdef float hw = w/2.0
        cdef list sq = [
            vec3(self.x-hl,self.y-hw,self.z),
            vec3(self.x+hl,self.y-hw,self.z),
            vec3(self.x+hl,self.y+hw,self.z),
            vec3(self.x-hl,self.y+hw,self.z)]
        return sq

    # translate self to the center of mass of a set of vectors
    cdef vec3 com_c(self,os):
        cdef int pcnt = len(os)
        cdef int px
        cdef float pcntf = float(pcnt)
        cdef vec3 n = vec3(0,0,0)
        for px in range(pcnt):n.trn_c(os[px])
        self.trn_c(n.uscl_c(1.0/pcntf))
        return self
        #return n.uscl_c(1.0/pcntf)

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

    # return an independent copy of the reciprocal of this point
    cpdef vec3 cpr(self):
        '''return an independent copy of the reciprocal of this point'''
        return self.cpr_c()

    # return an independent flipped copy of this point
    cpdef vec3 cpf(self):
        '''return an independent flipped copy of this point'''
        return self.cpf_c()

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

    # return the signed angle between self and vec3 o given z-up vec3 n
    cpdef float sang(self,vec3 o,vec3 n):
        '''determine the signed angle between this point and another'''
        return self.sang_c(o,n)

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

    # project a set of vectors onto this one, returning the min/max 
    cpdef tuple prjps(self,ps):
        '''project a set of vectors onto this one, returning the min/max'''
        return self.prjps_c(ps)

    # return the u,v barycentric coordinates of self given 3 corners of a triangle
    # everything is assumed to be in the xy plane
    cpdef tuple baryxy(self,vec3 a,vec3 b,vec3 c): 
        '''find the barycentric coords of this point in a triangle in the xy plane'''
        return self.baryxy_c(a,b,c)

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

    # is self within the edges between any two adjacent points in ps
    cpdef bint inbxy(self,ps):
        '''determine if self is within the edges between any two adjacent points in ps'''
        return self.inbxy_c(ps)

    # determine if the point pt is inside the triangle abc
    # assume all points are in the xy plane 
    cpdef bint intrixy(self,vec3 a,vec3 b,vec3 c):
        '''
        determine if the point pt is inside the triangle abc
        assume all points are in the xy plane
        '''
        return self.intrixy_c(a,b,c)

    # is self on an edge between any two points
    cpdef bint onsxy(self,vec3 s1,vec3 s2,bint ie = 0):
        '''is self on an edge between any two points'''
        return self.onsxy_c(s1,s2,ie)

    # is self on an edge between any two adjacent points in ps
    cpdef bint onbxy(self,ps):
        '''determine if self on an edge between any two adjacent points in ps'''
        return self.onbxy_c(ps)

    # is self on the boundary or any holes of a polygon
    cpdef bint onpxy(self,py):
        '''determine if self on an edge between any two adjacent points in ps'''
        return self.onpxy_c(py)

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

    # translate a set of vec3s by self
    cpdef vec3 trnos(self,os):
        '''translate a set of vec3s by self'''
        return self.trnos_c(os)

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

    # 1-1 multiplication by vec3 o
    cpdef vec3 scl(self,vec3 o):
        '''1-1 multiplication by another vector'''
        return self.scl_c(o)

    # multiply each component by a scalar of and return self
    cpdef vec3 uscl(self,float s):
        '''multiply components of this point by a scalar'''
        return self.uscl_c(s)

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

    # rotate a set of points around self
    cpdef vec3 fulc(self,quat q,pts):
        '''rotate a set of points around self'''
        return self.fulc_c(q,pts)

    # return a clockwise ordered version of a set of vectors
    # fst : a vector whose angle is considered to be minimized
    cpdef list cwoxy(self,vec3 fst,ps):
        '''
        return a clockwise ordered version of a set of vectors
        fst : a vector whose angle is considered to be minimized
        '''
        return self.cwoxy_c(fst,ps)

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

    # return a vector point from self to vec3 o in the xy plane
    cpdef vec3 tovxy(self,vec3 o):
        '''return a vector point from self to vec3 o in the xy plane'''
        return self.tovxy_c(o)

    # return a vector at the midpoint between self and vec3 o
    cpdef vec3 mid(self,vec3 o):
        '''return the midpoint between this point and another'''
        return self.mid_c(o)

    # linearly interpolate between self and vec3 o proportionally to ds
    cpdef vec3 lerp(self,vec3 o,float ds):
        '''create a new point linearly interpolated between this point and another'''
        return self.lerp_c(o,ds)

    # generate a polyline between seelf and vec3 o with n points between ends
    # self and o are not modified nor contained in the result
    cpdef list pline(self,vec3 o,int n):
        '''create a polyline between this and another'''
        return self.pline_c(o,n)

    # return the n points appearing on a spline between self and o
    #   st is the tangent vector as self
    #   ot is the tangent vector at o (the end point)
    cpdef list spline(self,vec3 o,vec3 st,vec3 ot,int n):
        '''generate the n points appearing on a spline between self and o'''
        return self.spline_c(o,st,ot,n)

    # return a ring of points of radius r with n corners
    cpdef list pring(self,float r,int n):
        '''return a ring of points of radius r with n corners centered at self'''
        return self.pring_c(r,n)

    # return a rectangle of dims l by w, centered at self
    cpdef list sq(self,float l,float w):
        '''return a rectangle of dims l by w, centered at self'''
        return self.sq_c(l,w)

    cpdef vec3 com(self,os):
        '''compute the center of mass for a set of vectors'''
        return self.com_c(os)

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

def pickle_vec3(v):
    print("pickling a vec3 instance...")
    return vec3,(v.x,v.y,v.z)

copyreg.pickle(vec3,pickle_vec3)

###########################################################################

