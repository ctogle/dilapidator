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



__doc__ = '''dilapidator\'s implementation of a curve in R3'''
# dilapidators implementation of a curve in R3
cdef class curve:

    ###########################################################################
    ### python object overrides ###############################################
    ###########################################################################

    def __str__(self):
        return 'curve:\n\t'+'\n\t'.join([x.__str__() for x in self.pts])
    #def __iter__(self):yield self.x;yield self.y;yield self.z
    #def __add__(self,o):return vec3(self.x+o.x,self.y+o.y,self.z+o.z)
    #def __sub__(self,o):return vec3(self.x-o.x,self.y-o.y,self.z-o.z)
    #def __mul__(self,o):return self.cp().scl(o)
    def __is_equal(self,o):
        return self.tail.isnear(o.tail) and self.tip.isnear(o.tip)
    def __richcmp__(x,y,op):
        if op == 2:return x.__is_equal(y)
        else:assert False

    ###########################################################################

    ###########################################################################
    ### c methods #############################################################
    ###########################################################################

    def __cinit__(self,vec3 tail,vec3 tip,int segs = 1,int meth = 0):
        self.meth = meth
        self.segs = segs
        self.tail = tail
        self.tip = tip
        self.pts = []
        self.nms = []
        self.tns = []

    # return an independent copy of this curve
    cdef curve cp_c(self):
        cdef curve n = curve(self.tail.cp_c(),self.tip.cp_c())
        return n

    # return an independent copy of the xy projection of this curve
    cdef curve cpxy_c(self):
        cdef curve n = curve(self.tail.cpxy_c(),self.tip.cpxy_c())
        return n

    # erase the current curve data so calc is required
    cdef curve cln_c(self):
        self.pts = []
        self.nms = []
        self.tns = []

    # add the position/normal/tangent of 
    # an edge between vec3 p1 and vec3 p2
    cdef void calcone_c(self,vec3 p1,vec3 p2):
        tn = p1.tov_c(p2).nrm_c()
        nm = p1.tov_c(p2).nrm_c()
        self.pts.append(p1)
        self.tns.append(tn)
        self.nms.append(nm)

    # calculate positions/normals/tangents of the polyline of the curve
    # there should be self.segs edge segments between self.tail and self.tip
    # this requires the calculation of n-1 new points/normals/tangents
    # NOTE: tail and tip are not included in the stored results
    # NOTE: tangents at tail and tip should not be calculated here
    cdef curve calc_c(self):
        cdef list npts
        if self.meth == 0:npts = self.tail.pline(self.tip,self.segs-1)
        else:raise ValueError
        cdef vec3 p1,p2,tn,nm
        cdef int npcnt = len(npts)
        cdef int x
        self.cln_c()
        for x in range(1,npcnt):
            p1,p2 = npts[x-1],npts[x]
            self.calcone_c(p1,p2)
            #tn = p1.tov_c(p2).nrm_c()
            #nm = p1.tov_c(p2).nrm_c()
            #self.pts.append(p1)
            #self.tns.append(tn)
            #self.nms.append(nm)
        self.calcone_c(npts[npcnt-1],self.tip)
        #p1 = npts[npcnt-1]
        #tn = p1.tov_c(self.tip).nrm_c()
        #nm = p1.tov_c(self.tip).nrm_c()
        #self.pts.append(p1)
        #self.tns.append(tn)
        #self.nms.append(nm)
        return self

    ###########################################################################

    ###########################################################################
    ### python wrappers for c methods #########################################
    ###########################################################################

    # return an independent copy of this curve
    cpdef curve cp(self):
        '''create an independent copy of this curve'''
        return self.cp_c()

    # return an independent copy of the xy projection of this curve
    cpdef curve cpxy(self):
        '''create an independent copy of the xy projection of this curve'''
        return self.cpxy_c()

    # erase the current curve data so calc is required
    cpdef curve cln(self):
        '''empty the curve data in order to recalculate it'''
        return self.cln_c()

    # add the position/normal/tangent of 
    # an edge between vec3 p1 and vec3 p2
    cpdef calcone(self,vec3 p1,vec3 p2):
        '''add the point/normal/tangent associated with an edge of the curve'''
        self.calcone_c(p1,p2)

    # calculate positions/normals/tangents of the polyline of the curve
    # there should be self.segs edge segments between self.tail and self.tip
    # this requires the calculation of n-1 new points/normals/tangents
    # NOTE: tail and tip are not included in the stored results
    # NOTE: tangents at tail and tip should not be calculated here
    cpdef curve calc(self):
        '''recalculate the positions/normals/tangents of the curve'''
        return self.calc_c()

    ###########################################################################

    








