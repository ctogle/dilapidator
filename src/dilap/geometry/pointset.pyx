# cython: profile=True
#cimport cython

#cimport dilap.core.tools as dpr

from dilap.geometry.vec3 cimport vec3
from dilap.geometry.quat cimport quat

stuff = 'hi'



__doc__ = '''dilapidator\'s implementation of a pointset'''
# dilapidators implementation of a pointset
cdef class pointset:

    ###########################################################################
    ### python object overrides ###############################################
    ###########################################################################

    def __str__(self):return 'pointset of size:'+str(self.pcnt)
    def __iter__(self):return self.ps.__iter__()

    #def __add__(self,o):return vec3(self.x+o.x,self.y+o.y,self.z+o.z)
    #def __sub__(self,o):return vec3(self.x-o.x,self.y-o.y,self.z-o.z)
    #def __mul__(self,o):return vec3(self.x*o.x,self.y*o.y,self.z*o.z)
    #def __is_equal(self,o):return self.isnear(o)
    #def __richcmp__(x,y,op):
    #    if op == 2:return x.__is_equal(y)
    #    else:assert False

    ###########################################################################

    ###########################################################################
    ### c methods #############################################################
    ###########################################################################

    def __cinit__(self):
        self.ps = []
        self.pcnt = 0

    # return a list of copies of each index specified point
    cdef list gpscp_c(self,rng):
        if rng is None:rng = range(self.pcnt)
        return [self.ps[x].cp() for x in rng]

    # return a list of the index specified points
    cdef list gps_c(self,rng):
        if rng is None:rng = range(self.pcnt)
        return [self.ps[x] for x in rng]

    # add a point and return its index
    cdef int ap_c(self,np):
        cdef int px = self.pcnt
        self.ps.append(np)
        self.pcnt += 1
        return px

    # add a sequence of points and return their indices
    cdef list aps_c(self,list nps):
        cdef int npcnt = len(nps)
        cdef int npx
        cdef list npxs = []
        for npx in range(npcnt):
            self.ps.append(nps[npx])
            npxs.append(self.pcnt)
            self.pcnt += 1
        return npxs

    # add a point to the pointset or find a point which 
    # is already present and near to the new point and return
    # the associated index
    cdef int np_c(self,np):
        x = self.fp_c(np)
        if x == -1:return self.ap_c(np)
        else:return x

    # for each point in a list of points, do what np does and 
    # return a list of the resulting indices
    cdef list nps_c(self,list nps):
        cdef int npcnt = len(nps)
        cdef int nx
        cdef list npxs = []
        for nx in range(npcnt):
            np = nps[nx]
            x = self.fp_c(np)
            if x == -1:npxs.append(self.ap_c(np))
            else:npxs.append(x)
        return npxs

    # given a point, return the index of a point that is near it
    # or return -1 if no such point exists in the pointset
    cdef int fp_c(self,p):
        cdef int px
        for px in range(self.pcnt):
            if self.ps[px].isnear(p):
                return px
        return -1

    # given a sequence of points, return whatever 
    # find_point would in a one to one sequence
    ##### NOTE: could be faster
    cdef list fps_c(self,list ps):
        return [self.fp_c(p) for p in ps]
        
    # is self disjoint from pointset o
    cdef bint disjoint_c(self,pointset o):
        cdef int x,y
        for x in range(self.pcnt):
            p = self.ps[x]
            for y in range(o.pcnt):
                q = o.ps[y]
                if p.isnear(q):return 0
        return 1

    # translate all points by vec3 v
    cdef pointset trn_c(self,vec3 v):
        cdef int px
        for px in range(self.pcnt):
            #self.ps[px].trn_c(v)
            self.ps[px].trn(v)
        return self

    # rotate all points by quat q
    cdef pointset rot_c(self,quat q):
        cdef int px
        for px in range(self.pcnt):
            #self.ps[px].rot_c(q)
            self.ps[px].rot(q)
        return self

    # scale all points by vec3 0
    cdef pointset scl_c(self,vec3 o):
        cdef int px
        for px in range(self.pcnt):
            self.ps[px].scl(o)
        return self

    # scale all points by float f
    cdef pointset uscl_c(self,float f):
        cdef int px
        for px in range(self.pcnt):
            self.ps[px].uscl(f)
        return self

    ###########################################################################

    ###########################################################################
    ### python wrappers for c methods #########################################
    ###########################################################################

    # return a list of copies of each index specified point
    cpdef list gpscp(self,rng):
        '''create a list of copies of every point in this pointset'''
        return self.gpscp_c(rng)

    # return a list of the index specified points
    cpdef list gps(self,rng):
        '''create a list of every point in this pointset'''
        return self.gps_c(rng)

    # add a point and return its index
    cpdef int ap(self,np):
        '''add a point to this pointset'''
        return self.ap_c(np)

    # add a sequence of points and return their indices
    cpdef list aps(self,list nps):
        '''add points to this pointset'''
        return self.aps_c(nps)

    # add a point to the pointset or find a point which 
    # is already present and near to the new point and return
    # the associated index
    cpdef int np(self,np):
        '''add point or return exist point in close proximity'''
        return self.np_c(np)

    # for each point in a list of points, do what np does and 
    # return a list of the resulting indices
    cpdef list nps(self,list nps):
        '''perform np on a list of points'''
        return self.nps_c(nps)

    # given a point, return the index of a point that is near it
    # or return -1 if no such point exists in the pointset
    cpdef int fp(self,p):
        '''return the index of a point in this pointset near another point,
        or return -1 if no such point exists'''
        return self.fp_c(p)

    # given a sequence of points, return whatever 
    # find_point would in a one to one sequence
    ##### NOTE: could be faster
    cpdef list fps(self,list ps):
        '''find a set of points'''
        return self.fps_c(ps)
        
    # is self disjoint from pointset o
    cpdef bint disjoint(self,pointset o):
        '''determine if this pointset is disjoint from another'''
        return self.disjoint_c(o)

    # translate all points by vec3 v
    cpdef pointset trn(self,vec3 v):
        '''translate the points in this pointset'''
        return self.trn_c(v)

    # rotate all points by quat q
    cpdef pointset rot(self,quat q):
        '''rotate the points in this pointset'''
        return self.rot_c(q)

    # scale all points by vec3 0
    cpdef pointset scl(self,vec3 o):
        '''scale (x,y,z scale) the points in this pointset'''
        return self.scl_c(o)

    # scale all points by float f
    cpdef pointset uscl(self,float f):
        '''scale (uniform scale) the points in this pointset'''
        return self.uscl_c(f)

    ###########################################################################




        



 

