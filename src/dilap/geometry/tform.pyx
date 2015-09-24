# cython: profile=True
#cimport cython

cimport dilap.core.tools as dpr

from dilap.geometry.vec3 cimport vec3
from dilap.geometry.quat cimport quat

stuff = 'hi'










__doc__ = '''dilapidator\'s implementation of a transform object in R3'''
# dilapidator\'s implementation of a transform object in R3'''
cdef class tform:

    ###########################################################################
    ### python object overrides ###############################################
    ###########################################################################

    #def __repr__(self):return self.__str__()
    def __str__(self):
        strr  = 'tform:'
        strr += '\n\t'+'pos:'+self.pos.__str__()
        strr += '\n\t'+'rot:'+self.rot.__str__()
        strr += '\n\t'+'scl:'+self.scl.__str__()
        return strr

    def __is_equal(self,o):
        if self.pos.isnear(o.pos):
            if self.rot.isnear(o.rot):
                if self.scl.isnear(o.scl):
                    return 1
        return 0
    def __richcmp__(x,y,op):
        if op == 2:return x.__is_equal(y)
        else:assert False

    ###########################################################################

    ###########################################################################
    ### c methods #############################################################
    ###########################################################################

    def __cinit__(self,pos,rot,scl):
        self.pos = pos
        self.rot = rot
        self.scl = scl

    # return an independent copy of this tform
    cdef tform cp_c(self):
        cdef tform n = tform(self.pos.cp_c(),self.rot.cp_c(),self.scl.cp_c())
        return n

    # should return a tform to world space, possibly given a parent 
    # tform from which to inherit position,rotation,scale
    # NOTE: the parent should be in world space already
    cdef tform true_c(self,parent):
        cdef vec3 np = self.pos.cp_c()
        cdef quat nr = self.rot.cp_c()
        cdef vec3 ns = self.scl.cp_c()
        if parent:
            np.rot_c(parent.rot)
            np.trn_c(parent.pos)
            #nr.rot_c(parent.rot)
            nr = parent.rot*nr
            ns.scl_c(parent.scl)
        cdef tform n = tform(np,nr,ns)
        return n

    ###########################################################################

    ###########################################################################
    ### python wrappers for c methods #########################################
    ###########################################################################

    # return an independent copy of this tform
    cpdef tform cp(self):
        '''create an independent copy of this tform'''
        return self.cp_c()

    # should return a tform to world space, possibly given a parent 
    # tform from which to inherit position,rotation,scale
    # NOTE: the parent should be in world space already
    cpdef tform true(self,parent):
        '''return a new tform representing this tform in world space'''
        return self.true_c(parent)

    ###########################################################################










