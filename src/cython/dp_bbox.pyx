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
        self._edge_data()

    # establishes data about the xy projection of this bbox
    def _edge_data(self):
        cs = [dpv.vector(self.x.x,self.y.x,0),dpv.vector(self.x.y,self.y.x,0),
            dpv.vector(self.x.y,self.y.y,0),dpv.vector(self.x.x,self.y.y,0)]
        self.corners = cs
        self.edgenorms = dpv.edge_normals_xy(self.corners)
        self.edgecount = len(self.edgenorms)
        self.center = dpv.com(self.corners)
        cvs = [dpv.v1_v2_c(self.center,v) for v in self.corners]
        cdists = [c.magnitude() for c in cvs]
        self.radius = max(cdists)

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

    cpdef bint point_inside(self,dpv.vector point):
        if not p_in_rng(point.x,self.x.x,self.x.y):return 0
        if not p_in_rng(point.y,self.y.x,self.y.y):return 0
        if not p_in_rng(point.z,self.z.x,self.z.y):return 0
        return 1

    cpdef bint intersect_tri(self,list tri):
        for p in tri:
            if self.point_inside(p):
                return 1
        return 0

cpdef list intersect_tri_filter(bbox bb,list tris,list tpts):
    cdef list isected = []
    cdef int tcnt = len(tris)
    cdef int tdx
    for tdx in range(tcnt):
        tri = [tpts[tris[tdx][x]] for x in range(3)]
        if bb.intersect_tri(tri):
            isected.append(tdx)
    return isected

cdef bint p_in_rng(float p,float x,float y):
    if p < x:return 0
    if p > y:return 0
    return 1

cdef bbox zero_c():
    cdef dpv.vector2d z = dpv.zero2d()
    cdef bbox new = bbox(z.copy(),z.copy(),z.copy())
    return new

cpdef bbox zero():
    return zero_c()

cdef bint overlap_c(dpv.vector2d rng1,dpv.vector2d rng2):
    if   rng1.y < rng2.x:return 0
    elif rng2.y < rng1.x:return 0
    else:return 1

# return 1 if a separating axis IS found
# return 0 if none are found
cdef bint separating_axis_c(bbox bb1,bbox bb2):
    cdef list ns1 = bb1.edgenorms
    cdef list ns2 = bb2.edgenorms
    cdef int egcnt1 = bb1.edgecount
    cdef int egcnt2 = bb2.edgecount
    cdef dpv.vector edgenorm
    cdef int egdx
    cdef dpv.vector2d proj1
    cdef dpv.vector2d proj2
    for egdx in range(egcnt1):
        edgenorm = <dpv.vector>ns1[egdx]
        proj1 = <dpv.vector2d>dpv.project_coords_c(bb1.corners,edgenorm)
        proj2 = <dpv.vector2d>dpv.project_coords_c(bb2.corners,edgenorm)
        if not <bint>overlap_c(proj1,proj2):
            return 1
    for egdx in range(egcnt2):
        edgenorm = <dpv.vector>ns2[egdx]
        proj1 = <dpv.vector2d>dpv.project_coords_c(bb1.corners,edgenorm)
        proj2 = <dpv.vector2d>dpv.project_coords_c(bb2.corners,edgenorm)
        if not <bint>overlap_c(proj1,proj2):
            return 1
    return 0

cpdef bint separating_axis(bbox bb1,bbox bb2):
    return separating_axis_c(bb1,bb2)


