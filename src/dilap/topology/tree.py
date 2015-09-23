import dilap.core.tools as dpr

import dilap.topology.vert as dvt
import dilap.topology.edge as deg

from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat










__doc__ = '''dilapidator\'s implementation of a topological tree'''
# dilapidators implementation of a topological tree
class tree:

    vertclass = dvt.vert
    edgeclass = deg.edge

    def __str__(self):return 'tree:'+str(self.vertcount)+str(edgecount)

    def __init__(self,*args,**kwargs):
        self.verts = []
        self.edges = []
        self.vertcount = 0
        self.edgecount = 0
        self.root = self.aroot(*args,**kwargs)

    # return the d-cells which are incident upon v and/or e
    def mask(self,d = 0,v = None,e = None):
        results = []
        if d == 0:
            if not v is None:
                for ve in self.mask(1,v,None):
                    for vv in self.mask(0,None,ve):
                        if not vv is v:results.append(vv)
            if not e is None:results.extend([e.one,e.two])
        elif d == 1:
            if not v is None:results.extend(v.ring)
            if not e is None:
                for vv in self.mask(0,None,e):
                    for ve in self.mask(1,vv,None):
                        if not ve is e:results.append(ve)
        return results

    # add a new root vertex in the tree (only needed once)
    def aroot(self,*args,**kwargs):
        vrt = self.vertclass(self.vertcount,*args,**kwargs)
        self.verts.append(vrt)
        self.vertcount += 1
        return vrt

    # create new vertex connected to par by edge ending at new vertex
    def avert(self,par,*args,**kwargs):
        vrt = self.vertclass(self.vertcount,*args,**kwargs)
        self.verts.append(vrt)
        self.vertcount += 1
        if par is None:par = self.root
        edg = self.aedge(par,vrt)
        return vrt

    # provided two vertices and possibly a geometric index
    # create a new edge from vrt1 to vrt2
    def aedge(self,vrt1,vrt2):
        edg = self.edgeclass(vrt1,vrt2,self.edgecount)
        self.edges.append(edg)
        self.edgecount += 1
        vrt1.connect(edg)
        vrt2.connect(edg)
        return edg

    # provided a vertex, return its parent which lies at the 
    # beginning of the one edge ending at vrt
    def above(self,vrt):
        es = self.mask(1,vrt,None)
        if es:return es[0].one

    # provided a vertex, return its children which lie at the 
    # ends of the edges starting at vrt
    def below(self,vrt):
        par = self.above(vrt)
        vs = self.mask(0,vrt,None)
        chn = [v for v in vs if not v is par]
        return chn






 



