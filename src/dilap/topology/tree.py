#import dilap.core.tools as dpr
import dilap.topology.vert as dvt
import dilap.topology.edge as deg

import dilap.geometry.tools as dpr
from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat

import pdb










__doc__ = '''dilapidator\'s implementation of a topological tree'''
# dilapidators implementation of a topological tree
class tree:

    vertclass = dvt.vert
    edgeclass = deg.edge

    def __str__(self):return 'tree:'+str(self.vertcount)+str(edgecount)
    def vcnt(self):return self.cnt(self.verts)
    def ecnt(self):return self.cnt(self.edges)
    def cnt(self,tobjs):
        x = 0
        for v in tobjs:
            if not v is None:
                x += 1
        return x

    def __init__(self,*args,**kwargs):
        self.verts = []
        self.edges = []
        self.vertcount = 0
        self.edgecount = 0
        self.vistack = []
        self.eistack = []
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

    # provided a vertex, return the edge to its parent 
    # which lies at the beginning of the one edge ending at vrt
    def ebove(self,vrt):
        es = self.mask(1,vrt,None)
        if es:return es[0]

    # provided a vertex, return its parent which lies at the 
    # beginning of the one edge ending at vrt
    def above(self,vrt):
        eb = self.ebove(vrt)
        if eb is None:return
        else:return eb.one

    # provided a vertex, return its children which lie at the 
    # ends of the edges starting at vrt
    def below(self,vrt):
        par = self.above(vrt)
        vs = self.mask(0,vrt,None)
        chn = [v for v in vs if not v is par]
        return chn

    # add a new root vertex in the tree (only needed once)
    def aroot(self,*args,**kwargs):
        vrt = self.vertclass(self.vertcount,*args,**kwargs)
        self.verts.append(vrt)
        self.vertcount += 1
        return vrt

    # create new vertex connected to par by edge ending at new vertex
    def avert(self,par,*args,**kwargs):
        vrt = self.vertclass(self.vertcount,*args,**kwargs)
        if self.vistack:
            vx = self.vistack.pop(0)
            self.verts[vx] = vrt
        else:
            self.verts.append(vrt)
            vx = self.vertcount
            self.vertcount += 1
        if par is None:par = self.root
        edg = self.aedge(par,vrt)
        return vrt

    # delete existing vertex
    # recursively delete all children as well
    def rvert(self,vrt,re = True):
        for b in self.below(vrt):self.rvert(b)
        if re:self.redge(self.ebove(vrt),False)
        self.verts[vrt.ix] = None
        self.vistack.append(vrt.ix)

    # provided two vertices and possibly a geometric index
    # create a new edge from vrt1 to vrt2
    def aedge(self,vrt1,vrt2):
        edg = self.edgeclass(vrt1,vrt2,self.edgecount)
        if self.eistack:
            ex = self.eistack.pop(0)
            self.edges[ex] = edg
        else:
            self.edges.append(edg)
            ex = self.edgecount
            self.edgecount += 1
        vrt1.connect(edg)
        vrt2.connect(edg)
        return edg

    # given edge, delete child vertex and its children
    # and the edge itself
    # NOTE: child is the vertex at the end of the edge
    def redge(self,edge,rv = True):
        if not edge in self.edges:return
        edge.one.disconnect(edge)
        edge.two.disconnect(edge)
        self.edges[edge.ix] = None
        self.eistack.append(edge.ix)
        if rv:self.rvert(edge.two,False)

    # provided an edge, create a new vert between 
    # edge.one and edge.two and connect
    def sedge(self,edge,*args,**kwargs):
        one,two = self.mask(0,None,edge)
        new = self.avert(one,*args,**kwargs)
        self.redge(self.ebove(two),False)
        self.aedge(new,two)
        return new

    # remove a vertex by connecting its children to its parent
    def medge(self,vrt):
        children = self.below(vrt)
        parent = self.above(vrt)
        es = self.mask(1,vrt,None)
        for e in es:self.redge(e,False)
        for ch in children:self.aedge(parent,ch)
        self.rvert(vrt,False)







 



