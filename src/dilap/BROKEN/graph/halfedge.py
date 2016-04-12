import dilap.core.base as db
import dilap.core.tools as dpr
import dilap.core.vector as dpv
import dilap.core.pointset as dps

import dilap.mesh.tools as dtl

import matplotlib.pyplot as plt
import pdb





class node(db.base):

    def __init__(self,index,px,**kwargs):
        self.index = index
        self.px = px
        self.halfedge = None

class edge(db.base):

    def __init__(self,index,**kwargs):
        self.index = index
        self.end = None
        self.fce = None
        self.nxt = None
        self.lst = None

class face(db.base):

    def __init__(self,index,**kwargs):
        self.index = index
        self.halfedge = None
        self.fnorm = None

# a twomanifold graph represents a topological shell with holes
class halfedge_graph(db.base):

    # returns all the d-cells incident on the 
    # (d-1)-cells and adjacent to the (d+1)-cells, after 
    # considering the restrictions imposed by non-None arguments
    def mask(self,d = 0,v = None,e = None,l = None):
        results = []
        if d == 0:
            if not v is None:
                for ve in self.mask(1,v,None,None):
                    for vv in self.mask(0,None,ve,None):
                        if not vv is v:results.append(vv)
            if not e is None:
                raise NotImplemented
                #results.extend([e.one,e.two])
            if not l is None:
                for ve in self.mask(1,None,None,l):
                    for vv in self.mask(0,None,ve,None):
                        if not vv in results:results.append(vv)
        elif d == 1:
            if not v is None:
                raise NotImplemented
                #results.extend(v.ring)
            if not e is None:
                for vv in self.mask(0,None,e,None):
                    for ve in self.mask(1,vv,None,None):
                        if not ve is e:results.append(ve)
            if not l is None:
                results.append(l.halfedge)
                nhe = results[-1].nxt
                while not nhe in results:
                    results.append(nhe)
                    nhe = results[-1].nxt
        elif d == 2:
            if not v is None:
                for ve in self.mask(1,v,None,None):
                    for vl in self.mask(2,None,ve,None):
                        if not vl in results:results.append(vl)
            if not e is None:
                results.append(e.fce)
                results.append(e.opposite.fce)
            if not l is None:
                for ve in self.mask(1,None,None,l):
                    for vl in self.mask(2,None,ve,None):
                        if not vl is l:results.append(vl)
        return results

    def __init__(self,**kwargs):
        self.nodes = []
        self.edges = []
        self.faces = []
        self.nodecount = 0
        self.edgecount = 0
        self.facecount = 0

    def add_node(self,npx,**kwargs):
        newnode = node(self.nodecount,npx,**kwargs)
        self.nodes.append(newnode)
        self.nodecount += 1
        return newnode

    def add_edge(self,nd1,nd2,**kwargs):
        nd1 = self.nodes[nd1]
        nd2 = self.nodes[nd2]
        newedge = edge(self.edgecount,**kwargs)
        newedge.end = nd2
        nd1.halfedge = newedge
        self.edges.append(newedge)
        self.edgecount += 1
        return newedge

    def add_face(self,*edges,**kwargs):
        newface = face(self.facecount,**kwargs)
        for eg in edges:eg.fce = newface
        newface.halfedge = edges[0]
        self.faces.append(newface)
        self.facecount += 1
        return newface

    def new_face(self,pxs,**kwargs):
        nxs = [self.add_node(x) for x in pxs]
        exs = [self.add_edge(nxs[x-1].index,nxs[x].index) 
                                for x in range(len(nxs))]
        return self.add_face(*exs,**kwargs)







 
