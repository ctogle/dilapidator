import dilap.core.base as db
import dilap.core.tools as dpr
import dilap.core.vector as dpv
import dilap.core.pointset as dps

import dilap.mesh.tools as dtl

import dilap.graph.graph as dgg

import matplotlib.pyplot as plt
import pdb





class node(db.base):

    def connect(self,other):
        if not other in self.ring:
            self.ring.append(other)

    def disconnect(self,other):
        if other in self.ring:
            self.ring.remove(other)

    def __init__(self,index,px,**kwargs):
        self.index = index
        self.px = px
        self.ring = []

class edge(db.base):

    def connect(self,other):
        if not other in self.ring:
            self.ring.append(other)

    def disconnect(self,other):
        if other in self.ring:
            self.ring.remove(other)

    def __init__(self,index,one,two,**kwargs):
        self.index = index
        self.one = one
        self.two = two
        self.ring = []

# a wire graph simply allows for nodes and edges
class wire_graph(db.base):

    # returns all the d-cells incident on the 
    # (d-1)-cells and adjacent to the (d+1)-cells, after 
    # considering the restrictions imposed by non-None arguments
    def mask(self,d = 0,v = None,e = None):
        results = []
        if d == 0:
            if not v is None:
                for ve in self.mask(1,v,None,None):
                    for vv in self.mask(0,None,ve,None):
                        if not vv is v:results.append(vv)
            if not e is None:results.extend([e.one,e.two])
        elif d == 1:
            if not v is None:results.extend(v.ring)
            if not e is None:
                for vv in self.mask(0,None,e,None):
                    for ve in self.mask(1,vv,None,None):
                        if not ve is e:results.append(ve)
        return results

    def __init__(self,**kwargs):
        self.nodes = []
        self.nlook = {}
        self.edges = []
        self.elook = {}
        self.nodecount = 0
        self.edgecount = 0

    def add_node(self,npx,**kwargs):
        if npx in self.nlook:return self.nlook[npx]
        newnode = node(self.nodecount,npx,**kwargs)
        self.nodes.append(newnode)
        self.nlook[newnode.px] = self.nodecount
        self.nodecount += 1
        return newnode

    def add_edge(self,nd1,nd2,**kwargs):
        nd1 = self.nodes[nd1]
        nd2 = self.nodes[nd2]
        ekey1 = (nd1.index,nd2.index)
        ekey2 = (nd2.index,nd1.index)
        if nd2 in nd1.ring:return self.elook[ekey1]
        newedge = edge(self.edgecount,nd1,nd2,**kwargs)
        nd1.connect(nd2)
        nd2.connect(nd1)
        self.edges.append(newedge)
        self.elook[ekey1] = newedge.index
        self.elook[ekey2] = newedge.index
        self.edgecount += 1
        return newedge









