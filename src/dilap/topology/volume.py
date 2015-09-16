import dilap.core.base as db

import dilap.topology.vertex as dvt
import dilap.topology.edges as dve
import dilap.topology.shell as dvs

import dilap.mesh.tools as dtl

import pdb

vix = 0
def index():
    global vix
    vix +=1 
    return vix

def mask(vol,d = 0,v = None,e = None,l = None):
    results = []
    if d == 0:
        if not v is None:
            for ve in mask(vol,1,v,None,None):
                for vv in mask(vol,0,None,ve,None):
                    if not vv is v:results.append(vv)
        if not e is None:results.extend([e.one,e.two])
        if not l is None:
            for ve in mask(vol,1,None,None,l):
                for vv in mask(vol,0,None,ve,None):
                    if not vv in results:results.append(vv)
    elif d == 1:
        if not v is None:results.extend(v.ring)
        if not e is None:
            for vv in mask(vol,0,None,e,None):
                for ve in mask(vol,1,vv,None,None):
                    if not ve is e:results.append(ve)
        if not l is None:results.extend(l.edges)
    elif d == 2:
        if not v is None:
            for ve in mask(vol,1,v,None,None):
                for vl in mask(vol,2,None,ve,None):
                    if not vl in results:results.append(vl)
        if not e is None:results.extend(e.ring)
        if not l is None:
            for ve in mask(vol,1,None,None,l):
                for vl in mask(vol,2,None,ve,None):
                    if not vl is l:results.append(vl)
    return results

class volume(db.base):

    def plot(self,geom,ax):
        for vdx in range(self.vertexcount):
            vt = self.vertices[vdx]
            if vt is None:continue
            vt.plot(geom,ax)
        for edx in range(self.edgecount):
            eg = self.edges[edx]
            if eg is None:continue
            eg.plot(geom,ax)
        return ax

    def new_vertex(self):
        nv = dvt.vertex()
        self.vertices.append(nv)
        self.vertexcount += 1
        return nv

    # needs to detect 3 cases!!
    # MEKL,MEFL,MEKBFL
    def new_edge(self,sv,ev):
        ne = dve.halfedge()
        ne.tip = ev
        ne.tail = sv
        self.edges.append(ne)
        self.edgecount += 1
        return ne

    # MEV
    def new_edge_vertex(self,v):
        nv = self.new_vertex()

        pdb.set_trace()

        ne = self.new_edge(v,nv)
        return nv,ne

    # MBFLV
    def new_shell(self):
        ns = dvs.shell()
        nv,nl,nf = ns.new_face()
        self.vertices.append(nv)
        self.vertexcount += 1
        self.faces.append(nf)
        self.facecount += 1
        self.loops.append(nl)
        self.loopcount += 1
        self.shells.append(ns)
        self.shellcount += 1
        return nv,nl,nf,ns

    def __init__(self):
        self.ix = index()
        self.vertices = []
        self.vertexcount = 0
        self.edges = []
        self.edgecount = 0
        self.loops = []
        self.loopcount = 0
        self.faces = []
        self.facecount = 0
        self.shells = []
        self.shellcount = 0



