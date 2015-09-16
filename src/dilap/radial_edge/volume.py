import dilap.core.base as db
import dilap.core.vector as dpv

import dilap.radial_edge.vertex as dvt
import dilap.radial_edge.edges as dve

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

    def plot(self,ax = None):
        if ax is None:ax = dtl.plot_axes()
        for vdx in range(self.vertexcount):
            vt = self.vertices[vdx]
            if vt is None:continue
            vt.plot(self.geom,ax)
        for edx in range(self.edgecount):
            eg = self.edges[edx]
            if eg is None:continue
            eg.plot(self.geom,ax)
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

    def __init__(self,geom):
        self.ix = index()
        self.geom = geom
        self.vertices = []
        self.vertexcount = 0
        self.edges = []
        self.edgecount = 0
        self.loops = []
        self.loopcount = 0
        self.faces = []
        self.facecount = 0

    def cube(self):
        v1 = self.new_vertex()
        self.geom.add_point(v1.ix,dpv.vector(0,0,1))
        return self






