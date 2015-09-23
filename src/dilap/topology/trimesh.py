import dilap.core.tools as dpr

import dilap.topology.vert as dvt
import dilap.topology.edge as deg
import dilap.topology.loop as dlp
import dilap.topology.face as dfc

from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat










__doc__ = '''dilapidator\'s implementation of a triangle mesh'''
# dilapidators implementation of a triangle mesh
# analogous to a single manifold surface
# NOTE: no explicit v,e,f classes are used
class trimesh:

    def __str__(self):return 'trimesh:'
    def vcnt(self):return self.cnt(self.verts)
    def ecnt(self):return self.cnt(self.edges)
    def fcnt(self):return self.cnt(self.faces)
    def cnt(self,tobjs):
        x = 0
        for v in tobjs:
            if not v is None:
                x += 1
        return x

    # vertices are tuples (pos_index,nrm_index,uv_index)
    # edges are tuples of vertex indices (vx1,vx2)
    #   edges are implicited directed (ccw about faces)
    # faces are tuples of vertex indices (vx1,vx2,vx3)
    # vertices form keys pointing to a list of other vertices
    #   with which that vertex forms an edge
    # edges form keys pointing to faces which they bound
    # when a vertex is deleted: del self.verts[v]
    # when an edge is deleted: del self.edges[e]
    # when a face is deleted: self.faces[fx] = None
    #   del the edges and verts (if no other edges connect to them)
    #
    # when creating geometry, it is not a consideration if an edge
    # exists, since it only exists if a face is attached to it
    #
    def __init__(self,*args,**kwargs):
        self.verts = []         # list of vertex tuples (px,nx,ux)
        self.vertcount = 0      # 
        self.vistack = []
        self.edges = []         # list of edge tuples (vx1,vx2)
        self.edgecount = 0      # 
        self.eistack = []
        self.faces = []         # list of face tuples (ex1,ex2,ex3)
        self.facecount = 0      # 
        self.fistack = []

        self.ve_rings = {}      # lookup of edge rings per vertex
        self.ef_rings = {}      # lookup of faces by edge tuples
        self.fs_mats = {}       # loopup of face material per face

    # return the d-cells which are incident upon any of v,e,f
    # v is a tuple of (px,nx,ux)
    # e is a tuple of (vx1,vx2)
    # f is a tuple of (ex1,ex2,ex3)
    def mask(self,d = 0,v = None,e = None,f = None):
        results = []
        if d == 0:
            if not v is None:
                for ve in self.mask(1,v,None,None):
                    for vv in self.mask(0,None,ve,None):
                        if not vv is v and not vv in results:
                            results.append(vv)
            if not e is None:
                results.append(self.verts[e[0]])
                results.append(self.verts[e[1]])
            if not f is None:
                results.append(self.verts[self.edges[f[0]][0]])
                results.append(self.verts[self.edges[f[1]][0]])
                results.append(self.verts[self.edges[f[2]][0]])
        elif d == 1:
            if not v is None:results.extend(self.ve_rings[v])
            if not e is None:
                for vv in self.mask(0,None,e,None):
                    for ve in self.mask(1,vv,None,None):
                        if not ve is e:results.append(ve)
            if not f is None:
                results.append(self.edges[f[0]])
                results.append(self.edges[f[1]])
                results.append(self.edges[f[2]])
        elif d == 2:
            if not v is None:
                for ve in self.mask(1,v,None,None):
                    vf = self.ef_rings[ve]
                    if not vf in results:results.append(vf)
            if not e is None:results.append(self.ef_rings[e])
            if not f is None:
                for ve in self.mask(1,None,None,f):
                    ofekey = (ve[1],ve[0])
                    if ofekey in self.ef_rings:
                        results.append(self.ef_rings[ofekey])
        return results

    def vonb(self,v):
        for ve in self.mask(1,v,None,None):
            if self.eonb(ve):return True
        return False

    def eonb(self,e):
        if (e[1],e[0]) in self.ef_rings:return False
        else:return True

    # create new vertex and return its index
    # a vert can be created disconnected from the mesh
    #   but should immediately be joined via two edges 
    #   and a face at minimum
    # a vertex normal and uv can be calculated later, 
    #   but should be unique
    def avert(self,px,nx,ux):
        vrt = (px,nx,ux)
        self.ve_rings[vrt] = []
        if self.vistack:
            vx = self.vistack.pop(0)
            self.verts[vx] = vrt
        else:
            self.verts.append(vrt)
            vx = self.vertcount
            self.vertcount += 1
        return vx

    # destory existing vertex
    # destroy any faces touching this vertex
    # destroy any edges touching this vertex
    def rvert(self,v):
        es = self.mask(1,v,None,None)
        fs = self.mask(2,v,None,None)
        for f in fs:self.rface(f)
        for e in es:self.redge(e)
        del self.ve_rings[v]
        vx = self.verts.index(v)
        self.verts[vx] = None
        self.vistack.append(vx)

    # create new edge and return its index
    # NOTE: call only from self.aface
    def aedge(self,u,v):
        edg = (u,v)
        self.ve_rings[self.verts[u]].append(edg)
        self.ve_rings[self.verts[v]].append(edg)
        self.ef_rings[edg] = None
        if self.eistack:
            ex = self.eistack.pop(0)
            self.edges[ex] = edg
        else:
            self.edges.append(edg)
            ex = self.edgecount
            self.edgecount += 1
        return ex

    # destroy existing edge
    # destroy faces connected to this edge
    # destroy verts disconnected by removal of this edge
    def redge(self,e):
        vs = self.mask(0,None,e,None)
        fs = self.mask(2,None,e,None)
        for f in fs:self.rface(f)
        for v in vs:
            self.ve_rings[v].remove(e)
            if not self.ve_rings[v]:
                self.rvert(v)
        del self.ef_rings[e]
        ex = self.edges.index(e)
        self.edges[ex] = None
        self.eistack.append(ex)

    # given the indices of three existing vertices,
    # create new face and return its index
    def aface(self,u,v,w):
        fac = (u,v,w)
        self.fs_mats[fac] = None
        self.aedge(u,v)
        self.aedge(v,w)
        self.aedge(w,u)
        self.ef_rings[(u,v)] = fac
        self.ef_rings[(v,w)] = fac
        self.ef_rings[(w,u)] = fac
        if self.fistack:
            fx = self.fistack.pop(0)
            self.faces[fx] = fac
        else:
            self.faces.append(fac)
            fx = self.facecount
            self.facecount += 1
        return fx

    # destroy existing face
    def rface(self,f):
        es = self.mask(1,None,None,f)
        for e in es:self.redge(e)
        fx = self.faces.index(f)
        self.faces[fx] = None
        self.fistack.append(fx)

    # create 3 new triangles from an existing triangle
    def sface(self,tri,*args,**kwargs):
        raise NotImplemented

    # NEED ALL METHODS FOR TRANSFORMATIONS!
    # translate, rotate, scale






 




