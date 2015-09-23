import dilap.core.tools as dpr

import dilap.topology.vert as dvt
import dilap.topology.edge as deg
import dilap.topology.loop as dlp
import dilap.topology.face as dfc

from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat










class pmeshvert(dvt.vert):

    def __init__(self,ix,p):
        dvt.vert.__init__(self,ix)
        self.p = p

# pmeshface is a concave polygon with holes
class pmeshface(dfc.face):

    def __init__(self,bound,holes,ix):
        dfc.face.__init__(self,ix)
        self.bound = bound
        self.holes = holes

__doc__ = '''dilapidator\'s implementation of a polygonal mesh'''
# dilapidators implementation of a polygonal mesh
# analogous to a single solid or topological shell
class polygonmesh:

    vertclass = pmeshvert
    edgeclass = deg.edge
    loopclass = dlp.loop
    faceclass = pmeshface

    # produce a trimesh representing all faces of this mesh
    # this requires geometric data; should it be done elsewhere??
    def trimesh(self):
        raise NotImplemented

    def __str__(self):return 'polygonmesh:'

    def vcnt(self):return self.cnt(self.verts)
    def ecnt(self):return self.cnt(self.edges)
    def lcnt(self):return self.cnt(self.loops)
    def fcnt(self):return self.cnt(self.faces)
    def cnt(self,tobjs):
        x = 0
        for v in tobjs:
            if not v is None:
                x += 1
        return x

    def __init__(self,*args,**kwargs):
        if args:pinit = args[0]
        else:pinit = vec3(0,0,0)
        self.verts = []
        self.edges = []
        self.loops = []
        self.faces = []
        self.vertcount = 0
        self.edgecount = 0
        self.loopcount = 0
        self.facecount = 0
        self.aroot(pinit,*args,**kwargs)

    # return the d-cells which are incident upon any of v,e,l,f
    def mask(self,d = 0,v = None,e = None,l = None,f = None):
        results = []
        if d == 0:
            if not v is None:
                for ve in self.mask(1,v,None):
                    for vv in self.mask(0,None,ve):
                        if not vv is v:results.append(vv)
            if not e is None:results.extend([e.one,e.two])
            if not l is None:pass
            if not f is None:pass
        elif d == 1:
            if not v is None:results.extend(v.ring)
            if not e is None:
                for vv in self.mask(0,None,e):
                    for ve in self.mask(1,vv,None):
                        if not ve is e:results.append(ve)
            if not l is None:pass
            if not f is None:pass
        elif d == 2:
            if not v is None:pass
            if not e is None:pass
            if not l is None:pass
            if not f is None:pass
        elif d == 3:
            if not v is None:pass
            if not e is None:pass
            if not l is None:pass
            if not f is None:pass
        return results

    # add a new root vertex in the mesh (only needed once)
    # MBFLV - analog
    # NOTE: allows passing to vertclass but not loopclass...
    def aroot(self,*args,**kwargs):
        vrt = self.vertclass(self.vertcount,*args,**kwargs)
        self.verts.append(vrt)
        self.vertcount += 1
        edg = self.aedge(vrt,vrt)
        lop = self.aloop(edg)
        fac = self.aface(lop,[])
        return vrt,edg,lop,fac

    # create new vertex connected to par by edge ending at new vertex
    # MEV - analog
    # split the edge by adding a new vertex between its endpoints
    # remove the current edge and add two move, connect each 
    #  of its endpoints to the newly added vertex
    def sedge(self,edg,*args,**kwargs):
        nvrt = self.avert(*args,**kwargs)
        vrt1,vrt2 = edg.one,edg.two
        self.redge(edg)
        eg1 = self.aedge(vrt1,nvrt)
        eg2 = self.aedge(nvrt,vrt2)
        return nvrt,eg1,eg2

    # create new vertex connected to par by edge ending at new vertex
    # MV - analog - non-euler?
    def avert(self,*args,**kwargs):
        vrt = self.vertclass(self.vertcount,*args,**kwargs)
        self.verts.append(vrt)
        self.vertcount += 1
        return vrt

    # provided two vertices, create a new edge from vrt1 to vrt2
    # ME - analog
    # subcases:
    #   MEKL - analog - when two loops on the same face are connected
    #   MEFL - analog - this is how faces are added incrementally?
    #   MEKBFL - analog - this is essentially merging with another mesh
    def aedge(self,vrt1,vrt2,*args,**kwargs):
        edg = self.edgeclass(vrt1,vrt2,self.edgecount,*args,**kwargs)
        self.edges.append(edg)
        self.edgecount += 1
        vrt1.connect(edg)
        vrt2.connect(edg)
        return edg

    # GLUE - analog - glue two faces together, possibly from disjoint meshes
    # can create handles through solids?
    # subcases:
    #   KFLEVMG
    #   KFLEVB
    # UNGLUE - analog - glue two faces together, possibly from disjoint meshes
    # can create handles through solids?
    # subcases:
    #   MFLEVKG
    #   MFLEVB

    # provided a sequence of ordered edges, create a new loop
    def aloop(self,edg,*args,**kwargs):
        lop = self.loopclass(edg,self.loopcount,*args,**kwargs)
        self.loops.append(lop)
        self.loopcount += 1
        return lop

    # provided a sequence of loops, create a new face
    def aface(self,elp,ilps,*args,**kwargs):
        fac = self.faceclass(elp,ilps,self.facecount,*args,**kwargs)
        self.faces.append(fac)
        self.facecount += 1
        return fac

    # provided an edge, remove it from the mesh
    # KE - analog
    def redge(self,edg):
        self.edges[edg.ix] = None
        vrt1,vrt2 = edg.one,edg.two
        vrt1.disconnect(edg)
        vrt2.disconnect(edg)

    ###################################################
    ###################################################
    ###################################################
    ###################################################
    ###################################################

    '''#

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
    '''#






 




