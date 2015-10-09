import dilap.core.tools as dpr

import dilap.topology.vert as dvt
import dilap.topology.edge as deg
import dilap.topology.loop as dlp
import dilap.topology.face as dfc

from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat

import pdb

# dilapidators implementation of a topological edge
class halfedge:

    def __init__(self,one,two,l):
        self.one = one # associated edge start
        self.two = two # associated edge end
        self.l = l # associated loop
        self.nxt = self # next clockwise edge
        self.lst = self # next counter-clockwise edge 
        self.ops = self # opposite halfedge

__doc__ = '''dilapidator\'s implementation of a polygonal mesh'''
# dilapidators implementation of a polygonal mesh
# analogous to a single solid or topological shell
class polygonmesh:

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

    # get an available index for a new vertex
    def vi(self):
        if self.vistack:return self.vistack.pop(0)
        else:
            self.verts.append(None)
            self.vertcount += 1
            return self.vertcount-1
    # get an available index for a new edge
    def ei(self):
        if self.eistack:return self.eistack.pop(0)
        else:
            self.edges.append(None)
            self.edgecount += 1
            return self.edgecount-1
    # get an available index for a new loop
    def li(self):
        if self.listack:return self.listack.pop(0)
        else:
            self.loops.append(None)
            self.loopcount += 1
            return self.loopcount-1
    # get an available index for a new face
    def fi(self):
        if self.fistack:return self.fistack.pop(0)
        else:
            self.faces.append(None)
            self.facecount += 1
            return self.facecount-1

    # return the d-cells which are incident upon any of v,e,l,f
    # with the exception that d = 4 indicates halfedges
    def mask(self,d = 0,v = None,e = None,l = None,f = None):
        results = []

        if d == 0:
            if not l is None:
                for he in self.mask(4,None,None,l,None):
                    if not he.one in results:results.append(he.one)
                    if not he.two in results:results.append(he.two)
          
            else:raise NotImplemented

        elif d == 1:
            if not v is None:
                for ve in self.vrings[v]:
                    if not ve in results:
                        results.append(ve)
          
            else:raise NotImplemented

        elif d == 2:
            if not v is None:
                for ve in self.mask(1,v,None,None,None):
                    for vl in self.mask(2,None,ve,None,None):
                        if not vl in results:
                            results.append(vl)
            elif not e is None:
                for he in self.erings[e]:
                    hel = self.loops[he.l]
                    if not hel in results:
                        results.append(hel)
            elif not f is None:
                for fl in self.frings[f]:
                    if not fl in results:
                        results.append(fl)
                    
            else:raise NotImplemented

        elif d == 3:raise NotImplemented

        elif d == 4:
            if l is None:raise ValueError
            else:
                hes = [l[1]]
                while not hes[-1].nxt in hes:
                    hes.append(hes[-1].nxt)
            if v is None:return hes
            else:return [h for h in hes if h.one == v or h.two == v]

        return results

    def __init__(self,*args,**kwargs):
        self.verts = []     # list of tuples (geo-id,)
        self.vrings = {}    # list of edges connected to each vertex
        self.vistack = []   # list of indices of currently deleted verts
        self.vertcount = 0

        self.edges = []     # list of tuples (v1-id,v2-id,geo-id)
        self.erings = {}    # list of half-edges connected to each each
        self.eistack = []   # list of indices of currently deleted edges
        self.edgecount = 0

        self.loops = []     # list of tuples (e1-id,)
        self.lrings = {}    # list faces which use this loop
        self.listack = []   # list of indices of currently deleted loops
        self.loopcount = 0

        self.faces = []     # list of tuples (eloop-id,geo-id)
        self.frings = {}    # list of loops in each face (e-bound first)
        self.fistack = []   # list of indices of currently deleted faces
        self.facecount = 0

    ###########################################################################
    ### non eulerian access to the datastructure
    ###########################################################################

    # place a new vertex in the datastructure
    # initialize a corresponding ring and return the v index
    def avert(self,g = 0):
        vx = self.vi()
        vrt,vring = (vx,g),[]
        self.verts[vx] = vrt
        self.vrings[vrt] = vring
        return vrt

    # place a new edge in the datastructure
    # initialize a corresponding ring and return the e index
    # u is the start vertex, v is the end vertex
    def aedge(self,u,v,g = 0):
        ex = self.ei()
        edg,ering = (u,v,ex,g),[]
        self.edges[ex] = edg
        self.vrings[u].append(edg)
        self.vrings[v].append(edg)
        self.erings[edg] = ering
        return edg

    # place a new loop in the datastructure
    # initialize a corresponding ring and return the l index
    # the edge e gets a half-edge with this loop on its left
    # b indicates whether e.one,e.two is backward wrt to the new loop
    def aloop(self,e,b = 0):
        eu,ev = e[0],e[1]
        if   b == 0:u,v = eu,ev
        elif b == 1:v,u = eu,ev
        lx = self.li()
        he = halfedge(u,v,lx)
        lop,lring = (e,he,lx),[]
        self.loops[lx] = lop
        self.erings[e].append(he)
        self.lrings[lop] = lring
        return lop
  
    # place a new face in the datastructure
    # initialize a corresponding ring and return the f index
    def aface(self,v,g = 0):
        fx = self.fi()
        fac,fring = (v,fx,g),[]
        self.faces[fx] = fac
        self.frings[fac] = fring
        return fac

    ###########################################################################
    ### basic euler operators
    ###########################################################################

    # MBFV
    # add a new root vertex in the mesh (only needed once)
    def mfv(self,vgx = 0,vfx = 0):
        rv = self.avert(vgx)
        rf = self.aface(rv,vfx)
        return rv,rf

    # MEV
    #   1 - new spur edge/vertex
    #   2 - split an edge 
    #   3 - split a vertex
    # add an edge and a vertex
    #   bv is the base vertex index
    #   oe is the orienting edge
    #   vgx is the geometric index of the new vertex
    #   egx is the geometric index of the new edge
    def mev(self,bv,oe = None,ol = None,vgx = 0,egx = 0):
        vr = self.vrings[bv]
        rc = len(vr)

        nv = self.avert(vgx)
        ne = self.aedge(bv,nv,egx)

        if rc == 0: # the vertex is the root of a face
            nl = self.aloop(ne)
            ops = halfedge(nl[1].two,nl[1].one,nl[2])
            ops.nxt = nl[1]
            ops.lst = nl[1]
            ops.ops = nl[1]
            nl[1].nxt = ops
            nl[1].lst = ops
            nl[1].ops = ops
            self.erings[nl[0]].append(ops)
            for f in self.faces:
                if f is None:continue
                if f[0] == bv:
                    self.frings[f].append(nl)

        elif rc == 1: # the vertex is a spur 
            se = vr[0]
            hel,her = self.erings[se]
            nhel = halfedge(ne[0],ne[1],hel.l)
            nher = halfedge(ne[1],ne[0],hel.l)
            self.erings[ne].append(nhel)
            self.erings[ne].append(nher)
            nhel.nxt = nher;nhel.lst = hel
            nhel.nxt.lst = nhel;nhel.lst.nxt = nhel
            nher.lst = nhel;nher.nxt = her
            nher.nxt.lst = nher;nher.lst.nxt = nher

        elif rc >= 2: # the vertex is in a 2+ edge loop
            
            pdb.set_trace()

        return nv,ne

    # MFE
    #   split a face into two
    #   slice an edge
    # add an edge and a face
    #   u,v are the start,end vertices
    #   e1,e2 specify the new loop boundaries
    #   egx,fgx are geometric indices
    #def mfe(self,u,v,e1 = None,e2 = None,egx = 0,fgx = 0):
    def mfe(self,u,v,f = None,egx = 0,fgx = 0):
        # find the common loop to u and v
        uls = self.mask(2,u,None,None,None)
        vls = self.mask(2,v,None,None,None)
        uvls = tuple(x for x in uls if x in vls)

        if f is None:
            if len(uvls) == 0:
                print('mfe u,v share no loops!')
                raise ValueError
            elif len(uvls) == 1:oll = uvls[0]
        else:
            oll = tuple(x for x in uvls if x in self.frings[f])[0]

        # find the edges of the common loop incident on u and v
        lue = self.mask(4,u,None,oll,None)
        lve = self.mask(4,v,None,oll,None)

        if not lue or not lve:
            print('wtffff')
            pdb.set_trace()

        ne = self.aedge(u,v,egx)
        nl = self.aloop(ne)
        nf = self.aface(u,fgx)
        self.frings[nf].append(nl)

        nlx = nl[2]
        olhe = halfedge(u,v,oll[2])
        nlhe = nl[1]
        olhe.ops = nlhe
        nlhe.ops = olhe

        olue = tuple(x for x in lue if x.one == u)[0]
        nlue = tuple(x for x in lue if x.two == u)[0]
        olve = tuple(x for x in lve if x.two == v)[0]
        nlve = tuple(x for x in lve if x.one == v)[0]
        olhe.nxt = olue;olhe.lst = olve
        nlhe.nxt = nlve;nlhe.lst = nlue
        olhe.nxt.lst = olhe;olhe.lst.nxt = olhe
        nlhe.nxt.lst = nlhe;nlhe.lst.nxt = nlhe

        # need to traverse each loop setting the loop reference appropriately
        c = nlhe.nxt
        while not c is nlhe:
            c.l = nl[2]
            c = c.nxt

        return ne,nf







    # return the d-cells which are incident upon any of v,e,l,f
    def mask_____(self,d = 0,v = None,e = None,l = None,f = None):
        results = []
        if d == 0:
            if not v is None:
                for ve in self.mask(1,v,None):
                    for vv in self.mask(0,None,ve):
                        if not vv is v:results.append(vv)
            if not e is None:results.extend([e.one,e.two])
            if not l is None:
                for e in self.mask(1,None,None,l,None):
                    results.append(e.one)
            if not f is None:
                el,ils = f.bound,f.holes
                ev,ivs = [],[]
                for evv in self.mask(0,None,None,el,None):
                    ev.append(evv)
                #for il in ils:
                #    ivs.append([])
                #    for ivv in self.mask(0,None,None,il,None):
                #        ivs[-1].append(ivv)
                results.extend(ev)
        elif d == 1:
            if not v is None:results.extend(v.ring)
            if not e is None:
                for vv in self.mask(0,None,e):
                    for ve in self.mask(1,vv,None):
                        if not ve is e:results.append(ve)
            if not l is None:
                first = l.edge
                results.append(first)
                nxt = first.nxt
                while not nxt in results:
                    results.append(nxt)
                    nxt = nxt.nxt
            if not f is None:pass
        elif d == 2:
            if not v is None:pass
            if not e is None:
                for l in self.loops:
                    if l is None:continue
                    if l.edge == e:results.append(l)
            if not l is None:pass
            if not f is None:pass
        elif d == 3:
            if not v is None:pass
            if not e is None:pass
            if not l is None:pass
            if not f is None:pass
        return results

    '''#
    # add a new root vertex in the mesh (only needed once)
    # MBFLV - analog
    # NOTE: allows passing to vertclass but not loopclass...
    def aroot___(self,*args,**kwargs):
        vrt = self.vertclass(self.vertcount,*args,**kwargs)
        self.verts.append(vrt)
        self.vertcount += 1
        edg = self.aedge(vrt,vrt)
        lop = self.aloop(edg)
        fac = self.aface(lop,[])
        return vrt,edg,lop,fac

    # create new vertex connected to par by edge ending at new vertex
    # MV - analog - non-euler?
    def avert___(self,*args,**kwargs):
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
    def aedge____(self,vrt1,vrt2,*args,**kwargs):
        edg = self.edgeclass(vrt1,vrt2,self.edgecount,*args,**kwargs)
        self.edges.append(edg)
        self.edgecount += 1
        vrt1.connect(edg)
        vrt2.connect(edg)
        return edg






    # create new vertex connected to par by edge ending at new vertex
    # MEV - analog
    # split the edge by adding a new vertex between its endpoints
    # remove the current edge and add two move, connect each 
    #  of its endpoints to the newly added vertex
    def sedge___(self,edg,nvrt,*args,**kwargs):
        vrt1,vrt2 = edg.one,edg.two
        erng = self.mask(1,None,edg,None,None)
        self.redge(edg)
        eg1 = self.aedge(vrt1,nvrt)
        eg2 = self.aedge(nvrt,vrt2)
        eg1.nxt = eg2
        eg2.lst = eg1
        if vrt1 is vrt2:
            eg2.nxt = eg1
            eg1.lst = eg2
            loop = edg.loop
            eg1.loop = loop
            loop.edge = eg1
        else:
            eg2.nxt = edg.nxt
            eg1.lst = edg.lst
            eg1.lst.nxt = eg1
            eg2.nxt.lst = eg2
        return eg1,eg2

    # MEFL - analog - this is how faces are added incrementally?
    def mefl(self,vrt1,vrt2,*args,**kwargs):
        edg1 = self.aedge(vrt1,vrt1,*args,**kwargs)
        lop = self.aloop(edg1)
        fac = self.aface(lop,[])
        edg1,edg2 = self.sedge(edg1,vrt2)
        return edg1,edg2,lop,fac

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
    def aloop____(self,edg,*args,**kwargs):
        lop = self.loopclass(edg,self.loopcount,*args,**kwargs)
        self.loops.append(lop)
        self.loopcount += 1
        return lop

    # provided a sequence of loops, create a new face
    def aface____(self,elp,ilps,*args,**kwargs):
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

    # given a list of vertex indices, return their positions
    def gvps(self,f):
        vs = self.mask(0,None,None,None,f)
        vps = []
        for vx in range(len(vs)):
            vps.append(vs[vx].p)
        return vps
    '''#

    ###################################################
    ###################################################





 




