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


# remove each item in one that is not in two
def isect(one,two):
    extras = 0
    for x in range(len(one)):
        if not one[x] in two:
            one[x] = None
            extras += 1
    while extras:
        one.remove(None)
        extras -= 1


# add any item in two to one if its not already present
def union(one,two):
    raise NotImplemented


__doc__ = '''
dilapidator\'s implementation of a polygonal mesh

a polygonmesh has lists for six topological element types:

    vertices - tuple(vertex_index,geometric_index)
        vertices have associated rings of edges

    edges - tuple(start_vertex,end_vertex,edge_index,geometric_index)
        edges have associated rings of halfedges

    loops - tuple(edge,halfedge,loop_index)
        loops have associated rings of faces which use them

    faces - tuple(vertex,face_index,geometric_index)
        faces have associated rings of loops they use

    facegroups - tuple(face,facegroup_index)
        facegroups have associated rings of included faces

    shells - tuple(vertex,shell_index)
        shells have associated rings of facegroups
'''
# dilapidators implementation of a polygonal mesh
# analogous to a single solid or topological shell
class polygonmesh:

    def __str__(self):return 'polygonmesh'
    def  vcnt(self):return self.cnt(self.verts)
    def  ecnt(self):return self.cnt(self.edges)
    def  lcnt(self):return self.cnt(self.loops)
    def  fcnt(self):return self.cnt(self.faces)
    def fgcnt(self):return self.cnt(self.fgroups)
    def  scnt(self):return self.cnt(self.shells)
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
    # get an available index for a new facegroup
    def fgi(self):
        if self.fgistack:return self.fgistack.pop(0)
        else:
            self.fgroups.append(None)
            self.fgroupcount += 1
            return self.fgroupcount-1
    # get an available index for a new shell
    def si(self):
        if self.sistack:return self.sistack.pop(0)
        else:
            self.shells.append(None)
            self.shellcount += 1
            return self.shellcount-1

    # return verts incident on v,e,l,f,fg,s
    def mask0(self,v = None,e = None,l = None,f = None,fg = None,s = None):
        results = []
        if not l is None:
            for he in self.mask(6,None,None,l,None):
                if not he.one in results:results.append(he.one)
                if not he.two in results:results.append(he.two)
        else:raise NotImplemented
        return results

    # return edges incident on v,e,l,f,fg,s
    def mask1(self,v = None,e = None,l = None,f = None,fg = None,s = None):
        results = []
        if not v is None:
            for ve in self.vrings[v]:
                if not ve in results:
                    results.append(ve)
        else:raise NotImplemented
        return results

    # return loops incident on v,e,l,f,fg,s
    def mask2(self,v = None,e = None,l = None,f = None,fg = None,s = None):
        results = []
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
        return results

    # return faces incident on v,e,l,f,fg,s
    def mask3(self,v = None,e = None,l = None,f = None,fg = None,s = None):
        raise NotImplemented

    # return facegroups incident on v,e,l,f,fg,s
    def mask4(self,v = None,e = None,l = None,f = None,fg = None,s = None):
        raise NotImplemented

    # return shells incident on v,e,l,f,fg,s
    def mask5(self,v = None,e = None,l = None,f = None,fg = None,s = None):
        raise NotImplemented

    # return halfedges incident on v,e,l,f,fg,s
    def mask6(self,v = None,e = None,l = None,f = None,fg = None,s = None):
        results = []
        if not v is None:
            for ve in self.mask(1,v,None,None,None):
                for he in self.mask(6,None,ve,None,None):
                    if not he in results:
                        results.append(he)
        if not e is None:
            subresults = []
            for he in self.erings[e]:
                if not he in subresults:
                    subresults.append(he)
            if results:isect(results,subresults)
            else:results.extend(subresults)
        if not l is None:
            subresults = []
            firsthe = l[1]
            subresults.append(firsthe)
            subresults.append(firsthe.nxt)
            while not subresults[-1] is firsthe:
                subresults.append(subresults[-1].nxt)
            if results:isect(results,subresults)
            else:results.extend(subresults)
        if not f is None:
            subresults = []
            for l in self.frings[f]:
                for he in self.mask(6,None,None,l,None):
                    if not he in subresults:
                        subresults.append(he)
            if results:isect(results,subresults)
            else:results.extend(subresults)
        return results

    # return the d-cells which are incident upon any of v,e,l,f,fg,s
    # with the exception that d = 6 indicates halfedges
    def mask(self,d = 0,v = None,e = None,l = None,f = None,fg = None,s = None):
        if   d == 0:return self.mask0(v,e,l,f,fg,s)
        elif d == 1:return self.mask1(v,e,l,f,fg,s)
        elif d == 2:return self.mask2(v,e,l,f,fg,s)
        elif d == 3:return self.mask3(v,e,l,f,fg,s)
        elif d == 4:return self.mask4(v,e,l,f,fg,s)
        elif d == 5:return self.mask5(v,e,l,f,fg,s)
        elif d == 6:return self.mask6(v,e,l,f,fg,s)
        else:raise ValueError

    # given a vertex and a face, 
    # return the halfedge pointing away from v
    def hefrom(self,v,f):
        hes = self.mask(6,v,None,None,f)
        oe = tuple(x for x in hes if x.one is v)[0]
        return oe

    def __init__(self,*args,**kwargs):
        self.verts = []     # list of tuples (vindex,geo-id)
        self.vrings = {}    # list of edges connected to each vertex
        self.vistack = []   # list of indices of currently deleted verts
        self.vertcount = 0

        self.edges = []     # list of tuples (v1-id,v2-id,eindex,geo-id)
        self.erings = {}    # list of half-edges connected to each each
        self.eistack = []   # list of indices of currently deleted edges
        self.edgecount = 0

        self.loops = []     # list of tuples (e1-id,lindex)
        self.lrings = {}    # list faces which use this loop
        self.listack = []   # list of indices of currently deleted loops
        self.loopcount = 0

        self.faces = []     # list of tuples (eloop-id,findex,geo-id)
        self.frings = {}    # list of loops in each face (e-bound first)
        self.fistack = []   # list of indices of currently deleted faces
        self.facecount = 0

        self.fgroups = []    # list of tuples (fgindex,)
        self.fgrings = {}    # list of faces in each facegroup
        self.fgistack = []   # list of indices of currently deleted facegroups
        self.fgroupcount = 0

        self.shells = []    # list of tuples (vertex-id,face-id,sindex)
        self.srings = {}    # list of facegroups in each shell
        self.sistack = []   # list of indices of currently deleted shells
        self.shellcount = 0

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

    # place a new facegroup in the datastructure
    # initialize a corresponding ring and return the fg index
    def afgroup(self,f):
        fgi = self.fgi()
        fg,fgring = (fgi,),[f]
        self.fgroups[fgi] = fg
        self.fgrings[fg] = fgring
        return fg

    # place a new shell in the datastructure
    # initialize a corresponding ring and return the s index
    def ashell(self,fg):
        si = self.si()
        sh,sring = (si,),[fg]
        self.shells[si] = sh
        self.srings[sh] = sring
        return sh

    ###########################################################################
    ### basic euler operators
    ###########################################################################

    # MBFV - add a new shell to the mesh
    # add a new root vertex in the mesh (only needed once)
    def mbfv(self,vgx = 0,fgx = 0):
        rv = self.avert(vgx)
        rf = self.aface(rv,fgx)
        fg = self.afgroup(rf)
        sh = self.ashell(fg)
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
            raise NotImplemented

        return nv,ne

    # MEV - split an edge 
    # add an edge and a vertex
    #   u,v are the start,end vertices
    #   vgx is the geometric index of the new vertex
    #   egx is the geometric index of the new edge
    def esplit(self,u,v,vgx = 0,egx = 0):
        
        raise NotImplemented

    # MEV - split a vertex
    # add an edge and a vertex
    #   sv is the vertex to split
    #   be1,be2 are edges bounding the set of edges at sv
    #   vgx is the geometric index of the new vertex
    #   egx is the geometric index of the new edge
    def vsplit(self,sv,be1 = None,be2 = None,vgx = 0,egx = 0):

        raise NotImplemented

    # MEFL - split a face into two
    # add an edge and a face
    #   u,v are the start,end vertices
    #   egx,fgx are geometric indices
    #   e1 is the halfedge following a new halfedge in the old loop
    #   e2 is the halfedge following a new halfedge in the new loop
    #   the new loop is to the right of the halfedge u,v
    def mefl(self,u,v,e1 = None,e2 = None,egx = 0,fgx = 0):
        if e1 is None or e2 is None:
            uls = self.mask(2,u,None,None,None)
            vls = self.mask(2,v,None,None,None)
            uvls = tuple(x for x in uls if x in vls)
            if len(uvls) == 1:ol = uvls[0]
            else:
                print('loop ambiguity!')
                raise ValueError
            if e1 is None:
                lve = self.mask(6,v,None,ol)
                e1 = tuple(x for x in lve if x.one is v)[0]
            if e2 is None:
                lue = self.mask(6,u,None,ol)
                e2 = tuple(x for x in lue if x.one is u)[0]
        else:ol = self.loops[e1.l]
        ne = self.aedge(u,v,egx)
        nl = self.aloop(ne,b = 1)
        nf = self.aface(u,fgx)
        self.frings[nf].append(nl)
        olhe = halfedge(u,v,ol[2])
        nlhe = nl[1]
        olhe.ops = nlhe
        nlhe.ops = olhe
        olhe.nxt = e1;e1.lst = olhe
        nlhe.nxt = e2;e2.lst = nlhe
        c = e2
        while not c.nxt is e1:
            c.l = nl[2]
            c = c.nxt
        c.nxt = nlhe;nlhe.lst = c
        c = e1
        while not c.nxt is e2:
            c.l = ol[2]
            c = c.nxt
        c.nxt = olhe;olhe.lst = c
        return ne,nf

    # MEFL - slice an edge 
    # add an edge and a face
    #   u,v are the start,end vertices of the edge to slice
    #   egx is the geometric index of the new edge
    #   fgx is the geometric index of the new face
    def eslice(self,u,v,egx = 0,fgx = 0):
        
        raise NotImplemented

    # MEKL - connect two loops in the same face
    # add an edge and kill a hole loop in in a face
    #   u,v are the start,end vertices
    #   egx,fgx are geometric indices
    #   e1 is the halfedge following a new halfedge in the old loop
    #   e2 is the halfedge following a new halfedge in the new loop
    #   the new loop is to the right of the halfedge u,v
    def mekl(self,u,v,e1 = None,e2 = None,egx = 0,fgx = 0):

        raise NotImplemented

    # MEKBFL - connect two faces from two disconnected bodies
    # add an edge and kill a body
    #   u,v are the start,end vertices
    #   egx,fgx are geometric indices
    #   e1 is the halfedge following a new halfedge in the old loop
    #   e2 is the halfedge following a new halfedge in the new loop
    #   the new loop is to the right of the halfedge u,v
    def mekbfl(self,u,v,e1 = None,e2 = None,egx = 0,fgx = 0):

        raise NotImplemented

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

    ###################################################
    ###################################################





 




