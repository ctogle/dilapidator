import dilap.core.base as db
import dilap.core.tools as dpr
import dilap.core.vector as dpv
import dilap.core.pointset as dps
#import dilap.core.graphnode as gnd
#import dilap.core.graphedge as geg

import dilap.mesh.tools as dtl

import matplotlib.pyplot as plt
import pdb



class geometry(db.base):

    def radius(self):
        pmags = []
        for x in range(self.points.pcnt):
            pmags.append(self.points.ps[x].magnitude())
        return max(pmags)
    
    def __init__(self,*args,**kwargs):
        self.points = dps.pointset()
        self.curves = []
        self.surfaces = []

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

class loop(db.base):

    def __init__(self,index,*edges,**kwargs):
        self.index = index
        self.edges = edges

class face(db.base):

    def __init__(self,index,*loops,**kwargs):
        self.index = index
        self.loops = loops
        self.fnorm = None

class body(db.base):

    def __init__(self,index,*faces,**kwargs):
        self.index = index
        self.faces = faces

class hole(db.base):

    def __init__(self,index,body,**kwargs):
        self.index = index
        self.body = body

class mesh(db.base):

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
            if not e is None:results.extend([e.one,e.two])
            if not l is None:
                for ve in self.mask(1,None,None,l):
                    for vv in self.mask(0,None,ve,None):
                        if not vv in results:results.append(vv)
        elif d == 1:
            if not v is None:results.extend(v.ring)
            if not e is None:
                for vv in self.mask(0,None,e,None):
                    for ve in self.mask(1,vv,None,None):
                        if not ve is e:results.append(ve)
            if not l is None:results.extend(l.edges)
        elif d == 2:
            if not v is None:
                for ve in self.mask(1,v,None,None):
                    for vl in self.mask(2,None,ve,None):
                        if not vl in results:results.append(vl)
            if not e is None:results.extend(e.ring)
            if not l is None:
                for ve in self.mask(1,None,None,l):
                    for vl in self.mask(2,None,ve,None):
                        if not vl is l:results.append(vl)
        return results

    def __init__(self,index,owner,**kwargs):
        self.index = index
        self.owner = owner
        self.nodes = []
        self.edges = []
        self.loops = []
        self.faces = []
        self.bdies = []
        self.holes = []

        self.nlook = {}
        self.elook = {}

        self.nodecount = 0
        self.edgecount = 0
        self.loopcount = 0
        self.facecount = 0
        self.bodycount = 0
        self.holecount = 0

    def topo(self):
        return (self.nodecount,self.edgecount,self.loopcount,
                self.facecount,self.bodycount,self.holecount)

    def eulerstate(self):
        v,e,l,f,s,g = self.topo()
        state = v-e+f-(l-f)-2*(s-g)
        return state

    # make-body-face-loop-vertex
    def mbflv(self,nx):
        newnode = self.add_node(nx)
        newloop = self.add_loop()
        newface = self.add_face(newloop)
        newbody = self.add_body(newface)
        return newface,newloop,newnode

    # make-edge-vertex
    # create an edge starting at e with angle a
    def mev(self,v,e = None,a = None):
        dp = dpv.vector(0,0,0)
        np = self.owner.geom.points.ps[v.index].copy().translate(dp)
        nx = self.owner.geom.new_point(np)
        newnode = self.add_node(nx)
        newedge = self.add_edge(v,newnode)
        return newedge,newnode

    # make-edge
    #
    def me(self,v1,v2,e1 = None,a1 = None,e2 = None,a2 = None):
        # call either mekl,mefl,or mekbfl
        #newedge = self.add_edge(v1,v2)
        #return newedge.index
        raise NotImplemented
    
    # make-edge-kill-loop
    def mekl(self):raise NotImplemented
    # make-edge-face-loop
    def mefl(self):raise NotImplemented
    # make-edge-kill-body-face-loop
    def mekbfl(self):raise NotImplemented


    def glue(self):raise NotImplemented

    def dkflev(self):raise NotImplemented
    def kev(self):raise NotImplemented
    def ke(self):raise NotImplemented
    def unglue(self):raise NotImplemented

    def mme(self):raise NotImplemented
    def esplit(self):raise NotImplemented
    def kve(self):raise NotImplemented

    def lmove(self):raise NotImplemented

    def add_node(self,npx,**kwargs):
        if npx in self.nlook:return self.nlook[npx]
        newnode = node(self.nodecount,npx,**kwargs)
        self.nodes.append(newnode)
        self.nlook[newnode.px] = self.nodecount
        self.nodecount += 1
        return newnode

    def add_edge(self,nd1,nd2,**kwargs):
        if type(nd1) == type(1):nd1 = self.nodes[nd1]
        if type(nd2) == type(1):nd2 = self.nodes[nd2]
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

    def add_loop(self,*edges,**kwargs):
        newloop = loop(self.loopcount,edges,**kwargs)
        for eg in edges:eg.connect(newloop)
        self.loops.append(newloop)
        self.loopcount += 1
        return newloop

    def add_face(self,*loops,**kwargs):
        newface = face(self.facecount,loops,**kwargs)
        #for lp in loops:lp.connect(newface)
        self.faces.append(newface)
        self.facecount += 1
        return newface

    def add_body(self,*faces,**kwargs):
        newbody = body(self.bodycount,faces,**kwargs)
        #for lp in loops:lp.connect(newface)
        self.bdies.append(newbody)
        self.bodycount += 1
        return newbody

    def delete_edge(self,ex):
        eg = self.edges[ex]
        nd1,nd2 = eg.one,eg.two
        nd1.disconnect(nd2)
        nd2.disconnect(nd1)
        ekey1 = (nd1.index,nd2.index)
        ekey2 = (nd2.index,nd1.index)
        self.edges[eg.index] = None
        del self.elook[ekey1]
        del self.elook[ekey2]

    def replace_edge(self,rx,*exs):
        ering = self.edges[rx].ring
        self.delete_edge(rx)
        for ex in exs:
            for er in ering:
                self.edges[ex].connect(er)

class brep(db.base):

    def plot(self,ax = None):
        if ax is None:ax = dtl.plot_axes()
        for n in self.mesh.nodes:
            dtl.plot_point(self.geom.points.ps[n.px],ax)
        for e in self.mesh.edges:
            if e is None:continue
            eps = self.geom.points.ps[e.one.px],self.geom.points.ps[e.two.px]
            dtl.plot_edges(eps,ax)
        r = self.geom.radius()
        ax.set_xlim([-r,r])
        ax.set_ylim([-r,r])
        ax.set_zlim([-r,r])
        return ax

    def __init__(self,*args,**kwargs):
        self.geom = geometry()
        self.mesh = mesh(0,self)

    def new_body(self,np):
        nx = self.geom.points.new_point(np)
        self.mesh.mbflv(nx)

    def add_node(self,p):
        x = self.points.new_point(p)
        self.mesh.add_node(x)

    def add_edge(self,p1,p2):
        x1,x2 = self.points.new_points(p1,p2)
        self.mesh.add_edge(x1,x2)

    def add_face(self,*ps):
        pxs = self.geom.points.new_points(*ps)
        pxs = [self.mesh.add_node(pxs[x]) for x in range(len(pxs))]
        exs = [self.mesh.add_edge(pxs[x-1],pxs[x]) for x in range(len(pxs))]
        self.mesh.add_face(*exs)

    def split_edge(self,ex):
        e = self.mesh.edges[ex]
        nds = self.mesh.mask(0,None,e,None)
        mp = dpv.midpoint(*self.geom.points.get_points(nds[0].px,nds[1].px))
        mx = self.mesh.add_node(self.geom.points.add_point(mp))
        e1 = self.mesh.add_edge(nds[0].index,mx)
        e2 = self.mesh.add_edge(nds[1].index,mx)
        self.mesh.replace_edge(e.index,e1,e2)

    def model(self):
        pdb.set_trace()

    def fun(self):


        self.add_face(
            dpv.vector( 0, 0,0),dpv.vector(10, 0,0),
            dpv.vector(10,10,0),dpv.vector( 0,10,0))
        self.add_face(
            dpv.vector(10, 0,0),dpv.vector(20, 0,0),
            dpv.vector(20,10,0),dpv.vector(10,10,0))
        self.geom.points.ps[1].translate_z(10)
        self.geom.points.ps[2].translate_z(10)
        self.split_edge(2)
        self.geom.points.ps[6].translate_z(2)
        return

def btest():
    br = brep()

    nbp = dpv.vector(1,1,0)
    br.new_body(nbp).mev(0)

    #br.fun()

    br.plot()
    plt.show()
























class graph(db.base):

    def plot(self,ax = None):
        if ax is None:ax = dtl.plot_axes()
        for n in self.nodes:
            if not n is None:n.plot(ax)
        for e in self.edges:
            if not e is None:e.plot(ax)
        for f in self.faces:
            if not f is None:f.plot(ax)
        r = self.radius()
        ax.set_xlim([-r,r])
        ax.set_ylim([-r,r])
        ax.set_zlim([-r,r])
        return ax

    def radius(self):
        pmags = []
        for nx in range(self.nodecount):
            nd = self.nodes[nx]
            if nd is None:continue
            pmags.append(nd.p.magnitude())
        return max(pmags)

    def euler_poincare(self):
        #V – E + F = H + 2 * (S – G)
        left = len(self.nodes)-len(self.edges)+len(self.faces)
        right = 0 + 2 * (1 - 0)
        return left == right

    nodeclass = node
    edgeclass = edge
    faceclass = face
    def __init__(self,**kwargs):
        self._def('index',None,**kwargs)
        self.nodes = []
        self.nodes_lookup = {}
        self.edges = []
        self.edges_lookup = {}
        self.faces = []
        self.faces_lookup = {}
        self.nodecount = 0
        self.edgecount = 0
        self.facecount = 0

    # add a new node to the graph if no existing node exists
    # ndkey is a tuple(x,y,z)
    def add_node(self,ndkey,**kwargs):
        if ndkey in self.nodes_lookup:
            nd = self.nodes[self.nodes_lookup[ndkey]]
            if not nd is None:return nd.index
        nx,ny,nz = ndkey
        newnode = self.nodeclass(self.nodecount,
            dpv.vector(nx,ny,nz),**kwargs)
        self.nodes.append(newnode)
        self.nodes_lookup[newnode.key()] = newnode.index
        self.nodecount += 1
        return newnode.index

    # delete an existing node from the graph
    def del_node(self,ndkey):
        if ndkey in self.nodes_lookup:
            nd = self.nodes[self.nodes_lookup[ndkey]]
            if nd is None:return
        for ekey in self.edges_lookup:
            if nd.index in ekey:
                self.del_edge(*ekey)
        self.nodes[nd.index] = None
        del self.nodes_lookup[nd.key()]

    # add a new edge to the graph, or return existing index
    def add_edge(self,ndkey1,ndkey2,**kwargs):
        if ndkey1 in self.nodes_lookup:
            nd1 = self.nodes[self.nodes_lookup[ndkey1]]
        else:nd1 = self.nodes[self.add_node(ndkey1)]
        if ndkey2 in self.nodes_lookup:
            nd2 = self.nodes[self.nodes_lookup[ndkey2]]
        else:nd2 = self.nodes[self.add_node(ndkey2)]
        egkey = (ndkey1,ndkey2)
        if egkey in self.edges_lookup:
            eg = self.edges[self.edges_lookup[egkey]]
            if not eg is None:return eg.index
        egkey = (ndkey2,ndkey1)
        if egkey in self.edges_lookup:
            eg = self.edges[self.edges_lookup[egkey]]
            if not eg is None:return eg.index
        newedge = self.edgeclass(self.edgecount,nd1,nd2,**kwargs)
        self.edges_lookup[(ndkey1,ndkey2)] = newedge.index
        self.edges_lookup[(ndkey2,ndkey1)] = newedge.index
        self.edges.append(newedge)
        self.edgecount += 1
        return newedge.index

    # delete an existing edge from the graph
    def del_edge(self,ndkey1,ndkey2):
        ekey = (ndkey1,ndkey2)
        if not ekey in self.edges_lookup:return
        edge = self.edges[self.edges_lookup[ekey]]
        del self.edges_lookup[ekey]
        del self.edges_lookup[ekey[::-1]]
        self.edges[edge.index] = None

    # add a new face to the graph if no existing face exists
    # fkey is a tuple(nkey1,...,nkeyN)
    def add_face(self,fkey,**kwargs):
        pdb.set_trace()

    # delete an existing face from the graph
    def del_face(self,fkey):
        pdb.set_trace()

def test():
    g = graph()
    g.add_edge(( 0, 0,0),(10, 10,0))
    g.add_edge(( 0, 0,0),(10,-10,0))
    g.add_edge((10,10,0),( 0,-10,0))

    g.plot()
    plt.show()









    '''#
    def get_nodes(self,nkeys):
        nds = []
        for nk in nkeys:
            nd = self.nodes[self.nodes_lookup[nk]]
            nds.append(nd)
        return nds

    def get_node_points(self,nkeys):
        nps = []
        for nk in nkeys:
            np = self.nodes[self.nodes_lookup[nk]].p
            nps.append(np.copy())
        return nps

    # add a new edges to the graph, return indicies
    def _add_edges(self,ndkeys,**kwargs):
        edgexs = []
        for kdx in range(1,len(ndkeys)):
            ndkey1,ndkey2 = ndkeys[kdx-1],ndkeys[kdx]
            edgexs.append(self._add_edge(ndkey1,ndkey2,**kwargs))
        return edgexs

    # remove existing edge from ndx1 to ndx2
    # add two new edges, connecting n to ndx1 and to ndx2
    def _split_edge(self,ndx1,ndx2,newndx,**kwargs):
        sekey = self.edges_lookup[(ndx1,ndx2)]
        if not sekey is None:
            sedge = self.edges[sekey]
            kwargs['interpolated'] = sedge.interpolated
        self._del_edge(ndx1,ndx2)
        nd1,nd2 = self.nodes[ndx1],self.nodes[newndx]

        if nd1.p.near(nd2.p):pdb.set_trace()

        newedge = edge(nd1,nd2,**kwargs)
        self._add_edge(newedge)
        nd1,nd2 = self.nodes[ndx2],self.nodes[newndx]

        if nd1.p.near(nd2.p):pdb.set_trace()

        newedge = edge(nd1,nd2,**kwargs)
        self._add_edge(newedge)
    '''#





