import dilap.core.base as db

import dilap.core.tools as dpr
import dilap.core.mesh.tools as dtl
import dilap.core.lsystem as dls

import dp_vector as dpv

import matplotlib.pyplot as plt
import random as rm
import pdb



class infralsystem(dls.lsystem):

    def _realize(self,p,d):
        self.edges = []
        self.nodes = []
        return dls.lsystem._realize(self,p,d)

    loadouts = []
    loadouts.append(('FX',[
        ('X','X/YF/'),
        ('Y','\FX\Y')]))

    def __init__(self,ldx,*args,**kwargs):
        axiom,rules = self.loadouts[ldx]
        self._def('axiom',axiom,**kwargs)
        self._def('rules',rules,**kwargs)

        self._def('seed',1,**kwargs) # seed used for random numbers
        self._def('iterations',5,**kwargs) # number of iterations

        self._def('rho',100.0,**kwargs) # the length of a new edge
        self._def('angle',dpr.PI/2.0,**kwargs) # angle used for rotations
        self._def('minangle',dpr.PI/12.0,**kwargs) # angle used for rotations
        self._def('maxangle',dpr.PI*2.0,**kwargs) # angle used for rotations

        self._def('branchdraw',self.draw_branch,**kwargs)
        self._def('leafdraw',self.draw_leaf,**kwargs)
        self._def('finaldraw',self.draw,**kwargs)
        dls.lsystem.__init__(self,*args,**kwargs)

    def draw_branch(self,p,n):
        self.edges.append((p.copy(),n.copy()))
        dls.draw_branch(p,n)

    def draw_leaf(self,p):
        self.nodes.append((p.copy()))
        dls.draw_leaf(p)

    def draw(self):
        es,ls = self.edges,self.nodes
        nps = []
        nes = []
        for ex in range(len(es)):
            e1,e2 = es[ex]
            e1x,e2x = None,None
            for pdx in range(len(nps)):
                if e1.near(nps[pdx]):e1x = pdx
                if e2.near(nps[pdx]):e2x = pdx
                if not e1x is None and not e2x is None:break
            if e1x is None:
                e1x = len(nps)
                nps.append(e1)
            if e2x is None:
                e2x = len(nps)
                nps.append(e2)
            nes.append((e1x,e2x))
        self.nodes,self.edges = nps,nes
        dls.draw()

class node(db.base):

    # d is the direction to o.p from self.p
    def connect(self,d,o):
        self.targetring[o.index] = d
        self.ring[o.index] = d

    def plot_xy(self,ax = None):
        if ax is None:ax = dtl.plot_axes_xy()
        dtl.plot_point_xy(self.p,ax)
        for d in self.ring.keys():
            f = self.p.copy().translate(self.spikes[d])
            dtl.plot_edges_xy([self.p,f],ax)
        return ax

    def __init__(self,p,**kwargs):
        self.index = None
        self.p = p
        self.targetring = {}
        self.ring = {}
        self.spikes = {}

    def _spikes_one(self):
        tkys = list(self.targetring.keys())
        k = tkys[0]
        self.ring[k] = self.targetring[k]
        spike = dpv.xhat.copy().rotate_z(dpr.rad(self.ring[k])).scale_u(5)
        self.spikes[k] = spike

    def _spikes_two(self):
        tkys = list(self.targetring.keys())
        rkys = (self.targetring[tkys[0]],self.targetring[tkys[1]])
        da,sda = adist(*rkys),signedadist(*rkys)
        err,mpt = (180.0-da)/2.0,(rkys[0]+rkys[1])/2.0
        mpt = dpr.clamp_periodic(mpt + 180,0,360) if sda > 180 else mpt
        nrk1 = apply_err(rkys[0],err,mpt)
        nrk2 = apply_err(rkys[1],err,mpt)
        self.ring[tkys[0]] = nrk1
        self.ring[tkys[1]] = nrk2
        spike1 = dpv.xhat.copy().rotate_z(dpr.rad(nrk1)).scale_u(5)
        spike2 = dpv.xhat.copy().rotate_z(dpr.rad(nrk2)).scale_u(5)
        self.spikes[tkys[0]] = spike1
        self.spikes[tkys[1]] = spike2

    def _spikes_three(self):
        tkys = list(self.targetring.keys())
        tk1,tk2,tk3 = tkys
        ta1 = self.targetring[tk1]
        ta2 = self.targetring[tk2]
        ta3 = self.targetring[tk3]
        da1,da2,da3 = adist(ta1,ta2),adist(ta2,ta3),adist(ta3,ta1)

        for x in range(10):ta1,ta2,ta3 = threeway_nudge(ta1,ta2,ta3,0.5)
        self.ring[tk1] = ta1
        self.ring[tk2] = ta2
        self.ring[tk3] = ta3
        for d in self.targetring:
            spike = dpv.xhat.copy().rotate_z(dpr.rad(self.ring[d])).scale_u(5)
            self.spikes[d] = spike

    def _spikes_four(self):
        tkys = list(self.targetring.keys())
        tk1,tk2,tk3,tk4 = tkys
        ta1 = self.targetring[tk1]
        ta2 = self.targetring[tk2]
        ta3 = self.targetring[tk3]
        ta4 = self.targetring[tk4]

        print('calc spikes:',tkys)
        for d in self.targetring:
            self.ring[d] = self.targetring[d]
            spike = dpv.xhat.copy().rotate_z(dpr.rad(self.ring[d])).scale_u(5)
            self.spikes[d] = spike

        #pdb.set_trace()

    # update the actual positions of spikes based on targetring
    def _spikes(self):
        self.ring,self.spikes = {},{}
        tkys = list(self.targetring.keys())
        rcnt = len(tkys)
        if rcnt == 0:return
        elif rcnt == 1:self._spikes_one()
        elif rcnt == 2:self._spikes_two()
        elif rcnt == 3:self._spikes_three()
        elif rcnt == 4:self._spikes_four()
        elif rcnt > 4:
            print("WHATT")
            for d in self.targetring:
                self.ring[d] = self.targetring[d]
                spike = dpv.xhat.copy().rotate_z(dpr.rad(self.ring[d])).scale_u(5)
                self.spikes[d] = spike

class edge(db.base):
    
    def plot_xy(self,ax = None):
        if ax is None:ax = dtl.plot_axes_xy()
        np1 = self.one.p
        np2 = self.two.p
        dtl.plot_edges_xy([np1,np2],ax)
        dtl.plot_edges_xy(self.rpts,ax)
        return ax

    def __init__(self,node1,node2,**kwargs):
        self.index = None
        self.one = node1
        self.two = node2
        self._elaborate()

    def _elaborate(self,snapping = None):
        tang = dpv.v1_v2(self.one.p,self.two.p)
        self.tangent = tang

        ndir = dpr.deg(dpv.angle_from_xaxis_xy(tang))
        ndirf = ndir+180 if ndir < 180 else ndir - 180

        self.one.connect(ndir ,self.two)
        self.two.connect(ndirf,self.one)

    # use vector spline to add road points to plot!!
    def _place_road(self):
        r1 = self.one.p.copy()
        r2 = self.one.p.copy().translate(self.one.spikes[self.two.index])
        r3 = self.two.p.copy().translate(self.two.spikes[self.one.index])
        r4 = self.two.p.copy()
        self.rpts = dpv.vector_spline(r1,r2,r3,r4,self.tangent.magnitude()/5.0)

class graph(db.base):
    
    def plot_xy(self,ax = None):
        if ax is None:ax = dtl.plot_axes_xy()
        for n in self.nodes:n.plot_xy(ax)
        for eg in self.edges:eg.plot_xy(ax)
        return ax

    def __init__(self,**kwargs):
        self.nodes = []
        self.edges = []
        self.edges_lookup = {}
        self.nodecount = 0
        self.edgecount = 0

    # add a new node to the graph
    def _node(self,n):
        n.index = self.nodecount
        self.nodes.append(n)
        self.nodecount += 1
        return self.nodecount-1


    # ADD THE ABILITY TO REMOVE NODES/EDGES
    # ADD THE ABILITY TO SPLIT A PAIR OF INTERSECTING EDGES
    # ADD THE ABILITY TO APLIT AN EDGE WITH A NODE


    # add a new edge to the graph, or return existing index
    def _edge(self,n1,n2,**kwargs):
        nd1 = self.nodes[n1]
        nd2 = self.nodes[n2]
        existing = self._find_edge(n1,n2)
        if not existing is None:return existing
        newedge = edge(nd1,nd2,**kwargs)
        newedge.index = self.edgecount
        self.edges_lookup[(n1,n2)] = self.edgecount
        self.edges_lookup[(n2,n1)] = self.edgecount
        self.edges.append(newedge)
        self.edgecount += 1
        return self.edgecount-1

    # return index of closest node to p within e, or None
    def _find_node(self,p,e):
        nps = [nd.p for nd in self.nodes]
        ndx = dpv.find_closest(p,nps,self.nodecount,1.0)
        if self.nodes[ndx].p.neighborhood(p,e):return ndx

    # return index of closest node within a cone,
    # cone is pointed from o towards p, 
    # has halfangle e, and ends at p, 
    # or return None if no node exists
    def _find_node_cone(self,o,p,e):
        ca = dpr.deg(dpv.angle_from_xaxis_xy(dpv.v1_v2(o,p).normalize()))
        #best,margin = None,dpv.distance(o,p)
        best,margin = None,100000000000000000
        for ndx in range(self.nodecount):
            nd = self.nodes[ndx]
            tn = dpv.v1_v2(o,nd.p).normalize()
            npa = dpr.deg(dpv.angle_from_xaxis_xy(tn))
            if adist(ca,npa) < e:
                ndd = dpv.distance(o,nd.p)
                if ndd < margin:
                    best = ndx
                    margin = ndd
        return best

    # return index of edge within connecting ndx1,ndx2, or None
    def _find_edge(self,ndx1,ndx2):
        if (ndx1,ndx2) in self.edges_lookup:
            return self.edges_lookup[(ndx1,ndx2)]

    # which is the index of a node, delta is the translation
    # if the new node is within epsilon of an existing node,
    # simply connect which and existing node
    def _extrude(self,which,delta,epsilon = 22.5):
        oldp = self.nodes[which].p
        newp = oldp.copy().translate(delta)
        #exnd = self._find_node(newp,epsilon)
        exnd = self._find_node_cone(oldp,newp,epsilon)
        if exnd:newx = exnd
        else:newx = self._node(node(newp))
        edgx = self._edge(which,newx)

    # verify graph is correct where possible
    def _update(self):
        for nd in self.nodes:nd._spikes()
        for eg in self.edges:eg._place_road()

def adist(a1,a2):
    da = dpr.clamp_periodic(a1-a2,0,360)
    return da if da < 180 else 360 - da

def signedadist(a1,a2):
    return a1-a2 if a1 > a2 else a2-a1

def apply_err(x,e,m):
    xpe = dpr.clamp_periodic(x+e,0,360)
    xme = dpr.clamp_periodic(x-e,0,360)
    da1,da2 = adist(xpe,m),adist(xme,m)
    if da1 > da2:return xpe
    else:return xme

def threeway_nudge(t1,t2,t3,f):
    nt1,nt2,nt3 = t1,t2,t3
    d1,d2,d3 = adist(t1,t2),adist(t2,t3),adist(t3,t1)
    e = f*dpr.clamp(90-d1,0,90)/2.0
    if e:t1,t2,t3 = apply_err(t1,e,t2),apply_err(t2,e,t1),t3
    e = f*dpr.clamp(90-d2,0,90)/2.0
    if e:t1,t2,t3 = t1,apply_err(t2,e,t3),apply_err(t3,e,t2)
    e = f*dpr.clamp(90-d3,0,90)/2.0
    if e:t1,t2,t3 = apply_err(t1,e,t3),t2,apply_err(t3,e,t1)
    return t1,t2,t3

def lsystem_graph():
    g = graph()

    p = dpv.zero()
    d = dpv.xhat.copy()
    lsys = infralsystem(0)._realize(p,d)
    ndps,edgs = lsys.nodes,lsys.edges

    ndxs = []
    nexs = []
    for ndp in ndps:ndxs.append(g._node(node(ndp)))
    for edg in edgs:nexs.append(g._edge(ndxs[edg[0]],ndxs[edg[1]]))
    return g

def spine_test():
    g = graph()

    sd = dpv.zero()
    g._node(node(sd))
    for x in range(5):
        nndel = dpv.vector(25+2*x,50-5*x,-10)
        g._extrude(x,nndel)

    for x in range(5):
        nm = g.edges[x].tangent.copy().rotate_z(dpr.rad(90))
        nm.normalize().scale_u(100)
        an1 = node(g.nodes[x].p.copy().translate(nm))
        an2 = node(g.nodes[x].p.copy().translate(nm.flip()))
        g._node(an1)
        g._node(an2)
        g._edge(x,g.nodecount-1)
        g._edge(x,g.nodecount-2)

    return g

def alt_test():
    g = graph()

    g._node(node(dpv.vector(-100, 0,0)))
    g._node(node(dpv.vector(100,  0,0)))
    g._node(node(dpv.vector(100,100,0)))
    g._node(node(dpv.vector(0  ,100,0)))
    g._edge(0,1)
    g._edge(1,2)
    g._edge(2,3)
    g._edge(3,0)
    g._node(node(dpv.vector(150,80,0)))
    g._edge(2,4)
    g._edge(4,1)
    g._extrude(0,dpv.vector(-100,-100,0))
    g._extrude(0,dpv.vector(-100,50,0))

    return g

def arc_test():
    g = graph()

    g._node(node(dpv.vector(0,0,0)))
    for x in range(5):
        nddel = dpv.vector(x*rm.choice(
                [-500,-250,250,500]),
                500+rm.choice([-100,100])*x,0)
        edx = g._extrude(x,nddel)

    return g

def smatter():
    g = graph()

    g._node(node(dpv.vector(0,0,0)))
    for x in range(100):
        ndx = rm.choice(range(g.nodecount))
        rcnt = len(g.nodes[ndx].ring)
        if rcnt == 0:
            ndd = dpv.xhat.copy()
            ndd.rotate_z(dpr.rad(rm.choice(range(360))))
            ndd.scale_u(rm.choice([100,200,300]))
        elif rcnt == 1:
            ndir = next(iter(g.nodes[ndx].ring.values()))
            nda = ndir+180 if ndir < 180 else ndir - 180
            ndd = dpv.xhat.copy().rotate_z(dpr.rad(nda))
            ndd.scale_u(rm.choice([100,200,300]))
        elif rcnt == 2:
            r1,r2 = tuple(g.nodes[ndx].ring.values())
            mpt = (r1+r2)/2.0
            nda = mpt+180 if mpt < 180 else mpt - 180
            ndd = dpv.xhat.copy().rotate_z(dpr.rad(nda))
            ndd.scale_u(rm.choice([100,200,300]))
        elif rcnt == 3:
            t1,t2,t3 = tuple(g.nodes[ndx].ring.values())
            d1,d2,d3 = adist(t1,t2),adist(t2,t3),adist(t3,t1)
            if   d1 > d2 and d1 > d3:nda = (t1+t2)/2.0
            elif d2 > d1 and d2 > d3:nda = (t2+t3)/2.0
            elif d3 > d1 and d3 > d2:nda = (t3+t1)/2.0
            ndd = dpv.xhat.copy().rotate_z(dpr.rad(nda))
            ndd.scale_u(rm.choice([100,200,300]))
        elif rcnt == 4:
            print('this node cannot be more connected!',ndx)
        g._extrude(ndx,ndd)

        #g._update()

        #ax = g.plot_xy()
        #ax.set_aspect('equal')
        #plt.show()

    return g

def generate_graph(plugs = []):
    #g = graph()

    g = smatter()
    
    #g = arc_test()
    #g = spine_test()
    #g = alt_test()
    #g = lsystem_graph()

    g._update()

    ax = g.plot_xy()
    ax.set_aspect('equal')
    plt.show()



