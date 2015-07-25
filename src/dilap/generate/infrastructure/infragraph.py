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

    # o is another node currently connected to self
    def disconnect(self,o):
        del self.targetring[o.index]
        del self.ring[o.index]

    def plot(self,ax = None):
        if ax is None:ax = dtl.plot_axes()
        plotp = self.p.copy()
        dtl.plot_point(plotp,ax)
        for d in self.ring.keys():
            f = plotp.copy().translate(self.spikes[d])
            dtl.plot_edges([plotp,f],ax)
        ax.set_xlim([-100,100])
        ax.set_ylim([-100,100])
        ax.set_zlim([-100,100])
        return ax

    def plot_xy(self,ax = None):
        if ax is None:ax = dtl.plot_axes_xy()
        dtl.plot_point_xy(self.p,ax)
        for d in self.ring.keys():
            f = self.p.copy().translate(self.spikes[d])
            dtl.plot_edges_xy([self.p,f],ax)
        return ax

    def __init__(self,p,**kwargs):
        self.index = None
        self._def('layer',0,**kwargs)
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

    def _spikes_two(self,graph,target = 180):
        ws = []
        angles = []
        tkeys = list(self.targetring.keys())
        for tk in tkeys:
            edx = graph.edges_lookup[(tk,self.index)]
            ws.append(1.0 if graph.edges[edx].interpolated else 0.0)
            angles.append(self.targetring[tk])

        iangles = nudge(angles,ws,target = target)
        for tk,na in zip(tkeys,iangles):
            self.ring[tk] = na
            spike = dpv.xhat.copy().rotate_z(dpr.rad(na)).scale_u(5)
            self.spikes[tk] = spike

        '''#
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
        '''#

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

        tas = nudge([ta1,ta2,ta3,ta4],[1.0,1.0,1.0,1.0])
        ta1,ta2,ta3,ta4 = tas

        print('calc spikes:',tkys)
        for d in self.targetring:
            self.ring[d] = self.targetring[d]
            spike = dpv.xhat.copy().rotate_z(dpr.rad(self.ring[d])).scale_u(5)
            self.spikes[d] = spike

        #pdb.set_trace()

    # update the actual positions of spikes based on targetring
    def _spikes(self,graph):
        self._retarget(graph)
        self.ring,self.spikes = {},{}
        tkys = list(self.targetring.keys())
        rcnt = len(tkys)
        if rcnt == 0:return
        elif rcnt == 1:self._spikes_one()
        elif rcnt == 2:self._spikes_two(graph,180)
        elif rcnt == 3:self._spikes_two(graph,90)
        #elif rcnt == 3:self._spikes_three()
        #elif rcnt == 4:self._spikes_four()
        elif rcnt == 4:self._spikes_two(graph,90)
        elif rcnt > 4:
            print("WHATT")
            for d in self.targetring:
                self.ring[d] = self.targetring[d]
                spike = dpv.xhat.copy().rotate_z(dpr.rad(self.ring[d])).scale_u(5)
                self.spikes[d] = spike

    # verify the directions in target ring and update 
    # tangents of any connecting edges
    def _retarget(self,graph):
        rgxs = list(self.ring.keys())
        for rgx in rgxs:

            rkey = (rgx,self.index)
            if not rkey in graph.edges_lookup:
                pdb.set_trace()

            regx = graph.edges_lookup[(rgx,self.index)]
            if regx is None:
                del self.targetring[rgx]
                del self.ring[rgx]
                continue

            reg = graph.edges[regx]
            ndir1,ndir2 = reg._elaborate()
            if not reg.one is self:ndir1,ndir2 = ndir2,ndir1
            self.targetring[rgx] = ndir1
            graph.nodes[rgx].targetring[self.index] = ndir2

class edge(db.base):
    
    def plot(self,ax = None):
        if ax is None:ax = dtl.plot_axes()
        np1 = self.one.p
        np2 = self.two.p
        dtl.plot_edges([np1,np2],ax)
        dtl.plot_edges(self.rpts,ax)
        return ax

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
        self._def('interpolated',True,**kwargs)

    def _elaborate(self):
        self.tangent = dpv.v1_v2(self.one.p,self.two.p)
        ndir = dpr.deg(dpv.angle_from_xaxis_xy(self.tangent))
        ndirf = ndir+180 if ndir < 180 else ndir - 180
        return ndir,ndirf

    # return list of layers which self and o have in common
    # otherwise return None
    def _layers_intersect(self,o):
        l1,l2,l3,l4 = self.one.layer,self.two.layer,o.one.layer,o.two.layer
        shared = []
        if not l1 in shared and (l1 == l3 or l1 == l4):shared.append(l1)
        if not l2 in shared and (l2 == l3 or l2 == l4):shared.append(l2)
        return shared

    # return pt of intersection if the tangents of self and o intersect
    # otherwise return None
    def _tangents_intersect(self,o):
        s1 = (self.one.p.copy().xy(),self.two.p.copy().xy())
        s2 = (o.one.p.copy().xy(),o.two.p.copy().xy())
        segisect = dtl.segments_intersect_at(s1,s2)
        return segisect 

    # use vector spline to add road points to plot!!
    def _place_road(self,graph):                       
        rcnt = int(self.tangent.magnitude()/5.0)
        if self.interpolated:
            r1 = self.one.p.copy()
            r2 = self.one.p.copy().translate(self.one.spikes[self.two.index])
            r3 = self.two.p.copy().translate(self.two.spikes[self.one.index])
            r4 = self.two.p.copy()
            self.rpts = dpv.vector_spline(r1,r2,r3,r4,rcnt)
        else:self.rpts = dpr.point_line(self.one.p,self.two.p,rcnt)

class graph(db.base):
    
    def plot(self,ax = None):
        if ax is None:ax = dtl.plot_axes()
        for n in self.nodes:
            if not n is None:
                n.plot(ax)
        for eg in self.edges:
            if not eg is None:
                eg.plot(ax)
        return ax

    def plot_xy(self,ax = None):
        if ax is None:ax = dtl.plot_axes_xy()
        for n in self.nodes:
            if not n is None:
                n.plot_xy(ax)
        for eg in self.edges:
            if not eg is None:
                eg.plot_xy(ax)
        return ax

    def __str__(self):
        st = '\tinfragraph with:\n\t'
        st += str(self._count_nodes())+'\tnodes\n\t'
        st += str(self._count_edges())+'\tedges\n\t'
        return st

    def __init__(self,**kwargs):
        self.nodes = []
        self.nodes_lookup = {}
        self.edges = []
        self.edges_lookup = {}
        self.nodecount = 0
        self.edgecount = 0

    def _count_nodes(self):
        ncnt = 0
        for x in range(self.nodecount):
            if self.nodes[x]:ncnt += 1
        return ncnt

    def _count_edges(self):
        ecnt = 0
        for x in range(self.edgecount):
            eg = self.edges[x]
            if self.edges[x]:ecnt += 1
        return ecnt

    # safely add a new node to the graph
    # if node n is within e of existing nodes, try to stack/merge
    def _node(self,n,e = 50):
        exnds = self._find_nodes(n.p,e)
        if len(exnds) == 0:return [self._add_node(n)]
        exnps = [self.nodes[x].p for x in exnds]
        exnls = [self.nodes[x].layer for x in exnds]
        exnps.append(n.p)
        exnls.append(n.layer)
        exnds.append(None)

        lnp = dpv.center_of_mass(exnps)
        exlys = {x:[] for x in exnls}
        for x in range(len(exnls)):
            if not exnds[x] is None:
                exlys[exnls[x]].append(exnds[x])

        newxs = [None for x in range(max(exlys)+1)]

        print('xlys',exlys)

        for l in exlys:
            ln = node(lnp.copy().xy().translate_z(20*l),layer = l)
            newx = self._add_node(ln)
            self._merge_nodes(exlys[l],newx)
            self.nodes_lookup[(lnp.x,lnp.y,l)] = newx
            newxs[l] = newx

        return newxs

    # add a new node to the graph
    def _add_node(self,n):
        n.index = self.nodecount
        self.nodes.append(n)
        self.nodes_lookup[(n.p.x,n.p.y,n.layer)] = n.index
        self.nodecount += 1
        return self.nodecount-1

    # delete an existing node from the graph
    def _del_node(self,ndx):
        n = self.nodes[ndx]
        for ekey in self.edges_lookup:
            if ndx in ekey:
                self._del_edge(*ekey)
        self.nodes[ndx] = None
        self.nodes_lookup[(n.p.x,n.p.y,n.layer)] = None

    # delete existing nodes from the graph and replace all
    # connectivity with new edges to a new node
    # ndxs is a list of existing nodes indices being merged
    # nndx is the index of the new node which replaces them
    def _merge_nodes(self,ndxs,nndx,**kwargs):
        mnds = [self.nodes[x] for x in ndxs]
        for ndx in ndxs:
            for ndrk in list(self.nodes[ndx].ring.keys()):
                ekey = (ndrk,ndx)
                if ekey in list(self.edges_lookup.keys()):
                    eg = self.edges_lookup[ekey]
                    if eg is None:continue
                    iterp = self.edges[eg].interpolated
                    self._del_edge(*ekey)
                    if not ndrk in ndxs:
                        nd1,nd2 = self.nodes[nndx],self.nodes[ndrk]
                        newedge = edge(nd1,nd2,interpolated = iterp)
                        #newedge = edge(nd1,nd2,**kwargs)
                        self._add_edge(newedge)
            self._del_node(ndx)

    # return index of closest node to p within e, or None
    def _find_node(self,p,e):
        nps = [nd.p for nd in self.nodes]
        ndx = dpv.find_closest(p,nps,self.nodecount,1.0)
        if self.nodes[ndx].p.neighborhood(p,e):return ndx

    # return indices of all nodes within e of p
    def _find_nodes(self,p,e):
        within = []
        for ndx in range(self.nodecount):
            nd = self.nodes[ndx]
            if nd is None:continue
            if nd.p.neighborhood(p,e):
                within.append(ndx)
        return within

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

    # cast a line from op to np + epsilon
    # if an edge is found, split
    # if a node is found within epsilon, return it
    def _find_node_cast(self,op,np,epsilon):

        pdb.set_trace()

    # THIS SHOULD NOT EXIST; IT DOES NOT MATCH WATCH OTHER FIND FUNCTIONS DO
    # merge all existing nodes within e of p into a new node at p
    #def _find_node_gather(self,p,e):
    #    merging = self._find_nodes(p,e)
    #    newx = self._node(node(p))
    #    self._merge_nodes(merging,newx)
    #    return newx

    #   CHECK   # ADD THE ABILITY TO REMOVE NODES/EDGES
    # ADD THE ABILITY TO SPLIT A PAIR OF INTERSECTING EDGES
    #   CHECK   # ADD THE ABILITY TO SPLIT AN EDGE WITH A NODE

# more than one type of intersection can occur:

# two edges can outright overlap, meaning that their straight 
# tangent vectors produce two overlapping segments

# two edges splined points might overlap whlie their
# tangent vectors do not produce overlapping segments

# can the former happen without the latter?

# a node should have a layer
# an edge is the layer that is the union of its nodes' layers
# two edges can overpass if they are in nonoverlapping layers
# two edges can meet at a new node if they each
#   exist in one layer and that layer is the same for both
#
# any other form of intersection should eliminate both edges??

    # add a new edge to the graph, or return existing index
    # this function should do this safely, so the resulting graph
    # does not carry improper intersections!
    # return None if the desired edge could not be created properly
    def _edge(self,nkey1,nkey2,**kwargs):
        n1 = self.nodes_lookup[nkey1]
        n2 = self.nodes_lookup[nkey2]
        existing = self._find_edge(n1,n2)
        if not existing is None:return existing
        else:
            nd1,nd2 = self.nodes[n1],self.nodes[n2]
            newedge = edge(nd1,nd2,**kwargs)
            ipts = []
            ilys = []
            for edx in range(self.edgecount):
                eg = self.edges[edx]
                if eg is None:continue
                ipt = eg._tangents_intersect(newedge)
                if not ipt is None:
                    ily = eg._layers_intersect(newedge)
                    ilycnt = len(ily)
                    if ilycnt > 1:raise ValueError

                    if type(ipt) is type(tuple()):
                        if ilycnt == 0:
                            continue
                        print('overlapping intersection!')
                        return None
                        #pdb.set_trace()

                    if ily:
                        l = ily[0]
                        ilys.append(l)
                    elif eg.one.layer == eg.two.layer:
                        l = eg.one.layer
                        if newedge.one.layer == newedge.two.layer:

                            iptndxxs = self._node(node(ipt,layer = newedge.one.layer))

                            print('shit: just like you said')
                            ilys.append(newedge.one.layer)
                        else:
                            print('shit: layer ambiguity')
                            pdb.set_trace()
                    else:
                        print('shit: layer ambiguity')
                        pdb.set_trace()

                    iptndxs = self._node(node(ipt,layer = l))
                    iptndx = iptndxs[l]

                    self._split_edge(eg.one.index,eg.two.index,iptndx)
                    ipts.append(ipt)

            if not ipts:return self._add_edge(newedge,**kwargs)

            newedgexs = []
            ipts.insert(0,nd1.p)
            ilys.insert(0,nd1.layer)
            ipts.append(nd2.p)
            ilys.append(nd2.layer)
            siptxs = dpv.proximity_order_xy(nd1.p,ipts)
            for iptx in range(1,len(ipts)):
                ipt1,ipt2 = ipts[siptxs[iptx-1]],ipts[siptxs[iptx]]
                ily1,ily2 = ilys[siptxs[iptx-1]],ilys[siptxs[iptx]]
                n1 = self.nodes_lookup[(ipt1.x,ipt1.y,ily1)]
                n2 = self.nodes_lookup[(ipt2.x,ipt2.y,ily2)]
                nd1,nd2 = self.nodes[n1],self.nodes[n2]

                print('illlys',ilys,ily1,ily2)
                #pdb.set_trace()

                newedge = edge(nd1,nd2,**kwargs)
                newedgexs.append(self._add_edge(newedge,**kwargs))
            return newedgexs

    def _connect_nodes(self,nx1,nx2,**kwargs):
        nd1,nd2 = self.nodes[nx1],self.nodes[nx2]
        nkey1 = (nd1.p.x,nd1.p.y,nd1.layer)
        nkey2 = (nd2.p.x,nd2.p.y,nd2.layer)
        return self._edge(nkey1,nkey2,**kwargs)

    # add a new edge to the graph, or return existing index
    def _add_edge(self,newedge,**kwargs):
        n1,n2 = newedge.one.index,newedge.two.index
        nd1,nd2 = self.nodes[n1],self.nodes[n2]
        newedge.index = self.edgecount
        ndir1,ndir2 = newedge._elaborate()
        newedge.one.connect(ndir1,newedge.two)
        newedge.two.connect(ndir2,newedge.one)
        self.edges_lookup[(n1,n2)] = self.edgecount
        self.edges_lookup[(n2,n1)] = self.edgecount
        self.edges.append(newedge)
        self.edgecount += 1
        return self.edgecount-1

    # delete an existing edge from the graph
    def _del_edge(self,n1,n2):
        edx = self.edges_lookup[(n1,n2)]
        if edx is None:return
        eg = self.edges[edx]
        eg.one.disconnect(eg.two)
        eg.two.disconnect(eg.one)
        self.edges_lookup[(n1,n2)] = None
        self.edges_lookup[(n2,n1)] = None
        self.edges[edx] = None

    # return index of edge within connecting ndx1,ndx2, or None
    def _find_edge(self,ndx1,ndx2):
        if (ndx1,ndx2) in self.edges_lookup:
            return self.edges_lookup[(ndx1,ndx2)]

    # remove existing edge from ndx1 to ndx2
    # add two new edges, connecting n to ndx1 and to ndx2
    def _split_edge(self,ndx1,ndx2,newndx,**kwargs):
        self._del_edge(ndx1,ndx2)


        nd1,nd2 = self.nodes[ndx1],self.nodes[newndx]

        if nd1.p.near(nd2.p):pdb.set_trace()

        newedge = edge(nd1,nd2,**kwargs)
        self._add_edge(newedge)


        nd1,nd2 = self.nodes[ndx2],self.nodes[newndx]

        if nd1.p.near(nd2.p):pdb.set_trace()

        newedge = edge(nd1,nd2,**kwargs)
        self._add_edge(newedge)

    # which is the index of a node, delta is the translation
    # if the new node is within epsilon of an existing node,
    # simply connect which and existing node
    def _extrude(self,which,delta,layer = 0):
        oldp = self.nodes[which].p
        newp = oldp.copy().translate(delta)
        newx = self._node(node(newp,layer = layer),e = 25)[layer]
        edgx = self._connect_nodes(which,newx)
        #edgx = self._edge(which,newx)


    def _extrude_safe(self,which,delta,epsilon = 22.5):
        oldp = self.nodes[which].p
        newp = oldp.copy().translate(delta)

        #exnd = self._find_node(newp,epsilon)
        exnd = self._find_node_cone(oldp,newp,epsilon)
        #exnd = self._find_node_cast(oldp,newp,epsilon)
        #####exnd = self._find_node_gather(newp,epsilon)

        if exnd:newx = exnd
        else:newx = self._node(node(newp))
        edgx = self._edge(which,newx)


    # verify graph is correct where possible
    def _update(self):
        for nd in self.nodes:
            if not nd is None:
                nd._spikes(self)
        for eg in self.edges:
            if not eg is None:
                eg._place_road(self)

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

def min_adist(alpha,angles,exempt = None):
    acnt = len(angles)
    if acnt < 2:return None,None
    else:
        minad = None
        minax = None
        for a2x in range(acnt):
            if a2x == exempt:continue
            ad = adist(angles[a2x],alpha)
            if minad is None or ad < minad:
                minad = ad
                minax = a2x
        return minad,minax

# given a list of angles, gradually move them to satisfy a condition
# the condition is that da for each angle be nearly equal to their avg
# da is the minimum distance to the nearest angle, which should be ~90
# ws are weights to affect the motion per nudge
def nudge(angles,ws,target = 90,error = 1):
    def acceptable(das):return abs(min(das)-target) < error
    def measure(angles):
        das = []
        dxs = []
        for x in range(acnt):
            xda,xdx = min_adist(angles[x],angles,x)
            das.append(xda)
            dxs.append(xdx)
        return das,dxs

    oas = angles[:]
    acnt = len(angles)
    das,dxs = measure(oas)

    tries = 0

    while not acceptable(das):

        tries += 1
        if tries > 100:
            pdb.set_trace()

        das,dxs = measure(oas)
        tdas = [dpr.clamp(target-da,0,target) for da in das]
        es = [tda*w*0.5 for tda,w in zip(tdas,ws)]
        oas = [apply_err(oas[x],es[x],oas[dxs[x]]) for x in range(acnt)]
    return oas

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

def del_test():
    g = graph()

    g._node(node(dpv.vector(0,0,0)))
    g._node(node(dpv.vector(100,50,0)))
    g._node(node(dpv.vector(100,100,0)))
    g._node(node(dpv.vector(50,100,0)))
    g._node(node(dpv.vector(50,50,0)))
    #g._edge(0,1)
    #g._edge(1,2)
    #g._edge(2,3)
    #g._edge(3,0)
    #g._edge(0,2)
    g._edge(0,4)
    g._edge(1,4)
    g._edge(2,4)
    g._edge(3,4)
    g._extrude(1,dpv.vector(0,-50,0))
    g._extrude(3,dpv.vector(-50,0,0))

    g._update()
    ax = g.plot_xy()
    ax.set_aspect('equal')
    plt.show()

    #g._del_node(1)
    #g._del_edge(0,2)
    #g._split_edge(0,2,4)

    #g._update()
    #ax = g.plot_xy()
    #ax.set_aspect('equal')
    #plt.show()

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
        #g._extrude_safe(ndx,ndd)
        g._extrude(ndx,ndd)

        #g._update()

        #ax = g.plot_xy()
        #ax.set_aspect('equal')
        #plt.show()

    return g

def ramp():
    g = graph()

    '''#
    p = dpv.vector(-100,-100,0)
    q = p.copy().translate_x(100)
    #q = p.copy().translate_x(100).translate_y(100)
    g._node(node(p))
    g._node(node(q))
    g._edge((p.x,p.y,0),(q.x,q.y,0),interpolated = False)
    #q = p.copy().translate_y(-100)

    p = q.copy()
    q = p.copy().translate_x(100).translate_y(100)
    g._node(node(p))
    g._node(node(q))
    g._edge((p.x,p.y,0),(q.x,q.y,0))

    p = q.copy()
    q = p.copy().translate_x(-100).translate_y(100)
    g._node(node(p))
    g._node(node(q,layer = 0))
    g._edge((p.x,p.y,0),(q.x,q.y,0))

    p = q.copy()
    q = p.copy().translate_x(-100).translate_y(-100)
    g._node(node(p,layer = 0))
    g._node(node(q,layer = 1))
    g._edge((p.x,p.y,0),(q.x,q.y,1))

    p = q.copy()
    q = p.copy().translate_x(100).translate_y(-100)
    g._node(node(p,layer = 1))
    g._node(node(q,layer = 1))
    g._edge((p.x,p.y,1),(q.x,q.y,1))

    p = q.copy()
    q = p.copy().translate_x(100)
    #q = p.copy().translate_x(100).translate_y(100)
    g._node(node(p,layer = 1))
    g._node(node(q,layer = 1))
    g._edge((p.x,p.y,1),(q.x,q.y,1),interpolated = False)

    p = q.copy()
    q = p.copy().translate_x(100)
    g._node(node(p,layer = 1))
    g._node(node(q,layer = 1))
    g._edge((p.x,p.y,1),(q.x,q.y,1))
    '''#



    #g._node(node(dpv.vector( 100,100,0)))
    g._node(node(dpv.vector(  0,100,0)))
    g._node(node(dpv.vector(200,100,0)))
    g._node(node(dpv.vector(100,  0,0),layer = 1))
    g._node(node(dpv.vector(100,200,0),layer = 1))
    g._edge((0,100,0),(200,100,0))
    g._edge((100,0,1),(100,200,1))
    g._edge((200,100,0),(100,200,1))
    g._node(node(dpv.vector(-100,300,0),layer = 1))
    g._edge((-100,300,1),(100,200,1))

    g._node(node(dpv.vector(   0,0,0)))
    g._node(node(dpv.vector(   0,0,0),layer = 1))
    g._node(node(dpv.vector(-100,0,0)))
    g._node(node(dpv.vector( 100,0,0)))
    g._node(node(dpv.vector(0,-100,0)))
    g._node(node(dpv.vector(0, 100,0)))

    g._edge((0,0,0),(-100,0,0))
    g._edge((0,0,0),( 100,0,0))
    g._edge((0,0,1),(0,-100,0))
    g._edge((0,0,1),(0, 100,0))

    #pdb.set_trace()

    g._update()

    print(g)

    ax = g.plot()
    plt.show()

    quit()
    return g

def mergetest():
    g = graph()
    
    g._node(node(dpv.vector(0,0,0)))
    g._extrude(0,dpv.vector(0,100,0))
    ndd = dpv.xhat.copy().scale_u(10)
    for x in range(4):
        ndd.rotate_z(dpr.rad(90))
        g._extrude(x+1,ndd,epsilon = 0.0)

    newp = dpv.vector(10,110,0)
    merging = g._find_nodes(newp,50)
    newx = g._node(node(newp))
    g._merge_nodes(merging,newx)

    return g




def generate_graph(plugs = []):
    #g = graph()

    #g = mergetest()
    #g = smatter()
    #g = del_test()
    
    #g = arc_test()
    #g = spine_test()
    #g = alt_test()
    #g = lsystem_graph()

    g = ramp()

    g._update()

    ax = g.plot()
    ax.set_aspect('equal')
    plt.show()



