import dilap.core.base as db
import dilap.core.vector as dpv
import dilap.core.tools as dpr
import dilap.core.lsystem as dls

import dilap.mesh.tools as dtl
import dilap.mesh.piecewisecomplex as pwc

import dilap.infrastructure.graphnode as gnd
import dilap.infrastructure.graphedge as geg
import dilap.infrastructure.graphregion as grg
import dilap.infrastructure.infralsystem as ifl

import matplotlib.pyplot as plt
import random as rm
import pdb



class graph(db.base):
    
    def plot_regions(self,ax = None):
        if ax is None:ax = dtl.plot_axes()

        pdb.set_trace()

    def plot(self,ax = None):
        if ax is None:ax = dtl.plot_axes()
        for n in self.nodes:
            if not n is None:
                n.plot(ax)
        for eg in self.edges:
            if not eg is None:
                eg.plot(ax)
        ax.set_xlim([-100,100])
        ax.set_ylim([-100,100])
        ax.set_zlim([-100,100])
        return ax

    def plot_xy(self,ax = None):
        if ax is None:ax = dtl.plot_axes_xy()
        for n in self.nodes:
            if not n is None:
                n.plot_xy(ax)
        for eg in self.edges:
            if not eg is None:
                eg.plot_xy(ax)
        ax.set_aspect('equal')
        return ax

    def __str__(self):
        st = '\tinfragraph with:\n\t'
        st += str(self._count_nodes())+'\tnodes\n\t'
        st += str(self._count_edges())+'\tedges\n\t'
        return st

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

    # verify graph is correct where possible
    def _update(self):
        for nd in self.nodes:
            if not nd is None:
                nd._spikes(self)
        for eg in self.edges:
            if not eg is None:
                eg._place_road(self)
        self._regions()

    def __init__(self,**kwargs):
        self.nodes = []
        self.nodes_lookup = {}
        self.edges = []
        self.edges_lookup = {}
        self.nodecount = 0
        self.edgecount = 0

    # given an edge e, direction 0 or 1, and cw or ccw
    # return the forward path of e
    def _loopwalk(self,e,d,w):
        if d:inpath = [e.one.key(),e.two.key()]
        else:inpath = [e.two.key(),e.one.key()]
        while True:
            ekey = (inpath[-2],inpath[-1])
            e = self.edges[self.edges_lookup[ekey]]
            nx = e._walk(inpath[-1],w)
            if nx is None:return inpath
            nxndkey = self.nodes[nx].key()
            if nxndkey in inpath:return inpath
            else:inpath.append(nxndkey)

    # return a collection of points outlining all edge loops in the graph
    def _edge_loops(self):
        edgelloops = []
        edgerloops = []
        edgestodo = self.edges[:]
        while edgestodo:
            e = edgestodo.pop(0)
            ewalkrcw  = self._loopwalk(e,1,1)
            ewalkrccw = self._loopwalk(e,0,0)
            ewalklccw = self._loopwalk(e,0,1)
            ewalklcw  = self._loopwalk(e,1,0)

            if set(ewalkrcw)  == set(ewalkrccw):
                print('closed loop!',len(edgestodo))
                rloop = tuple(ewalkrcw)
            else:
                print('unclosed loop!',len(edgestodo))
                pdb.set_trace()

            if set(ewalklccw) ==  set(ewalklcw):
                print('closed loop!',len(edgestodo))
                lloop = tuple(ewalklccw)
            else:
                print('unclosed loop!',len(edgestodo))
                pdb.set_trace()

            edgelloops.append(lloop)
            edgerloops.append(rloop)
        return edgelloops,edgerloops

    # eloop is a list of node keys which are connected in a loop by edges
    # side is either 0 (right) or 1 (left) relative to the first edge
    # in the loop - other edges must be handled carefully
    def _edge_loop_points(self,eloop,side):
        elcnt = len(eloop)
        looppts = []

        ne = self.edges[self.edges_lookup[eloop[0],eloop[1]]]
        if   side == 0:
            looppts.extend(ne.rbpts)
            lnkey = ne.two.key()
        elif side == 1:
            looppts.extend(ne.lbpts)
            lnkey = ne.one.key()

        for elx in range(2,elcnt+1):
            elx1,elx2 = elx-1,elx if elx < elcnt else 0
            ne = self.edges[self.edges_lookup[eloop[elx1],eloop[elx2]]]
            if   ne.one.key() == lnkey:
                if   side == 0:rbpts,lbpts = ne.rbpts,ne.lbpts
                elif side == 1:rbpts,lbpts = ne.lbpts,ne.rbpts
                lnkey = ne.two.key()
            elif ne.two.key() == lnkey:
                if   side == 0:rbpts,lbpts = ne.lbpts,ne.rbpts
                elif side == 1:rbpts,lbpts = ne.rbpts,ne.lbpts
                lnkey = ne.one.key()
            if   side == 0:looppts.extend(rbpts)
            elif side == 1:looppts.extend(lbpts)

        return looppts

    # return a collection of points outlining all nodes/edges in the graph
    def _edge_loop_boundaries(self):
        # a stem is a set of edges which are not part of a loop
        # every step exists within some loop (possibly implicit exterior loop)
        # a loop should be represented by a set of tuples which correspond 
        # to the edges which comprise the loop
        edgelloops,edgerloops = self._edge_loops()
        bedgeloops = {}
        for ex in range(len(edgelloops)):
            lloop,rloop = edgelloops[ex],edgerloops[ex]

            isperm = False
            for rps in bedgeloops:
                if dpr.cyclic_permutation(rloop,rps):isperm = True;break
            if not isperm:bedgeloops[rloop] = self._edge_loop_points(rloop,0)

            isperm = False
            for lps in bedgeloops:
                if dpr.cyclic_permutation(lloop,lps):isperm = True;break
            if not isperm:bedgeloops[lloop] = self._edge_loop_points(lloop,1)
        return bedgeloops

    # calculate polygons representing regions to place terrain
    def _regions(self):
        rpts = []
        for eg in self.edges:rpts.extend([x.copy() for x in eg.rpts])

        convexbnd = dpr.pts_to_convex_xy(rpts)
        convexbnd = dpr.inflate(convexbnd,100)
        eloops = self._edge_loop_boundaries()

        # rank the loops based on containment hierarchy to then 
        # describe polygons with holes using the loops
        '''#
        ax = dtl.plot_axes_xy()
        ax = self.plot_xy(ax)
        for bedge in eloops:
            ax = dtl.plot_edges_xy(eloops[bedge],ax)
        plt.show()
        '''#

        eloopkeys = eloops.keys()
        eloops1 = eloops[[x for x in eloopkeys][0]]
        eloops2 = eloops[[x for x in eloopkeys][1]]

        print('rpts',len(rpts))
        ax = dtl.plot_axes_xy()
        dtl.plot_edges_xy(convexbnd,ax)
        dtl.plot_edges_xy(eloops1,ax)
        dtl.plot_edges_xy(eloops2,ax)
        plt.show()

        eloops3 = [p.copy() for p in eloops2]
        #dpv.translate_coords_x(eloops3,200)
        #polygons = [(convexbnd,(eloops1,)),(eloops1,(eloops2,)),(eloops2,())]
        #polygons = [(convexbnd,(eloops1,))]
        #polygons = [(eloops2,())]

        print('amen')

        tpolygons = [(convexbnd,(eloops1,)),(eloops3,())]
        #tpolygons = [(convexbnd,(eloops1,))]
        #tpolygons = [(convexbnd,())]
        #tpolygons = [(eloops3,())]
        tplc = pwc.piecewise_linear_complex()
        tplc.add_polygons(*tpolygons)
        #tplc.triangulate_xy()
        #tplc.triangulate()
        self.tplc = tplc

        #another = dtl.icosphere()
        #another = dtl.box(10,20,50)
        print('amen2')

        rpolygons = [(eloops1,(eloops2,))]
        #rpolygons = [(eloops1,())]
        #rpolygons = [(eloops2,())]
        rplc = pwc.piecewise_linear_complex()
        rplc.add_polygons(*rpolygons)
        #rplc.triangulate()

        self.tplc = tplc
        self.rplc = rplc
        #self.rplc = another
        ax = dtl.plot_axes_xy()
        #ax = self.plot()
        ax = dtl.plot_polygon_xy(convexbnd,ax,True)
        ax = dtl.plot_polygon_xy(eloops2,ax,True)
        ax = dtl.plot_polygon_xy(eloops1,ax,True)
        ax = self.tplc.plot_xy(ax)
        ax = self.rplc.plot_xy(ax)
        #ax = another.plot_xy(ax)
        plt.show()

    # add a new node to the graph or existing node index
    # ndkey is a tuple(x,y,layer)
    def _add_node(self,ndkey):
        if ndkey in self.nodes_lookup:
            nd = self.nodes[self.nodes_lookup[ndkey]]
            if not nd is None:return nd.index
        nx,ny,nl = ndkey
        newnode = gnd.node(dpv.vector(nx,ny,20*nl),layer = nl)
        newnode.index = self.nodecount
        self.nodes.append(newnode)
        self.nodes_lookup[ndkey] = newnode.index
        self.nodecount += 1
        return newnode.index

    # delete an existing node from the graph
    def _del_node(self,ndkey):
        if ndkey in self.nodes_lookup:
            nd = self.nodes[self.nodes_lookup[ndkey]]
            if nd is None:return
        for ekey in self.edges_lookup:
            if nd.index in ekey:
                self._del_edge(*ekey)
        self.nodes[nd.index] = None
        del self.nodes_lookup[nd.key()]

    # add a new edge to the graph, or return existing index
    def _add_edge(self,ndkey1,ndkey2,**kwargs):
        if ndkey1 in self.nodes_lookup:
            nd1 = self.nodes[self.nodes_lookup[ndkey1]]
        else:nd1 = self.nodes[self._add_node(ndkey1)]
        if ndkey2 in self.nodes_lookup:
            nd2 = self.nodes[self.nodes_lookup[ndkey2]]
        else:nd2 = self.nodes[self._add_node(ndkey2)]
        newedge = geg.edge(nd1,nd2,**kwargs)
        newedge.index = self.edgecount
        ndir1,ndir2 = newedge._directions()
        newedge.one.connect(ndir1,newedge.two)
        newedge.two.connect(ndir2,newedge.one)
        self.edges_lookup[(ndkey1,ndkey2)] = newedge.index
        self.edges_lookup[(ndkey2,ndkey1)] = newedge.index
        self.edges.append(newedge)
        self.edgecount += 1
        return newedge.index

    # add a new edges to the graph, return indicies
    def _add_edges(self,ndkeys,**kwargs):
        edgexs = []
        for kdx in range(1,len(ndkeys)):
            ndkey1,ndkey2 = ndkeys[kdx-1],ndkeys[kdx]
            edgexs.append(self._add_edge(ndkey1,ndkey2,**kwargs))
        return edgexs

    # delete an existing edge from the graph
    def _del_edge(self,ndkey1,ndkey2):
        ekey = (ndkey1,ndkey2)
        if not ekey in self.edges_lookup:return
        edge = self.edges[self.edges_lookup[ekey]]
        edge.one.disconnect(edge.two)
        edge.two.disconnect(edge.one)
        del self.edges_lookup[ekey]
        del self.edges_lookup[ekey[::-1]]
        self.edges[edge.index] = None

    '''#
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
            lpp = lnp.copy().xy().translate_z(20*l)
            lky = (lpp.x,lpp.y,l)
            if lky in self.nodes_lookup:
                newx = self.nodes_lookup[lky]
                if self.nodes[newx] is None:
                    ln = node(lpp,layer = l)
                    newx = self._add_node(ln)
            else:
                ln = node(lpp,layer = l)
                newx = self._add_node(ln)
            self._merge_nodes(exlys[l],newx)
            self.nodes_lookup[(lnp.x,lnp.y,l)] = newx
            newxs[l] = newx

        return newxs
    '''#

    # delete existing nodes from the graph and replace all
    # connectivity with new edges to a new node
    # ndxs is a list of existing nodes indices being merged
    # nndx is the index of the new node which replaces them
    def _merge_nodes(self,ndxs,nndx,**kwargs):
        mnds = [self.nodes[x] for x in ndxs if not x == nndx]
        for ndx in ndxs:
            if ndx == nndx:continue
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

    # return index of edge within connecting ndx1,ndx2, or None
    def _find_edge(self,ndx1,ndx2):
        if (ndx1,ndx2) in self.edges_lookup:
            return self.edges_lookup[(ndx1,ndx2)]

    # remove existing edge from ndx1 to ndx2
    # add two new edges, connecting n to ndx1 and to ndx2
    def _split_edge(self,ndx1,ndx2,newndx,**kwargs):
        sekey = self.edges_lookup[(ndx1,ndx2)]
        if not sekey is None:
            sedge = self.edges[sekey]
            kwargs['interpolated'] = sedge.interpolated
        else:
            print('IM BULLSHITTING OVER HERE')
            return
        self._del_edge(ndx1,ndx2)


        nd1,nd2 = self.nodes[ndx1],self.nodes[newndx]

        if nd1.p.near(nd2.p):pdb.set_trace()

        newedge = edge(nd1,nd2,**kwargs)
        self._add_edge(newedge)


        nd1,nd2 = self.nodes[ndx2],self.nodes[newndx]

        if nd1.p.near(nd2.p):pdb.set_trace()

        newedge = edge(nd1,nd2,**kwargs)
        self._add_edge(newedge)





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

    g._add_edge((-100,-100,0),(0,-100,0),interpolated = False)
    g._add_edge((0,-100,0),(100,0,0))
    g._add_edge((100,0,0),(0,100,0))
    g._add_edge((0,100,0),(-100,0,1))
    g._add_edge((-100,0,1),(0,-100,1))
    g._add_edge((0,-100,1),(100,-100,1),interpolated = False)
    g._add_edge((100,-100,1),(200,-100,1))
    g._add_edge((200,-100,1),(300,0,0))
    g._add_edge((300,0,0),(300,100,0))

    return g

def hairpin():
    g = graph()

    g._add_edge((0,0,0),(100,50,1))
    g._add_edge((100,50,1),(0,100,2))
    g._add_edge((0,100,2),(100,150,3))
    #g._add_edge((0,0,0),(-50,100,1))
    #g._add_edge((-50,100,2),(100,150,3))

    return g

def circle():
    g = graph()

    g._add_edge((0,0,0),(50,50,0),interpolated = False)
    #g._add_edge((0,0,0),(50,50,0))
    g._add_edge((50,50,0),(0,100,0),interpolated = False)
    #g._add_edge((50,50,0),(0,100,0))
    g._add_edge((0,100,0),(-50,50,0),interpolated = True)
    #g._add_edge((0,100,0),(-50,50,0))
    g._add_edge((-50,50,0),(0,0,0),interpolated = True)
    #g._add_edge((-50,50,0),(0,0,0))

    return g

def opass():
    g = graph()

    bnd = dpr.square(100,100)
    oprgn1 = grg.overpass(bnd)
    oprgn1._graph(g)

    bnd = dpr.square(100,100,dpv.vector(500,0,0),dpr.rad(75))
    oprgn2 = grg.overpass(bnd)
    oprgn2._graph(g)

    bnd = dpr.square(100,100,dpv.vector(250,200,0),dpr.rad(35))
    oprgn3 = grg.overpass(bnd)
    oprgn3._graph(g)

    oprgn1._connect(oprgn2,g)
    oprgn2._connect(oprgn3,g)
    oprgn3._connect(oprgn1,g)

    g._update()

    ax = dtl.plot_axes()
    ax = oprgn1.plot(ax)
    ax = oprgn2.plot(ax)
    ax = oprgn3.plot(ax)
    ax = g.plot(ax)
    ax.set_zlim([0,40])

    plt.show()

    return g

def generate_graph():
    #g = graph()

    #g = smatter()
    #g = lsystem_graph()
    g = ramp()
    #g = opass()
    #g = hairpin()

    g._update()

    ax = g.plot()
    plt.show()

    return g





