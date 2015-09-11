import dilap.core.base as db
import dilap.core.vector as dpv
import dilap.core.tools as dpr
import dilap.core.lsystem as dls

import dilap.mesh.tools as dtl
import dilap.mesh.pointset as dps

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
    def _loopwalk(self,ie,d,w):
        def complete(inp):
            i1,i2 = inp[0],inp[1]
            cnt = len(inp)
            for x in range(1,cnt-1):
                if inp[x] == i1:
                    if inp[x+1] == i2:
                        return inp[:x+1]

        if d:inpath = [ie.one.key(),ie.two.key()]
        else:inpath = [ie.two.key(),ie.one.key()]
        while True:
            ekey = (inpath[-2],inpath[-1])
            e = self.edges[self.edges_lookup[ekey]]
            nx = e._walk(inpath[-1],w)
            if nx is None:
                inpath.append(inpath[-2])
                nx = e._walk(inpath[-3],w)
            nxndkey = self.nodes[nx].key()
            #if ie.two.key() == (0.0,100.0,0.0):
            #    print('going',nxndkey,inpath)
            #    pdb.set_trace()
            res = complete(inpath)
            if not res is None:return res

            #if inpath[-1] == inpath[0] and nxndkey == inpath[1]:return inpath
            #if inpath.count(ie.one.key()) > 1 and nxndkey == inpath[1]:
            #if inpath.count(ie.one.key()) > 1 and inpath.count(ie.two.key()) > 1:
            #    return inpath
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
                #print('closed loop!',len(edgestodo))
                rloop = tuple(ewalkrcw)
            else:
                print('unclosed loop!',len(edgestodo))
                pdb.set_trace()
                rloop = tuple(ewalkrccw[::-1][:-1]+ewalkrcw[1:])

            if set(ewalklccw) ==  set(ewalklcw):
                #print('closed loop!',len(edgestodo))
                lloop = tuple(ewalklccw)
            else:
                print('unclosed loop!',len(edgestodo))
                pdb.set_trace()
                lloop = tuple(ewalklccw[::-1][:-1]+ewalklcw[1:])

            rlloop = lloop[::-1]
            if not dpr.cyclic_permutation(rlloop,rloop):
                edgelloops.append(lloop)
                edgerloops.append(rloop)
        
        #pdb.set_trace()
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
        le = ne
        for elx in range(2,elcnt+1):
            elx1,elx2 = elx-1,elx if elx < elcnt else 0
            nekey = (eloop[elx1],eloop[elx2])
            if nekey[0] == nekey[1]:return looppts
            ne = self.edges[self.edges_lookup[nekey]]
            nelooppts = self._find_road_points(looppts[-1],le,ne)
            looppts.extend(nelooppts)
            le = ne
        return looppts

    # given the last and next edge, and the last point in a loop
    # properly return the set of road points which connects
    def _find_road_points(self,tip,le,ne):
        # create the shortest line segment from el1 to el2
        def tiptail(el1,el2):
            d1 = dpv.distance(el1[ 0],el2[ 0])
            d2 = dpv.distance(el1[ 0],el2[-1])
            d3 = dpv.distance(el1[-1],el2[ 0])
            d4 = dpv.distance(el1[-1],el2[-1])
            md = min(d1,d2,d3,d4)
            if   md == d1:return el1[ 0],el2[ 0]
            elif md == d2:return el1[ 0],el2[-1]
            elif md == d3:return el1[-1],el2[ 0]
            elif md == d4:return el1[-1],el2[-1]
        def closer(p,r,l):
            if dpv.distance(p,r) < dpv.distance(p,l):return r
            else:return l

        '''#
        ax = dtl.plot_axes_xy()
        ax = dtl.plot_point_xy(tip,ax)
        ax = dtl.plot_edges_xy(le.rpts,ax)
        ax = dtl.plot_edges_xy(ne.rpts,ax)
        s1,s2 = tiptail(le.rpts,ne.rpts)
        ax = dtl.plot_edges_xy([s1,s2],ax,lw = 5.0)
        ax = dtl.plot_edges_xy([tip,closer(tip,ne.rbpts[0],ne.rbpts[-1])],ax)
        ax = dtl.plot_edges_xy([tip,closer(tip,ne.lbpts[0],ne.lbpts[-1])],ax)
        plt.show()
        '''#

        '''#
        this function is verrrrry sloppy.... rewrite it....
        '''#
        def same_side(lp):
            lp0d = dpv.distance(lp[ 0],tip)
            lp1d = dpv.distance(lp[-1],tip)
            lpt = lp[0] if lp0d < lp1d else lp[-1]
            s1,s2 = tiptail(le.rpts,ne.rpts)
            segsect = dpr.segments_intersect_noncolinear(s1,s2,lpt,tip)
            if not segsect:return lpt
        def connect_end(lp,lpt):
            d1,d2 = dpv.distance(lp[0],lpt),dpv.distance(lp[-1],lpt)
            if d1 < d2:return lp[:]
            else:return lp[::-1]

        if le is ne:
            if tip in le.rbpts:return connect_end(ne.lbpts,tip)
            else:return connect_end(ne.rbpts,tip)
        else:
            lrpt = same_side(ne.rbpts)
            llpt = same_side(ne.lbpts)
            if lrpt is None and llpt is None:
                lsd = dpv.distance(tip,ne.lbpts[ 0])
                led = dpv.distance(tip,ne.lbpts[-1])
                rsd = dpv.distance(tip,ne.rbpts[ 0])
                red = dpv.distance(tip,ne.rbpts[-1])
                sxs = dpr.order_ascending([lsd,led,rsd,red]) 
                nelooppts = None
                for sx in sxs:
                    if   sx == 0 and not tip in ne.lbpts:nelooppts = ne.lbpts[:]
                    elif sx == 1 and not tip in ne.lbpts:nelooppts = ne.lbpts[:-1]
                    elif sx == 2 and not tip in ne.rbpts:nelooppts = ne.rbpts[:]
                    elif sx == 3 and not tip in ne.rbpts:nelooppts = ne.rbpts[:-1]
                    if not nelooppts is None:break
                return nelooppts
            if not lrpt is None:return connect_end(ne.rbpts,lrpt)
            else:return connect_end(ne.lbpts,llpt)

    # return a collection of points outlining all nodes/edges in the graph
    def _edge_loop_boundaries(self):
        def uniq_loop(eloops,elp):
            uniq = True
            for elps in eloops:
                for x in range(len(elps)):
                    p = elps[x]
                    for y in range(len(elp)):
                        q = elp[y]
                        if p.near(q):return False
            return True

        edgelloops,edgerloops = self._edge_loops()
        rperms = {}
        lperms = {}
        for ex in range(len(edgelloops)):
            lloop,rloop = edgelloops[ex],edgerloops[ex]
            rkey = rloop[:-1]
            isperm = False
            for rps in rperms:
                if dpr.cyclic_permutation(rkey,rps):isperm = True;break
            if not isperm:rperms[rkey] = self._edge_loop_points(rloop,0) 
            lkey = lloop[:-1]
            isperm = False
            for lps in lperms:
                if dpr.cyclic_permutation(lkey,lps):isperm = True;break
            if not isperm:lperms[lkey] = self._edge_loop_points(lloop,1) 
        eloops = []
        for el in lperms:
            elp = [v for v in lperms[el]]
            if uniq_loop(eloops,elp):eloops.append(elp)
        for el in rperms:
            elp = [v for v in rperms[el]]
            if uniq_loop(eloops,elp):eloops.append(elp)
        return self._rank_edge_loops(eloops)

    # determine how the loops are arranged based on containment
    # so that they can be properly triangulated
    def _rank_edge_loops(self,eloops):
        bedgeloops = {}

        #ax = dtl.plot_axes_xy()
        #ax = self.plot_xy(ax)
        #for bedge in eloops:ax = dtl.plot_edges_xy(bedge,ax)
        #plt.show()

        containments = [[] for el in eloops]
        for elx in range(len(eloops)):
            elp = tuple(eloops[elx])
            for elxo in range(len(eloops)):
                if elxo == elx:continue
                elpo = tuple(eloops[elxo])
                isect = dpr.concaves_intersect(elp,elpo)
                elins = dpr.inconcave_xy(elpo[0],elp)
                if isect:raise ValueError
                elif elins:containments[elx].append(elxo)
        looplook = {'king':[],'interiors':[]}
        for elx in range(len(eloops)):
            cont = containments[elx]
            if cont:looplook['king'].append(eloops[elx])
            else:looplook['interiors'].append(eloops[elx])
        return looplook
        # provide a polygon for the terrain
        #
        # provide a polygon for the road
        #
        # the terrain runs from convex bound to the loop that contains 
        # all other loops
        # the terrain also contains the interiors of all loops of road
        # 
        # the road extends from the loop that contains all others to the
        # collection of all other loops of road
        #
        # assume the graph is connected? fix if not?

    # calculate polygons representing regions to place terrain
    def _regions(self):
        rpts = []
        for eg in self.edges:rpts.extend([x.copy() for x in eg.rpts])
        convexbnd = dpr.pts_to_convex_xy(rpts)
        convexbnd = dpr.inflate(convexbnd,50)
        eloops = self._edge_loop_boundaries()

        #ax = dtl.plot_axes_xy()
        #dtl.plot_edges_xy(eloops['king'][0],ax)
        #for loop in eloops['interiors']:dtl.plot_edges_xy(loop,ax)
        #plt.show()

        self.tpolygons = [(tuple(convexbnd),(tuple(eloops['king'][0]),))]+\
                              [(tuple(i),()) for i in eloops['interiors']]
        self.rpolygons = [(eloops['king'][0],tuple(eloops['interiors']))]

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

    g._add_edge((0,0,0),(100,50,0))
    g._add_edge((100,50,0),(0,100,0))
    g._add_edge((0,100,0),(100,150,0))
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

def newcastle():
    g = graph()

    l,w = 200,200
    g._add_edge((0,0,0),(l*0.5,w*0.5,0),interpolated = True)
    g._add_edge((l*0.5,w*0.5,0),(0,w,0),interpolated = True)
    g._add_edge((0,w,0),(-l*0.5,l*0.5,0),interpolated = True)
    g._add_edge((-l*0.5,w*0.5,0),(0,0,0),interpolated = True)
    g._add_edge((0,0,0),(0,w*0.5,0),interpolated = True)

    return g

def eight():
    g = graph()

    r = 100
    g._add_edge((0,0,0),(r,0,0),interpolated = True)
    g._add_edge((r,0,0),(2*r,0,0),interpolated = True)
    g._add_edge((2*r,0,0),(2*r,r,0),interpolated = True)
    g._add_edge((2*r,r,0),(r,r,0),interpolated = True)
    g._add_edge((r,r,0),(0,r,0),interpolated = True)
    g._add_edge((0,r,0),(0,0,0),interpolated = True)
    g._add_edge((r,r,0),(r,0,0),interpolated = True)

    g._update()

    #ax = dtl.plot_axes()
    #ax = g.plot(ax)
    #ax.set_zlim([0,40])
    #plt.show()

    return g

def clover():
    g = graph()

    r = 100
    g._add_edge((0,0,0),( r,0,0),interpolated = True)
    g._add_edge((0,0,0),(-r,0,0),interpolated = True)
    g._add_edge((0,0,0),(0, r,0),interpolated = True)
    g._add_edge((0,0,0),(0,-r,0),interpolated = True)
    g._add_edge(( r,0,0),(2*r,-r,0),interpolated = True)
    g._add_edge((2*r,-r,0),(3*r,0,0),interpolated = True)
    g._add_edge((3*r,0,0),(2*r,r,0),interpolated = True)
    g._add_edge((2*r,r,0),(r,0,0),interpolated = True)

    g._update()

    #ax = dtl.plot_axes()
    #ax = g.plot(ax)
    #ax.set_zlim([0,40])
    #plt.show()

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





