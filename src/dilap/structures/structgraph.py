import dilap.core.base as db
import dilap.core.vector as dpv
import dilap.core.tools as dpr
import dilap.core.graph as dgr

import dilap.mesh.tools as dtl
import dilap.mesh.pointset as dps
import dilap.mesh.piecewisecomplex as pwc

import dilap.structures.tools as dstl
import dilap.structures.graphnode as gnd
import dilap.structures.graphedge as geg

import matplotlib.pyplot as plt
import random as rm
import math
import pdb



class graph(dgr.graph):

    nodeclass = gnd.node
    edgeclass = geg.edge

    # given a collection of relevant walls
    # produce polygons for one side of a wall
    def model_interior_wall(self,ekeys):
        pass

    def model_rooms(self):
        ww = 0.75      # THIS IS THE NOT THE RIGHT VALUE...
        mpolys = []
        for r in self.rooms:
            outlines = r
            #bbnd = []
            #for x in outlines[0]:
            #    nd = self.nodes[self.nodes_lookup[x]]
            #    bbnd.append(nd.p.copy())
            bbnd = self.get_node_points(outlines[0])
            #tbnd = []
            tbnd = self.get_node_points(outlines[-1])

            for x in range(len(outlines[-1])):
                nd = self.nodes[self.nodes_lookup[outlines[-1][x]]]
                tbnd[x].translate_z(nd.height)

            dpr.inflate(bbnd,-ww/math.sqrt(2))
            dpr.inflate(tbnd,-ww/math.sqrt(2))
            mpolys.append((tuple(bbnd),()))
            mpolys.append((tuple(tbnd),()))
        return mpolys

    def model_walls(self):
        mpolys = []
        for e in self.edges:
            if e is None:continue

            print('lookup',e.key(),self.rooms_lookup[e.key()])
            if len(self.rooms_lookup[e.key()]) == 1:
                #e.cut_window(3,2,1,0.5)
                pass

                '''#
                wh1,wh2 = e.one.height+e.one.gap,e.two.height+e.two.gap
                wargs = (e.one.p,e.two.p,wh1,wh2,e.width)
                wkwargs = {'doors':e.doors,'windows':e.windows}
                mpolys.extend(dstl.wall(*wargs,**wkwargs))
                '''#

            wh1,wh2 = e.one.height+e.one.gap,e.two.height+e.two.gap
            wargs = (e.one.p,e.two.p,wh1,wh2,e.width)
            wkwargs = {'doors':e.doors,'windows':e.windows}
            mpolys.extend(dstl.wall(*wargs,**wkwargs))
        return mpolys

    def model_corners(self):
        ww = 0.75      # THIS IS THE NOT THE RIGHT VALUE...
        mpolys = []
        for n in self.nodes:
            if n is None:continue

            print('node',n.key(),n.ring)

            mpolys.extend(dstl.post(n.p,4,ww,n.height+n.gap))
        return mpolys

    def model_roof(self):
        ww = 0.75      # THIS IS THE NOT THE RIGHT VALUE...
        mpolys = []

        for r in self.rooms:
            bnd = []
            outlines = r
            for nd in self.get_nodes(outlines[-1]):
                np = nd.p.copy().translate_z(nd.height+nd.gap)
                bnd.append(np)
            mpolys.append((tuple(bnd),()))

        return mpolys

    def model(self):
        mpolys = []

        plc1 = dtl.box(5,5,5)
        plc2 = dtl.box(5,4,4).translate(dpv.vector(4,0,0))

        print('union input')
        ax = dtl.plot_axes()
        ax = plc1.plot(ax)
        ax = plc2.plot(ax)
        plt.show()

        #plc3 = pwc.union(plc1,plc2)
        plc3 = pwc.difference(plc1,plc2)
        #plc3 = pwc.intersection(plc1,plc2)

        print('union output')
        ax = dtl.plot_axes()
        ax = plc3.plot(ax)
        plt.show()

        pys = []
        for px in range(plc3.polygoncount):
            pys.append(plc3.get_polygon_points(px))
        mpolys = pys

        #mpolys.extend(self.model_rooms())
        #mpolys.extend(self.model_walls())
        #mpolys.extend(self.model_corners())
        #mpolys.extend(self.model_roof())
        #mpolys = dtl.merge_polygons(mpolys)
        return mpolys

    def plot(self,ax = None):
        ax = dgr.graph.plot(self,ax)
        dtl.plot_polygon(list(self.boundary),ax)
        return ax

    def __init__(self,boundary,**kwargs):
        dgr.graph.__init__(self,**kwargs)

        self.boundary = boundary
        self.rooms_lookup = {}
        self.rooms = []
        self.roomcount = 0

    # given the index of a node, apply the effects of its layer
    def _apply_node_layer(self,ndx):
        nd = self.nodes[ndx]
        if nd is None or nd.layer == 0:return
        zkey = (nd.p.x,nd.p.y,nd.layer-1)
        below = self.nodes[self.nodes_lookup[zkey]]
        zoff = below.p.z + below.height + below.gap
        nd.p.translate_z(zoff)

    def _add_node(self,ndkey,**kwargs):
        kwargs['height'] = 8
        kwargs['gap'] = 1
        return dgr.graph._add_node(self,ndkey,**kwargs)

    # add a new edge to the graph, or return existing index
    def _add_edge(self,ndkey1,ndkey2,**kwargs):
        nex = dgr.graph._add_edge(self,ndkey1,ndkey2,**kwargs)
        if nex == self.edgecount-1:
            ekey1,ekey2 = (ndkey1,ndkey2),(ndkey2,ndkey1)
            self.rooms_lookup[ekey1] = []
            self.rooms_lookup[ekey2] = []
        return nex

    # delete an existing edge from the graph
    def _del_edge(self,ndkey1,ndkey2):
        dgr.graph._del_edge(self,ndkey1,ndkey2)
        ekey1,ekey2 = (ndkey1,ndkey2),(ndkey2,ndkey1)
        del self.rooms_lookup[ekey1]
        del self.rooms_lookup[ekey2]

    # find edges which bound both rooms u and v
    def _find_edges(self,u,v):
        found = []
        foundkeys = []
        for ekey in self.rooms_lookup:
            if ekey in foundkeys:continue
            foundkeys.append(ekey)
            foundkeys.append(ekey[::-1])
            ering = self.rooms_lookup[ekey]
            if u in ering and v in ering:
                found.append(self.edges_lookup[ekey])
            elif v == -1 and len(ering) == 1:
                found.append(self.edges_lookup[ekey])
        return found

    # add door to an edge which bound both rooms u and v
    def _connect_rooms(self,u,v):
        bwn = self._find_edges(u,v)
        if len(bwn) == 0:print('rooms are nonadjacent...',u,v)
        elif len(bwn) == 1:self.edges[bwn[0]].cut_door(2,3,0.5)
        else:
            elys = [self.edges[x].layer() for x in bwn]
            elym = min(elys)
            bwn = [bwn[x] for x in range(len(bwn)) if elys[x] == elym]
            self.edges[bwn[0]].cut_door(3,4,0.25)

    def _add_room(self,rmnds):
        self.rooms.append(rmnds)
        rl,rdex = self.rooms_lookup,self.roomcount
        outlines = rmnds
        for ox in range(len(outlines)):
            out = outlines[ox]
            for x in range(len(out)):
                ndkey1,ndkey2 = out[x-1],out[x]
                nex = self._add_edge(ndkey1,ndkey2)
                ekey1,ekey2 = (ndkey1,ndkey2),(ndkey2,ndkey1)
                if not rdex in rl[ekey1]:rl[ekey1].append(rdex)
                if not rdex in rl[ekey2]:rl[ekey2].append(rdex)
        self.roomcount += 1
        return rdex












