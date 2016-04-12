import dilap.core.base as db
import dilap.core.tools as dpr
import dilap.core.vector as dpv
import dilap.core.graphnode as gnd
import dilap.core.graphedge as geg

import dilap.mesh.tools as dtl

import matplotlib.pyplot as plt
import pdb



class graph(db.base):

    nodeclass = gnd.node
    edgeclass = geg.edge

    def plot(self,ax = None):
        if ax is None:ax = dtl.plot_axes()
        for n in self.nodes:
            if not n is None:n.plot(ax)
        for e in self.edges:
            if not e is None:e.plot(ax)
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

    def __init__(self,**kwargs):
        self.nodes = []
        self.nodes_lookup = {}
        self.edges = []
        self.edges_lookup = {}
        self.nodecount = 0
        self.edgecount = 0

    # given the index of a node, apply the effects of its layer
    def _apply_node_layer(self,ndx):
        nd = self.nodes[ndx]
        if nd is None:return
        nd.p.translate_z(20*nd.layer)

    # add a new node to the graph or existing node index
    # ndkey is a tuple(x,y,layer)
    def _add_node(self,ndkey,**kwargs):
        if ndkey in self.nodes_lookup:
            nd = self.nodes[self.nodes_lookup[ndkey]]
            if not nd is None:return nd.index
        nx,ny,nl = ndkey
        nz = 0.0
        newnode = self.nodeclass(dpv.vector(nx,ny,nz),layer = nl,**kwargs)
        newnode.index = self.nodecount
        self.nodes.append(newnode)
        self.nodes_lookup[ndkey] = newnode.index
        self.nodecount += 1
        self._apply_node_layer(newnode.index)
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
        egkey = (ndkey1,ndkey2)
        if egkey in self.edges_lookup:
            eg = self.edges[self.edges_lookup[egkey]]
            if not eg is None:return eg.index
        egkey = (ndkey2,ndkey1)
        if egkey in self.edges_lookup:
            eg = self.edges[self.edges_lookup[egkey]]
            if not eg is None:return eg.index
        newedge = self.edgeclass(nd1,nd2,**kwargs)
        newedge.index = self.edgecount
        ndir1,ndir2 = newedge.directions()
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





