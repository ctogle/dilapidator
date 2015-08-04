import dilap.core.base as db

import dilap.core.tools as dpr
import dilap.core.mesh.tools as dtl
import dilap.core.lsystem as dls

import dilap.generate.infrastructure.tools as itl

import dp_vector as dpv

import matplotlib.pyplot as plt
import random as rm
import pdb



def graph_dec(graphfunc):
    def graphdec(self,g):
        print('decorated!',g)
    return graphdec

# a region is represents a selfcontained area that can be systematically
# connected to other regions
#
# a region has a boundary in the xy-plane
# regions have plugs which encode how roads can enter the region from the boundary
# if the region is bounded, its boundary has a road along it
# a region has a sealevel, where standing water will be found (if visible)
# a region can only intersect another region along a boundary of both
# a region can contain another region which it will connect its plugs to
#
# roads on the interior of regions follow the terrain
# roads on the boundary of regions shape the terrain
class region(db.base):

    def plot_boundary(self,layer,ax = None):
        if ax is None:ax = dtl.plot_axes()
        bnd = [b.copy().translate_z(20*layer) for b in self.boundary]
        dtl.plot_polygon(bnd,ax)
        return ax

    def plot(self,ax = None,limit = 250):
        if ax is None:ax = dtl.plot_axes()
        for pg in self.plugs:
            edx,epos,elay = pg
            e0,e1 = edx,edx+1 if edx+1<len(self.boundary) else 0
            s,e = self.boundary[e0],self.boundary[e1]
            pgp = s.linterpolate(e,epos).translate_z(20*elay)
            ndkey = (pgp.x,pgp.y,elay)
            dtl.plot_point(pgp,ax)
        for l in self._layers():self.plot_boundary(l,ax)
        ax.set_xlim([-limit,limit])
        ax.set_ylim([-limit,limit])
        ax.set_zlim([0,40])
        for ch in self.children:ax = ch.plot(ax,limit = limit)
        return ax

    # produce a set of node keys for self.plugs
    def _plugkeys(self):
        plugkeys = []
        for pg in self.plugs:
            edx,epos,elay = pg
            e0,e1 = edx,edx+1 if edx+1<len(self.boundary) else 0
            s,e = self.boundary[e0],self.boundary[e1]
            pgp = s.linterpolate(e,epos)
            plugkeys.append((pgp.x,pgp.y,elay))
        return plugkeys

    # return the subset of self.plugs which reside on edge edgex
    def _plugsonedge(self,edgex,layer = None):
        if not layer is None:lays = [layer]
        else:lays = self._layers()
        onedge = []
        for pgx in range(len(self.plugs)):
            pg = self.plugs[pgx]
            if pg[0] == edgex and pg[2] in lays:
                onedge.append(pgx)
        return onedge

    def _graph_boundary(self,g):
        pgkeys = self._plugkeys()
        bpkeys = []
        bpcnt = len(self.boundary)
        for bpx in range(bpcnt):
            bp = self.boundary[bpx]
            bpkeys.append((bp.x,bp.y,0))
            epxs = self._plugsonedge(bpx,layer = 0)
            bplugs = [dpv.vector(*pgkeys[x]) for x in epxs]
            if bplugs:
                poepxs = dpv.proximity_order_xy(bp,bplugs)
                for poepx in poepxs:bpkeys.append(pgkeys[epxs[poepx]])
        bpkeys.append(bpkeys[0])
        g._add_edges(bpkeys,interpolated = False)
        return g

    # using the provided graph g, add nodes/edges to represent
    # this region and its plugs
    @graph_dec
    def _graph(self,g):
        pkeys = self._plugkeys()
        for pk in pkeys:g._add_node(pk)
        if self.bounded:g = self._graph_boundary(g)
        for ch in self.children:g = ch._graph(g)
        return g

    # given another region o and graph g
    # add nodes and edges which connect the regions properly
    def _connect_edges(self,o,sedgex,oedgex,g,**kwargs):
        splugs,oplugs = self._plugsonedge(sedgex),o._plugsonedge(oedgex)
        spgcnt,opgcnt = len(splugs),len(oplugs)

        spkeys = self._plugkeys()
        subspkeys = [spkeys[splugs[x]] for x in range(spgcnt)]
        opkeys = o._plugkeys()
        subopkeys = [opkeys[oplugs[x]] for x in range(opgcnt)]
        b1,b2 = self._tangent(o).rotate_z(dpr.rad(90)).normalize(),dpv.zhat.copy()
        sps = dpv.projected_order([dpv.vector(*spk) for spk in subspkeys],b1,b2)
        ops = dpv.projected_order([dpv.vector(*opk) for opk in subopkeys],b1,b2)

        if spgcnt == opgcnt:pgpairs = [(x,y) for x,y in zip(sps,ops)]
        else:
            print('plug ambiguity!!!')
            pdb.set_trace()

        for pgpair in pgpairs:
            pkey1,pkey2 = spkeys[splugs[pgpair[0]]],opkeys[oplugs[pgpair[1]]]
            g._add_edge(pkey1,pkey2,**kwargs)
            print('added edge',pkey1,pkey2)

    # given another region o and graph g
    # add nodes and edges which connect the regions properly
    def _connect(self,o,g,**kwargs):
        sedgex,oedgex = self._facing(o)
        self._connect_edges(o,sedgex,oedgex,g)

    # given another region o, whose boundary is contained by self.boundary
    # give the proper parent/child relationship to each
    def _embed(self,o,**kwargs):
        self.children.append(o)
        o.parent = self

    def __init__(self,boundary,**kwargs):
        self._def('bounded',False,**kwargs)
        self._def('sealevel',0.0,**kwargs)
        self._def('plugs',[],**kwargs)
        self.boundary = boundary
        self._def('children',[],**kwargs)
        self._def('parent',None,**kwargs)

    # return a list of layers for which this region
    # has connecting plugs
    def _layers(self):
        lays = []
        for pg in self.plugs:
            edx,epos,elay = pg
            if not elay in lays:
                lays.append(elay)
        if not lays:lays.append(0)
        return lays

    # return a vector pointing from com of self.boundary to com of o.boundary
    def _tangent(self,o):
        scom = dpv.center_of_mass(self.boundary)
        ocom = dpv.center_of_mass(o.boundary)
        sotang = dpv.v1_v2(scom,ocom)
        return sotang

    # given another region o, return the edge indices
    # for self,o for the best choice of facing angles
    def _facing(self,o):
        sfaces = dpv.edge_normals_xy(self.boundary)
        ofaces = dpv.edge_normals_xy(o.boundary)
        #scom = dpv.center_of_mass(self.boundary)
        #ocom = dpv.center_of_mass(o.boundary)
        #sotang = dpv.v1_v2(scom,ocom).normalize()
        sotang = self._tangent(o).normalize()
        sdiffs = [dpr.deg(dpv.angle_between_xy(sotang,ev)) for ev in sfaces]
        sdiffs = [dpr.clamp_periodic(sd,0,360) for sd in sdiffs]
        sdiffs = [sd if sd < 180 else 360 - sd for sd in sdiffs]
        sfacex = sdiffs.index(min(sdiffs))
        sotang.flip()
        odiffs = [dpr.deg(dpv.angle_between_xy(sotang,ev)) for ev in ofaces]
        odiffs = [dpr.clamp_periodic(od,0,360) for od in odiffs]
        odiffs = [od if od < 180 else 360 - od for od in odiffs]
        ofacex = odiffs.index(min(odiffs))
        #print('facing',sotang.flip(),sdiffs,sfacex)
        #print('facing',sotang.flip(),odiffs,ofacex)
        return sfacex,ofacex

class overpass(region):

    # using the provided graph g, add nodes/edges to represent
    # this region and its plugs
    def _graph(self,g):
        pkeys = self._plugkeys()

        g._add_edge(pkeys[0],pkeys[1],interpolated = False)
        g._add_edge(pkeys[2],pkeys[3],interpolated = False)
        g._add_edge(pkeys[4],pkeys[5],interpolated = False)
        g._add_edge(pkeys[6],pkeys[7],interpolated = False)

        region._graph(self,g)
        return g

    def __init__(self,boundary,**kwargs):
        opassplugs = [
            (0,0.25,0),(2,0.75,0),(0,0.75,0),(2,0.25,0),
            (1,0.25,1),(3,0.75,1),(1,0.75,1),(3,0.25,1)]
        self._def('plugs',opassplugs,**kwargs)
        region.__init__(self,boundary,**kwargs)
    
def neighborhood(*args,**kwargs):
    kwargs['plugs'] = [(1,0.5,0)]
    nbhd = region(*args,**kwargs)

    #@nbhd._graph
    #def nbhd_graph(self,g):
        #pkeys = self._plugkeys()
        #pkey0 = pkeys[0]
        #pkey1 = (pkey0[0],pkey0[1]-25,pkey0[2])
        #g._add_edge(pkey0,pkey1,interpolated = False)
    #    print('decorated!',g)
    #    return g

    return nbhd






