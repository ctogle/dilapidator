import dilap.core.base as db
import dilap.core.tools as dpr
import dilap.core.vector as dpv

import dilap.graph.graph as dgg
import dilap.graph.twomanifold as tmg
#import dilap.graph.halfedge as the
import dilap.graph.wire as twr

import dilap.mesh.tools as dtl
import matplotlib.pyplot as plt

import pdb



# a room is a topological shell with holes
class room(db.base):

    def plot(self,ax):
        for fx in range(self.shell.facecount):
            sf = self.shell.faces[fx]
            if sf is None:continue
            sfvs = self.shell.mask(0,None,None,sf)
            sfps = [self.geom.points.ps[v.px] for v in sfvs]
            dtl.plot_polygon(sfps,ax)
        return ax

    def __init__(self,*args,**kwargs):
        self._def('boundary',None,**kwargs)                                         

        self.geom = dgg.geometry()
        self.shell = tmg.twomanifold_graph()
        #self.shell = the.halfedge_graph()

        h = 8
        bbnd  = [b.copy() for b in self.boundary[0]]
        tbnd  = [b.copy().translate_z(h) for b in bbnd]
        #front = bbnd[0],bbnd[1],tbnd[1],tbnd[0]
        #back  = bbnd[2],bbnd[3],tbnd[3],tbnd[2]
        #left  = bbnd[1],bbnd[2],tbnd[2],tbnd[1]
        #right = bbnd[3],bbnd[0],tbnd[0],tbnd[3]
        self.shell.new_face(self.geom.points.new_points( *bbnd))
        self.shell.new_face(self.geom.points.new_points( *tbnd))
        #self.shell.new_face(self.geom.points.new_points(*front))
        #self.shell.new_face(self.geom.points.new_points( *back))
        #self.shell.new_face(self.geom.points.new_points( *left))
        #self.shell.new_face(self.geom.points.new_points(*right))

# sgraph will represent the topology of a building
# and create a mesh representing this topology
# a room is a node, a door/window is an edge, possibly
# to the exterior, which is implicitly a room
class sgraph(db.base):

    def plot_rgraph(self,ax):
        for n in self.rgraph.nodes:
            if not n is None:
                n.room.plot(ax)
                dtl.plot_point(self.geom.points.ps[n.px],ax)
        for e in self.rgraph.edges:
            if not e is None:
                e1 = self.geom.points.ps[e.one.px]
                e2 = self.geom.points.ps[e.two.px]
                dtl.plot_edges([e1,e2],ax)

    def plot(self,ax = None):
        if ax is None:ax = dtl.plot_axes()
        self.plot_rgraph(ax)
        self.plot_rgraph(ax)
        dtl.plot_polygon(self.boundary,ax)
        r = self.radius()
        ax.set_xlim([-r,r])
        ax.set_ylim([-r,r])
        ax.set_zlim([-r,r])
        return ax

    def radius(self):
        pmags = []
        for nx in range(self.rgraph.nodecount):
            nd = self.rgraph.nodes[nx]
            if nd is None:continue
            pmags.append(self.geom.points.ps[nd.px].magnitude())
        return max(pmags)

    def __init__(self,*args,**kwargs):
        self.geom = dgg.geometry()
        self.rgraph = twr.wire_graph()
        self._def('boundary',None,**kwargs)

    def block(self,bnd):
        bl = (tuple(dpr.square(5,10,dpv.vector(2.5,5,0))),())
        nbnd = dtl.polygon_difference(bnd,bl)
        rmprints = [nbnd,bl]
        return rmprints

    def blueprint(self):
        root = (tuple(b.copy() for b in self.boundary),())

        rmprints = self.block(root)

        self.geom.points.add_point(dpv.vector(0,0,0))
        self.geom.points.add_point(dpv.vector(0,0,10))

        entry = self.rgraph.add_node(0)
        room1 = self.rgraph.add_node(1)
        self.rgraph.add_edge(0,1)

        entry.room = room(boundary = rmprints[0])
        room1.room = room(boundary = rmprints[1])

    def fun(self):
        self.blueprint()

        #self.geom.points.add_point(dpv.vector(0,0,0))
        #self.geom.points.add_point(dpv.vector(0,0,10))

        #entry = self.rgraph.add_node(0)
        #room1 = self.rgraph.add_node(1)
        #self.rgraph.add_edge(0,1)

        #entry.room = room()
        #room1.room = room()

        return self

def test():
    bnd = dpr.square(10,20,dpv.vector(5,10,0))
    g = sgraph(boundary = bnd).fun()
    g.plot()
    plt.show()



