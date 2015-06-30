import dilap.core.mesh.tools as dtl
import dilap.core.mesh.pointset as dps
import dilap.core.mesh.tetrahedralization as dth
import dilap.core.mesh.triangulation as dtg

import dp_vector as dpv

import matplotlib.pyplot as plt
import pdb

class piecewise_linear_complex:

    def plot_xy(self,ax = None):
        ax = dtl.plot_points_xy(self.points.ps,ax)
        for edx in range(len(self.edges)):
            eg = self.edges[edx]
            if eg is None:continue
            veg = self.points.get_points(*eg)
            ax = dtl.plot_edges_xy(veg,ax)

        for key in self.covers:self.covers[key].plot_xy(ax)
        return ax

    def plot(self,ax = None):
        ax = plt.figure().add_subplot(111,projection = '3d')
        for pdx in range(self.points.pcnt):
            p = self.points.ps[pdx]
            ax.plot([p.x],[p.y],[p.z],marker = '+')
        for edx in range(len(self.edges)):
            up,vp = self.points.get_points(*self.edges[edx])
            ax.plot([up.x,vp.x],[up.y,vp.y],zs = [up.z,vp.z])
        return ax

    def __init__(self):
        self.points = dps.pointset()
        self.edges = []
        self.edgecount = 0
        self.polygons = []
        self.polygoncount = 0
        self.polyhedra = []
        self.polyhedroncount = 0
        self.covers = {}

        self.eg_lookup = {}
        self.eg_poly_lookup = {}

    def add_points(self,*nps):
        pst = self.points.pcnt
        self.points.add_points(*nps)
        return [x for x in range(pst,self.points.pcnt)]

    # u and v are indices of vertices; 
    # create a new edge which connects them
    def add_edge(self,u,v):
        ekey = (v,u)
        if ekey in self.eg_lookup:return self.eg_lookup[ekey]
        ekey = (u,v)
        if ekey in self.eg_lookup:return self.eg_lookup[ekey]
        edex = self.edgecount
        self.edges.append(ekey)
        self.edgecount += 1
        self.eg_lookup[ekey] = edex
        return edex

    def add_edges(self,*epts):
        exs = []
        for x in range(len(epts)):
            exs.append(self.add_edge(epts[x-1],epts[x]))
        return exs

    # given u,v which form an edge,
    # set all data to do with this edge to None
    # return the index of the deleted edge
    def delete_edge(self,u,v):
        edex = self.eg_lookup[(u,v)]
        self.edges[edex] = None
        self.eg_lookup[(u,v)] = None
        return edex

    # u,v likely form an existing edge;
    # replace this edge with two new edges
    def split_edge(self,u,v):

        def replace_edge(loop,which,new1,new2):
            loop.pop(which)
            loop.insert(which,new2)
            loop.insert(which,new1)

        def add_edge_lookup(poly,new):
            if not new in self.eg_poly_lookup:
                self.eg_poly_lookup[new] = []
            self.eg_poly_lookup[new].append(poly)

        newp = dpv.midpoint(*self.points.get_points(u,v))
        edex = self.delete_edge(u,v)
        w = self.points.add_point(newp)
        ne1dex = self.add_edge(u,w)
        ne2dex = self.add_edge(w,v)
        if edex in self.eg_poly_lookup:
            poly = self.eg_poly_lookup[edex]
            for pygn in poly:
                eb = pygn[0]
                if edex in eb:
                    epdex = eb.index(edex)
                    replace_edge(eb,epdex,ne1dex,ne2dex)
                else:
                    for ib in pygn[1]:
                        if edex in ib:
                            ipdex = ib.index(edex)
                            replace_edge(ib,ipdex,ne1dex,ne2dex)
                            break
                add_edge_lookup(pygn,ne1dex)
                add_edge_lookup(pygn,ne2dex)
        return (u,w),(w,v)

    def add_polygons(self,*polygons):
        pxs = []
        for px in range(len(polygons)):
            polygon = polygons[px]
            polyeb,polyibs = polygon
            pxs.append(self.polygon_frompoints(polyeb,*polyibs))
        return pxs

    def polygon_frompoints(self,ebnd,*ibnds):
        plcxs = self.add_points(*ebnd)
        ebnddexes = self.add_edges(*plcxs)
        ebnddexes.append(self.add_edge(plcxs[-1],plcxs[0]))
        ibnddexes = []
        for ibnd in ibnds:
            plcxs = self.add_points(*ibnd)
            hedgedexes = self.add_edges(*plcxs)
            hedgedexes.append(self.add_edge(plcxs[-1],plcxs[0]))
            ibnddexes.append(hedgedexes)
        return self.polygon_fromedges(ebnddexes,*ibnddexes)

    def polygon_fromedges(self,ebnd,*ibnds):
        pdex = self.polygoncount
        polygon = (ebnd,ibnds)
        self.polygons.append(polygon)
        self.polygoncount += 1
        for e in ebnd:
            if not e in self.eg_poly_lookup:
                self.eg_poly_lookup[e] = []
            self.eg_poly_lookup[e].append(polygon)
        for ib in ibnds:
            for i in ib:
                if not i in self.eg_poly_lookup:
                    self.eg_poly_lookup[i] = []
                self.eg_poly_lookup[i].append(polygon)
        return pdex

    def add_polyhedra(self,*polyhedra):
        return
        raise NotImplemented

    def tetrahedralize(self):
        tetra = dth.tetrahedralization(self)
        self.covers['tetra'] = tetra

    def triangulate_xy(self):
        tri = dtg.triangulation(self)
        self.covers['tri'] = tri

def model_plc(points = None,edges = None,polygons = None,polyhedra = None):
    plc = piecewise_linear_complex()
    if not points is None:plc.add_points(*points)
    if not edges is None:plc.add_edges(*edges)
    if not polygons is None:plc.add_polygons(*polygons)
    if not polyhedra is None:plc.add_polyhedra(*polyhedra)
    plc.triangulate_xy()

    ax = plc.plot_xy()
    plt.show()

    pelt = plc.covers['tri'].pelt()
    return pelt
    










