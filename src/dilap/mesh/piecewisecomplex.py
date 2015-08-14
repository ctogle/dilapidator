import dilap.core.tools as dpr
import dilap.core.vector as dpv
import dilap.core.model as dmo
import dilap.mesh.tools as dtl
import dilap.mesh.pointset as dps
import dilap.mesh.tetrahedralization as dth
import dilap.mesh.triangulation as dtg
import dilap.mesh.triangulate as dtg2

import matplotlib.pyplot as plt
import pdb
import math

sqrt3 = math.sqrt(3)

class piecewise_linear_complex:

    def plot(self,ax = None):
        ax = dtl.plot_points(self.points.ps,ax)
        for edx in range(len(self.edges)):
            e = self.edges[edx]
            if e is None:continue
            dtl.plot_edges(self.points.get_points(*e),ax)
        for smp in self.simplices:dtl.plot_polygon(list(smp),ax)
        for gst in self.ghostbnds:dtl.plot_edges(gst,ax,lw = 5.0)
        return ax

    def plot_xy(self,ax = None):
        if ax is None:ax = dtl.plot_axes_xy()
        ax = dtl.plot_points_xy(self.points.ps,ax)
        for edx in range(len(self.edges)):
            eg = self.edges[edx]
            if eg is None:continue
            veg = self.points.get_points(*eg)
            ax = dtl.plot_edges_xy(veg,ax)
        for smp in self.simplices:dtl.plot_polygon_xy(list(smp),ax)
        for gst in self.ghostbnds:dtl.plot_edges_xy(gst,ax,lw = 5.0)
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
    def split_edge(self,u,v,newp = None):

        def replace_edge(loop,which,new1,new2):
            loop.pop(which)
            loop.insert(which,new2)
            loop.insert(which,new1)

        def add_edge_lookup(poly,new):
            if not new in self.eg_poly_lookup:
                self.eg_poly_lookup[new] = []
            self.eg_poly_lookup[new].append(poly)

        if newp is None:
            newp = dpv.midpoint(*self.points.get_points(u,v))
        edex = self.delete_edge(u,v)
        w = self.points.add_point(newp)
        ne1dex = self.add_edge(u,w)
        ne2dex = self.add_edge(w,v)
        if edex in self.eg_poly_lookup:
            poly = self.eg_poly_lookup[edex]
            for pygn in poly:
                eb,ibs = pygn
                if edex in eb:
                    epdex = eb.index(edex)
                    replace_edge(eb,epdex,ne1dex,ne2dex)
                else:
                    for ib in ibs:
                        if edex in ib:
                            ipdex = ib.index(edex)
                            replace_edge(ib,ipdex,ne1dex,ne2dex)
                            break
                add_edge_lookup(pygn,ne1dex)
                add_edge_lookup(pygn,ne2dex)
        return (u,w),(w,v)

    # is the edge formed by uv (or vu) a segment in the domain?
    def segment(self,u,v):
        pkey = (u,v)
        if pkey in self.eg_lookup:
            px = self.eg_lookup[pkey]
            if not px is None:
                if not self.edges[px] is None:
                    return True
        pkey = (v,u)
        if pkey in self.eg_lookup:
            px = self.eg_lookup[pkey]
            if not px is None:
                if not self.edges[px] is None:
                    return True
        return False

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
        ibnddexes = []
        for ibnd in ibnds:
            plcxs = self.add_points(*ibnd)
            hedgedexes = self.add_edges(*plcxs)
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

    # given a point p, return the index of a polygon 
    # which contains it or None if none exists
    def find_polygon(self,p):
        for px in range(self.polygoncount):
            poly = self.polygons[px]
            if poly is None:continue
            ebnd,ibnds = poly
            ebnd = [self.edges[e][0] for e in ebnd]
            ebnd = self.points.get_points(*ebnd)
            if dpv.inside(p,ebnd):
                isin = True
                for ib in ibnds:
                    ib = [self.edges[e][0] for e in ib]
                    ib = self.points.get_points(*ib)
                    if dpv.inside(p,ib):
                        isin = False
                        break
                if isin:return px

    def add_polyhedra(self,*polyhedra):
        return
        raise NotImplemented

    # return a dictionary of the length of 
    # every edge currently in the mesh
    def edge_lengths(self):
        elengths = {}
        for edx in range(self.edgecount):
            e = self.edges[edx]
            if e is None or e in elengths:continue
            ep1,ep2 = self.points.get_points(*e)
            d = dpv.distance(ep1,ep2)
            elengths[e] = d
            elengths[e[::-1]] = d
        return elengths

    # given v1,v2, the positions of the endpoints of an edge, 
    # return True if locally delaunay
    def locally_delaunay_edge(self,v1,v2):
        cc = dpv.midpoint(v1,v2)
        cr = dpv.distance(cc,v1)
        for p in self.points:
            if p.near(v1) or p.near(v2):continue
            if dpr.inside_circle(p,cc,cr,(dpv.zero(),dpv.zhat)):
                return False
        return True

    # force edges to be locally delaunay by splitting as needed
    def subdivide_edges(self):
        unfinished = [e for e in self.edges]
        while unfinished:
            unfin = unfinished.pop(0)
            if unfin is None:continue
            v1,v2 = self.points.get_points(*unfin)
            if not self.locally_delaunay_edge(v1,v2):
                ne1,ne2 = self.split_edge(*unfin)
                unfinished.append(ne1)
                unfinished.append(ne2)

    def chew1_subdivide_polygon(self,px):
        elengs = self.edge_lengths()

        poly = self.polygons[px]
        polyes = poly[0][:]
        for ib in poly[1]:polyes += ib[:]

        es,erng = self.edges,range(self.edgecount)
        unfinished = [es[x] for x in erng if x in polyes]
        hmin = min([elengs[x] for x in elengs if x in unfinished])
        print('hmin:',hmin)

        while unfinished:
            unfin = unfinished.pop(0)
            if unfin is None:continue
            ex1,ex2 = unfin
            ep1,ep2 = self.points.get_points(ex1,ex2)
            eleng = elengs[unfin]
            m = 1
            while eleng/m > sqrt3*hmin:m += 1
            divpts = dpr.point_line(ep1,ep2,m)[1:-1]
            curr = ex1
            for dpt in divpts:
                ne1,ne2 = self.split_edge(curr,ex2,dpt)
                curr = self.points.find_point(dpt)
        return hmin

    # force edges to abide by chew1 precondition on edge lengths
    def chew1_subdivide_edges(self):
        elengs = self.edge_lengths()
        hmin = min([elengs[x] for x in elengs])

        print('hmin:',hmin)

        unfinished = [e for e in self.edges]
        while unfinished:
            unfin = unfinished.pop(0)
            if unfin is None:continue
            ex1,ex2 = unfin
            ep1,ep2 = self.points.get_points(ex1,ex2)
            eleng = elengs[unfin]
            m = 1
            while eleng/m > sqrt3*hmin:m += 1
            divpts = dpr.point_line(ep1,ep2,m)[1:-1]
            curr = ex1
            for dpt in divpts:
                ne1,ne2 = self.split_edge(curr,ex2,dpt)
                curr = self.points.find_point(dpt)
        elengs = self.edge_lengths()
        hmin = min([elengs[x] for x in elengs])
        return hmin

    def tetrahedralize(self):
        tetra = dth.tetrahedralization(self)
        self.covers['tetra'] = tetra

# BEGIN HERE
# MUST TRIANGULATE POLYGONS ONE BY ONE, 
# PERMITTING NON XY TRI.TION AS WELL!!

    # given the index of an edge, return vectors for its endpoints
    # or return None if the edge is missing
    def get_edge_points(self,ex):
        edge = self.edges[ex]
        if edge is None:return
        return self.points.get_points(*edge)

    # given a loop of edge indices, return a loop of vectors
    # which represent that edgeloop
    def get_edgeloop_points(self,eloop):
        ps,es = self.points.ps,self.edges
        ploop = tuple((ps[es[elx][0]] for elx in eloop))
        return ploop

    # given the index of a polygon, return its representation
    # using a set of vectors, or None if the polygon is missing
    def get_polygon_points(self,px):
        poly = self.polygons[px]
        if poly is None:return
        eb,ibs = poly
        ebps  = self.get_edgeloop_points(eb)
        ibsps = tuple((self.get_edgeloop_points(ib) for ib in ibs))
        return (ebps,ibsps)

    # iterate over each polygon, adding simplices which properly cover
    # a simplex is a tuple of indices pointing to self.points
    def triangulate(self):
        self.subdivide_edges()
        smps,bnds = [],[]
        for x in range(self.polygoncount):
            hmin = self.chew1_subdivide_polygon(x)
            polypts = self.get_polygon_points(x)
            #polysmp = dtg2.triangulate(*polypts,hmin = hmin)
            polysmp,polybnd = dtg2.triangulate(*polypts,hmin = hmin)
            smps.extend(polysmp)
            bnds.extend(polybnd)
        self.simplices = smps
        self.ghostbnds = bnds

    def pelt(self):
        s = dmo.model()
        for smp in self.simplices:
            s._triangle(*smp)
        return s

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
    










