import dilap.core.base as db
import dilap.core.tools as dpr
import dilap.core.vector as dpv
import dilap.core.ray as dry
import dilap.core.model as dmo
import dilap.core.pointset as dps
import dilap.mesh.tools as dtl
#import dilap.mesh.pointset as dps
#import dilap.mesh.tetrahedralization as dth
#import dilap.mesh.triangulation as dtg
import dilap.mesh.triangulate as dtg2

import matplotlib.pyplot as plt
import pdb
import math

sqrt3 = math.sqrt(3)



def inpolyhedron(point,triangles):
    pray = dry.ray(point,dpv.x())
    isects = []

    isects = dry.intersect_hits(pray,triangles)

    print('inininin',len(triangles),isects)

    '''#
    for x in range(len(triangles)):
        t1,t2,t3 = triangles[x]
        isct = pray.intersect_tri(t1,t2,t3)
        if isct == 1:
            t = pray.cast.x
            print('t',t)
            isects.append(t)
        else:print('noinsct',x)
    '''#
    ins = len(isects) % 2 > 0

    pray = dry.ray(point,dpv.y())
    isects = []

    isects = dry.intersect_hits(pray,triangles)

    ins = ins or len(isects) % 2 > 0

    print('insideeee',ins,len(isects))
    '''#
    ax = dtl.plot_axes()
    for x in range(len(triangles)):
        ax = dtl.plot_polygon(list(triangles[x]),ax)  
    ax = dtl.plot_edges([
        pray.origin,pray.origin.copy().translate(
            pray.direction.copy().scale_u(100))],ax)
    plt.show()
    '''#

    return ins

# return the subset of p1s that does not intersect the interior of plc2
def polygons_outcomplex(p1s,plc2):
    plc2.refine,plc2.smooth = False,False
    plc2.triangulate()
    tris = plc2.simplices
    outpolyhedron = []
    for x in range(len(p1s)):

        pcom = dpv.center_of_mass(list(p1s[x][0]))
        #point = p1s[x][0][0].copy()
        #point.translate(dpv.v1_v2(point,pcom).scale_u(0.1))
        point = pcom

        if not inpolyhedron(point,tris):outpolyhedron.append(p1s[x])
    return outpolyhedron

# return the subset of p1s that does not intersect the interior of plc2
def polygons_incomplex(p1s,plc2):
    plc2.refine,plc2.smooth = False,False
    plc2.triangulate()
    tris = plc2.simplices
    outpolyhedron = []
    for x in range(len(p1s)):

        #point = dpv.center_of_mass(list(p1s[x][0]))
        pcom = dpv.center_of_mass(list(p1s[x][0]))
        #point = p1s[x][0][0].copy()
        #point.translate(dpv.v1_v2(point,pcom).scale_u(0.1))
        point = pcom

        if inpolyhedron(point,tris):outpolyhedron.append(p1s[x])
    return outpolyhedron

# given a poly and a set of other polygons, return pieces of the input 
# polygon upon consideration their intersections
def break_polygon(py,p2s,subop = 'union'):
    breakers = [x for x in range(len(p2s))]
    broken = [dpr.copy_polygon(py)]
    pieces = []

    while broken:
        print('prepop',len(broken))
        piece = broken.pop(0)

        '''#
        print('PEACE!!!!!')
        ax = dtl.plot_axes()
        dtl.plot_polygon_full(piece,ax,lw = 4.0)
        plt.show()
        '''#

        eb1,ibs1 = piece
        ebn1 = dpr.polygon_normal(eb1)

        brokeone = False
        for x in breakers:
            breaker = p2s[x]
            eb2,ibs2 = breaker
            ebn2 = dpr.polygon_normal(eb2)
            pj1 = dpv.project_coords(list(eb1),ebn1)
            pj2 = dpv.project_coords(list(eb2),ebn1)

            if dpr.isnear(pj2.x,pj2.y):
                if not dpr.isnear(pj1.x,pj2.x):
                    print('disjoint polygons!')
                else:
                    if subop == 'union':
                        pu = dtl.polygon_union(piece,breaker)
                    elif subop == 'intersection':
                        pu = dtl.polygon_intersection(piece,breaker)
                    elif subop == 'difference':
                        pu = dtl.polygon_difference(piece,breaker)
                    else:print('unknown subop!',subop)
                    print('known subop!',subop)
                    print('coplanar polygons!',pu)

                    if not pu is None:

                        ax = dtl.plot_axes()
                        #dtl.plot_polygon_full(piece,ax)
                        #dtl.plot_polygon_full(breaker,ax)
                        dtl.plot_polygon(list(pu[0]),ax,lw = 5.0)
                        dtl.plot_polygon(list(pu[1][0]),ax,lw = 1.0)
                        plt.show()

                        #broken.append(pu)
                        pieces.append(pu)
                        print("breaker will continue to intersect this guy...")
                        brokeone = True
                        break
            else:
                print('skew polygons')
                plintersect = dtl.planes_intersection(ebn1,eb1[0],ebn2,eb2[0])
                if not plintersect is None:
                    pli1,pli2 = plintersect

                    intins = subop == 'difference'
                    #intins = dpr.inconcave_xy(dpv.midpoint(pli1,pli2),eb2)

                    print('inis',subop,intins)
                    breakersect = dtl.line_intersects_polygon_at(
                                        pli1,pli2,breaker,intins)
                    if not breakersect is None:

                        if not len(breakersect) == 2:
                            print('breakersect is confusing!!')

                            ax = dtl.plot_axes()
                            dtl.plot_polygon_full(breaker,ax)
                            dtl.plot_line(pli1,pli2,25,ax,lw = 4.0)
                            for x in range(len(breakersect)):
                                dtl.plot_point(breakersect[x],ax)
                            plt.show()

                            pdb.set_trace()
                        
                        b1,b2 = breakersect
                        lsp = dtl.segment_split_polygon(b1,b2,piece)
                        if not lsp is None:
                            
                            '''#
                            print('lsssspppp',len(lsp))
                            ax = dtl.plot_axes()
                            dtl.plot_line(b1,b2,25,ax,lw = 8.0)
                            for ls in lsp:
                                ax = dtl.plot_polygon_full(ls,ax,lw = 1.0)
                            ax = dtl.plot_polygon_full(breaker,ax,lw = 3.0)
                            plt.show()
                            '''#

                            broken.extend(lsp)
                            brokeone = True
                            #pdb.set_trace()
                            break
        if not brokeone:
            pieces.append(piece)
            print('nice piece',len(pieces))
    
    '''#
    if len(pieces) > 1:
      print('broke some shit??')
      ax = dtl.plot_axes()
      for p in pieces:dtl.plot_polygon_full(p,ax)
      plt.show()
    '''#
  
    return pieces

# given the bounding polygons of two polyhedrons, break p1s based 
# upon the intersection with p2s into properly nonintersecting polygons
def break_complex(p1s,p2s,subop = 'union'):
    print('entry to break complex',len(p1s),len(p2s))
    broken = []
    for x in range(len(p1s)):
        py = p1s[x]
        pieces = break_polygon(py,p2s,subop)
        broken.extend(pieces)
    return broken

# given two piecewise linear complexes, make a broken set of polygons
# representing nonintersecting proper polygons to construct a plc
# representing the union,intersection,difference of the two input plcs
def break_complexes(plc1,plc2,subop = 'union'):
    p1s = []
    for px in range(plc1.polygoncount):

        #if not px == 5:continue
        #if not px == 3:continue
        #if px == 0 or px == 1:continue

        p1s.append(plc1.get_polygon_points(px))
    p2s = []
    for px in range(plc2.polygoncount):

        #if px == 0 or px == 1:continue
        #if px in [0,2,3,4,5]:continue
        #if not px == 3:continue

        p2s.append(plc2.get_polygon_points(px))
    p1s = [x for x in p1s if not x is None]
    p2s = [x for x in p2s if not x is None]

    print('before breaking something!!!!!')
    print('before breaking something!!!!!')
    print('before breaking something!!!!!')
    ax = dtl.plot_axes()
    for s in p1s:dtl.plot_polygon_full(s,ax)
    for s in p2s:dtl.plot_polygon_full(s,ax)
    plt.show()

    p1seg = break_complex(p1s,p2s,subop)

    print('after breaking something!!!!!')
    print('after breaking something!!!!!')
    print('after breaking something!!!!!')
    ax = dtl.plot_axes()
    for s in p1seg:dtl.plot_polygon_full(s,ax,lw = 4.0)
    for s in p1s:dtl.plot_polygon_full(s,ax)
    for s in p2s:dtl.plot_polygon_full(s,ax)
    plt.show()

    p2seg = break_complex(p2s,p1s,subop)

    print('after breaking another!!!!!')
    print('after breaking another!!!!!')
    print('after breaking another!!!!!')
    ax = dtl.plot_axes()
    for s in p1seg:dtl.plot_polygon_full(s,ax,lw = 4.0)
    for s in p1s:dtl.plot_polygon_full(s,ax)
    for s in p2s:dtl.plot_polygon_full(s,ax)
    plt.show()

    return p1seg,p2seg

# return a plc representing the union of two plcs
def union(plc1,plc2):
    broken = break_complexes(plc1,plc2,subop = 'union')
    if broken is None:return
    else:p1seg,p2seg = broken
    p1inp2 = polygons_outcomplex(p1seg,plc2)
    p2inp1 = polygons_outcomplex(p2seg,plc1)
    union = piecewise_linear_complex()
    union.add_polygons(*(p1inp2+p2inp1))
    return union

# return a plc representing the intersection of two plcs
def intersection(plc1,plc2):
    broken = break_complexes(plc1,plc2,subop = 'intersection')
    if broken is None:return
    else:p1seg,p2seg = broken
    p1inp2 = polygons_incomplex(p1seg,plc2)
    p2inp1 = polygons_incomplex(p2seg,plc1)
    union = piecewise_linear_complex()
    union.add_polygons(*(p1inp2+p2inp1))
    return union

# return a plc representing the difference of two plcs
def difference(plc1,plc2):
    broken = break_complexes(plc1,plc2,subop = 'difference')
    if broken is None:return
    else:p1seg,p2seg = broken
    p1inp2 = polygons_outcomplex(p1seg,plc2)
    p2inp1 = polygons_incomplex(p2seg,plc1)
    #p1inp2 = p1seg
    #p2inp1 = p2seg
    diffr = piecewise_linear_complex()
    diffr.add_polygons(*(p1inp2+p2inp1))
    return diffr

class piecewise_linear_complex(db.base):

    def plot(self,ax = None):
        if ax is None:
            l = self.radius()
            ax = dtl.plot_axes(x = l)
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

    def radius(self):
        rs = [p.magnitude() for p in self.points.ps]
        return max(rs)

    def __init__(self,*args,**kwargs):
        self.points = dps.pointset()
        self.edges = []
        self.edgecount = 0
        self.polygons = []
        self.polygoncount = 0
        self.polyhedra = []
        self.polyhedroncount = 0
        self.covers = {}
        self._def('refine',False,**kwargs)
        self._def('smooth',False,**kwargs)

        self.eg_lookup = {}
        self.eg_poly_lookup = {}

        self.simplices = []
        self.ghostbnds = []

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
        ekey = (u,v)
        if not ekey in self.eg_lookup:return
        edex = self.eg_lookup[ekey]
        self.edges[edex] = None
        self.eg_lookup[ekey] = None
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
            for pygnx in poly:
                pygn = self.polygons[pygnx]
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
                add_edge_lookup(pygnx,ne1dex)
                add_edge_lookup(pygnx,ne2dex)
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

    def delete_polygon(self,px):
        poly = self.polygons[px]
        if poly is None:return
        eb,ibs = poly
        for ex in eb:self.delete_edge(*self.edges[ex])
        for ib in ibs:
            for ex in ib:self.delete_edge(*self.edges[ex])
        self.polygons[px] = None

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
        for e in ebnd:
            if not e in self.eg_poly_lookup:
                self.eg_poly_lookup[e] = []
            #self.eg_poly_lookup[e].append(polygon)
            self.eg_poly_lookup[e].append(self.polygoncount)
        for ib in ibnds:
            for i in ib:
                if not i in self.eg_poly_lookup:
                    self.eg_poly_lookup[i] = []
                #self.eg_poly_lookup[i].append(polygon)
                self.eg_poly_lookup[i].append(self.polygoncount)
        self.polygoncount += 1
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
    def locally_delaunay_edge(self,u,v):
        v1,v2 = self.points.get_points(u,v)
        cc = dpv.midpoint(v1,v2)
        cr = dpv.distance(cc,v1)
        polyxs = self.eg_poly_lookup[self.eg_lookup[(u,v)]]
        polyps = [self.get_polygon_points(px) for px in polyxs]
        for polyp in polyps:
            ebps,ibs = polyp
            for ep in ebps:
                if ep.near(v1) or ep.near(v2):continue
                if dpr.inside_circle(ep,cc,cr):return False
            for ibps in ibs:
                for ip in ibps:
                    if ip.near(v1) or ip.near(v2):continue
                    if dpr.inside_circle(ip,cc,cr):return False
        return True

    # force edges to be locally delaunay by splitting as needed
    def subdivide_edges(self):
        unfinished = [e for e in self.edges]
        while unfinished:
            unfin = unfinished.pop(0)
            if unfin is None:continue
            if not self.locally_delaunay_edge(*unfin):
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
                curr = ne1[1]
        return hmin

    # given a vector, move all points by the vector tn
    def translate(self,tn):
        for p in self.points.ps:
            p.translate(tn)
        return self

    # given a vector, rotate all points by the quaternion qn
    def rotate(self,qn):
        for p in self.points.ps:
            p.rotate(qn)
        return self

    # given the index of a polygon and a vector, move
    # the associated points by the vector tn
    def translate_polygon(self,px,tn):
        polyp = self.get_polygon_points(px)
        for ep in polyp[0]:ep.translate(tn)
        return self

    # given the index of a polygon and a vector, extrude the 
    # polygon adding other polygons to maintain connectivity
    def extrude_polygon(self,px,tn):
        pre = self.get_polygon_points_copy(px)
        self.translate_polygon(px,tn)
        post = self.get_polygon_points_copy(px)
        for ex in range(len(pre[0])):
            p1,p2,p3,p4 = pre[0][ex-1],pre[0][ex],post[0][ex],post[0][ex-1]
            self.add_polygons(((p1,p2,p3,p4),()))

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

    # the same as "get_polygon_points" except return copies of the points
    def get_polygon_points_copy(self,px):
        poly = self.get_polygon_points(px)
        copy = (tuple(p.copy() for p in poly[0]),
            tuple(tuple(p.copy() for p in ib) for ib in poly[1]))
        return copy

    # verify/correct any unacceptable polygons
    def clean_polygons(self):
        def repair_boundary(eb,ib):
            raise NotImplemented

        for px in range(self.polygoncount):
            polybs = self.get_polygon_points(px)
            eb,ibs = polybs
            neb = eb[:]
            nibs = []
            for ib in ibs:
                isect = dpr.concaves_intersect(eb,ib)
                ins = dpr.inconcave_xy(ib[0],eb)
                if isect:
                    print('polygon hole intersects boundary',px)
                    #neb = repair_boundary(eb,ib)
                elif ins:nibs.append(ib)
                else:print('polygon hole found outside of boundary',px)

            #r = self.radius()
            #ax = dtl.plot_axes(x = r)
            #ax = dtl.plot_polygon(list(eb),ax)
            #for ib in ibs:ax = dtl.plot_polygon(list(ib),ax)
            #plt.show()

            npoly = (neb,tuple(nibs))
            if not polybs == npoly:
                self.delete_polygon(px)
                self.add_polygons(npoly)

                #r = self.radius()
                #ax = dtl.plot_axes(x = r)
                #ax = dtl.plot_polygon(list(neb),ax)
                #for ib in nibs:ax = dtl.plot_polygon(list(ib),ax)
                #plt.show()

    # iterate over each polygon, adding simplices which properly cover
    # a simplex is a tuple of indices pointing to self.points
    def triangulate(self):
        #self.clean_polygons()
        self.subdivide_edges()
        smps,bnds = [],[]
        ref,smo = self.refine,self.smooth
        for x in range(self.polygoncount):
            if self.polygons[x] is None:continue
            #hmin = self.chew1_subdivide_polygon(x)
            if ref:hmin = self.chew1_subdivide_polygon(x)
            else:hmin = 1.0
            polypts = self.get_polygon_points(x)

            #eb,ibs = polypts
            #olysmp,polybnd = dtg2.triangulate(eb,ibs,hmin,ref,smo)

            eb,ibs = polypts
            polysmp,polybnd =\
                dtg2.triangulate_nonplanar(eb,ibs,hmin,ref,smo,
                           dpv.vector(0,0,1),dpv.vector(0,0,0))

            smps.extend(polysmp)
            bnds.extend(polybnd)
        self.simplices = smps
        self.ghostbnds = bnds

    # 
    def tetrahedralize(self):
        self.simplices = smps
        self.ghostbnds = bnds

    def pelt(self):
        s = dmo.model()
        for smp in self.simplices:
            t1,t2,t3 = smp
            s._triangle(t3,t2,t1)
        return s

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
    










