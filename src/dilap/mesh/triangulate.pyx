import dilap.core.tools as dpr
import dilap.mesh.tools as dtl

cimport dilap.core.pointset as dps
import dilap.core.pointset as dps
cimport dilap.core.vector as dpv
import dilap.core.vector as dpv

import math
#stuff = 'hi'

cdef class triangulation:

    # given the index of a point, return the index of a ghost
    # whose real edge intersects the point, or None
    def point_on_boundary(self,u):
        up = self.points.ps[u]
        for gdx in range(self.ghostcnt):
            gst = self.ghosts[gdx]
            if gst is None:continue
            g1,g2 = self.points.get_points(gst[0],gst[1])
            dx = g2.x - g1.x
            dy = g2.y - g1.y
            dv = math.sqrt(dx**2 + dy**2)
            norm = dpv.vector(dy/dv,-dx/dv,0)
            nrmd = dpv.distance_to_edge(up,g1,g2,norm)
            if dtl.isnear(nrmd,0):
                linkd = dpv.distance(up,g1)+dpv.distance(up,g2)
                if dtl.isnear(linkd,dpv.distance(g1,g2)):
                    return gdx

    # return a vertex x such that uv
    # is a positively oriented edge
    def adjacent(self,u,v):
        ekey = (u,v)
        if ekey in self.eg_tri_lookup:
            tri = self.eg_tri_lookup[(u,v)]
            if not tri is None:
                triv = [x for x in self.triangles[tri] if not x in ekey][0]
                return triv
        if ekey in self.eg_ghost_lookup:
            tri = self.eg_ghost_lookup[(u,v)]
            if not tri is None:return self.ghosts[tri][2]

    cdef bint locally_delaunay(self,int u,int v):
        o1 = self.adjacent(u,v)
        o2 = self.adjacent(v,u)
        if o1 is None or o2 is None:return 1
        if o1 == 'g' or o2 == 'g':return 1
        up,vp,op1,op2 = self.points.get_points(u,v,o1,o2)
        if dtl.segments_intersect((up,vp),(op1,op2)):
            if dtl.incircle(up,vp,op1,op2) > 0:return 0
            if dtl.incircle(vp,up,op2,op1) > 0:return 0
        return 1

    # if triangle uvw is skinny
    # return the center of its circumcircle otherwise return None
    def skinny_triangle(self,u,v,w,h):
        vs = self.points.get_points(u,v,w)
        tcp,tcr = dpr.circumscribe_tri(*vs)
        if tcr/h > 1.0:return tcp 

    def __cinit__(self):
        self.points = dps.pointset()
        self.triangles = []
        self.tricnt = 0
        self.eg_tri_lookup = {}
        self.ghosts = []
        self.ghostcnt = 0
        self.eg_ghost_lookup = {}

    # given the edge u,v which bounds two non-ghost triangles
    # remove those triangles and replace with the alternative two that are
    # bounded by the same 4 vertices
    cdef tuple flip_edge(self,int u,int v):
        print('flipping an edge')
        o1 = self.adjacent(u,v)
        o2 = self.adjacent(v,u)
        vs = self.points.get_points(u,v,o1,o2)
        tcp1,tcr1 = dpr.circumscribe_tri(vs[0],vs[1],vs[2])
        tcp2,tcr2 = dpr.circumscribe_tri(vs[1],vs[0],vs[3])
        if tcp1.near(tcp2) and dtl.isnear(tcr1,tcr2):
            print('4way!',tcp1,tcp2,tcr1,tcr2,u,v,o1,o2)
            return ()
        self.delete_triangle(u,v,o1)
        self.delete_triangle(v,u,o2)
        self.add_triangle(o1,o2,v)
        self.add_triangle(o2,o1,u)
        return (o2,v),(v,o1),(o1,u),(u,o2)

    # add a positively oriented triangle u,v,w
    cdef void add_triangle(self,int u,int v,int w):
        self.triangles.append((u,v,w))
        self.eg_tri_lookup[(u,v)] = self.tricnt
        self.eg_tri_lookup[(v,w)] = self.tricnt
        self.eg_tri_lookup[(w,u)] = self.tricnt
        self.tricnt += 1

    # delete a positively oriented triangle u,v,w
    cdef void delete_triangle(self,int u,int v,int w):
        tri = self.eg_tri_lookup[(u,v)]
        if not tri is None:self.triangles[tri] = None
        self.eg_tri_lookup[(u,v)] = None
        self.eg_tri_lookup[(v,w)] = None
        self.eg_tri_lookup[(w,u)] = None

    # add a positively oriented ghost triangle u,v,g
    cdef void add_ghost(self,int u,int v):
        self.ghosts.append((u,v,'g'))
        self.eg_ghost_lookup[(u,v)] = self.ghostcnt
        self.ghostcnt += 1

    # delete a positively oriented ghost triangle u,v,g
    cdef void delete_ghost(self,int u,int v):
        ghost = self.eg_ghost_lookup[(u,v)]
        if not ghost is None:self.ghosts[ghost] = None
        self.eg_ghost_lookup[(u,v)] = None

    # u is a new vertex; is the oriented triangle u,v,w delaunay?
    cdef void dig_cavity(self,int u,int v,int w):
        # find triangle wvx opposite the facet vw from u
        x = self.adjacent(w,v)
        if x is None:return
        elif x == 'g':self.add_triangle(u,v,w)
        else:
            vu,vv,vw,vx = self.points.get_points(u,v,w,x)
            # uvw is not delaunay, dig the adjacent two triangles
            if dtl.incircle(vu,vv,vw,vx) > 0:
                self.delete_triangle(w,v,x)
                self.dig_cavity(u,v,x)
                self.dig_cavity(u,x,w)
            # w,v is a facet of the cavity and uvw is delaunay
            else:self.add_triangle(u,v,w)

    # u is the vertex to insert. vwx is a positively oriented triangle whose
    # circumcircle encloses u
    cdef void insert_vertex(self,int u,int v,int w,int x):
        self.delete_triangle(v,w,x)
        self.dig_cavity(u,v,w)
        self.dig_cavity(u,w,x)
        self.dig_cavity(u,x,v)

    # u is the vertex to insert. vwg is a positively oriented ghost triangle whose
    # circumcircle encloses u
    cdef void insert_ghost_vertex(self,int u,int v,int w,str x):
        if not x == 'g':raise ValueError
        self.delete_ghost(v,w)
        self.add_ghost(v,u)
        self.add_ghost(u,w)
        onb = self.point_on_boundary(u)
        if onb is None:self.add_triangle(u,v,w)

    # delete all ghosts and add new ghosts according to where they
    # should be after covering the plc edges
    cdef void ghost_border(self,list edges):
        for gst in self.ghosts:
            if gst is None:continue
            g1,g2,g = gst
            self.delete_ghost(g1,g2)
        for plce in edges:
            if plce is None:continue
            plce1,plce2 = plce
            e1,e2 = self.points.find_points(plce1,plce2)

            if e1 is None or e2 is None:
                print('ggggunit',e1,e2,plce)
                print('im off the rails',edges,plce)
                raise ValueError

            eadj = self.adjacent(e1,e2)
            if eadj is None:self.add_ghost(e1,e2)
            eadj = self.adjacent(e2,e1)
            if eadj is None:self.add_ghost(e2,e1)

# initialize a triangulation data structure 
# based on a list of vectors bnd
cdef void initialize(triangulation data,list bnd):
    convexcom = dpv.com(bnd)
    convexrad = max([dpv.distance(cx,convexcom) for cx in bnd])+1000
    c01delta = dpv.vector(-1,-1,0).normalize().scale_u(convexrad)
    c02delta = dpv.vector( 1,-1,0).normalize().scale_u(convexrad)
    c03delta = dpv.vector( 0, 1,0).normalize().scale_u(convexrad)
    c01 = convexcom.copy().translate(c01delta)
    c02 = convexcom.copy().translate(c02delta)
    c03 = convexcom.copy().translate(c03delta)
    c01x,c02x,c03x = data.points.add_points(c01,c02,c03)
    data.add_triangle(c01x,c02x,c03x)
    data.add_ghost(c03x,c02x)
    data.add_ghost(c02x,c01x)
    data.add_ghost(c01x,c03x)

# insert the vertex y into the triangulation data structure
# return a list of indices pointing to newly created triangles 
cdef list point_location(triangulation data,dpv.vector y):
    pretricnt = data.tricnt
    if data.points.find_point(y):return []
    nv = data.points.add_point(y)
    onb = data.point_on_boundary(nv)
    if onb:
        v,w,x = data.ghosts[onb]
        data.delete_ghost(v,w)
        data.add_ghost(v,nv)
        data.add_ghost(nv,w)
        tu = data.adjacent(w,v)
        data.delete_triangle(w,v,tu)
        data.add_triangle(tu,nv,v)
        data.add_triangle(tu,w,nv)
        return [x for x in range(pretricnt,data.tricnt)]

    for tdx in range(data.tricnt):
        tri = data.triangles[tdx]
        if tri is None:continue
        else:u,v,w = tri
        vu,vv,vw = data.points.get_points(u,v,w)
        if dpv.inside(y,[vu,vv,vw]):
            data.insert_vertex(nv,u,v,w)
            return [x for x in range(pretricnt,data.tricnt)]

    for gdx in range(data.ghostcnt):
        ghost = data.ghosts[gdx]
        if ghost is None:continue
        else:u,v,w = ghost
        vu,vv = data.points.get_points(u,v)
        if not dtl.orient2d(vu,vv,y) < 0:
            data.insert_ghost_vertex(nv,u,v,w)
            return [x for x in range(pretricnt,data.tricnt)]

# given the exterior bound and interior bounds (holes) of a concave polygon
# remove triangles which are not part of the cover of the polygon
cdef void cover_polygon(triangulation data,tuple ebnd,tuple ibnds):
    extras = []
    for tdx in range(data.tricnt):
        tri = data.triangles[tdx]
        if tri is None:continue
        else:u,v,w = tri
        ptri = data.points.get_points(u,v,w)
        extras.append(tdx)
        if dtl.concaves_contains(ebnd,ptri):
            extras.remove(tdx)
            for ibnd in ibnds:
                if dtl.concaves_contains(ibnd,ptri):
                    extras.append(tdx)
                    break
    for x in extras:
        xtri = data.triangles[x]
        if xtri is None:continue
        x1,x2,x3 = xtri
        data.delete_triangle(x1,x2,x3)

# apply the flip algorithm until all edges are locally delaunay
cdef void constrain_delaunay(triangulation data):
    unfinished = [e for e in data.eg_tri_lookup]
    while unfinished:
        u,v = unfinished.pop(0)

        #plcu,plcv = self.plc.points.find_points(*self.points.get_points(u,v))
        #if self.plc.segment(plcu,plcv):return True

        if not data.locally_delaunay(u,v):
            nedges = data.flip_edge(u,v)
            unfinished.extend(nedges)

cdef void refine_chews_first(triangulation data,float h):
    unfinished = [t for t in data.triangles]
    while unfinished:
        unfin = unfinished.pop(0)
        if unfin is None:continue
        if not unfin in data.triangles:continue
        ufx1,ufx2,ufx3 = unfin
        tcp = data.skinny_triangle(ufx1,ufx2,ufx3,h)
        if not tcp is None:
            ntxs = point_location(data,tcp)
            for ntx in ntxs:
                unfinished.append(data.triangles[ntx])

# given poly, a tuple containing vectors representing a polygon
# provide a list of simplices which triangulates the polygon
# poly contains an exterior bound and da tuple of interior bounds
# bounds are ordered loops of points with no duplicates
# bounds are possibly concave; interior bounds represent holes
def triangulate(ebnd,ibnds,hmin = None):
    data = triangulation()
    initialize(data,list(ebnd))
    plcedges = []
    for px in range(len(ebnd)):
        point_location(data,ebnd[px-1].copy())
        plcedges.append((ebnd[px-1],ebnd[px]))
    for ib in ibnds:
        for px in range(len(ib)):
            point_location(data,ib[px-1].copy())
            plcedges.append((ib[px-1],ib[px]))
    cover_polygon(data,ebnd,ibnds)
    data.ghost_border(plcedges)
    constrain_delaunay(data)
    refine_chews_first(data,hmin)

    smps = []
    for tx in range(data.tricnt):
        tri = data.triangles[tx]
        if tri is None:continue
        smp = tuple((data.points.ps[px] for px in tri))
        smps.append(smp)
    gsts = []
    for gdx in range(data.ghostcnt):
        gst = data.ghosts[gdx]
        if gst is None:continue
        gpair = data.points.get_points(gst[0],gst[1])
        gsts.append(gpair)
    return smps,gsts



