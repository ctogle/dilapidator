#imports
# cython: profile=True
#cimport cython
cimport dilap.geometry.tools as gtl
import dilap.geometry.tools as gtl
cimport dilap.geometry.pointset as dps

cimport dilap.geometry.vec3 as dpv
import dilap.geometry.vec3 as dpv
from dilap.geometry.vec3 cimport vec3
from dilap.geometry.quat cimport quat

import math,numpy

cdef class triangulation:

    # given the index of a point, return the index of a ghost
    # whose real edge intersects the point, or -1
    cdef int point_on_boundary(self,int u):
        cdef vec3 up = self.points.ps[u]
        cdef vec3 g1,g2
        cdef int gdx,gx1,gx2,gx3
        cdef float dx,dy,dv,nx,ny,prj1,prj2,prj3
        for gdx in range(self.ghostcnt):
            gst = self.ghosts[gdx]
            if gst is None:continue
            gx1,gx2,gx3 = gst
            #g1,g2 = self.points.get_points(gx1,gx2)
            g1,g2 = self.points.gps_c((gx1,gx2))
            if gtl.inseg_xy_c(up,g1,g2):
                return gdx
        return -1

    # given the indices of the endpoints of an edge, 
    # return 1 if locally delaunay, 0 otherwise
    cdef bint locally_delaunay(self,int u,int v):
        o1 = self.adjacent(u,v)
        o2 = self.adjacent(v,u)
        if o1 == -1 or o2 == -1:return 1
        if o1 == -2 or o2 == -2:return 1
        #up,vp,op1,op2 = self.points.get_points(u,v,o1,o2)
        up,vp,op1,op2 = self.points.gps_c((u,v,o1,o2))
        if gtl.segs_isect_perp_c(up,vp,op1,op2):
            if gtl.incircle_c(up,vp,op1,op2) > 0:return 0
            if gtl.incircle_c(vp,up,op2,op1) > 0:return 0
        return 1

    # return a vertex x such that uv
    # is a positively oriented edge
    cdef int adjacent(self,int u,int v):
        ekey = (u,v)
        if ekey in self.eg_tri_lookup:
            tri = self.eg_tri_lookup[(u,v)]
            if not tri is None:
                triv = [x for x in self.triangles[tri] if not x in ekey][0]
                return triv
        if ekey in self.eg_ghost_lookup:
            tri = self.eg_ghost_lookup[ekey]
            if not tri is None:return -2
        return -1

    def __cinit__(self,vec3 p0,vec3 pn):
        self.p0,self.pn = p0,pn
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
        #print('flipping an edge')
        o1 = self.adjacent(u,v)
        o2 = self.adjacent(v,u)
        #vs = self.points.get_points(u,v,o1,o2)
        vs = self.points.gps_c((u,v,o1,o2))
        #tcp1,tcr1 = dpr.circumscribe_tri_c(vs[0],vs[1],vs[2])
        #tcp2,tcr2 = dpr.circumscribe_tri_c(vs[1],vs[0],vs[3])
        #if tcp1.near(tcp2) and dpr.isnear_c(tcr1,tcr2):
        #    print('4way!',tcp1,tcp2,tcr1,tcr2,u,v,o1,o2)
        #    #return ()
        self.delete_triangle(u,v,o1)
        self.delete_triangle(v,u,o2)
        self.add_triangle(o1,o2,v)
        self.add_triangle(o2,o1,u)
        return (o2,v),(v,o1),(o1,u),(u,o2)

    # u is the vertex to insert. vwx is a positively oriented triangle whose
    # circumcircle encloses u
    cdef void insert_vertex(self,int u,int v,int w,int x):
        #vu,vv,vw,vx = self.points.get_points(u,v,w,x)
        vu,vv,vw,vx = self.points.gps_c((u,v,w,x))
        self.delete_triangle(v,w,x)
        self.dig_cavity(u,v,w)
        self.dig_cavity(u,w,x)
        self.dig_cavity(u,x,v)

    # u is the vertex to insert. vwg is a positively oriented ghost triangle whose
    # circumcircle encloses u
    cdef void insert_ghost_vertex(self,int u,int v,int w,int x):
        if not x == -2:raise ValueError
        self.delete_ghost(v,w)
        self.add_ghost(v,u)
        self.add_ghost(u,w)
        onb = self.point_on_boundary(u)
        if not onb == -1:self.add_triangle(u,v,w)

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
        self.ghosts.append((u,v,-2))
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
        if x == -1:return
        elif x == -2:self.add_triangle(u,v,w)
        else:
            #vu,vv,vw,vx = self.points.get_points(u,v,w,x)
            vu,vv,vw,vx = self.points.gps_c((u,v,w,x))
            # uvw is not delaunay, dig the adjacent two triangles
            if gtl.incircle_c(vu,vv,vw,vx) > 0:
                self.delete_triangle(w,v,x)
                self.dig_cavity(u,v,x)
                self.dig_cavity(u,x,w)
            # w,v is a facet of the cavity and uvw is delaunay
            else:self.add_triangle(u,v,w)

# initialize a triangulation data structure 
# based on a list of vectors bnd
cdef void initialize(triangulation data,list bnd):
    b1,b2,b3 = bnd[0],bnd[1],bnd[2]
    #bnorm,btang = dpr.norm_c(b1,b2,b3),dpr.tangent_c(b1,b2,b3)
    bnorm = gtl.nrm_c(b1,b2,b3)
    btang = b1.tov(b2).nrm()
    convexcom = vec3(0,0,0).com_c(bnd)
    convexrad = max([cx.d(convexcom) for cx in bnd])+1000
    c01delta = vec3(-1,-1,0).nrm().uscl(convexrad)
    c02delta = vec3( 1,-1,0).nrm().uscl(convexrad)
    c03delta = vec3( 0, 1,0).nrm().uscl(convexrad)
    c01 = convexcom.cp().trn(c01delta)
    c02 = convexcom.cp().trn(c02delta)
    c03 = convexcom.cp().trn(c03delta)
    c01x = data.points.ap_c(c01)
    c02x = data.points.ap_c(c02)
    c03x = data.points.ap_c(c03)
    #c01x,c02x,c03x = data.points.aps_c(c01,c02,c03)
    data.add_triangle(c01x,c02x,c03x)
    data.add_ghost(c03x,c02x)
    data.add_ghost(c02x,c01x)
    data.add_ghost(c01x,c03x)

# insert the vertex y into the triangulation data structure
# return a list of indices pointing to newly created triangles 
#cdef list point_location(triangulation data,vector y):
cdef list point_location(triangulation data,vec3 y):
    cdef int pretricnt,nv,onb,tdx,gdx,u,v,w,x,tu
    cdef vec3 vu,vv,vw
    pretricnt = data.tricnt
    #if data.points.find_point(y):return []
    #if data.points.fp_c(y):
    if not data.points.fp_c(y) == -1:return []
    #nv = data.points.add_point(y)
    nv = data.points.ap_c(y)
    onb = data.point_on_boundary(nv)
    if not onb == -1:
        v,w,x = data.ghosts[onb]
        data.delete_ghost(v,w)
        data.add_ghost(v,nv)
        data.add_ghost(nv,w)
        tu = data.adjacent(w,v)
        data.delete_triangle(w,v,tu)
        data.add_triangle(tu,nv,v)
        data.add_triangle(tu,w,nv)
        return [x for x in range(pretricnt,data.tricnt)]

    for tdx in range(data.tricnt-1,-1,-1):
        tri = data.triangles[tdx]
        if tri is None:continue
        else:u,v,w = tri
        #vu,vv,vw = data.points.get_points(u,v,w)
        vu,vv,vw = data.points.gps_c((u,v,w))
        if gtl.intri_xy_c(y,vu,vv,vw):
            data.insert_vertex(nv,u,v,w)
            return [x for x in range(pretricnt,data.tricnt)]

    for gdx in range(data.ghostcnt-1,-1,-1):
        ghost = data.ghosts[gdx]
        if ghost is None:continue
        else:u,v,w = ghost
        #vu,vv = data.points.get_points(u,v)
        vu,vv = data.points.gps_c((u,v))
        if not gtl.orient2d_c(vu,vv,y) < 0:
            data.insert_ghost_vertex(nv,u,v,w)
            return [x for x in range(pretricnt,data.tricnt)]

    print('never should have come to this')

# given a loop of points, add them to a triangulation 
# and return the line segments which bound the loop
cdef list loop_location(triangulation data,tuple loop):
    cdef list bnd = []
    cdef int lcnt = len(loop)
    cdef int x
    for x in range(lcnt):
        p1,p2 = loop[x-1],loop[x]
        #point_location(data,p1.xy())
        #bnd.append((p1.xy(),p2.xy()))
        #point_location(data,p1.cp_c())
        #bnd.append((p1.cp_c(),p2.cp_c()))
        point_location(data,p1.cp())
        bnd.append((p1.cp(),p2.cp()))
    return bnd

# locate each loop associated with polygon
cdef list polygon_location(triangulation data,tuple ebnd,tuple ibnds):
    cdef list bnd = []
    cdef int icnt = len(ibnds)
    cdef int ix
    bnd.extend(loop_location(data,ebnd))
    for ix in range(icnt):
        bnd.extend(loop_location(data,ibnds[ix]))
    return bnd

# given the exterior bound and interior bounds (holes) of a concave polygon
# remove triangles which are not part of the cover of the polygon
cdef void cover_polygon(triangulation data,tuple ebnd,tuple ibnds):
    extras = []
    for tdx in range(data.tricnt):
        tri = data.triangles[tdx]
        if tri is None:continue
        else:u,v,w = tri
        #ptri = tuple(data.points.get_points(u,v,w))
        ptri = tuple(data.points.gps_c((u,v,w)))
        extras.append(tdx)
        if gtl.polyinpoly_c(ebnd,ptri):
            extras.remove(tdx)
            for ibnd in ibnds:
                if gtl.polyinpoly_c(ibnd,ptri):
                    extras.append(tdx)
                    break
    for x in extras:
        xtri = data.triangles[x]
        if xtri is None:continue
        x1,x2,x3 = xtri
        data.delete_triangle(x1,x2,x3)

# delete all ghosts and add new ghosts according to where they
# should be after covering the plc edges
cdef void ghost_border(triangulation data,list edges):
    for gst in data.ghosts:
        if gst is None:continue
        g1,g2,g = gst
        data.delete_ghost(g1,g2)
    for plce in edges:
        if plce is None:continue
        plce1,plce2 = plce
        #e1,e2 = data.points.find_points(plce1,plce2)
        e1 = data.points.fp_c(plce1)
        e2 = data.points.fp_c(plce2)
        eadj = data.adjacent(e1,e2)
        if eadj == -1:data.add_ghost(e1,e2)
        eadj = data.adjacent(e2,e1)
        if eadj == -1:data.add_ghost(e2,e1)

# apply the flip algorithm until all edges are locally delaunay
cdef void constrain_delaunay(triangulation data):
    unfinished = [e for e in data.eg_tri_lookup]
    while unfinished:
        u,v = unfinished.pop(0)
        if not data.locally_delaunay(u,v):
            nedges = data.flip_edge(u,v)
            unfinished.extend(nedges)

# perform chews first refinement algorithm on a triangulation
cdef void refine_chews_first(triangulation data,float h):
    cdef list unfinished = [t for t in data.triangles]
    cdef int ucnt = len(unfinished)
    cdef int ufx1,ufx2,ufx3
    cdef vec3 v1,v2,v3
    cdef list ntxs
    cdef int ntx
    while ucnt > 0:
        unfin = unfinished.pop(0)
        ucnt -= 1
        if unfin is None:continue
        if not unfin in data.triangles:continue
        ufx1,ufx2,ufx3 = unfin
        #v1,v2,v3 = data.points.get_points(ufx1,ufx2,ufx3)
        try:
          v1,v2,v3 = data.points.gps_c((ufx1,ufx2,ufx3))
          tcp,tcr = gtl.circumscribe_tri_c(v1,v2,v3)
        except TypeError:
          print('fuckoff')
        if tcr/h > 1.0:
            ntxs = point_location(data,tcp)
            for ntx in ntxs:
                unfinished.append(data.triangles[ntx])
                ucnt += 1

# given triangulation data, construct a dictionary encoding vertex 1-rings
cdef dict vertex_rings(triangulation data):
    vrings = {}
    for x in range(data.points.pcnt):
        vrings[x] = []
        for egl in data.eg_tri_lookup:
            if egl is None:continue
            if data.eg_tri_lookup[egl] is None:continue
            if not x in egl:continue
            u,v = egl
            if   x == u and not v in vrings[x]:vrings[x].append(v)
            elif x == v and not u in vrings[x]:vrings[x].append(u)
    return vrings

# perform laplacian smoothing on interior points of a triangulation
cdef void smooth_laplacian(triangulation data,dict vrings,int smooths,float d):
    for s in range(smooths):
        dels = {}
        for x in range(data.points.pcnt):
            ghost = data.point_on_boundary(x)
            ring = vrings[x]
            if ghost == -1 and ring:
                rcom = vec3(0,0,0).com_c([data.points.ps[j] for j in ring])
                #rdel = data.points.ps[x].tov_c(rcom).uscl_c(d)
                rdel = data.points.ps[x].tov(rcom).uscl_c(d)
            else:rdel = vec3(0,0,0)
            dels[x] = rdel
        for x in range(data.points.pcnt):
            data.points.ps[x].trn_c(dels[x])

# plot the triangles currently found in a triangulation
cdef void plot_triangulation(triangulation data,ax = None):
    import dilap.mesh.tools as dtl
    if ax is None:ax = dtl.plot_axes_xy()
    for t in data.triangles:
        if t is None:continue
        else:t1,t2,t3 = t
        #tps = data.points.get_points(t1,t2,t3)
        tps = data.points.gps_c((t1,t2,t3))
        ax = dtl.plot_polygon_xy(tps,ax,center = True)

# given poly, a tuple containing vectors representing a polygon
# provide a list of simplices which triangulates the polygon
# poly contains an exterior bound and da tuple of interior bounds
# bounds are ordered loops of points with no duplicates
# bounds are possibly concave; interior bounds represent holes
cdef tuple triangulate_c(tuple ebnd,tuple ibnds,float hmin,bint refine,bint smooth):
    cdef vec3 p0,pn
    cdef quat prot

    #if refine:
    #    ebnd,ibnds = subdivide_edges(ebnd,ibnds)
    #    hminnnn,ebnd,ibnds = chew1_subdivide_polygon(ebnd,ibnds)


    #p0 = ebnd[0].cp_c()
    p0 = ebnd[0].cp()
    #pn = gtl.poly_nrm_c(ebnd)
    pn = gtl.poly_nrm(ebnd)
    prot = gtl.q_to_xy_c(pn)
    gtl.rot_poly_c((ebnd,ibnds),prot)
    data = triangulation(p0,pn)

    initialize(data,list(ebnd))
    plcedges = polygon_location(data,ebnd,ibnds)

    cover_polygon(data,ebnd,ibnds)
    ghost_border(data,plcedges)              
    constrain_delaunay(data)
    if refine:refine_chews_first(data,hmin)
    if smooth:smooth_laplacian(data,vertex_rings(data),100,0.5)

    prot.flp_c()
    for p in data.points.ps:p.rot(prot)
    gtl.rot_poly_c((ebnd,ibnds),prot)

    smps = []
    for tx in range(data.tricnt):
        tri = data.triangles[tx]
        if tri is None:continue
        smp = tuple([data.points.ps[px] for px in tri])
        smps.append(smp)
    gsts = []
    for gdx in range(data.ghostcnt):
        gst = data.ghosts[gdx]
        if gst is None:continue
        #gpair = data.points.get_points(gst[0],gst[1])
        gpair = data.points.gps_c((gst[0],gst[1]))
        gsts.append(gpair)
    return smps,gsts

# given poly, a tuple containing vectors representing a polygon
# provide a list of simplices which triangulates the polygon
# poly contains an exterior bound and da tuple of interior bounds
# bounds are ordered loops of points with no duplicates
# bounds are possibly concave; interior bounds represent holes
cpdef tuple triangulate(tuple ebnd,tuple ibnds,float hmin,bint refine,bint smooth):
    return triangulate_c(ebnd,ibnds,hmin,refine,smooth)

def split_nondelauney_edges(eb,ibs):
    aps = list(eb)
    for ib in ibs:aps.extend(list(ib))

    def split_loop(loop):
        x = 0
        while x < len(loop):
            p1,p2 = loop[x-1],loop[x]
            found = False
            for p in aps:
                if p.isnear(p1) or p.isnear(p2):continue
                cc = p1.mid(p2)
                cr = cc.d(p1)
                if p.inneighborhood(cc,cr):
                    loop.insert(x,cc)
                    found = True
                    break
            if found:continue
            else:x += 1
        return loop

    neb = split_loop(list(eb))
    nibs = [split_loop(list(ib)) for ib in ibs]
    return tuple(neb),tuple(tuple(nib) for nib in nibs)

def split_nondelauney_edges_chew1(eb,ibs):
    els = [eb[x-1].d(eb[x]) for x in range(len(eb))]
    for ib in ibs:els.extend([ib[x-1].d(ib[x]) for x in range(len(ib))])
    hmin = min(els)*math.sqrt(3)

    def split_loop(loop):
        oloop = [loop[0]]
        for x in range(1,len(loop)+1):
            if x == len(loop):x = 0
            p1,p2 = oloop[-1],loop[x]
            el = p1.d(p2)
            m = 1
            while el/m > hmin:m += 1
            divpts = p1.pline(p2,m-1)
            if divpts:
                for dvp in divpts:
                    oloop.append(dvp)
            if not x == 0:oloop.append(p2)

        return oloop

    neb = split_loop(list(eb))
    nibs = [split_loop(list(ib)) for ib in ibs]
    return hmin,tuple(neb),tuple(tuple(nib) for nib in nibs)




