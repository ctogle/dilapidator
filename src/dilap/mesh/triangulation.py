import dilap.core.vector as dpv
import dilap.core.model as dmo
import dilap.core.tools as dpr

import dilap.mesh.pointset as dps
import dilap.mesh.tools as dtl
import dilap.mesh.pointset as dps

import matplotlib.pyplot as plt
import math
import pdb

sqrt3 = math.sqrt(3)

class triangulation:
    
    def plot(self,ax = None):
        if ax is None:ax = dtl.plot_axes()
        for tdx in range(self.tricnt):
            tri = self.triangles[tdx]
            if tri is None:continue
            vtri = self.points.get_points(*tri)
            dtl.plot_polygon(vtri,ax,center = True)
        return ax

    def plot_xy(self,ax = None):
        if ax is None:ax = dtl.plot_axes_xy()
        for tdx in range(self.tricnt):
            tri = self.triangles[tdx]
            if tri is None:continue
            vtri = self.points.get_points(*tri)
            dtl.plot_polygon_xy(vtri,ax,center = True)
        for gdx in range(self.ghostcnt):
            gst = self.ghosts[gdx]
            if gst is None:continue
            gpair = self.points.get_points(gst[0],gst[1])
            dtl.plot_edges_xy(gpair,ax,lw = 5.0)
        return ax

    # add a positively oriented ghost triangle u,v,g
    def add_ghost(self,u,v):
        self.ghosts.append((u,v,'g'))
        self.eg_ghost_lookup[(u,v)] = self.ghostcnt
        self.ghostcnt += 1

    # delete a positively oriented ghost triangle u,v,g
    def delete_ghost(self,u,v):
        ghost = self.eg_ghost_lookup[(u,v)]
        if not ghost is None:self.ghosts[ghost] = None
        self.eg_ghost_lookup[(u,v)] = None

    # add a positively oriented triangle u,v,w
    def add_triangle(self,u,v,w):
        self.triangles.append((u,v,w))
        self.eg_tri_lookup[(u,v)] = self.tricnt
        self.eg_tri_lookup[(v,w)] = self.tricnt
        self.eg_tri_lookup[(w,u)] = self.tricnt
        self.tricnt += 1

    # delete a positively oriented triangle u,v,w
    def delete_triangle(self,u,v,w):
        tri = self.eg_tri_lookup[(u,v)]
        if not tri is None:self.triangles[tri] = None
        self.eg_tri_lookup[(u,v)] = None
        self.eg_tri_lookup[(v,w)] = None
        self.eg_tri_lookup[(w,u)] = None

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

    # return vertices v,w such that uvw
    # is a positively oriented triangle
    def adjacent_one(self,u):raise NotImplemented

    # plc is a piecewise linear complex to be tetrahedralized
    def __init__(self,plc):
        self.points = dps.pointset()
        self.triangles = []
        self.tricnt = 0
        self.eg_tri_lookup = {}
        self.ghosts = []
        self.ghostcnt = 0
        self.eg_ghost_lookup = {}
        self.plc = plc
        self.cover(plc)

    def initial_cover(self,plc):
        convexbnd = dpr.pts_to_convex_xy(plc.points.get_points())
        convexcom = dpv.center_of_mass(convexbnd)
        convexrad = max([dpv.distance(cx,convexcom) for cx in convexbnd])+1000.0
        #convexrad = max([dpv.distance(cx,convexcom) for cx in convexbnd])+10
        c01delta = dpv.vector(-1,-1,0).normalize().scale_u(convexrad)
        c02delta = dpv.vector( 1,-1,0).normalize().scale_u(convexrad)
        c03delta = dpv.vector( 0, 1,0).normalize().scale_u(convexrad)
        c01 = convexcom.copy().translate(c01delta)
        c02 = convexcom.copy().translate(c02delta)
        c03 = convexcom.copy().translate(c03delta)
        c0psx = self.points.add_points(c01,c02,c03)
        self.add_triangle(*c0psx)
        ghost1 = (c0psx[2],c0psx[1])
        self.add_ghost(*ghost1)
        ghost2 = (c0psx[1],c0psx[0])
        self.add_ghost(*ghost2)
        ghost3 = (c0psx[0],c0psx[2])
        self.add_ghost(*ghost3)
        self.initial_cover_extras = c0psx
        #ax = self.plc.plot_xy()
        #self.plot_xy(ax)
        #plt.show()

    # generate tetrahedralization of the plc
    def cover(self,plc):
        hmin = plc.chew1_subdivide_edges()
        plc.subdivide_edges()
        self.initial_cover(plc)
        self.cover_points(plc)
        #self.cover_edges(plc)
        self.cover_polygons(plc)
        self.cover_edges(plc)
        self.chew1_refine(plc,hmin)
        #self.ruppert_refine(plc)

    def cover_points(self,plc):
        for plcx in range(plc.points.pcnt):
            plcp = plc.points.ps[plcx].copy()
            self.point_location(plcp)
        #ax = self.plc.plot_xy()
        #self.plot_xy(ax)
        #plt.show()

    # given v1,v2, the positions of the endpoints of an edge, 
    # return True if p encroaches upon the edge
    def encroaches_edge(self,v1,v2,p):
        cc = dpv.midpoint(v1,v2)
        cr = dpv.distance(cc,v1)
        if p.near(v1) or p.near(v2):return False
        if dpr.inside_circle(p,cc,cr,(dpv.zero(),dpv.zhat)):return True
        else:return False

    # given the edge u,v which bounds two non-ghost triangles
    # remove those triangles and replace with the alternative two that are
    # bounded by the same 4 vertices
    def flip_edge(self,u,v):
        print('flipping an edge')
        o1 = self.adjacent(u,v)
        o2 = self.adjacent(v,u)

        vs = self.points.get_points(u,v,o1,o2)
        tcp1,tcr1 = dpr.circumscribe_tri(vs[0],vs[1],vs[2])
        tcp2,tcr2 = dpr.circumscribe_tri(vs[1],vs[0],vs[3])
        if tcp1.near(tcp2) and dtl.isnear(tcr1,tcr2):
            print('4way!',tcp1,tcp2,tcr1,tcr2,u,v,o1,o2)
            return ()

        '''#
        ax = dtl.plot_axes_xy()
        dtl.plot_polygon_xy(self.points.get_points(u,v,o1),ax)
        vs = self.points.get_points(u,v,o1)
        tcp,tcr = dpr.circumscribe_tri(*vs)
        dtl.plot_circle_xy(tcp,tcr,ax,True)
        dtl.plot_polygon_xy(self.points.get_points(v,u,o2),ax)
        vs = self.points.get_points(v,u,o2)
        tcp,tcr = dpr.circumscribe_tri(*vs)
        dtl.plot_circle_xy(tcp,tcr,ax,True)

        self.plc.plot_xy(ax)
        self.plot_xy(ax)

        #plt.show()
        '''#

        self.delete_triangle(u,v,o1)
        self.delete_triangle(v,u,o2)
        self.add_triangle(o1,o2,v)
        self.add_triangle(o2,o1,u)
        #self.dig_cavity(o1,o2,v)
        #self.dig_cavity(o2,o1,u)


        ax = dtl.plot_axes_xy()
        dtl.plot_polygon_xy(self.points.get_points(o1,o2,v),ax,lw = 4.0)
        vs = self.points.get_points(o1,o2,v)
        tcp,tcr = dpr.circumscribe_tri(*vs)
        dtl.plot_circle_xy(tcp,tcr,ax,True)
        dtl.plot_polygon_xy(self.points.get_points(o2,o1,u),ax,lw = 4.0)
        vs = self.points.get_points(o2,o1,u)
        tcp,tcr = dpr.circumscribe_tri(*vs)
        dtl.plot_circle_xy(tcp,tcr,ax,True)

        self.plc.plot_xy(ax)
        self.plot_xy(ax)

        plt.show()

        return (o2,v),(v,o1),(o1,u),(u,o2)

    # given v1,v2, the positions of the endpoints of an edge, 
    # return True if locally delaunay
    def locally_delaunay_edge(self,u,v):
        plcu,plcv = self.plc.points.find_points(*self.points.get_points(u,v))
        if self.plc.segment(plcu,plcv):return True
        o1 = self.adjacent(u,v)
        o2 = self.adjacent(v,u)
        if o1 is None or o2 is None:return True
        if o1 == 'g' or o2 == 'g':return True
        up,vp,op1,op2 = self.points.get_points(u,v,o1,o2)
        if dtl.segments_intersect((up,vp),(op1,op2)):
            if dtl.incircle(up,vp,op1,op2) > 0:return False
            if dtl.incircle(vp,up,op2,op1) > 0:return False
        return True

    # apply the flip algorithm until all edges are locally delaunay
    def cover_edges(self,plc):
        unfinished = [e for e in self.eg_tri_lookup]
        while unfinished:
            pre = unfinished[:]
            unfin = unfinished.pop(0)
            if not self.locally_delaunay_edge(*unfin):
                nedges = self.flip_edge(*unfin)
                unfinished.extend(nedges)

    def cover_polygons(self,plc):
        extras = []
        for tdx in range(self.tricnt):
            tri = self.triangles[tdx]
            if tri is None:continue
            else:u,v,w = tri
            ptri = self.points.get_points(u,v,w)
            extras.append(tdx)
            for p in plc.polygons:
                eb,ibs = p
                ebnd = plc.points.get_points(*[plc.edges[x][0] for x in eb])
                if dtl.concaves_contains(ebnd,ptri):
                    extras.remove(tdx)
                    for ib in ibs:
                        ibndxs = [plc.edges[x][0] for x in ib]
                        ibnd = plc.points.get_points(*ibndxs)
                        if dtl.concaves_contains(ibnd,ptri):
                            extras.append(tdx)
                            break
        for x in extras:
            xtri = self.triangles[x]
            if xtri is None:continue
            x1,x2,x3 = xtri
            self.delete_triangle(x1,x2,x3)
        self.ghost_border(plc)

    # delete all ghosts and add new ghosts according to where they
    # should be after covering the plc
    def ghost_border(self,plc):
        for gst in self.ghosts:
            if gst is None:continue
            g1,g2,g = gst
            self.delete_ghost(g1,g2)
        for plce in plc.edges:
            if plce is None:continue
            plce1,plce2 = plc.points.get_points(*plce)
            e1,e2 = self.points.find_points(plce1,plce2)
            eadj = self.adjacent(e1,e2)
            if eadj is None:self.add_ghost(e1,e2)
            eadj = self.adjacent(e2,e1)
            if eadj is None:self.add_ghost(e2,e1)

    # return a dictionary of the length of 
    # every edge currently in the mesh
    def edge_lengths(self):
        elengths = {}
        for tdx in range(self.tricnt):
            tri = self.triangles[tdx]
            if tri is None:continue
            t1,t2,t3 = tri
            tp1,tp2,tp3 = self.points.get_points(t1,t2,t3)
            if not (t1,t2) in elengths:
                d12 = dpv.distance(tp1,tp2)
                elengths[(t1,t2)] = d12
                elengths[(t2,t1)] = d12
            if not (t2,t3) in elengths:
                d23 = dpv.distance(tp2,tp3)
                elengths[(t2,t3)] = d23
                elengths[(t3,t2)] = d23
            if not (t3,t1) in elengths:
                d31 = dpv.distance(tp3,tp1)
                elengths[(t3,t1)] = d31
                elengths[(t1,t3)] = d31
        return elengths

    def chew1_refine(self,plc,h):
        unfinished = [t for t in self.triangles]
        while unfinished:
            unfin = unfinished.pop(0)
            if unfin is None:continue
            if not unfin in self.triangles:continue
            ufx1,ufx2,ufx3 = unfin
            tcp = self.chew1_skinny_triangle(ufx1,ufx2,ufx3,h)
            if not tcp is None:
                trip = self.points.get_points(ufx1)
                tcppoly = plc.find_polygon(tcp)
                tripoly = plc.find_polygon(*trip)
                #print('chewref',tcppoly,tripoly)
                if not tcppoly is None and tcppoly == tripoly:
                    ntxs = self.point_location(tcp)
                    for ntx in ntxs:unfinished.append(self.triangles[ntx])
                '''#
                if ntxs:
                    ax = self.plot_xy()
                    for ntx in ntxs:
                        tvs = self.points.get_points(*self.triangles[ntx])
                        dtl.plot_polygon_xy(tvs,ax,False,10.0)
                    tvs = self.points.get_points(ufx1,ufx2,ufx3)
                    tcp,tcr = dpr.circumscribe_tri(*tvs)
                    dtl.plot_circle_xy(tcp,tcr,ax,True)
                    dtl.plot_polygon_xy(tvs,ax,True,5.0)
                    plt.show()
                '''#

    # if triangle uvw is skinny by chew1 standards, 
    # return the center of its circumcircle otherwise return None
    def chew1_skinny_triangle(self,u,v,w,h):
        vs = self.points.get_points(u,v,w)
        tcp,tcr = dpr.circumscribe_tri(*vs)
        if tcr/h > 1.0:return tcp 

    def chew2_refine(self,plc,b = 2.0):
        raise NotImplemented

    # if triangle uvw is skinny by ruppert standards, 
    # return the center of its circumcircle otherwise return None
    def ruppert_skinny_triangle(self,u,v,w,b = 2.0):
        vs = self.points.get_points(u,v,w)
        tcp,tcr = dpr.circumscribe_tri(*vs)
        if tcr/dtl.shortest_edge_tri(*vs) > b:return tcp

    def ruppert_refine(self,plc,b = 2.0):
        unfinished = [e for e in plc.edges]
        while unfinished:
            unfin = unfinished.pop(0)
            if unfin is None:continue
            v1,v2 = plc.points.get_points(*unfin)
            isld = self.locally_delaunay_edge(v1,v2)
            if not isld:
                ne1,ne2 = plc.split_edge(*unfin)
                unfinished.append(ne1)
                unfinished.append(ne2)
                newp = plc.points.ps[ne1[1]]
                self.point_location(newp)
        unfinished = [t for t in self.triangles]
        while unfinished:
            unfin = unfinished.pop(0)
            if unfin is None:continue
            tcp = self.ruppert_skinny_triangle(*unfin)
            if not tcp is None:
                print('refinnnning skinny guy!',tcp)
                dodig = True
                for e in plc.edges:
                    if e is None:continue
                    v1,v2 = plc.points.get_points(*e)
                    ench = self.encroaches_edge(v1,v2,tcp)
                    if ench:
                        dodig = False
                        ne1,ne2 = plc.split_edge(*e)
                        newp = plc.points.ps[ne1[1]]
                        self.point_location(newp)
                if dodig:
                    ntxs = self.point_location(tcp)
                    for ntx in ntxs:unfinished.append(self.triangles[ntx])

        print('ruppert refined the mesh!')
        ax = self.plot_xy()
        plt.show()

    def point_location(self,y):
        pretricnt = self.tricnt
        for pdx in range(self.points.pcnt):
            if self.points.ps[pdx].near(y):
                return ()
        nv = self.points.add_point(y)
        onb = self.point_on_boundary(nv)
        if onb:
            v,w,x = self.ghosts[onb]
            self.delete_ghost(v,w)
            self.add_ghost(v,nv)
            self.add_ghost(nv,w)
            tu = self.adjacent(w,v)
            self.delete_triangle(w,v,tu)
            self.add_triangle(tu,nv,v)
            self.add_triangle(tu,w,nv)
            #self.dig_cavity(tu,nv,v)
            #self.dig_cavity(tu,w,nv)
            return [x for x in range(pretricnt,self.tricnt)]
        for tdx in range(self.tricnt):
            tri = self.triangles[tdx]
            if tri is None:continue
            else:u,v,w = tri
            vu,vv,vw = self.points.get_points(u,v,w)
            if dpv.inside(y,[vu,vv,vw]):
                self.insert_vertex(nv,*self.triangles[tdx])
                return [x for x in range(pretricnt,self.tricnt)]
        for gdx in range(self.ghostcnt):
            ghost = self.ghosts[gdx]
            if ghost is None:continue
            else:u,v,w = ghost
            vu,vv = self.points.get_points(u,v)
            if not dtl.orient2d(vu,vv,y) < 0:
                self.insert_ghost_vertex(nv,u,v,w)
                return [x for x in range(pretricnt,self.tricnt)]

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

    # u is the vertex to insert. vwx is a positively oriented triangle whose
    # circumcircle encloses u
    def insert_vertex(self,u,v,w,x):
        self.delete_triangle(v,w,x)
        self.dig_cavity(u,v,w)
        self.dig_cavity(u,w,x)
        self.dig_cavity(u,x,v)

    # u is a new vertex; is the oriented triangle u,v,w delaunay?
    def dig_cavity(self,u,v,w):
        # find triangle wvx opposite the facet vw from u
        x = self.adjacent(w,v)
        if x is None:return
        elif x == 'g':self.add_triangle(u,v,w)
        else:
            vu,vv,vw,vx = self.points.get_points(u,v,w,x)
            if dtl.incircle(vu,vv,vw,vx) > 0:
                self.delete_triangle(w,v,x)
                self.dig_cavity(u,v,x)
                self.dig_cavity(u,x,w)
            else:
                # w,v is a facet of the cavity and uvw is delaunay
                self.add_triangle(u,v,w)

    # u is the vertex to insert. vwg is a positively oriented ghost triangle whose
    # circumcircle encloses u
    def insert_ghost_vertex(self,u,v,w,x):
        if not x == 'g':raise ValueError
        self.delete_ghost(v,w)
        self.add_ghost(v,u)
        self.add_ghost(u,w)
        onb = self.point_on_boundary(u)
        if onb is None:self.add_triangle(u,v,w)

    def pelt(self):
        s = dmo.model()
        for f in self.triangles:
            if f is None:continue
            v1,v2,v3 = self.points.get_points(*f)
            s._triangle(v1,v2,v3)
        return s









