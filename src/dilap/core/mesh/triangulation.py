import dilap.core.tools as dpr

import dilap.core.mesh.pointset as dps
import dilap.core.mesh.tools as dtl
import dilap.core.mesh.pointset as dps

import dp_vector as dpv

import matplotlib.pyplot as plt
import pdb

class triangulation:
    
    def plot(self,ax = None):
        if ax is None:ax = plt.figure().add_subplot(111)
        for tdx in range(self.tricnt):
            tri = self.triangles[tdx]
            if tri is None:continue
            vtri = self.points.get_points(*tri)
            dtl.plot_polygon_xy(vtri,ax,center = True)
        return ax

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

    # generate tetrahedralization of the plc
    def cover(self,plc):
        c01 = dpv.vector(-70,-50,0)
        c02 = dpv.vector( 70,-50,0)
        c03 = dpv.vector(  0, 50,0)
        c0psx = self.points.add_points(c01,c02,c03)
        self.add_triangle(*c0psx)

        ax = plc.plot_xy()
        self.plot(ax)
        plt.show()

        ghost1 = (c0psx[2],c0psx[1])
        self.add_ghost(*ghost1)
        ghost2 = (c0psx[1],c0psx[0])
        self.add_ghost(*ghost2)
        ghost3 = (c0psx[0],c0psx[2])
        self.add_ghost(*ghost3)

        self.cover_points(plc)
        self.cover_edges(plc)
        
        for cx in c0psx:
            for tx in range(self.tricnt):
                tri = self.triangles[tx]
                if tri is None:continue
                if cx in tri:self.delete_triangle(*tri)

        ax = plc.plot_xy()
        self.plot(ax)
        plt.show()

        self.cover_polygons(plc)

        ax = plc.plot_xy()
        self.plot(ax)
        plt.show()

        pdb.set_trace()

    def cover_points(self,plc):
        for plcx in range(plc.points.pcnt):
            plcp = plc.points.ps[plcx].copy()
            nv = self.points.add_point(plcp)
            tloc = self.point_location(plcp)
            self.insert_vertex(nv,*self.triangles[tloc])

            #ax = plc.plot_xy()
            #self.plot(ax)
            #plt.show()

    def cover_edges(self,plc):

        def locally_delaunay(eg):
            v1,v2 = plc.points.get_points(*eg)
            cc = dpv.vector((v1.x+v2.x)/2.0,(v1.y+v2.y)/2.0,(v1.z+v2.z)/2.0)
            cr = dpv.distance(cc,v1)
            for p in plc.points:
                if p is v1 or p is v2:continue
                if dpr.inside_circle(p,cc,cr,(dpv.zero(),dpv.zhat)):
                    #ax = plc.plot_xy()
                    #dtl.plot_circle_xy(cc,cr,ax)
                    #dtl.plot_points_xy([p,v1,v2],ax)
                    #plt.show()
                    return False,cc
            #ax = plc.plot_xy()
            #dtl.plot_circle_xy(cc,cr,ax)
            #dtl.plot_points_xy([v1,v2],ax)
            #plt.show()
            return True,None

        unfinished = [e for e in plc.edges]
        new = []
        while unfinished:
            unfin = unfinished.pop(0)
            isld,newp = locally_delaunay(unfin)
            if not isld:
                ne1,ne2 = plc.split_edge(*unfin)
                unfinished.append(ne1)
                unfinished.append(ne2)
                new.append(ne1)
                new.append(ne2)
                newp = plc.points.ps[ne1[1]]
                nv = self.points.add_point(newp)
                tloc = self.point_location(newp)
                self.insert_vertex(nv,*self.triangles[tloc])

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
                if dtl.concaves_intersect(ebnd,ptri):
                    extras.remove(tdx)
                    for ib in ibs:
                        ibndxs = [plc.edges[x][0] for x in ib]
                        ibnd = plc.points.get_points(*ibndxs)
                        if dtl.concaves_intersect(ibnd,ptri):
                            extras.append(tdx)
                            break
        for x in extras:
            xtri = self.triangles[x]
            if xtri is None:continue
            self.delete_triangle(*xtri)

    def point_location(self,y):
        for tdx in range(self.tricnt):
            tri = self.triangles[tdx]
            if tri is None:continue
            else:u,v,w = tri
            vu,vv,vw = self.points.get_points(u,v,w)
            inc = dtl.incircle(vu,vv,vw,y)
            #if dtl.incircle(vu,vv,vw,y) > 0:
            if dpv.inside(y,[vu,vv,vw]):
                return tdx

    # add a positively oriented ghost triangle u,v,g
    def add_ghost(self,u,v):
        self.ghosts.append((u,v,'g'))
        self.eg_ghost_lookup[(u,v)] = self.ghostcnt
        self.ghostcnt += 1

    # delete a positively oriented ghost triangle u,v,g
    def delete_ghost(self,u,v):
        ghost = self.eg_ghost_lookup[(u,v)]
        self.ghosts[ghost] = None
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
            if tri is None:return
            else:
                triv = [x for x in self.triangles[tri] if not x in ekey][0]
                return triv
        elif ekey in self.eg_ghost_lookup:
            tri = self.eg_ghost_lookup[(u,v)]
            if tri is None:return
            else:return self.ghosts[tri][2]
        else:return

    # return vertices v,w such that uvw
    # is a positively oriented triangle
    def adjacent_one(self,u):raise NotImplemented

    # to add a vertex:
    # 1) find one tetrahedon whose circumsphere encloses v (point location)
    # 2) a depth-first search in triangulation finds all other tetrahedra whose
    #       circumspheres enclose v
    # 3) delete these tetrahedra; the union of the deleted tetrahedra 
    #       is a star-shaped polyhedral cavity

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










