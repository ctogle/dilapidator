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
            dtl.plot_polygon_xy(vtri,ax)
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
        self.cover(plc)

    # generate tetrahedralization of the plc
    def cover(self,plc):
        c01 = dpv.vector(-50,-50,0)
        c02 = dpv.vector(  0, 50,0)
        c03 = dpv.vector( 50,-50,0)
        c0psx = self.points.add_points(c01,c02,c03)
        self.add_triangle(*c0psx)

        ax = plc.plot_xy()
        self.plot(ax)
        plt.show()

        #ghost1 = (c0psx[2],c0psx[1],c0psx[0])
        #self.add_ghost_tetrahedron(*ghost1)
        #ghost2 = (c0psx[3],c0psx[2],c0psx[0])
        #self.add_ghost_tetrahedron(*ghost2)
        #ghost3 = (c0psx[1],c0psx[2],c0psx[3])
        #self.add_ghost_tetrahedron(*ghost3)
        #ghost4 = (c0psx[0],c0psx[1],c0psx[3])
        #self.add_ghost_tetrahedron(*ghost4)

        for plcx in range(plc.points.pcnt):
            plcp = plc.points.ps[plcx].copy()
            tloc = self.point_location(plcp)
            if tloc is None:
                print('shit')
                pdb.set_trace()
            else:
                nv = self.points.add_point(plcp)
                self.insert_vertex(nv,*self.triangles[tloc])
                print('inserted vertex',nv-3)

            ax = plc.plot_xy()
            self.plot(ax)
            plt.show()

        for cx in c0psx:
            for tx in range(self.tricnt):
                tri = self.triangles[tx]
                if tri is None:continue
                if cx in tri:self.delete_triangle(*tri)

        self.plot()
        plt.show()

        pdb.set_trace()

    # for a point y, find a tetrahedron whose circumsphere encloses y
    # return the index of the first found acceptable tetrahedron
    # return None if no such tetrahedron is found
    def point_location(self,y):
        for tdx in range(self.tricnt):
            tri = self.triangles[tdx]
            if tri is None:continue
            else:u,v,w = tri
            vu,vv,vw = self.points.get_points(u,v,w)
            if dtl.incircle(vu,vv,vw,y):return tdx
        # find a ghost triangle instead

    # add a positively oriented ghost triangle u,v,g
    def add_ghost_triangle(self,u,v):
        self.ghosts.append((u,v,'g'))
        self.eg_ghost_lookup[(u,v)] = self.ghostcnt
        self.ghostcnt += 1

    # delete a positively oriented ghost triangle u,v,g
    def delete_ghost_triangle(self,u,v):
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
        self.triangles[tri] = None
        self.eg_tri_lookup[(u,v)] = None
        self.eg_tri_lookup[(v,w)] = None
        self.eg_tri_lookup[(w,u)] = None

    # return a vertex x such that uv
    # is a positively oriented edge
    def adjacent(self,u,v):
        ekey = (u,v)
        if not ekey in self.eg_tri_lookup:return
        tri = self.eg_tri_lookup[(u,v)]
        if tri is None:return
        else:return self.triangles[tri][2]

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
        # find triangle wvx opposite the facet vwx from u
        x = self.adjacent(w,v)
        print('adjacent',w,v,'is',x)
        if x is None:#pass # do nothing if the tetrahedron was deleted
            self.add_triangle(u,v,w)
        else:
            vu,vv,vw,vx = self.points.get_points(u,v,w,x)
            if dtl.incircle(vu,vv,vw,vx):
                # tetrahedra uvwx and vwxy are not delaunay
                self.delete_triangle(w,v,x)
                print('i just deleted',w,v,x)
                print('i will consider',u,v,x)
                print('i will consider',u,x,w)
                self.dig_cavity(u,v,x)
                self.dig_cavity(u,x,w)
            else:
                # w,v is a facet of the cavity and uvw is delaunay
                self.add_triangle(u,v,w)


