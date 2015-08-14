import dilap.core.vector as dpv
import dilap.mesh.tools as dtl
import dilap.mesh.pointset as dps
                      
import matplotlib.pyplot as plt
import pdb

class tetrahedralization:

    def plot(self,ax = None):
        if ax is None:ax = plt.figure().add_subplot(111,projection = '3d')
        for tdx in range(self.tetracnt):
            tet = self.tetrahedrons[tdx]
            if tet is None:continue
            vu,vv,vw,vx = self.points.get_points(*tet)
            ax.plot([vu.x,vv.x,vw.x,vx.x],[vu.y,vv.y,vw.y,vx.y],
                        zs = [vu.z,vv.z,vw.z,vx.z],marker = 'o')
            ax.plot([vw.x,vu.x,vx.x,vv.x],[vw.y,vu.y,vx.y,vv.y],
                        zs = [vw.z,vu.z,vx.z,vv.z],marker = 'o')
        return ax

    # plc is a piecewise linear complex to be tetrahedralized
    def __init__(self,plc):
        self.points = dps.pointset()
        self.tetrahedrons = []
        self.tetracnt = 0
        self.tri_tetra_lookup = {}
        self.ghosts = []
        self.ghostcnt = 0
        self.tri_ghost_lookup = {}
        self.cover(plc)

    # generate tetrahedralization of the plc
    def cover(self,plc):
        c01 = dpv.vector(-50,-50,-50)
        c02 = dpv.vector(  0, 50,-50)
        c03 = dpv.vector( 50,-50,-50)
        c04 = dpv.vector(  0,  0, 50)
        c0psx = self.points.add_points(c01,c02,c03,c04)
        self.add_tetrahedron(*c0psx)

        ax = plc.plot()
        self.plot(ax)
        plt.show()

        ghost1 = (c0psx[2],c0psx[1],c0psx[0])
        self.add_ghost_tetrahedron(*ghost1)
        ghost2 = (c0psx[3],c0psx[2],c0psx[0])
        self.add_ghost_tetrahedron(*ghost2)
        ghost3 = (c0psx[1],c0psx[2],c0psx[3])
        self.add_ghost_tetrahedron(*ghost3)
        ghost4 = (c0psx[0],c0psx[1],c0psx[3])
        self.add_ghost_tetrahedron(*ghost4)

        for plcx in range(plc.points.pcnt):
            plcp = plc.points.ps[plcx].copy()
            tloc = self.point_location(plcp)
            if tloc is None:
                print('shit')
                pdb.set_trace()
            else:
                nv = self.points.add_point(plcp)
                self.insert_vertex(nv,*self.tetrahedrons[tloc])
                print('inserted vertex',nv-4)

            ax = plc.plot()
            self.plot(ax)
            plt.show()

        for cx in c0psx:
            for tx in range(self.tetracnt):
                tet = self.tetrahedrons[tx]
                if tet is None:continue
                if cx in tet:self.delete_tetrahedron(*tet)

        self.plot()
        plt.show()

        pdb.set_trace()

    # for a point y, find a tetrahedron whose circumsphere encloses y
    # return the index of the first found acceptable tetrahedron
    # return None if no such tetrahedron is found
    def point_location(self,y):
        for tdx in range(self.tetracnt):
            tetra = self.tetrahedrons[tdx]
            if tetra is None:continue
            else:u,v,w,x = tetra
            vu,vv,vw,vx = self.points.get_points(u,v,w,x)
            if dtl.insphere(vu,vv,vw,vx,y) > 0:return tdx
        # find a ghost tetrahedron instead

    # add a positively oriented ghost tetrahedron u,v,w,x
    def add_ghost_tetrahedron(self,u,v,w):
        self.ghosts.append((u,v,w,'g'))
        self.tri_ghost_lookup[(u,v,w)] = self.ghostcnt
        self.ghostcnt += 1

    # delete a positively oriented ghost tetrahedron u,v,w,x
    def delete_ghost_tetrahedron(self,u,v,w):
        ghost = self.tri_ghost_lookup[(u,v,w)]
        self.ghosts[ghost] = None
        self.tri_ghost_lookup[(u,v,w)] = None

    # add a positively oriented tetrahedron u,v,w,x
    def add_tetrahedron(self,u,v,w,x):
        self.tetrahedrons.append((u,v,w,x))
        self.tri_tetra_lookup[(u,v,w)] = self.tetracnt
        self.tri_tetra_lookup[(u,x,v)] = self.tetracnt
        self.tri_tetra_lookup[(u,w,x)] = self.tetracnt
        self.tri_tetra_lookup[(v,x,w)] = self.tetracnt
        self.tetracnt += 1

    # delete a positively oriented tetrahedron u,v,w,x
    def delete_tetrahedron(self,u,v,w,x):
        tetra = self.tri_tetra_lookup[(u,v,w)]
        self.tetrahedrons[tetra] = None
        self.tri_tetra_lookup[(u,v,w)] = None
        self.tri_tetra_lookup[(u,x,v)] = None
        self.tri_tetra_lookup[(u,w,x)] = None
        self.tri_tetra_lookup[(v,x,w)] = None

    # return a vertex x such that uvwx 
    # is a positively oriented tetrahedron
    def adjacent(self,u,v,w):
        tkey = (u,v,w)
        if not tkey in self.tri_tetra_lookup:return
        tetra = self.tri_tetra_lookup[(u,v,w)]
        if tetra is None:return
        else:
            try:return self.tetrahedrons[tetra][3]
            except:pdb.set_trace()

    # return vertices v,w,x such that uvwx 
    # is a positively oriented tetrahedron
    def adjacent_one(self,u):raise NotImplemented

    # to add a vertex:
    # 1) find one tetrahedon whose circumsphere encloses v (point location)
    # 2) a depth-first search in triangulation finds all other tetrahedra whose
    #       circumspheres enclose v
    # 3) delete these tetrahedra; the union of the deleted tetrahedra 
    #       is a star-shaped polyhedral cavity

    # u is the vertex to insert. vwxy is a positively oriented tetrahedon whose
    # circumsphere encloses u
    def insert_vertex(self,u,v,w,x,y):
        self.delete_tetrahedron(v,w,x,y)
        self.consider_tetrahedron(u,x,w,v)
        self.consider_tetrahedron(u,y,v,w)
        self.consider_tetrahedron(u,v,y,x)
        self.consider_tetrahedron(u,w,x,y)

    # u is a new vertex; is the oriented tetrahedron u,v,w,x delaunay?
    def consider_tetrahedron(self,u,v,w,x):
        # find tetrahedon vwxy opposite the facet vwx from u
        y = self.adjacent(v,w,x)
        print('adjacent',v,w,x,'is',y)
        if y is None:#pass # do nothing if the tetrahedron was deleted
            self.add_tetrahedron(u,v,w,x)
        else:
            vu,vv,vw,vx,vy = self.points.get_points(u,v,w,x,y)
            if dtl.insphere(vu,vv,vw,vx,vy) > 0:
                # tetrahedra uvwx and vwxy are not delaunay
                self.delete_tetrahedron(v,w,x,y)
                print('i just deleted',v,w,x,y)
                print('i will consider',u,v,w,y)
                print('i will consider',u,w,x,y)
                print('i will consider',u,x,v,y)

                self.consider_tetrahedron(u,v,w,y)
                self.consider_tetrahedron(u,w,x,y)
                self.consider_tetrahedron(u,x,v,y)
            else:
                # vwx is a facet of the cavity and uvwx is delaunay
                self.add_tetrahedron(u,v,w,x)


