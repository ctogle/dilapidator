#imports
# cython: profile=True
#cimport cython
cimport dilap.core.tools as dpr
cimport dilap.core.pointset as dps
cimport dilap.core.vector as dpv
cimport dilap.core.quaternion as dpq

import math,numpy

cdef class tetrahedralization:

    '''#
    # given the index of a point, return the index of a ghost
    # whose real edge intersects the point, or -1
    cdef int point_on_boundary(self,int u):
        cdef dpv.vector up = self.points.ps[u]
        cdef dpv.vector g1,g2
        cdef int gdx,gx1,gx2,gx3
        cdef float dx,dy,dv,nx,ny,prj1,prj2,prj3
        for gdx in range(self.ghostcnt):
            gst = self.ghosts[gdx]
            if gst is None:continue
            gx1,gx2,gx3 = gst
            g1,g2 = self.points.get_points(gx1,gx2)
            dx = g2.x - g1.x
            dy = g2.y - g1.y
            dv = math.sqrt(dx**2 + dy**2)
            nx =  dy/dv
            ny = -dx/dv
            prj1 = g1.x*nx + g1.y*ny
            prj3 = up.x*nx + up.y*ny
            if dpr.isnear_c(prj1,prj3):
                prj1 = -g1.x*ny + g1.y*nx
                prj3 = -up.x*ny + up.y*nx
                if dpr.near_c(prj3,prj1) >= prj1:
                    prj2 = -g2.x*ny + g2.y*nx
                    if dpr.near_c(prj3,prj2) <= prj2:
                        return gdx
        return -1

    # given the indices of the endpoints of an edge, 
    # return 1 if locally delaunay, 0 otherwise
    cdef bint locally_delaunay(self,int u,int v):
        o1 = self.adjacent(u,v)
        o2 = self.adjacent(v,u)
        if o1 == -1 or o2 == -1:return 1
        if o1 == -2 or o2 == -2:return 1
        up,vp,op1,op2 = self.points.get_points(u,v,o1,o2)
        if dpr.segments_intersect_c(up,vp,op1,op2):
            if dpr.incircle_c(up,vp,op1,op2) > 0:return 0
            if dpr.incircle_c(vp,up,op2,op1) > 0:return 0
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
    '''#

    # return a vertex x such that uvwx 
    # is a positively oriented tetrahedron
    def adjacent(self,u,v,w):
        tkey = (u,v,w)
        if not tkey in self.tri_tetra_lookup:return
        tetra = self.tri_tetra_lookup[(u,v,w)]
        if tetra is None:return
        else:
            return self.tetrahedrons[tetra][3]

    def __cinit__(self,dpv.vector p0,dpv.vector pn):
        self.p0,self.pn = p0,pn
        self.points = dps.pointset()
        self.tetrahedrons = []
        self.tetracnt = 0
        self.tri_tetra_lookup = {}
        self.ghosts = []
        self.ghostcnt = 0
        self.tri_ghost_lookup = {}

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
            if dpr.insphere(vu,vv,vw,vx,vy) > 0:
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

#
cdef tuple tetrahedralize_c():
    cdef list smps = []
    cdef list gsts = []
    return smps,gsts

#
cpdef tuple tetrahedralize():
    return tetrahedralize_c()










