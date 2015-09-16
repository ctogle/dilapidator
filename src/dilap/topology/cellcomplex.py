import dilap.core.base as db
import dilap.core.vector as dpv

import dilap.mesh.tools as dtl

import dilap.topology.volume as dtv

cix = 0
def index():
    global cix
    cix +=1 
    return cix

class cellcomplex(db.base):

    def plot(self,ax = None):
        if ax is None:ax = dtl.plot_axes()
        for vdx in range(self.volcount):
            vol = self.volumes[vdx]
            if vol is None:continue
            vol.plot(self.geom,ax)
        return ax

    def new_volume(self):
        nvol = dtv.volume()
        nv,nl,nf,ns = nvol.new_shell()
        self.volumes.append(nvol)
        self.volcount += 1
        return nv,nl,nf,ns,nvol

    def __init__(self,geom):
        self.ix = index()
        self.geom = geom
        self.volumes = []
        self.volcount = 0

    def cube(self):
        nv1,nl,nf,ns,nvl = self.new_volume()
        nv2,ev1 = nvl.new_edge_vertex(nv1)
        #vx1 = nvl.new_vertex()
        #vx2,ex1 = nvl.new_edge_vertex(vx1)
        #vx3,ex2 = nvl.new_edge_vertex(vx2)
        #vx4,ex3 = nvl.new_edge_vertex(vx3)
        #ex4 = nvl.new_edge(vx4,vx1)
        #lx1 = nvl.new_loop()

        p1 = dpv.vector(0,0,0)
        p2 = dpv.vector(1,0,0)
        p3 = dpv.vector(1,1,0)
        p4 = dpv.vector(0,1,0)
        self.geom.add_point(nv1.ix,p1)
        self.geom.add_point(nv2.ix,p2)
        #self.geom.add_point(vx2.ix,p2)
        #self.geom.add_point(vx3.ix,p3)
        #self.geom.add_point(vx4.ix,p4)

        return self



