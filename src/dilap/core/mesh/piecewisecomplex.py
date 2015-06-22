import dilap.core.mesh.tools as dtl
import dilap.core.mesh.pointset as dps
import dilap.core.mesh.tetrahedralization as dth
import dilap.core.mesh.triangulation as dtg

import matplotlib.pyplot as plt
import pdb

class piecewise_linear_complex:

    def plot_xy(self,ax = None):
        ax = dtl.plot_points_xy(self.points.ps,ax)
        for edx in range(len(self.edges)):
            veg = self.points.get_points(*self.edges[edx])
            ax = dtl.plot_edges_xy(veg,ax)
        return ax

    def plot(self,ax = None):
        ax = plt.figure().add_subplot(111,projection = '3d')
        for pdx in range(self.points.pcnt):
            p = self.points.ps[pdx]
            ax.plot([p.x],[p.y],[p.z],marker = '+')
        for edx in range(len(self.edges)):
            up,vp = self.points.get_points(*self.edges[edx])
            ax.plot([up.x,vp.x],[up.y,vp.y],zs = [up.z,vp.z])
        return ax

    def __init__(self):
        self.points = dps.pointset()
        self.edges = []
        self.polygons = []
        self.polyhedra = []
        self.covers = {}

    def add_points(self,*nps):
        pst = self.points.pcnt
        self.points.add_points(*nps)
        return range(pst,self.points.pcnt)

    def connect(self,u,v):
        self.edges.append((u,v))

    def polygon(self,ebnd,*ibnds):
        self.polygons.append((ebnd,ibnds))

    def tetrahedralize(self):
        tetra = dth.tetrahedralization(self)
        self.covers['tetra'] = tetra

    def triangulate_xy(self):
        tri = dtg.triangulation(self)
        self.covers['tri'] = tri



