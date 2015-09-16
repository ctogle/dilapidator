import dilap.core.base as db
import dilap.core.pointset as dps

import dilap.topology.cellcomplex as dtc

import matplotlib.pyplot as plt

class geometry(db.base):

    def radius(self):
        pmags = []
        for x in range(self.points.pcnt):
            pmags.append(self.points.ps[x].magnitude())
        return max(pmags)
    
    # given the index of a vertex, return its geometric point
    def point(self,vx):
        return self.points.ps[self.plook[vx]]

    def add_point(self,vx,p):
        px = self.points.add_point(p)
        self.plook[vx] = px

    def __init__(self,*args,**kwargs):
        self.points = dps.pointset()
        self.plook = {}

        self.curves = []
        self.surfaces = []

def test():
    geom = geometry()
    cplex = dtc.cellcomplex(geom).cube()
    print('cplex!',cplex)
    cplex.plot()
    plt.show()

test()

