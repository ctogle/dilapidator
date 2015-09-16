import dilap.core.base as db

import dilap.mesh.tools as dtl

vix = 0
def index():
    global vix
    vix +=1 
    return vix

class vertex(db.base):

    def plot(self,geom,ax):
        p = geom.point(self.ix)
        dtl.plot_point(p,ax)
        return ax

    def __init__(self):
        self.ix = index()
        self.hedge = None

        self.loop = None
        self.shell = None
        self.cplex = None



