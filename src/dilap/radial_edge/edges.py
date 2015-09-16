import dilap.core.base as db

import dilap.mesh.tools as dtl

eix = 0
def eindex():
    global eix
    eix +=1 
    return eix

heix = 0
def heindex():
    global heix
    heix +=1 
    return heix

class edge(db.base):

    def __init__(self):
        self.ix = eindex()
        self.one = None
        self.two = None

class halfedge(db.base):

    def plot(self,geom,ax):
        p1 = geom.point(self.tail.ix)
        p2 = geom.point(self.tip.ix)
        dtl.plot_edges([p1,p2],ax)
        return ax

    def __init__(self):
        self.ix = heindex()
        self.tip = None
        self.tail = None
        self.face = None
        self.next_hedge = None

        self.shell = None
        self.cplex = None



