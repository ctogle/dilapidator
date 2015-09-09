import dilap.core.base as db
import dilap.core.vector as dpv
import dilap.core.model as dmo
import dilap.core.lsystem as dls

import dilap.mesh.tools as dtl

import dilap.graph.graph as dgg

import matplotlib.pyplot as plt

import pdb

# a plant combines an lsystem with a tetrahedralization
# to represent an organic object such as a tree
class plant(db.base):

    def branchdraw(self,p,q):
        self.graph.add_edge(p.to_tuple(),q.to_tuple())
    def leafdraw(self,q):
        self.graph.add_node(q.to_tuple())
    def finaldraw(self):
        ax = self.graph.plot()
        
        # create a plc with a cloudlike pointset containing
        # the boundary of the solid
        # tetrahedralize the pointset
        # retain the triangulated bounding surface of the tetrahedralization

        m = dtl.box(5,5,5)
        m.triangulate()
        smod = m.pelt()
        self.model = smod
        m.plot(ax)

        plt.show()

    def __init__(self,*args,**kwargs):
        kwargs['branchdraw'] = self.branchdraw
        kwargs['leafdraw'] = self.leafdraw
        kwargs['finaldraw'] = self.finaldraw
        self.ltree = dls.ltree(-1,**kwargs)
        self.graph = dgg.graph()

    def model(self):
        p = dpv.zero()
        d = dpv.z()
        #ltree(-1)._realize(p,d)
        total = dmo.model()
        #tps = [(x,y) for x in range(3) for y in range(3)]
        tps = [(x,y) for x in range(1) for y in range(1)]
        for i,xy in enumerate(tps):
            x,y = xy
            kws = {'seed':i}
            self.ltree._realize(p,d)
            lmod = self.model
            lmod.translate(dpv.vector(10*x,10*y,0))
            total._consume(lmod)
        return total





