import dilap.core.base as db
import dilap.core.tools as dpr
import dilap.core.vector as dpv
import dilap.core.model as dmo

import dilap.mesh.piecewisecomplex as pwc
import dilap.mesh.tools as dtl

import dilap.structures.tools as dstl
import dilap.structures.structgraph as dst

import matplotlib.pyplot as plt
import pdb



class structure(db.base):

    def __init__(self,boundary,**kwargs):
        self.boundary = boundary
        self.graph = dst.graph(boundary)
        self.plc = pwc.piecewise_linear_complex()
        #self._def('wallwidth',0.75,**kwargs)

    # return a local space model representing this structure
    def model(self,t,q):
        self.plc.add_polygons(*self.graph.model())

        ax = self.graph.plot()
        ax = self.plc.plot(ax)
        plt.show()

        if not self.plc.polygons:self.plc.add_polygons((self.boundary,()))
        self.plc.triangulate()
        pelt = self.plc.pelt()
        pelt._project_uv_flat()
        pelt.rotate(q).translate(t)
        return pelt






