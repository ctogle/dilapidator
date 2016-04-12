import dilap.core.vector as dpv
import dilap.core.quaternion as dpq
import dilap.core.tools as dpr
import dilap.core.context as dgc
import dilap.core.profiler as dprf
import dilap.mesh.piecewisecomplex as pwc
import dilap.mesh.tools as dtl
import dilap.infrastructure.graphregion as grg
import dilap.infrastructure.infragraph as ifg
import dilap.generate.landscape as dls
import dilap.generate.city as dcy
import dilap.generate.area as dar

import dilap.generate.lot as dlt
import dilap.primitive.cylinder as dcyl

import matplotlib.pyplot as plt

# a continent should be a top level context for dilap
# it creates a full world (island), with a true boundary (ocean)
# given a collection of regions, it provides the earth and sea
class continent(dgc.context):

    def __init__(self,*args,**kwargs):
        dgc.context.__init__(self,*args,**kwargs)
        self._def('sealevel',-0.5,**kwargs)

    def layout_infrastructure(self):
        #g = ifg.graph()
        #g = ifg.hairpin()
        #g = ifg.circle()
        #g = ifg.newcastle()
        #g = ifg.eight()
        g = ifg.clover()
        #g = ifg.ramp()

        g._update()
        self.igraph = g

    # first generate the graph of infrastructure
    # this also yields area polygons and for the infrastructure area
    # this covers the entire ground of the continent
    #
    # construction of the road happens within the bounds of the infrastructure polygons
    # construction of each area happens within their respective bounds
    def generate(self,worn = 0):
        self.layout_infrastructure()
        for tpoly in self.igraph.tpolygons:
            tarea = dar.area(boundary = tpoly)
            self._consume(tarea.generate())

        rplc = pwc.piecewise_linear_complex(refine = True,smooth = True)
        rplc.add_polygons(*self.igraph.rpolygons)
        rplc.triangulate()
        rpelt = rplc.pelt()
        rnode = self._node_wrap(rpelt)
        self._nodes_to_graph(rnode)

        # add water models to scenegraph
        '''#
        wl = lscape.landbb.x.y - lscape.landbb.x.x + 100.0
        ww = lscape.landbb.y.y - lscape.landbb.y.x + 100.0
        wr = max([wl,ww])
        water = dcyl.cylinder(n = 8).translate_z(-0.5)
        water.scale_x(wr).scale_y(wr).scale_z(20)
        water.translate_z(self.sealevel)
        water.translate(lscape.landbb._center().xy())
        wnode = self._node_wrap(water)
        self._nodes_to_graph(wnode)
        '''#
        return self



