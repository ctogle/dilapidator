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
        self.define()

    def define(self):
        #g = ifg.graph()
        #g = ifg.hairpin()
        g = ifg.circle()
        #g = ifg.newcastle()
        #g = ifg.eight()
        #g = ifg.clover()
        #g = ifg.ramp()

        g._update()
        self.igraph = g

        cityseed = dpv.vector(0,0,150)
        self.cityseeds = [cityseed]
        #self.cityseeds.append(dpv.vector(250,250,100))

    def generate(self,worn = 0):
        '''#
        cities = []
        tpts = []
        hpts = []
        rpts = []

        for cd in self.cityseeds:
            cities.append(dcy.city().generate(cd,self.igraph,worn))
        for cy in cities:
            self._consume(cy)
            tpts.extend(cy._terrain_points())
            hpts.extend(cy._hole_points())
            rpts.extend(cy._region_points())
        '''#
        
        #lscape = dls.landscape(controls = tpts,holes = hpts,regions = rpts)
        #lscape.generate(worn)
        #self._consume(lscape.generate(worn))
        #tplc = pwc.piecewise_linear_complex()
        #tplc.add_polygons(*self.igraph.tpolygons)
        #dprf.profile_function(tplc.triangulate)

        for tpoly in self.igraph.tpolygons[1:]:
            tarea = dar.area(boundary = tpoly)
            self._consume(tarea.generate())

        rplc = pwc.piecewise_linear_complex()
        rplc.add_polygons(*self.igraph.rpolygons)
        #rplc.triangulate()
        #rpelt = rplc.pelt()
        #rnode = self._node_wrap(rpelt)
        #self._nodes_to_graph(rnode)

        #opelt = dtl.facade().pelt()
        #onode = self._node_wrap(opelt)
        #self._nodes_to_graph(onode)

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



