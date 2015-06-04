import dilap.core.tools as dpr
import dilap.core.context as dgc
import dilap.generate.landscape as dls
import dilap.generate.city as dcy

import dilap.generate.lot as dlt
import dilap.primitive.cube as dcu
import dilap.primitive.road as dr

import dp_vector as dpv
import dp_quaternion as dpq

# a continent should be a top level context for dilap
# it creates a full world (island), with a true boundary (ocean)
class continent(dgc.context):

    def __init__(self,*args,**kwargs):
        dgc.context.__init__(self,*args,**kwargs)
        self._def('boundary',dpr.point_ring(250,6),**kwargs)
        self._def('sealevel',-0.5,**kwargs)
        self.define()

    def define(self):
        cityseed = dpv.vector(0,0,50)
        self.cityseeds = [cityseed]

    def generate(self,worn = 0):
        cities = []
        tpts = []
        hpts = []
        rpts = []
        for cd in self.cityseeds:cities.append(dcy.city().generate(cd,worn))
        for cy in cities:
            self._consume(cy)
            tpts.extend(cy._terrain_points())
            hpts.extend(cy._hole_points())
            rpts.extend(cy._region_points())

        lscape = dls.landscape(controls = tpts,holes = hpts,regions = rpts)
        self._consume(lscape.generate(worn))

        # add water models to scenegraph
        water = dcu.cube().translate_z(-0.5)
        water.scale_x(2000).scale_y(2000).scale_z(20)
        water.translate_z(self.sealevel)
        wnode = self._node_wrap(water)
        self._nodes_to_graph(wnode)
        return self


