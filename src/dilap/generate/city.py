import dilap.core.context as dgc
import dilap.core.tools as dpr
import dilap.generate.landscape as dls

import dilap.generate.lot as dlt
import dilap.primitive.road as dr

import dp_vector as dpv
import dp_quaternion as dpq

# a continent should be a top level context for dilap
# it creates a full world (island), with a true boundary (ocean)
class city(dgc.context):

    def _terrain_points(self):
        return self.tpts

    def _hole_points(self):
        return self.hpts

    def _region_points(self):
        return self.rpts

    def generate(self,seed,worn = 0):
        self.tpts = [seed]
        self.hpts = []
        self.rpts = dpr.point_ring(100,6)

        start = dpv.vector(-100,-300, 20)
        end   = dpv.vector( 100, 300, 40)
        tip  = dpv.vector(0,1,0)
        tail = dpv.vector(1,1,0)
        cs = [dpv.vector(-100,-100, 30),dpv.vector( 100, 100, 40)]
        rd = dr.road(start,end,tip,tail,controls = cs)
        self.tpts.extend(rd._terrain_points())
        self._nodes_to_graph(self._node_wrap(rd))

        return self

