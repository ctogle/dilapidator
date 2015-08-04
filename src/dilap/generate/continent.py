import dilap.core.tools as dpr
import dilap.core.context as dgc
import dilap.core.mesh.tools as dtl
import dilap.generate.infrastructure as pif
import dilap.generate.infrastructure.graphregion as grg
import dilap.generate.infrastructure.infragraph as ifg
import dilap.generate.landscape as dls
import dilap.generate.city as dcy

import dilap.generate.lot as dlt
import dilap.primitive.cylinder as dcyl
import dilap.primitive.road as dr

import dp_vector as dpv
import dp_quaternion as dpq

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
        g = ifg.graph()

        contbnd = dpr.point_ring(250,8)
        gregion = grg.region(contbnd,sealevel = self.sealevel)

        nbhdbnd = dpr.point_ring(150,6)
        nbhdregion = grg.neighborhood(nbhdbnd)
        gregion._embed(nbhdregion)

        gregion._graph(g)
        g._update()

        ax = g.plot()
        #ax = gregion.plot()
        plt.show()

        print('im off the rails!!!')
        quit()



        cityseed = dpv.vector(0,0,150)
        self.cityseeds = [cityseed]
        self.cityseeds.append(dpv.vector(250,250,100))

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
        wl = lscape.landbb.x.y - lscape.landbb.x.x + 100.0
        ww = lscape.landbb.y.y - lscape.landbb.y.x + 100.0
        wr = max([wl,ww])
        water = dcyl.cylinder(n = 8).translate_z(-0.5)
        water.scale_x(wr).scale_y(wr).scale_z(20)
        water.translate_z(self.sealevel)
        water.translate(lscape.landbb._center().xy())
        wnode = self._node_wrap(water)
        self._nodes_to_graph(wnode)
        return self


