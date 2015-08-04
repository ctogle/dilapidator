import dilap.core.context as dgc
import dilap.generate.landscape as dls
import dilap.generate.lot as dlt
import dilap.primitive.road as dr

import dp_vector as dpv
import dp_quaternion as dpq

class street(dgc.context):

    def generate(self,worn = 0):
        start = dpv.vector(-100,-300, 20)
        end   = dpv.vector( 100, 300,-10)
        tip  = dpv.vector(0,1,0)
        tail = dpv.vector(1,1,0)
        cs = [dpv.vector(-100,-100, 10),dpv.vector( 100, 100,-10)]
        rd = dr.road(start,end,tip,tail,controls = cs)
        self._nodes_to_graph(self._node_wrap(rd))

        #bbs = []
        #lotspace = rd._lotspace(bbs)
        #dlot = dlt.lot(lotspace[0],lotspace[1]).generate(worn)
        #lsppos,lsprot = lotspace[2],lotspace[3]
        #dlot._transform(lsppos,lsprot,dpv.one())
        #self._consume(dlot)
        #lotspace = rd._lotspace(bbs)
        #dlot = dlt.lot(lotspace[0],lotspace[1]).generate(worn)
        #lsppos,lsprot = lotspace[2],lotspace[3]
        #dlot._transform(lsppos,lsprot,dpv.one())
        #self._consume(dlot)

        tpts = []
        #tpts.extend(dlot.terrain_points)
        tpts.extend(rd._terrain_points())

        lscape = dls.landscape(controls = tpts)
        self._consume(lscape.generate(worn))

   

