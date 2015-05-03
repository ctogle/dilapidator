import dilap.core.base as db
import dilap.core.sgraph as dsg
import dilap.core.context as dgc
import dilap.io.io as dio
import dilap.primitive.cube as dcu
import dilap.primitive.cone as dco
import dilap.primitive.cylinder as dcyl
import dilap.primitive.wall as dw
import dilap.primitive.floor as df

import dp_vector as dpv

import random

class landscape(dgc.context):

    def generate(self,other,worn = 0):
        dcube = dcu.cube().scale_x(10).scale_y(10).translate_z(-2)
        dnode = self._node_wrap(dcube)
        self._nodes_to_graph(dnode)


