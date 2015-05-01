import dilap.core.model as dmo
import dilap.primitive.tools as dpr

import dp_vector as dpv

import pdb

class cube(dmo.model):

    def __init__(self,*args,**kwargs):
        dmo.model.__init__(self,*args,**kwargs)
        self._geo()

    def _geo(self):
        #bottom = dpr.point_ring(0.5,4)
        og = dpv.zero()
        bottom = [og.copy(),
            og.copy().translate_x(1),og.copy().translate_y(1),
            og.copy().translate_x(1).translate_y(1)]
        top = [b.copy().translate_z(1) for b in bottom]

        pdb.set_trace()

        print 'cube geometry!'


