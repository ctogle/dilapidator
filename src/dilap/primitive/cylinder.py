import dilap.core.model as dmo
import dilap.core.tools as dpr

import dp_vector as dpv

import pdb

class cylinder(dmo.model):

    def __init__(self,*args,**kwargs):
        dmo.model.__init__(self,*args,**kwargs)
        self._geo(**kwargs)

    def _geo(self,n = 8):
        bottom = dpr.point_ring(1,n)
        dpv.translate_coords_z(bottom,-0.5)
        top = [b.copy().translate_z(1) for b in bottom]
        bottom.reverse()
        self._tripie(bottom)
        self._tripie(top)
        top.reverse()
        bottom.append(bottom[0])
        top.append(top[0])
        self._bridge(bottom,top)


