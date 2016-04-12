import dilap.core.model as dmo
import dilap.core.tools as dpr
import dilap.core.vector as dpv

import pdb

class cone(dmo.model):

    def __init__(self,*args,**kwargs):
        dmo.model.__init__(self,*args,**kwargs)
        self._geo(**kwargs)

    def _geo(self,n = 8):
        bottom = dpr.point_ring(1,n)
        dpv.translate_coords_z(bottom,-0.5)
        top = dpv.center_of_mass(bottom).translate_z(1)
        bottom.reverse()
        self._tripie(bottom)
        bottom.append(bottom[0])
        bottom.reverse()
        self._trifan(top,bottom)


