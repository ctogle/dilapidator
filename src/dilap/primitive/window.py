import dilap.core.model as dmo
import dilap.core.tools as dpr
import dilap.primitive.cube as dcu

import dp_vector as dpv
import dp_quaternion as dpq

import pdb

class window(dmo.model):

    def __init__(self,*args,**kwargs):
        dmo.model.__init__(self,*args,**kwargs)
        self._def('wall',None,**kwargs)
        self._def('z',0.25,**kwargs)
        self._def('w',1.5,**kwargs)
        self._def('h',2.0,**kwargs)
        self._geo()

    def _geo(self):
        w,ww,bh = self.w,self.wall.w,self.wall.h-self.h-self.z
        bottom = dcu.cube().translate_z(0.5)
        bottom.scale_x(w).scale_y(ww).scale_z(bh)
        bottom.rotate(dpq.q_from_uu(dpv.xhat,self.wall.tangent))
        bottom._project_uv_flat()
        th = self.wall.h-self.h-bh
        top = dcu.cube().translate_z(0.5)
        top.scale_x(w).scale_y(ww).scale_z(th)
        top.translate_z(bh+self.h)
        top.rotate(dpq.q_from_uu(dpv.xhat,self.wall.tangent))
        top._project_uv_flat()
        self._consume(bottom)._consume(top)


