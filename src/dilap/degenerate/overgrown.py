import dilap.core.dilapidor as dd
import dilap.core.model as dmo
import dilap.primitive.vine as dv

import dp_vector as dpv
import dp_ray as dr

import pdb

class ivy(dd.dilapidor):

    def __init__(self,*args,**kwargs):
        dd.dilapidor.__init__(self,*args,**kwargs)
        self._def('z_max',10,**kwargs)
        self.withers.append('ivy')

    # desire a model in world space of ALL nodes/children of the context
    def ivy(self,model,years):
        growth = dmo.model()

        mfaces = model._face_positions(model.faces)
        zray = dr.ray(dpv.zhat.copy().translate_z(20),dpv.zhat.copy().flip())
        rayfaces = dr.intersect_filter(zray,mfaces)
        if rayfaces:print('rayfaces!',len(rayfaces))
        ivy = dv.vine().grow(years).translate(zray.origin)

        growth._consume(ivy)
        return growth


