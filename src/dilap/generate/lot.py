import dilap.core.context as dgc
import dilap.generate.house as dlh

import dp_vector as dpv
import dp_quaternion as dpq

class lot(dgc.context):

    def __init__(self,l,w,*args,**kwargs):
        dgc.context.__init__(self,*args,**kwargs)
        self.l = l
        self.w = w

    def _transform(self,t,q,s):
        dgc.context._transform(self,t,q,s)
        dpv.scale_coords(self.terrain_points,s)
        dpv.rotate_coords(self.terrain_points,q)
        dpv.translate_coords(self.terrain_points,t)

    def _terrain_points(self):
        tpts = []
        [tpts.extend(s._terrain_points()) for s in self.structures]
        self.terrain_points = tpts

    def generate(self,worn = 0):
        house = dlh.house(self.l,self.w).generate(worn)
        self.structures = [house]
        for s in self.structures:self._consume(s)
        self._terrain_points()
        return self

   
