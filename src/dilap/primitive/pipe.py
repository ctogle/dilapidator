import dilap.core.model as dmo
import dilap.core.tools as dpr

import dilap.primitive.cube as dcu

import dp_vector as dpv
import dp_quaternion as dq

import pdb,random,numpy

class pipe(dmo.model):

    def __init__(self,*args,**kwargs):
        dmo.model.__init__(self,*args,**kwargs)
        self._loop(**kwargs)
        self._curve(**kwargs)
        self._geo()
    
    def _loop(self,**kwargs):
        if 'loop' in kwargs.keys():self.loop = kwargs['loop']
        else:
            self.loop = dpr.point_ring(1,16)
            self.loop.append(self.loop[0].copy())

    def _curve(self,**kwargs):
        if 'curve' in kwargs.keys():self.curve = kwargs['curve']
        else:
            b = dpv.zero()
            bottom = [b,b.copy().translate_z(10)]
            top = [b.copy().translate_z(11).translate_x(2) for b in bottom]
            top.append(top[-1].copy().translate_y(10))
            self.curve = bottom + top

    # extrude a loop along self.curve
    def _geo(self):
        self._extrude(self.loop,self.curve)


