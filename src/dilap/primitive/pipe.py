import dilap.core.model as dmo
import dilap.core.tools as dpr

import dilap.primitive.cube as dcu

import dp_vector as dpv
import dp_quaternion as dq

import pdb,random,numpy

class pipe(dmo.model):

    def __init__(self,*args,**kwargs):
        dmo.model.__init__(self,*args,**kwargs)
        self._def('m','generic',**kwargs)
        self._loop(**kwargs)
        self._curve(**kwargs)
        self._geo()
    
    def _loop(self,**kwargs):
        if 'loop' in kwargs.keys() and kwargs['loop']:
            self.loop = kwargs['loop']
        else:
            self.loop = dpr.point_ring(1,16)
            #self.loop = dpr.corners(4,4)
            self.loop.append(self.loop[0].copy())

    def _curve(self,**kwargs):
        if 'curve' in kwargs.keys() and kwargs['curve']:
            self.curve = kwargs['curve']
        else:
            self.curve = [dpv.zero(),
                dpv.zero().translate_z(5),
                dpv.zero().translate_z(5).translate_y(5),
                dpv.zero().translate_z(5).translate_y(5).translate(dpv.vector(1,1,2)),
                dpv.zero().translate_z(5).translate_y(5).translate(dpv.vector(1,1,2)).translate_x(4),
                dpv.zero().translate_z(5).translate_y(5).translate(dpv.vector(1,1,2)).translate_x(4).translate_z(-2)]

    # extrude a loop along self.curve
    def _geo(self):
        control = dpv.zero()
        self._extrude(self.loop,self.curve,control,m = self.m)


