import dilap.core.base as db
import dilap.core.tools as dpr
import dilap.core.vector as dpv

import dilap.mesh.tools as dtl

import matplotlib.pyplot as plt
import pdb



class node(db.base):

    def plot(self,ax = None):
        if ax is None:ax = dtl.plot_axes()
        plotp = self.p.copy()
        dtl.plot_point(plotp,ax)
        return ax

    # d is the direction to o.p from self.p
    def connect(self,d,o):
        self.ring[o.index] = d

    # o is another node currently connected to self
    def disconnect(self,o):
        del self.ring[o.index]

    def key(self):return (self.p.x,self.p.y,self.layer)
    def __init__(self,p,**kwargs):
        self._def('index',None,**kwargs)
        self._def('layer',0,**kwargs)
        self.ring = {}
        self.p = p






