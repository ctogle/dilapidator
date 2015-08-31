import dilap.core.base as db
import dilap.core.tools as dpr
import dilap.core.vector as dpv

import dilap.mesh.tools as dtl

import matplotlib.pyplot as plt
import pdb



class edge(db.base):
    
    def plot(self,ax = None):
        if ax is None:ax = dtl.plot_axes()
        np1 = self.one.p
        np2 = self.two.p
        dtl.plot_edges([np1,np2],ax)
        return ax

    def directions(self):
        self.tangent = dpv.v1_v2(self.one.p,self.two.p)
        ndir = dpr.angle_from_xaxis_xy(self.tangent)
        ndirf = ndir+dpr.PI if ndir < dpr.PI else ndir - dpr.PI
        return ndir,ndirf

    def key(self):return (self.one.key(),self.two.key())
    def __init__(self,node1,node2,**kwargs):
        self._def('index',None,**kwargs)
        self.one = node1
        self.two = node2

            

      

