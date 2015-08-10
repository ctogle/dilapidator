import dilap.core.base as db

import dilap.core.tools as dpr
import dilap.core.mesh.tools as dtl
import dilap.core.lsystem as dls

import dilap.generate.infrastructure.tools as itl

import dp_vector as dpv

import matplotlib.pyplot as plt
import random as rm
import pdb



class node(db.base):

    # d is the direction to o.p from self.p
    def connect(self,d,o):
        self.targetring[o.index] = d
        self.ring[o.index] = d

    # o is another node currently connected to self
    def disconnect(self,o):
        del self.targetring[o.index]
        del self.ring[o.index]

    def plot(self,ax = None):
        if ax is None:ax = dtl.plot_axes()
        plotp = self.p.copy()
        dtl.plot_point(plotp,ax)
        for d in self.ring.keys():
            f = plotp.copy().translate(self.spikes[d])
            dtl.plot_edges([plotp,f],ax)
        return ax

    def plot_xy(self,ax = None):
        if ax is None:ax = dtl.plot_axes_xy()
        dtl.plot_point_xy(self.p,ax)
        for d in self.ring.keys():
            f = self.p.copy().translate(self.spikes[d])
            dtl.plot_edges_xy([self.p,f],ax)
        return ax

    def key(self):return (self.p.x,self.p.y,self.layer)
    def __init__(self,p,**kwargs):
        self._def('index',None,**kwargs)
        self._def('layer',0,**kwargs)
        self.p = p
        self.targetring = {}
        self.ring = {}
        self.spikes = {}

    def _spikes_nudge(self,graph,target):
        ws = []
        angles = []
        tkeys = list(self.targetring.keys())
        for tk in tkeys:
            ekey = (self.key(),graph.nodes[tk].key())
            edx = graph.edges_lookup[ekey]
            ws.append(1.0 if graph.edges[edx].interpolated else 0.0)
            angles.append(self.targetring[tk])
        iangles = itl.nudge(angles,ws,target = target)
        for tk,na in zip(tkeys,iangles):
            self.ring[tk] = na
            spike = dpv.xhat.copy().rotate_z(dpr.rad(na)).scale_u(8)
            self.spikes[tk] = spike

    # update the actual positions of spikes based on targetring
    def _spikes(self,graph):
        self._retarget(graph)
        self.ring,self.spikes = {},{}
        tkys = list(self.targetring.keys())
        rcnt = len(tkys)
        if rcnt == 0:return
        else:
            if rcnt == 2:t = 180
            else:t = 90
            self._spikes_nudge(graph,t)

    # verify the directions in target ring and update 
    # tangents of any connecting edges
    def _retarget(self,graph):
        rgxs = list(self.ring.keys())
        for rgx in rgxs:
            rgnd = graph.nodes[rgx]

            rkey = (rgnd.key(),self.key())
            regx = graph.edges_lookup[rkey]
            if regx is None:
                del self.targetring[rgx]
                del self.ring[rgx]
                continue

            reg = graph.edges[regx]
            ndir1,ndir2 = reg._directions()
            if not reg.one is self:ndir1,ndir2 = ndir2,ndir1
            self.targetring[rgx] = ndir1
            graph.nodes[rgx].targetring[self.index] = ndir2










