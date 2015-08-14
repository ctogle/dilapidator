import dilap.core.base as db
import dilap.core.vector as dpv
import dilap.core.tools as dpr
import dilap.core.lsystem as dls
import dilap.mesh.tools as dtl
import dilap.infrastructure.graphnode as gnd
import dilap.infrastructure.infralsystem as ifl

import matplotlib.pyplot as plt
import random as rm
import pdb



class edge(db.base):
    
    def plot(self,ax = None):
        if ax is None:ax = dtl.plot_axes()
        np1 = self.one.p
        np2 = self.two.p
        dtl.plot_edges([np1,np2],ax)
        dtl.plot_edges(self.rpts,ax)
        return ax

    def plot_xy(self,ax = None):
        if ax is None:ax = dtl.plot_axes_xy()
        np1 = self.one.p
        np2 = self.two.p
        dtl.plot_edges_xy([np1,np2],ax)
        dtl.plot_edges_xy(self.rpts,ax)
        return ax

    def __init__(self,node1,node2,**kwargs):
        self._def('index',None,**kwargs)
        self._def('interpolated',True,**kwargs)
        self.one = node1
        self.two = node2
        self._def('width',10,**kwargs)

    def _directions(self):
        self.tangent = dpv.v1_v2(self.one.p,self.two.p)
        ndir = dpr.deg(dpv.angle_from_xaxis_xy(self.tangent))
        ndirf = ndir+180 if ndir < 180 else ndir - 180
        return ndir,ndirf

    # return list of layers which self and o have in common
    # otherwise return None
    def _layers_intersect(self,o):
        l1,l2,l3,l4 = self.one.layer,self.two.layer,o.one.layer,o.two.layer
        shared = []
        if not l1 in shared and (l1 == l3 or l1 == l4):shared.append(l1)
        if not l2 in shared and (l2 == l3 or l2 == l4):shared.append(l2)
        return shared

    # return pt of intersection if the tangents of self and o intersect
    # otherwise return None
    def _tangents_intersect(self,o):
        s1 = (self.one.p.copy().xy(),self.two.p.copy().xy())
        s2 = (o.one.p.copy().xy(),o.two.p.copy().xy())
        segisect = dtl.segments_intersect_at(s1,s2)
        return segisect 

    # use vector spline to add road points to plot!!
    def _place_road(self,graph):                       
        rcnt = int(self.tangent.magnitude()/8.0)
        if self.interpolated:
            r1 = self.one.p.copy()
            r2 = self.one.p.copy().translate(self.one.spikes[self.two.index])
            r3 = self.two.p.copy().translate(self.two.spikes[self.one.index])
            r4 = self.two.p.copy()
            rpts = dpv.vector_spline(r1,r2,r3,r4,rcnt)
        else:rpts = dpr.point_line(self.one.p,self.two.p,rcnt)
        self.rpts = rpts
        self.rtangents = []
        self.rnormals = []
        self.lbpts = []
        self.rbpts = []
        for rpx in range(1,len(rpts)):
            rtangent = dpv.v1_v2(rpts[rpx-1],rpts[rpx]).xy().normalize()
            rnormal = rtangent.copy().rotate_z(dpr.rad(90)).normalize()
            self.rtangents.append(rtangent)
            self.rnormals.append(rnormal)
            rwv = rnormal.copy().scale_u(self.width)             
            self.lbpts.append(rpts[rpx-1].copy().translate(rwv))
            self.rbpts.append(rpts[rpx-1].copy().translate(rwv.flip()))
        self.rtangents.append(self.rtangents[-1])
        self.rnormals.append(self.rnormals[-1])
        rwv = self.rnormals[-1].copy().scale_u(self.width)
        self.lbpts.append(rpts[-1].copy().translate(rwv))
        self.rbpts.append(rpts[-1].copy().translate(rwv.flip()))
        self.lbpts.reverse()

    # given a wise (cw or ccw), return a node key from an edge
    # which minimizes the angle between self and that edge
    # where the edges are joined by their shared node
    # sndkey is the key of the node which is to be shared
    def _walk(self,sndkey,wise):
        if   sndkey == self.one.key():share,unshare = self.one,self.two
        elif sndkey == self.two.key():share,unshare = self.two,self.one
        turns = [k for k in share.ring if not k == unshare.index]
        tncnt = len(turns)
        if tncnt == 0:return None
        elif tncnt == 1:return turns[0]
        else:
            print('turnnns',turns)
            pdb.set_trace()
            










