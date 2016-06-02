import dilap.core.base as db
from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat
import dilap.geometry.polymath as pym

import dilap.topology.tree as dtr

import dilap.core.plotting as dtl
import matplotlib.pyplot as plt

import math,numpy,random,pdb



###############################################################################

class tmesh(db.base):

    # tmesh represents a tree of nested loops of which have 
    #   equal z position in the terrain

    # add a new loop to the tree
    def al(self,loop,par):
        nv = self.looptree.avert(par)
        nv.loop = loop
        return nv

    # locate a parent with the tree and add a loop
    def ll(self,loop):
        print('IMPLEMENT LOCATE LOOP!?')
        par = self.root
        return self.al(loop,par)

    ###########################################################################

    # determine the one or two loops which locate this point topologically
    # find the vertex such that p is in the loop but non of the children
    def loc(self,p):
        unfn = [self.root]
        while unfn:
            v = unfn.pop(0)
            if p.onbxy(v.loop):
                return v,[]
            elif p.inbxy(v.loop):
                chs = self.looptree.below(v)
                unfn.extend(chs)
                lst = v,chs
        return lst

    # given a parent loop, a zoffset, and a radial offset
    # create a new loop which is a child of the parent 
    # vertex by contracting and lifting the parent loop
    def grow(self,par,z,r):
        inner = pym.contract(par.loop,r)
        for p in inner:p.ztrn(z)
        iv = self.al(inner,par)
        return iv

    ###########################################################################

    def plot(self,ax = None):
        if ax is None:ax = dtl.plot_axes(100)
        def pv(v,ax):
            ax = dtl.plot_polygon(v.loop,ax,lw = 2)
            for c in self.looptree.below(v):pv(c,ax)
        pv(self.root,ax)
        return ax
        
    def plotxy(self,ax = None):
        if ax is None:ax = dtl.plot_axes_xy(100)
        def pv(v,ax):
            ax = dtl.plot_polygon_xy(v.loop,ax,lw = 2)
            for c in self.looptree.below(v):pv(c,ax)
        pv(self.root,ax)
        return ax

    ###########################################################################

    def __init__(self):pass

    ###########################################################################

    def __call__(self,x,y):
        tp = vec3(x,y,0)
        v,chs = self.loc(tp)
        ccnt = len(chs)
        if ccnt == 0:z = v.loop[0].z
        elif ccnt == 1:
            ov = chs[0]
            olx,nearod = pym.bnearpxy(ov.loop,tp)
            lx ,nearvd = pym.bnearpxy( v.loop,tp)
            totd = nearvd+nearod
            z = (nearod*v.loop[lx].z+nearvd*ov.loop[olx].z)/totd
        else:raise ValueError
        return z

###############################################################################

def checkseq(fp,h,seq,show = False):
    print('check-tseq:',seq)
    t = tmesh()
    t.looptree = dtr.tree()
    t.root = t.al(fp,None)

    tip = t.root
    tip = t.grow(tip,1,2)
    tip = t.grow(tip,1,8)
    tip = t.grow(tip,1,2)
    tip = t.grow(tip,2,4)

    if show:
        t.plot()
        t.plotxy()
        plt.show()
    return t

###############################################################################



