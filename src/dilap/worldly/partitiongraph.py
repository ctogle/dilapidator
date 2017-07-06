import dilap.core.base as db

from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat
import dilap.geometry.polymath as pym

import dilap.worldly.topography as dtp

import dilap.worldly.blockletters as dbl

import dilap.core.plotting as dtl
import matplotlib.pyplot as plt

import pdb

def test():
    sq = vec3(0,0,0).sq(256,256)
    peg = vec3(128,0,0).pring(64,8)
    p = partition(loop = sq,holes = [],style = 'A')
    l,r = p.split(p.root,peg)
    p.plot()
    plt.show()

    sq = vec3(0,0,0).sq(256,256)
    peg = vec3(0,0,0).pring(64,8)
    p = partition(loop = sq,holes = [],style = 'A')
    l,r = p.split(p.root,peg)
    p.plot()
    plt.show()

    sq = vec3(0,0,0).sq(256,256)
    peg = vec3(0,0,0).pring(300,8)
    p = partition(loop = sq,holes = [],style = 'A')
    l,r = p.split(p.root,peg)
    p.plot()
    plt.show()

    sq = vec3(0,0,0).sq(256,256)
    peg = vec3(0,0,0).pring(168,8)
    p = partition(loop = sq,holes = [],style = 'A')
    l,r = p.split(p.root,peg)
    p.plot()
    plt.show()

    sq = vec3(0,0,0).sq(256,256)
    peg = vec3(128,0,0).sq(256,256)
    p = partition(loop = sq,holes = [],style = 'A')
    l,r = p.split(p.root,peg)
    p.plot()
    plt.show()

    sq = vec3(0,0,0).sq(256,256)
    peg = vec3(256,0,0).sq(256,256)
    p = partition(loop = sq,holes = [],style = 'A')
    l,r = p.split(p.root,peg)
    p.plot()
    plt.show()

class partition(dtp.topography):

    def plot(self,ax = None):
        if ax is None:ax = dtl.plot_axes(300)
        for v,depth in self.enum(False):
            if self.below(v):col = 'g'
            else:col = 'b'
            vloop = [p.ztrn(depth*10) for p in pym.contract(v.loop,2)]
            #vloop = [p.cp().ztrn(depth*10) for p in v.loop]
            vanchor = vloop[0]
            if depth > 0:
                panchor = self.above(v).loop[0].cp().ztrn((depth-1)*10)
            else:panchor = vec3(0,0,0)
            ax = dtl.plot_polygon(vloop,ax,lw = 3,col = col)
            ax = dtl.plot_edges((vanchor,panchor),ax,col = col,lw = 5)
            for hole in v.holes:
                hole = [p.ztrn(depth*10) for p in pym.contract(hole,-2)]
                #hole = [p.cp().ztrn(depth*10) for p in hole]
                ax = dtl.plot_polygon(hole,ax,lw = 2,col = 'r')
        return ax

    def split(self,sv,nvloop,svloop = None,**kws):
        '''add two child loops of sv'''
        print('former sv method... NEEDS CODE FROM DILAP.TOPOLOGY.PARTITIONGRAPH')
        print('par.loop becomes ebdxy of par.loop and nv.loop')
        nsvs = []
        nv1holes,nv2holes = sv.holes[:],[]
        if svloop is None:
            svloop = sv.loop[:]
            if pym.binbxy(nvloop,svloop):
                nv1holes.append(nvloop)
                nsvs.append(dtp.topography.avert(
                    self,sv,loop = svloop,holes = nv1holes))
            elif pym.binbxy(svloop,nvloop):
                nv2holes.append(svloop)
                nsvs.append(dtp.topography.avert(
                    self,sv,loop = svloop,holes = nv1holes))
            elif pym.bintbxy(nvloop,svloop,ie = False):
                svloop = pym.ebdxy(svloop,nvloop)
                if len(svloop) > 1:
                    for svl in svloop[:-1]:
                        nsvs.append(dtp.topography.avert(
                            self,sv,loop = svl,holes = nv1holes))
                svloop = svloop[-1]
                nsvs.append(dtp.topography.avert(
                    self,sv,loop = svloop,holes = nv1holes))
        nv = dtp.topography.avert(self,sv,loop = nvloop,holes = nv2holes)
        for k in kws:
            for nsv in nsvs:
                nsv.__setattr__(k,sv.__getattribute__(k))
            nv.__setattr__(k,kws[k])

        #self.plot()
        #plt.show()

        return nsvs,nv
