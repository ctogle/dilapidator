import dilap.core.base as db
import dilap.core.context as cx
import dilap.modeling.model as dmo

from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat

import dilap.worldly.treeskin as ltr
import dilap.worldly.building as blg
import dilap.worldly.blgsequencing as bseq

import dilap.core.plotting as dtl
import matplotlib.pyplot as plt

import random,pdb

###############################################################################
###
###############################################################################

class regiongraph(db.base):

    def plotxy(self,fp = None,ax = None):
        if ax is None:
            ax = dtl.plot_axes_xy(vec3(1,0,0).prjps(self.boundary)[1])

        for rmv in self.vs:
            if rmv is None:continue

            '''#
            ax = dtl.plot_polygon_xy(rmv[2]['bound'],ax,lw = 4)
            rc = vec3(0,0,0).com(rmv[2]['bound'])
            rs = str(rmv[0])+','+str(rmv[1])
            ax = dtl.plot_point_xy(rc,dtl.plot_point_xy_annotate(rc,ax,rs))
            if rmv[1]:
                for re in rmv[1]:
                    if self.vs[re] is None:continue
                    rt = (rc,vec3(0,0,0).com(self.vs[re][2]['bound']))
                    ax = dtl.plot_edges_xy(rt,ax,col = 'g')
            for exit in rmv[2]['exits']:
                ax = dtl.plot_point_xy(exit,dtl.plot_point_xy_annotate(exit,ax,'exit'))
            '''#

        ax = dtl.plot_polygon_xy(self.boundary,ax,lw = 2,col = 'r')
        return ax

    def av(self,os,kws):
        vx = self.vcnt
        self.vcnt += 1
        self.vs.append([vx,os,kws])
        return vx

    def __init__(self,boundary,*ags,**kws):
        self.vs = []
        self.vcnt = 0
        self.boundary = boundary

        self.prepare()

    # begin with self.boundary - the shoreline of a continent
    # establish a port region - all roads can trace back to this region
    # each segment of road specified potentially splits the boundary
    #   into additional subregions with distinct characteristics
    # regions grow iteratively
    #   region growth means additional roads/buildings/decorative meshes
    # when roads/buildings/etc are generated, 
    #   cover remaining void with terrain meshes
    def prepare(self):
        self.av([],{'b':[p.cp() for p in self.boundary]})

###############################################################################







###############################################################################
###############################################################################

def ushape(l,w,z):
    fp = [vec3(0,0,z),
        vec3(l/6,0,z),vec3(l/6,-w/3,z),vec3(l/3,-w/2,z),vec3(l/2,-w/2,z),
        #vec3(l/2,2*w/3,z),vec3(-l/2,2*w/3,z),vec3(-l/2,-w/2,z),vec3(-l/3,-w/2,z),
        vec3(l/2,w/2,z),vec3(-l/2,w/2,z),vec3(-l/2,-w/2,z),vec3(-l/3,-w/2,z),
        vec3(-l/6,-w/3,z),vec3(-l/6,0,z)]
    stack = [
        (0,(vec3(-5*l/6,w/6,z),vec3(5*l/6,w/6,z))),
        (1,(vec3(-l/4,-5*w/3,z),vec3(-l/4,5*w/3,z))),
        (2,(vec3(l/4,-5*w/3,z),vec3(l/4,5*w/3,z))),
        (1,(vec3(-5*l/6,0,z),vec3(5*l/6,0,z))),
        (3,(vec3(-5*l/6,0,z),vec3(5*l/6,0,z))),
        (0,(vec3(-l/24,-5*w/3,z),vec3(-l/24,5*w/3,z))),
        (6,(vec3(l/24,-5*w/3,z),vec3(l/24,5*w/3,z))),
        (0,(vec3(-l/4,-5*w/3,z),vec3(-l/4,5*w/3,z))),
        (7,(vec3(l/4,-5*w/3,z),vec3(l/4,5*w/3,z))),
            ]
    exs = [vec3(0,0,0),vec3(-l/2,w/6,0)]
    rtstack = [(1,4),(3,5),(0,8),(8,6),(6,7),(7,9)]
    #atstack = [(1,4),(3,5),(0,8),(8,6),(6,7),(7,9)]
    atstack = []
    shafts = [6]
    return fp,stack,exs,rtstack,atstack,shafts

###############################################################################
###############################################################################

class world(cx.context):

    def __init__(self,*args,**kwargs):
        self._def('name','buildingcontext',**kwargs)
        self._def('bfa',blg.blgfactory(),**kwargs)
        cx.context.__init__(self,*args,**kwargs)

    def populate_region(self,rg,rx,worn = 0):

        print('populate region',rg,rx)
        # several open ended aspects of the problem
        #   how to split arbitrary boundary shapes into rectangular plots

        p,q,s = vec3(0,0,0),quat(1,0,0,0),vec3(1,1,1)
        for m in range(3):

            lvls = 4
            l = 120.0
            w = 60.0
            z = 0

            print('lzod',l,w,z)

            fp = vec3(0,0,0).sq(l,w)
            seq = bseq.simplebuilding(lvls)
            fh = 7

            bkws = {
                'footprint':fp,'sequence':seq,
                'floorheight':fh,
                    }
            blgcx = self.bfa.new(p,q,s,**bkws)
            blgcx.generate(worn)
            self.achild(blgcx)

            p,q,s = p.cp().trn(vec3(l+10 if m % 2 == 0 else 0,w+10 if m % 2 == 1 else 0,0)),q.cp(),s.cp()

    def generate(self,worn = 0):
        bound = vec3(0,0,0).pring(250,8)

        rg = regiongraph(bound)
        self.populate_region(rg,0,worn)

        #rg.plotxy()
        #plt.show()

        return self





