import dilap.core.lsystem as lsy
from dilap.geometry import *
import dilap.geometry.tools as gtl
import dilap.geometry.polymath as pym
from dilap.topology import *
import dilap.topology.planargraph as pgr
import dilap.worldly.polygen as pyg
import dilap.core.plotting as dtl
import matplotlib.pyplot as plt
import random
import numpy
import math
import pdb

# tmesh represents a tree of nested loops of which have 
#   equal z position in the terrain
class topography(tree):
    '''
    a topography is a tree of loops with the ability 
    of locating an x,y location within the loop hierarchy
    '''

    # determine the one or two loops which locate this point topologically
    # find the vertex such that p is in the loop but none of its childrens'
    def loc(self,p):
        lst = None
        unfn = [self.root]
        while unfn:
            v = unfn.pop(0)
            if p.onbxy(v.loop):return v,[]
            elif p.inbxy(v.loop):
                chs = self.below(v)
                unfn.extend(chs)
                lst = v,chs
        if lst is None:lst = self.root,[]
        return lst

    # determine the one or two loops which locate this loop topologically
    # find the vertex such that l is in the loop but none of the children
    def locloop(self,l,override = 0,mingrad = 1.0):

        # modify the loop of v (and its children) to respect new loop l
        def correct_loops(v,l):
            lswollen = [p.cp().ztrn(v.loop[0].z-p.z) for p in l]
            #lswollen = pym.contract(lswollen,abs(l[0].z-v.loop[0].z)/mingrad)

            #ax = dtl.plot_axes(300)
            #ax = dtl.plot_polygon(lswollen,ax,lw = 3,col = 'b')
            #ax = dtl.plot_polygon(v.loop,ax,lw = 3,col = 'r')
            #plt.show()

            newloops = pym.ebdxy(v.loop,lswollen)
            v.loop = newloops[0]
            if len(newloops) > 1:
                for nl in newloops[1:]:
                    nv = self.avert(self.above(v),loop = nl)
            for ch in self.below(v):
                if pym.bintbxy(lswollen,ch.loop):
                    correct_loops(ch,lswollen)

        lst = None
        unfn = [self.root]
        while unfn:
            v = unfn.pop(0)
            if pym.binbxy(l,v.loop):
                chs = self.below(v)
                unfn.extend(chs)
                lst = v,chs
            elif pym.bintbxy(l,v.loop):
                # modify the existing loop to accomodate the new loop
                if override == 1:
                    correct_loops(v,l)
                    chs = self.below(v)
                    unfn.extend(chs)
                    lst = self.above(v),chs
                # modify the new loop to accomodate the existing loop
                elif override == 2:
                    raise NotImplemented
                    return v,[]
                # modify neither loop (likely causes poor results)
                else:return v,[]
        if lst is None:lst = self.root,[]
        return lst

    # given a parent loop, a zoffset, and a radial offset
    # create a new loop which is a child of the parent 
    # vertex by contracting and lifting the parent loop
    def grow(self,par,z,r,epsilon,base = None):
        if base is None:base = [p.cp() for p in par.loop]
        i = base[:]
        i = pym.contract(i,r,epsilon)
        if i is None:return None
        #i = pym.smoothxy(i,0.5,epsilon)
        #i = pym.smoothxy(i,0.5,epsilon)
        #i = pym.smoothxy(i,0.5,epsilon)
        i = pym.aggregate(i,r)
        #i = pym.cleanbxy(i,2.0)

        bv = pym.bvalidxy(i,epsilon)
        if   bv == -1:i.reverse()
        elif bv == -2:
            return None
            raise ValueError
        if len(i) > 3:lout = i
        else:return None

        for p in lout:p.ztrn(z)
        iv = self.avert(par,loop = lout)
        return iv

    def split(self,v,s,z,rad,epsilon):
        sp = vec3(0,0,0).com(v.loop)
        sd = vec3(0,1,0) if random.random() > 0.5 else vec3(1,0,0)
        s1 = sp.cp().trn(sd.cp().uscl( 10000))
        s2 = sp.cp().trn(sd.cp().uscl(-10000))
        loops = pym.bsegsxy(v.loop,s1,s2,0.1)

        lvs = []
        for el in loops:
            if not pym.bccw(el):el.reverse()
            el = pym.blimithmin(el,2,50)
            el = pym.aggregate(el,2)
            #el = pym.smoothxy(el,0.5,epsilon)
            nv = self.grow(v,z,rad,epsilon,el)
            lvs.append(nv)

        #ax = self.plot()
        #for el in loops:
        #    ax = dtl.plot_polygon(el,ax,lw = 2,col = 'r')
        #plt.show()

        return lvs

    ###########################################################################

    def plot(self,ax = None,l = 300,s = (vec3(-10000,0,0),vec3(10000,0,0))):
        if ax is None:ax = dtl.plot_axes(l)
        def pv(v,ax):
            dtl.plot_polygon(v.loop,ax,lw = 2)
            vp = v.loop[0]
            for c in self.below(v):
                cp = c.loop[0]
                dtl.plot_points((vp,cp),ax)
                dtl.plot_edges((vp,cp),ax,lw = 3,col = 'b')
                pv(c,ax)
        pv(self.root,ax)
        return ax
        
    def plotxy(self,ax = None,l =300):
        if ax is None:ax = dtl.plot_axes_xy(l)
        def pv(v,ax):
            ax = dtl.plot_polygon_xy(v.loop,ax,lw = 2)
            #for c in self.looptree.below(v):pv(c,ax)
            for c in self.below(v):pv(c,ax)
        pv(self.root,ax)
        return ax

    ###########################################################################

    def call_many(self,tp,v,chs):
        vex,vtpd = pym.bnearpxy(v.loop,tp)
        z1 = v.loop[vex].z
        if chs:
            chexs,chpds = zip(*[pym.bnearpxy(c.loop,tp) for c in chs])

            if min(chpds) < 0.1 and len(chs) > 1:
                print('chpds < 1 and len(chs) > 1')
                #pdb.set_trace()
                #ax = dtl.plot_axes(200)
                #ax = dtl.plot_point(tp,ax)
                #for c in chs:
                #    ax = dtl.plot_polygon(c.loop,ax,lw = 2,col = 'g')
                #ax = dtl.plot_polygon(v.loop,ax,lw = 2,col = 'r')
                #plt.show()

            chx = chpds.index(min(chpds))
            z2 = chs[chx].loop[chexs[chx]].z

            dmax = max([vtpd]+list(chpds))
            dtot = (dmax-vtpd)+sum([dmax-chpd for chpd in chpds])
            if dtot == 0:dtot = 1
            z = ((dmax-vtpd)*z1+sum([(dmax-chpd)*z2 for chpd in chpds]))/dtot
        else:
            z2 = z1
            z = z1
        return z

