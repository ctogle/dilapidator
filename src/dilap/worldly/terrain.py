import dilap.core.base as db
from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat
import dilap.geometry.polymath as pym

import dilap.topology.tree as dtr

import dilap.worldly.polygen as pyg

import dilap.core.plotting as dtl
import matplotlib.pyplot as plt

import math,numpy,random,pdb



###############################################################################

# tmesh represents a tree of nested loops of which have 
#   equal z position in the terrain
class tmesh(db.base):

    # add a new loop to the tree
    def al(self,loop,par):
        # DOES THIS NEED TO SPLIT EDGES AND NOT JUST ADD VERTICES!?!?!
        # DOES THIS NEED TO SPLIT EDGES AND NOT JUST ADD VERTICES!?!?!
        # DOES THIS NEED TO SPLIT EDGES AND NOT JUST ADD VERTICES!?!?!
        nv = self.looptree.avert(par)
        nv.loop = loop
        return nv

    ###########################################################################

    # determine the one or two loops which locate this point topologically
    # find the vertex such that p is in the loop but none of its childrens'
    def loc(self,p):
        lst = None
        unfn = [self.root]
        while unfn:
            v = unfn.pop(0)
            if p.onbxy(v.loop):return v,[]
            elif p.inbxy(v.loop):
                chs = self.looptree.below(v)
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
            lswollen = pym.contract(lswollen,abs(l[0].z-v.loop[0].z)/mingrad)
            newloops = pym.ebdxy(v.loop,lswollen)
            v.loop = newloops[0]
            if len(newloops) > 1:
                for nl in newloops[1:]:
                    nv = self.al(nl,self.looptree.above(v))
            for ch in self.looptree.below(v):
                correct_loops(ch,lswollen)

        lst = None
        unfn = [self.root]
        while unfn:
            v = unfn.pop(0)
            if pym.binbxy(l,v.loop):
                chs = self.looptree.below(v)
                unfn.extend(chs)
                lst = v,chs
            elif pym.bintbxy(l,v.loop):
                # modify the existing loop to accomodate the new loop
                if override == 1:
                    correct_loops(v,l)
                    chs = self.looptree.below(v)
                    unfn.extend(chs)
                    lst = self.looptree.above(v),chs
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
        iv = self.al(lout,par)
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

        ax = self.plot()
        for el in loops:
            ax = dtl.plot_polygon(el,ax,lw = 2,col = 'r')
        plt.show()

        return lvs

    ###########################################################################

    def plot(self,ax = None,l = 300,s = (vec3(-10000,0,0),vec3(10000,0,0))):
        if ax is None:ax = dtl.plot_axes(l)
        def pv(v,ax):
            dtl.plot_polygon(v.loop,ax,lw = 2)
            vp = v.loop[0]
            for c in self.looptree.below(v):
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
            for c in self.looptree.below(v):pv(c,ax)
        pv(self.root,ax)
        return ax

    ###########################################################################

    def __init__(self):
        self.looptree = dtr.tree()

    ###########################################################################

    def call_many(self,tp,v,chs):
        vex,vtpd = pym.bnearpxy(v.loop,tp)
        z1 = v.loop[vex].z
        if chs:
            chexs,chpds = zip(*[pym.bnearpxy(c.loop,tp) for c in chs])

            if min(chpds) < 0.1 and len(chs) > 1:
                #pdb.set_trace()
                ax = dtl.plot_axes(200)
                ax = dtl.plot_point(tp,ax)
                for c in chs:
                    ax = dtl.plot_polygon(c.loop,ax,lw = 2,col = 'g')
                ax = dtl.plot_polygon(v.loop,ax,lw = 2,col = 'r')
                plt.show()

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

    def __call__(self,x,y):
        tp = vec3(x,y,0)
        v,chs = self.loc(tp)
        z = self.call_many(tp,v,chs)
        return z

###############################################################################

# generate landmasses within a boundary polygon
def sow_earth(t,b,e):
    lm = pym.contract(b,e*50.0,e)
    for x in range(10):lm = pyg.ajagged(lm,e)
    lm = pym.bisectb(lm)
    lm = pym.smoothxy(lm,0.5,e)
    lm = pym.smoothxy(lm,0.5,e)
    lm = pym.smoothxy(lm,0.5,e)
    lm = pym.aggregate(lm,2)
    lm = pym.blimithmin(lm,2,50)
    lms = [lm]
    return lms

# generate the topographical decisions of a set of topo loops
def raise_earth(t,tips,e,mingrad = 1.0,mindelz = 5.0):
    newtips = None
    for tip in tips:
        if tip is None:continue

        newoffset = vec3(-30,-30,mindelz)
        newloop = [p.cp().trn(newoffset) for p in tip.loop]

        oldloop = [p.cp().ztrn(mindelz) for p in tip.loop]
        oldloop = pym.aggregate(pym.contract(oldloop,mindelz/mingrad),20)

        uloops = pym.ebixy(oldloop,newloop)
        if len(uloops) == 1:
            newtip = t.al(uloops[0],tip)
        else:
            newtip = None
            print('ebicnt',len(uloops))

            ax = dtl.plot_axes(200)
            ax = dtl.plot_polygon(tip.loop,ax,ls = '--',lw = 2,col = 'k')
            ax = dtl.plot_polygon(oldloop,ax,lw = 3,col = 'b')
            ax = dtl.plot_polygon(newloop,ax,lw = 3,col = 'g')
            for uloop in uloops:
                ax = dtl.plot_polygon(uloop,ax,lw = 3,col = 'r')
            plt.show()

            #pdb.set_trace()

        if newtip is None:continue
        else:newtips = [newtip]

        '''#
        ax = dtl.plot_axes(200)
        ax = dtl.plot_polygon(tip.loop,ax,ls = '--',lw = 2,col = 'k')
        ax = dtl.plot_polygon(oldloop,ax,lw = 3,col = 'b')
        ax = dtl.plot_polygon(newloop,ax,lw = 3,col = 'g')
        ax = dtl.plot_polygon(newtip.loop,ax,lw = 3,col = 'r')
        plt.show()
        '''#

        '''#
        try:tiparea = pym.bareaxy(tip.loop)
        except ValueError:continue
        if tiparea > (e*50.0)**2:
            d = 1.0
            grad = 0.5
            z = e*5.0

            newtip = t.grow(tip,d*z,abs(z/grad),e)

            if newtip is None:continue
            else:newtips = [newtip]

            newtiparea = pym.bareaxy(newtip.loop)
            print('EARTHRAISED',newtiparea,pym.bccw(newtip.loop))

            if newtiparea > (e*40.0)**2:
                print('break topo')
                newtips = t.split(newtip,None,d*z,abs(z/grad),e)
        '''#


    if not newtips is None:raise_earth(t,newtips,e)

    return



    tipsarea = [pym.bareaxy(tip.loop) for tip in tips]
    z = e*5.0
    while max(tipsarea) > (e*50.0)**2:
        for tip in tips:
            d = 1.0
            grad = 0.5

            tip = t.grow(tip,d*z,abs(z/grad),e)
            if tip is None:break
            tiparea = pym.bareaxy(tip.loop)
            print('EARTHRAISED',tiparea,pym.bccw(tip.loop))

            if tiparea > (e*10.0)**2:
                print('break topo')
                tips = t.split(tip,None,d*z,abs(z/grad),e)

        #ax = t.plot()
        #plt.show()

# create a terrain mesh for a continent
def continent(b,landmasses = None,epsilon = None):
    t = tmesh()
    t.root = t.al(b,None)
    if epsilon is None:
        xpj = vec3(1,0,0).prjps(b)
        epsilon = (xpj[1]-xpj[0])/1000.0
    if landmasses is None:
        landmasses = sow_earth(t,b,epsilon)
    for l in landmasses:
        tips = [t.al(l,t.root)]
        raise_earth(t,tips,epsilon)
    print('RAISED TERRAIN WITH EPSILON:',epsilon)
    return t

###############################################################################





