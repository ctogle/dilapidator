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
        nv = self.looptree.avert(par)
        nv.loop = loop
        return nv

    ###########################################################################

    # determine the one or two loops which locate this point topologically
    # find the vertex such that p is in the loop but non of the children
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

    # given a parent loop, a zoffset, and a radial offset
    # create a new loop which is a child of the parent 
    # vertex by contracting and lifting the parent loop
    def grow(self,par,z,r,epsilon):
        base = [p.cp() for p in par.loop]

        lin = base[:]

        #inner = lin[:]
        inner = pym.contract(lin,r,epsilon)
        inner = pym.smoothxy(inner,0.5,epsilon)
        #inner = pym.aggregate(inner,2)
        inner = pym.aggregate(inner,r)
        #inner = pym.cleanbxy(inner,2.0)

        bv = pym.bvalidxy(inner)
        if bv == -1:inner.reverse()

        if len(inner) > 3:lout = inner
        else:return None

        for p in lout:p.ztrn(z)
        iv = self.al(lout,par)
        return iv

    ###########################################################################

    def plot(self,ax = None):
        if ax is None:ax = dtl.plot_axes(700)
        def pv(v,ax):
            ax = dtl.plot_polygon(v.loop,ax,lw = 2)
            for c in self.looptree.below(v):pv(c,ax)
        pv(self.root,ax)
        return ax
        
    def plotxy(self,ax = None):
        if ax is None:ax = dtl.plot_axes_xy(700)
        def pv(v,ax):
            ax = dtl.plot_polygon_xy(v.loop,ax,lw = 2)
            for c in self.looptree.below(v):pv(c,ax)
        pv(self.root,ax)
        return ax

    ###########################################################################

    def __init__(self):
        self.looptree = dtr.tree()

    ###########################################################################

    def __call__(self,x,y):
        tp = vec3(x,y,0)
        v,chs = self.loc(tp)
        ccnt = len(chs)
        if ccnt == 0:z = v.loop[0].z
        elif ccnt == 1:
            ov = chs[0]
            olx1,nearod1 = pym.bnearpxy(ov.loop,tp)
            olx2,nearod2 = pym.bnearpxy(
                [ov.loop[x] for x in range(len(ov.loop)) if not x == olx1],tp)
            lx1 ,nearvd1 = pym.bnearpxy( v.loop,tp)
            lx2 ,nearvd2 = pym.bnearpxy(
                [v.loop[x] for x in range(len(v.loop)) if not x == lx1],tp)
            olx = olx1
            lx = lx1
            nearod = nearod1
            nearvd = nearvd1
            totd = nearvd+nearod
            z1 = (nearod*v.loop[lx].z+nearvd*ov.loop[olx].z)/totd
            olx = olx2
            lx = lx2
            nearod = nearod2
            nearvd = nearvd2
            totd = nearvd+nearod
            z2 = (nearod*v.loop[lx].z+nearvd*ov.loop[olx].z)/totd
            z = (z1+z2)/2.0
        else:raise ValueError
        return z

###############################################################################

'''#
def growlandmass(self,b,lm):
    ex = random.randint(0,len(lm)-1)
    p1,p2 = lm[ex-1].cp(),lm[ex].cp()
    etn = p1.tov(p2)
    enm = vec3(0,0,1).crs(etn).flp()
    p4,p3 = p1.cp().trn(enm).pline(p2.cp().trn(enm),2)
    nlm = [p1,p2,p3,p4]
    return nlm
'''#

def continent(b,epsilon = None):
    t = tmesh()
    t.root = t.al(b,None)

    xpj = vec3(1,0,0).prjps(b)
    xsc = xpj[1]-xpj[0]
    if epsilon is None:epsilon = xsc/1000.0
    print('approximated epsilon:',epsilon)

    lm = pym.contract(b,xsc/20.0,epsilon)

    for x in range(10):lm = pyg.ajagged(lm,epsilon)

    lm = pym.bisectb(lm)
    lm = pym.smoothxy(lm,0.5,epsilon)
    lm = pym.smoothxy(lm,0.5,epsilon)
    lm = pym.smoothxy(lm,0.5,epsilon)
    lm = pym.aggregate(lm,2)
    lm = pym.blimithmin(lm,2)
    lms = [lm]

    print('GROWING TERRAIN')

    z = xsc/200.0
    for l in lms:
        tip = t.al(l,t.root)
        tiparea = pym.bareaxy(tip.loop)
        while tiparea > (xsc/20.0)**2:
            #grad = random.uniform(-0.1,0.3)
            #d = random.choice([1,1,-1])
            d = 1.0
            #grad = d*random.uniform(0.2,0.8)
            grad = 0.5
            tip = t.grow(tip,d*z,abs(z/grad),epsilon)

            if tip is None:break
            tiparea = pym.bareaxy(tip.loop)

            print('GROW',grad,z,tiparea,pym.bccw(tip.loop))

    print('GREW TERRAIN')
    ax = t.plot()
    plt.show()

    return t

###############################################################################





