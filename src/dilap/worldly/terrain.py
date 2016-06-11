import dilap.core.base as db
from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat
import dilap.geometry.polymath as pym

import dilap.topology.tree as dtr

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
    def grow(self,par,z,r):
        base = [p.cp() for p in par.loop]

        '''#
        if random.random() > 1.9:
            s1,s2 = vec3(0,-1000,0),vec3(-100,1000,0)
            base = pym.bsegsxy(par.loop,s1,s2)

            print('rollll..',len(base))

            base = base[-1]

            def pl(bs):
                ax = dtl.plot_axes_xy(10)
                ax = dtl.plot_polygon_xy(par.loop,ax,col = 'g',lw = 2)
                ax = dtl.plot_edges_xy((s1,s2),ax,col = 'r',lw = 1)
                for bpy in bs:
                    bpy = pym.contract(bpy,0.1)
                    ax = dtl.plot_polygon_xy(bpy,ax,col = 'b',lw = 4)
                plt.show()
            pl([base])

        else:base = par.loop
        '''#

        lin = base[:]

        #inner = lin[:]
        inner = pym.contract(lin,r,0.1)
        inner = pym.smoothxy(inner,0.1)
        inner = pym.aggregate(inner,10)
        #inner = pym.cleanbxy(inner,2.0)

        if len(inner) > 3:lout = inner
        else:lout = lin

        for p in lout:p.ztrn(z)
        iv = self.al(lout,par)
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

# create a polygon which intersects b and return the difference
def jagged(b):
    xpj = vec3(1,0,0).prjps(b)
    xsc = xpj[1]-xpj[0]

    t = random.uniform(0,2.0*numpy.pi)
    r = random.uniform(xsc/8.0-100,xsc/8.0+100)
    n = random.randint(3,8)

    stamp = vec3(xsc/2.0*math.cos(t),xsc/2.0*math.sin(t),0).pring(r,n)
    b = pym.ebdxy(b,stamp)[0]
    b = pym.aggregate(b,10)

    print('ISVALID',pym.bvalidxy(b))
    ax = dtl.plot_axes_xy(700)
    ax = dtl.plot_polygon_xy(b,ax,lw = 4,col = 'b')
    ax = dtl.plot_polygon_xy(stamp,ax,lw = 2,col = 'r')
    plt.show()
    if not pym.bvalidxy(b) > 0:
        raise ValueError
    
    return b

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

def continent(b):
    t = tmesh()
    t.root = t.al(b,None)

    xpj = vec3(1,0,0).prjps(b)
    xsc = xpj[1]-xpj[0]

    lm = pym.contract(b,50.0)

    for x in range(8):lm = jagged(lm)

    #lm = pym.bisectb(lm)
    lm = pym.smoothxy(lm,0.5)
    lm = pym.aggregate(lm,10)

    lms = [lm]



    r = 5
    for l in lms:
        tip = t.al(l,t.root)

        #tip = t.root
        tip = t.grow(tip,5,r)
        r += 2
        tip = t.grow(tip,5,r)
        tip = t.grow(tip,5,r)
        r += 2
        r += 2
        tip = t.grow(tip,10,r)
        tip = t.grow(tip,15,r)
        tip = t.grow(tip,10,r)
        r += 2
        r += 2
        r += 2
        tip = t.grow(tip,5,r)
        r += 2
        r += 2
        r += 2
        r += 2
        tip = t.grow(tip,5,r)

    ax = t.plot()
    plt.show()

    return t

###############################################################################

def checkseq(fp,h,seq,show = False):
    print('check-tseq:',seq)
    t = tmesh()
    #t.looptree = dtr.tree()
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





