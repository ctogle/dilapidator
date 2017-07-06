import dilap.core.lsystem as lsy
from dilap.core.model import model
from dilap.geometry import *
import dilap.geometry.tools as gtl
import dilap.geometry.polymath as pym

import dilap.topology.tree as dtr
import dilap.topology.planargraph as pgr

import dilap.worldly.polygen as pyg
import dilap.worldly.topography as dtp

import dilap.core.plotting as dtl
import matplotlib.pyplot as plt

import math,numpy,random,pdb

class terrain(object):

    @staticmethod
    def sow_earth(b,e):
        '''generate landmasses within a boundary polygon'''
        lm = pym.contract(b,e*50.0,e)
        return [lm]
        for x in range(10):lm = pyg.ajagged(lm,e)
        lm = pym.bisectb(lm)
        lm = pym.smoothxy(lm,0.5,e)
        lm = pym.smoothxy(lm,0.5,e)
        lm = pym.smoothxy(lm,0.5,e)
        lm = pym.aggregate(lm,2)
        lm = pym.blimithmin(lm,2,50)
        lms = [lm]
        return lms

    @staticmethod
    def raise_earth(topo,tips,e,mingrad = 1.0,mindelz = -5.0,depth = 0):
        if depth == 0:
            for j,tip in enumerate(tips):
                tips[j].loop = pym.bisectb(tips[j].loop)
        #if depth == 10:return
        print('... raising earth ... (depth: %d)' % depth)
        newtips = []
        for tip in tips:
            if tip is None:continue
            newloop = [p.cp().ztrn(mindelz) for p in tip.loop]
            newloop = pym.aggregate(pym.contract(newloop,abs(mindelz/mingrad)),5)
            #uloops = [pym.smoothxyi(newloop,0.8,e,10,1)]
            #uloops = [pym.pinchb(u,10) for u in uloops if u]
            #uloops = [pym.pinchb(newloop,10)]
            uloops = pym.pinchb(newloop,10)
            uloops = [u for u in uloops if u]
            for u in uloops:
                if not pym.bccw(u):
                    u.reverse()
            '''
            if len(uloops) == 1:
                print('... 1 loop to raise ...')
                theloop = uloops[0]
                newtiparea = pym.bareaxy(theloop,True)
            elif uloops:
                print('... 2+ loops to raise ...')
                bareas = [pym.bareaxy(u,True) for u in uloops]
                theloop = uloops[bareas.index(max(bareas))]
                newtiparea = max(bareas)
            else:
                print('... no loops to raise ...')
                theloop = []
            '''
            for theloop in uloops:
                newtiparea = pym.bareaxy(theloop,True)
                newtip = None
                if theloop and newtiparea > 10:#2*(e)**2:
                    newtip = topo.avert(tip,loop = theloop)
                else:print('... abortive tip! ...')
                if newtip is None:continue
                else:newtips.append(newtip)
        if newtips:
            terrain.raise_earth(topo,newtips,e,mingrad,mindelz,depth+1)
        else:print('... exhausted tip! ...')

    @staticmethod
    def mountains(topo,tip,e = 2):
        '''topographical loop corresponding to the base of a mountain range'''
        p,d,i,axiom,rules = ((vec3(0,0,0),vec3(0,1,0),6,
            'X',dict([('X','{[[X}{]X}F]X'),('F','FA'),('A','F')])))
        params = dict(dazimuthal = gtl.rad(25.7),drho = 20)

        pg = pgr.planargraph()
        for piece in lsy.lgen(p,d,axiom,rules,i,**params):
            if isinstance(piece,tuple):
                p1,p2 = piece
                v1,v2 = pg.fp(p1,10),pg.fp(p2,10)
                e12 = pg.fe(v1,v2)
            elif isinstance(piece,vec3):pass
        py = pym.pgtopy(pg,5)[0]
        py = pym.smoothxy(py,0.5,2,0)
        py = pym.pinchb(py,5)[0]
  
        pyprjx = vec3(1,0,0).prjps(py)
        pyprjy = vec3(0,1,0).prjps(py)
        tipprjx = vec3(1,0,0).prjps(tip.loop)
        tipprjy = vec3(0,1,0).prjps(tip.loop)
        recenter = vec3(
            ((pyprjx[1]+pyprjx[0])-(tipprjx[1]+tipprjx[0]))/-2.0,
            ((pyprjy[1]+pyprjy[0])-(tipprjy[1]+tipprjy[0]))/-2.0,0)
        for p in py:p.trn(recenter)
        pyprjx = vec3(1,0,0).prjps(py)
        pyprjy = vec3(0,1,0).prjps(py)
        tipprjx = vec3(1,0,0).prjps(tip.loop)
        tipprjy = vec3(0,1,0).prjps(tip.loop)
        scale = 1.*(tipprjx[1]-tipprjx[0]+tipprjy[1]-tipprjy[0])/\
                    ( pyprjx[1]- pyprjx[0]+ pyprjy[1]- pyprjy[0])
        for p in py:p.scl(vec3(scale,scale,0))

        #com = vec3(0,0,0).com(tip.loop).tov(vec3(0,0,0).com(py)).uscl(-1)
        py = pym.ebixy(tip.loop,py)[0]
        #py = pym.smoothxy(py,0.5,2)
        py = pym.pinchb(py,5)[0]

        newtips = [topo.avert(tip,loop = py)]

        #print('mountains')
        #ax = dtl.plot_axes_xy(400)
        #ax = dtl.plot_polygon_xy(tip.loop,ax,lw = 3,col = 'b')
        #ax = dtl.plot_polygon_xy(py,ax,lw = 3,col = 'g')
        #plt.show()

        return newtips

    def __init__(self,b,e,**kws):
        self.boundary = b
        self.e = e
        self.natural = dtp.topography(loop = b)
        self.landmasses = []
        for lm in self.sow_earth(b,e):
            lmv = self.natural.avert(self.natural.root,loop = lm)
            tips = self.mountains(self.natural,lmv,e)
            self.raise_earth(self.natural,tips,e)
            self.landmasses.append(lmv)
        self.layers = {'natural':self.natural}

    def __call__(self,x,y):
        '''properly interpolate among various topographies'''
        '''provide a deterministic mapping (x,y)->z consistent with terrain'''
        tp = vec3(x,y,0)
        v,chs = self.natural.loc(tp)
        z = self.natural.call_many(tp,v,chs)
        return z

    



def pepper(t,tip,e = 2):
    print('pepper')

    c = vec3(0,0,0).com(tip.loop)
    d = vec3(1,0,0).uscl(1000)
    s1,s2 = c.cp().trn(d.flp()),c.cp().trn(d.flp())
    ch = lambda b : pyg.chunk(b,e,random.random()+0.2,edge = random.randint(0,1))
    bs = [ch(tip.loop) for x in range(10)]
    bs = pym.bsuxy(bs,e)

    #ax = dtl.plot_axes_xy(200)
    #for b in bs:
    #    dtl.plot_polygon_xy(b,ax,lw = 2)
    #plt.show()

    newtips = []
    for sb in bs:
        if not pym.bccw(sb):sb.reverse()

        sb = pyg.ajagged(sb,2)
        sb = pyg.ajagged(sb,2)
        sb = pyg.ajagged(sb,2)
        #sb = pyg.ajagged(sb,2)

        #for x in range(3):b = pyg.ajagged(b,e)
        sb = pym.bisectb(sb)
        #sb = pym.smoothxy(sb,0.5,e)
        sb = pym.smoothxy(sb,0.5,e)
        sb = pym.aggregate(sb,4)
        #sb = pym.blimithmin(sb,2,50)

        for p in sb:p.ztrn(2)
        newtips.append(t.al(sb,tip))

    #ax = dtl.plot_axes(200)
    #ax = dtl.plot_polygon(t.root.loop,ax,ls = '--',lw = 2,col = 'k')
    #ax = dtl.plot_polygon(tip.loop,ax,lw = 3,col = 'r')
    #for b in bs:
    #    ax = dtl.plot_polygon(b,ax,lw = 3,col = 'g')
    #plt.show()

    return newtips


