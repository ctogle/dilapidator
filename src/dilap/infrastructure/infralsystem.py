import dilap.core.base as db

import dilap.core.tools as dpr
import dilap.core.lsystem as dls

import dp_vector as dpv

import matplotlib.pyplot as plt
import random as rm
import pdb



class infralsystem(dls.lsystem):

    def _realize(self,p,d):
        self.edges = []
        self.nodes = []
        return dls.lsystem._realize(self,p,d)

    loadouts = []
    loadouts.append(('FX',[
        ('X','X/YF/'),
        ('Y','\FX\Y')]))

    def __init__(self,ldx,*args,**kwargs):
        axiom,rules = self.loadouts[ldx]
        self._def('axiom',axiom,**kwargs)
        self._def('rules',rules,**kwargs)

        self._def('seed',1,**kwargs) # seed used for random numbers
        self._def('iterations',5,**kwargs) # number of iterations

        self._def('rho',100.0,**kwargs) # the length of a new edge
        self._def('angle',dpr.PI/2.0,**kwargs) # angle used for rotations
        self._def('minangle',dpr.PI/12.0,**kwargs) # angle used for rotations
        self._def('maxangle',dpr.PI*2.0,**kwargs) # angle used for rotations

        self._def('branchdraw',self.draw_branch,**kwargs)
        self._def('leafdraw',self.draw_leaf,**kwargs)
        self._def('finaldraw',self.draw,**kwargs)
        dls.lsystem.__init__(self,*args,**kwargs)

    def draw_branch(self,p,n):
        self.edges.append((p.copy(),n.copy()))
        dls.draw_branch(p,n)

    def draw_leaf(self,p):
        self.nodes.append((p.copy()))
        dls.draw_leaf(p)

    def draw(self):
        es,ls = self.edges,self.nodes
        nps = []
        nes = []
        for ex in range(len(es)):
            e1,e2 = es[ex]
            e1x,e2x = None,None
            for pdx in range(len(nps)):
                if e1.near(nps[pdx]):e1x = pdx
                if e2.near(nps[pdx]):e2x = pdx
                if not e1x is None and not e2x is None:break
            if e1x is None:
                e1x = len(nps)
                nps.append(e1)
            if e2x is None:
                e2x = len(nps)
                nps.append(e2)
            nes.append((e1x,e2x))
        self.nodes,self.edges = nps,nes
        dls.draw()

def lsystem_graph():
    g = graph()

    l = 0
    p = dpv.zero()
    d = dpv.xhat.copy()
    lsys = infralsystem(0)._realize(p,d)
    ndps,edgs = lsys.nodes,lsys.edges

    ndxs = []
    nexs = []
    for ndp in ndps:ndxs.append(g._node(node(ndp,layer = l))[l])
    for edg in edgs:nexs.append(g._connect_nodes(ndxs[edg[0]],ndxs[edg[1]]))
    return g











