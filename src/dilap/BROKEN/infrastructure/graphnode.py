import dilap.core.base as db
import dilap.core.vector as dpv
import dilap.core.tools as dpr
import dilap.core.lsystem as dls
import dilap.mesh.tools as dtl
#import dilap.infrastructure.tools as itl

import matplotlib.pyplot as plt
import random as rm
import pdb



class node(db.base):

    # d is the direction to o.p from self.p
    def connect(self,d,o):
        self.targetring[o.index] = d
        self.ring[o.index] = d

    # o is another node currently connected to self
    def disconnect(self,o):
        del self.targetring[o.index]
        del self.ring[o.index]

    def plot(self,ax = None):
        if ax is None:ax = dtl.plot_axes()
        plotp = self.p.copy()
        dtl.plot_point(plotp,ax)
        for d in self.ring.keys():
            f = plotp.copy().translate(self.spikes[d])
            dtl.plot_edges([plotp,f],ax)
        return ax

    def plot_xy(self,ax = None):
        if ax is None:ax = dtl.plot_axes_xy()
        dtl.plot_point_xy(self.p,ax)
        for d in self.ring.keys():
            f = self.p.copy().translate(self.spikes[d])
            dtl.plot_edges_xy([self.p,f],ax)
        return ax

    def key(self):return (self.p.x,self.p.y,self.layer)
    def __init__(self,p,**kwargs):
        self._def('index',None,**kwargs)
        self._def('layer',0,**kwargs)
        self.p = p
        self.targetring = {}
        self.ring = {}
        self.spikes = {}

    def _node_radius(self,graph):
        rwidths = []
        for oi in self.targetring:
            skey,okey = self.key(),graph.nodes[oi].key()
            eg = graph.edges[graph._find_edge(skey,okey)]
            rwidths.append(eg.width)
        return 0.5*max(rwidths)+1

    def _spikes_nudge(self,graph,target):
        ws = []
        angles = []
        tkeys = list(self.targetring.keys())
        for tk in tkeys:
            ekey = (self.key(),graph.nodes[tk].key())
            edx = graph.edges_lookup[ekey]
            ws.append(1.0 if graph.edges[edx].interpolated else 0.0)
            angles.append(self.targetring[tk])

        iangles = nudge(angles,ws,target = target)
        nr = self._node_radius(graph)
        for tk,na in zip(tkeys,iangles):
            self.ring[tk] = na
            #spike = dpv.xhat.copy().rotate_z(dpr.rad(na)).scale_u(8)
            spike = dpv.vector(1,0,0).rotate_z(dpr.rad(na)).scale_u(nr)
            self.spikes[tk] = spike

    # update the actual positions of spikes based on targetring
    def _spikes(self,graph):
        self._retarget(graph)
        self.ring,self.spikes = {},{}
        tkys = list(self.targetring.keys())
        rcnt = len(tkys)
        if rcnt == 0:return
        else:
            if rcnt == 2:t = 180
            else:t = 90
            self._spikes_nudge(graph,t)

    # verify the directions in target ring and update 
    # tangents of any connecting edges
    def _retarget(self,graph):
        rgxs = list(self.ring.keys())
        for rgx in rgxs:
            rgnd = graph.nodes[rgx]

            rkey = (rgnd.key(),self.key())
            regx = graph.edges_lookup[rkey]
            if regx is None:
                del self.targetring[rgx]
                del self.ring[rgx]
                continue

            reg = graph.edges[regx]
            ndir1,ndir2 = reg._directions()
            if not reg.one is self:ndir1,ndir2 = ndir2,ndir1
            self.targetring[rgx] = ndir1
            graph.nodes[rgx].targetring[self.index] = ndir2

###############################################################################
### FUNCTIONS TO BE RELOCATED EVENTUALLY...
###############################################################################

def adist(a1,a2):
    da = dpr.clamp_periodic(a1-a2,0.0,360.0)
    return da if da < 180 else 360 - da

def signedadist(a1,a2):
    return a1-a2 if a1 > a2 else a2-a1

def apply_err(x,e,m):
    xpe = dpr.clamp_periodic(x+e,0.0,360.0)
    xme = dpr.clamp_periodic(x-e,0.0,360.0)
    da1,da2 = adist(xpe,m),adist(xme,m)
    if da1 > da2:return xpe
    else:return xme

def min_adist(alpha,angles,exempt = None):
    acnt = len(angles)
    if acnt < 2:return None,None
    else:
        minad = None
        minax = None
        for a2x in range(acnt):
            if a2x == exempt:continue
            ad = adist(angles[a2x],alpha)
            if minad is None or ad < minad:
                minad = ad
                minax = a2x
        return minad,minax
# given a list of angles, gradually move them to satisfy a condition
# the condition is that da for each angle be nearly equal to their avg
# da is the minimum distance to the nearest angle, which should be ~90
# ws are weights to affect the motion per nudge
def nudge(angles,ws,target = 90,error = 1):
    def acceptable(das):return abs(min(das)-target) < error
    def measure(angles):
        das = []
        dxs = []
        for x in range(acnt):
            xda,xdx = min_adist(angles[x],angles,x)
            das.append(xda)
            dxs.append(xdx)
        return das,dxs

    oas = angles[:]
    acnt = len(angles)
    if acnt == 1:return oas
    das,dxs = measure(oas)
    tries = 0
    maxtries = 100
    while not acceptable(das) and tries < maxtries:
        tries += 1
        das,dxs = measure(oas)
        tdas = [dpr.clamp(target-da,0,target) for da in das]
        es = [tda*w*0.5 for tda,w in zip(tdas,ws)]
        oas = [apply_err(oas[x],es[x],oas[dxs[x]]) for x in range(acnt)]
    if tries == maxtries:print('tries exceeded while nuding!')
    return oas

###############################################################################
### FUNCTIONS TO BE RELOCATED EVENTUALLY...
###############################################################################








