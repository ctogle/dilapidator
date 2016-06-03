import dilap.core.base as db

from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat

import dilap.geometry.tools as gtl
import dilap.geometry.polymath as pym

import dilap.core.plotting as dtl
import matplotlib.pyplot as plt

import math,numpy,random,pdb



###############################################################################
###############################################################################

# planar wire graph class (purely topological)
class graph(db.base):

    ###################################
    ### topological methods
    ###################################

    # add a new vertex and edge to another vertex
    def mev(self,ov,vkws,ekws):
        nv = self.av(**vkws)
        ne = self.ae(ov,nv,**ekws)
        return nv,ne

    # add a new intersection to the graph
    def av(self,**vkws):
        j = self.vcnt
        i = (j,vkws)
        self.vs.append(i)
        self.rings[j] = {}
        self.orings[j] = []
        self.vcnt += 1
        return j

    # remove vertex u
    def rv(self,u):
        i = self.vs[u]
        ur = self.rings[u]
        for v in list(ur.keys()):self.re(u,v)
        self.vs[u] = None
        del self.rings[u]
        del self.orings[u]
        return i

    # and a new road to the graph
    def ae(self,u,v,**ekws):
        m = self.ecnt
        ur,vr = self.rings[u],self.rings[v]
        uor,vor = self.orings[u],self.orings[v]
        r = (m,ekws)
        self.es.append(r)
        if not v in ur:ur[v] = r
        if not u in vr:vr[u] = r
        urcnt = len(uor)
        if urcnt < 2:uor.append(v)
        else:
            w = uor[0]
            nea = self.ea(u,w,v)
            f = False
            for eax in range(1,urcnt):
                if nea < uor[eax]:
                    uor.insert(eax,v)
                    f = True
                    break
            if not f:uor.append(v)
        vrcnt = len(vor)
        if vrcnt < 2:vor.append(u)
        else:
            w = vor[0]
            nea = self.ea(v,w,u)
            f = False
            for eax in range(1,vrcnt):
                if nea < vor[eax]:
                    vor.insert(eax,u)
                    f = True
                    break
            if not f:vor.append(u)
        self.ecnt += 1
        return m

    # remove an edge between u and v
    def re(self,u,v):
        r = self.rings[u][v]
        if r is None:return r
        self.es[r[0]] = None
        del self.rings[u][v]
        del self.rings[v][u]
        if v in self.orings[u]:self.orings[u].remove(v)
        if u in self.orings[v]:self.orings[v].remove(u)
        return r

    # split an edge/road into two edges/roads
    def se(self,u,v,w):
        ruv = self.re(u,v)
        ruw = self.ae(u,w,**ruv[1])
        rwv = self.ae(w,v,**ruv[1])
        return ruv,rwv

    # split an edge/road into two edges/roads
    # make the third vertex by interpolating between the other two
    def be(self,u,v,f = 0.5,**kws):
        if not 'p' in kws:
            up,vp = self.vs[u][1]['p'],self.vs[v][1]['p']
            kws['p'] = up.lerp(vp,f)
        nv = self.av(**kws)
        nr1,nr2 = self.se(u,v,nv)
        return nv,nr1,nr2

    ###################################

    # compute the counterclockwise angle between two edges
    def ea(self,u,v,w):
        up = self.vs[u][1]['p']
        vp = self.vs[v][1]['p']
        wp = self.vs[w][1]['p']
        e1 = up.tov(vp)
        e2 = up.tov(wp)
        sa = e1.sang(e2,vec3(0,0,1))
        if sa < 0:sa = 2*numpy.pi-sa
        return sa

    # return a list of vertex indices which form a loop
    #   the first edge will be from u to v, turns of direction 
    #   d (clockwise or counterclockwise) form the loop
    def loop(self,u,v,d = 'cw'):
        if not v in self.rings[u]:
            raise ValueError
        lp = [u,v]
        while True:
            uor = self.orings[lp[-2]]
            vor = self.orings[lp[-1]]
            uori = vor.index(lp[-2])
            ror = vor[uori+1:]+vor[:uori]
            if ror:
                if d == 'cw':tip = ror[0]
                elif d == 'ccw':tip = ror[-1]
                else:raise ValueError
            else:tip = lp[-2]
            lp.append(tip)
            if lp[-1] == lp[1] and lp[-2] == lp[0]:
                lp.pop(-1)
                lp.pop(-1)
                return lp

    # return a list of all unique loops of the graph
    def uloops(self,d = 'cw'):
        loops = {}
        unfn = [x for x in range(self.ecnt)]
        for vx in range(self.vcnt):
            v = self.vs[vx]
            if v is None:continue
            for ox in self.orings[vx]:
                r = self.rings[vx][ox]
                rx,rkws = r
                if rx in unfn:
                    lp = self.loop(vx,ox,'ccw')
                    lpk = tuple(set(lp))
                    if not lpk in loops:loops[lpk] = lp
                    unfn.remove(rx)
                else:continue
            if not unfn:break
        return [loops[lpk] for lpk in loops]

    # return a full polygon representing the partitioning of 
    # the xy plane affected by the edges of this graph
    def polygon(self,r = 2,d = 'cw'):
        loops,seams,eseam = self.uloops(d),[],None
        for lp in loops:
            seam = []
            seams.append(seam)
            for lpx in range(len(lp)):
                lp0,lp1,lp2 = lp[lpx-2],lp[lpx-1],lp[lpx]
                lpp0 = self.vs[lp0][1]['p']
                lpp1 = self.vs[lp1][1]['p']
                lpp2 = self.vs[lp2][1]['p']
                lptn1 = lpp0.tov(lpp1).nrm()
                lptn2 = lpp1.tov(lpp2).nrm()
                lpnm1 = vec3(0,0,1).crs(lptn1).uscl(r)
                lpnm2 = vec3(0,0,1).crs(lptn2).uscl(r)
                if lp0 == lp2:
                    stemoffset = lptn1.cp().uscl(r)
                    seam.append(lpp1.cp().trn(lpnm1+stemoffset))
                    seam.append(lpp1.cp().trn(lpnm2+stemoffset))
                else:
                    if lpnm1.isnear(lpnm2):cnm = lpnm1
                    else:
                        s1 = lpp0.cp().trn(lpnm1).trn(lptn1.cp().uscl(-1000))
                        s2 = lpp1.cp().trn(lpnm1).trn(lptn1.cp().uscl( 1000))
                        s3 = lpp1.cp().trn(lpnm2).trn(lptn2.cp().uscl(-1000))
                        s4 = lpp2.cp().trn(lpnm2).trn(lptn2.cp().uscl( 1000))
                        ip = pym.sintsxyp(s1,s2,s3,s4,col = 0)
                        cnm = lpp1.tov(ip)
                    seam.append(lpp1.cp().trn(cnm))
            if eseam is None:eseam = 0
            else:
                if pym.binbxy(seams[eseam],seam):
                    seams.append(eseam)
                    eseam = len(seams)-1
        py = (tuple(seams.pop(eseam)),tuple(tuple(s) for s in seams))
        return py

    ###################################

    # plot the vertices and edges of the graph
    def plot(self,ax = None):
        if ax is None:ax = dtl.plot_axes(50)
        for j in range(self.vcnt):
            i = self.vs[j]
            if i is None:continue
            ip = i[1]['p']
            ax = dtl.plot_point(ip,ax,col = 'r')
        for k in self.rings:
            vr = self.rings[k]
            for ov in vr:
                if vr[ov] is None:continue
                re = (self.vs[k][1]['p'].cp(),self.vs[ov][1]['p'].cp())
                rn = vec3(0,0,1).crs(re[0].tov(re[1])).nrm()
                re = (re[0].trn(rn),re[1].trn(rn))
                ax = dtl.plot_edges(re,ax,lw = 2,col = 'g')
        return ax

    ###################################

    def __init__(self):
        self.vs = []
        self.vcnt = 0
        self.es = []
        self.ecnt = 0
        self.rings = {}
        self.orings = {}

###############################################################################




 
