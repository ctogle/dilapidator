import dilap.core.base as db

import dilap.core.plotting as dtl
import matplotlib.pyplot as plt

import pdb



###############################################################################

# wire graph class (purely topological)
#
#   wire graph is a simple data structure 
#     for vertices and edges connecting them
#   it also provides a means for ordering edges
#
class wiregraph(db.base):

    ###################################
    ### fundamental topological methods
    ###################################

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
        ur,vr = self.rings[u],self.rings[v]
        uor,vor = self.orings[u],self.orings[v]
        m = self.ecnt
        r = (m,ekws)
        self.es.append(r)
        self.elook[(u,v)] = m
        self.elook[(v,u)] = m
        if not v in ur:ur[v] = r
        if not u in vr:vr[u] = r

        urcnt = len(uor)
        uor.append(v)
        '''#
        if urcnt < 2:uor.append(v)
        else:
            w = uor[0]
            nea = self.ea(u,w,v)
            f = False
            for eax in range(1,urcnt):
                if nea < self.ea(u,w,uor[eax]):
                    uor.insert(eax,v)
                    f = True
                    break
            if not f:uor.append(v)
        '''#
              
        vrcnt = len(vor)
        vor.append(u)
        '''#
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
        '''#

        self.ecnt += 1

        #self.plotxy(l = 200)
        #plt.show()

        return m

    # remove an edge between u and v
    def re(self,u,v):
        r = self.rings[u][v]
        if r is None:return r
        self.es[r[0]] = None
        del self.elook[(u,v)]
        del self.elook[(v,u)]
        del self.rings[u][v]
        del self.rings[v][u]
        if v in self.orings[u]:self.orings[u].remove(v)
        if u in self.orings[v]:self.orings[v].remove(u)
        return r

    ###################################
    ### additional topological methods
    ###################################

    # add a new vertex and edge to another vertex
    def mev(self,ov,vkws,ekws):
        nv = self.av(**vkws)
        ne = self.ae(ov,nv,**ekws)
        return nv,ne

    # split an edge/road into two edges/roads
    def se(self,u,v,w):
        ruv = self.re(u,v)
        ruw = self.ae(u,w,**ruv[1])
        rwv = self.ae(w,v,**ruv[1])
        return ruv,rwv

    ###################################
    ### edge ordering mechanisms
    ###################################

    # provide an ordering mechanism for edges in the graph
    #  the default behavoir is the order of edge addition
    def ea(self,u,w,v):
        return 0

    def cw(self,u,v):
        #uor = self.orings[u]
        vor = self.orings[v]
        uori = vor.index(u)
        ror = vor[uori+1:]+vor[:uori]
        if ror:tip = ror[0]
        else:tip = u
        return tip

    def ccw(self,u,v):
        #uor = self.orings[u]
        vor = self.orings[v]
        uori = vor.index(u)
        ror = vor[uori+1:]+vor[:uori]
        if ror:tip = ror[-1]
        else:tip = u
        return tip

    # return a list of vertex indices which form a loop
    #   the first edge will be from u to v, turns of direction 
    #   d (clockwise or counterclockwise) form the loop
    def loop(self,u,v,d = 'cw'):
        if not v in self.rings[u]:raise ValueError
        lp = [u,v]

        c = 0

        while True:

            c += 1
            if c > 250:
                print('shit',d,u,v,len(lp))
                for j in range(self.vcnt):
                    print(self.vs[j][1]['p'])
                ax = dtl.plot_axes_xy(400)
                ax = self.plotxy(ax)
                ax = dtl.plot_polygon_xy([self.vs[lx][1]['p'] for lx in lp],ax)
                #ax = dtl.plot_points_xy([self.vs[lx][1]['p'] for lx in lp],ax,number = True)
                plt.show()
                pdb.set_trace()

            if   d == 'cw': tip = self.cw( lp[-2],lp[-1])
            elif d == 'ccw':tip = self.ccw(lp[-2],lp[-1])
            else:raise ValueError
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
                    if not lpk in loops:
                        loops[lpk] = lp
                    unfn.remove(rx)
                else:continue
            if not unfn:break
        lps = [loops[lpk] for lpk in loops]
        return lps

    ###################################

    def __init__(self):
        self.vs = []
        self.vcnt = 0
        self.es = []
        self.ecnt = 0
        self.elook = {}
        self.rings = {}
        self.orings = {}

    ###################################




 
