import dilap.core.base as db
import dilap.core.context as cx
import dilap.modeling.model as dmo
import dilap.modeling.factory as dfa

from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat
import dilap.geometry.tools as gtl
import dilap.geometry.polymath as pym

import dilap.worldly.treeskin as ltr
import dilap.worldly.building as blg
import dilap.worldly.blgsequencing as bseq
import dilap.worldly.partitiongraph as ptg

import dilap.core.plotting as dtl
import matplotlib.pyplot as plt

import math,numpy,random,pdb



###############################################################################
###############################################################################

# wire graph class (purely topological)
class wgraph(db.base):

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

    ###################################

    # compute the counterclockwise angle between two edges
    def ea(self,u,v,w):
        up = self.vs[u][1]['p']
        vp = self.vs[v][1]['p']
        wp = self.vs[w][1]['p']
        e1 = up.tov(vp)
        e2 = up.tov(wp)
        #sa = gtl.wrap(e1.sang(e2,vec3(0,0,1)),0,2*numpy.pi)
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

    ###################################

###############################################################################
###############################################################################

# return a set of loops representing the boundary of all intersections and roads
def gloops(g):

    pdb.set_trace()











###############################################################################

def plot(fp,ixs,rds):
    print('make me a city')
    ax = dtl.plot_axes_xy(100)
    ax = dtl.plot_polygon_xy(fp,ax,col = 'b',lw = 5)
    for e in rds:
        r = rds[e]
        rp1,rp2 = ixs[r['s']][1],ixs[r['e']][1]
        ax = dtl.plot_edges_xy((rp1,rp2),ax,lw = 2,center = True)
    for i in ixs:
        ax = dtl.plot_point_xy(i[1],ax,mk = 's',col = 'r')

#
def placeroad(fp,ixs,rds,p1,p2):
    print('check with existing roads')
    for rk in rds:
        r = rds[rk]
        rp1,rp2 = ixs[r['s']][1],ixs[r['e']][1]

        ip = pym.sintsxyp(rp1,rp2,p1,p2,ie = False)
        if not ip is None:
            print('road hits road!')
            pdb.set_trace()

    # may have returned if ended at another road instead

    print('check with boundary')
    tfp = pym.contract(fp,5)
    ips = pym.sintbxyp(p1,p2,tfp)
    if len(ips) == 0:return p2
    elif len(ips) == 1:
        if p1.inbxy(tfp):
            print('hit boundary!')
            return p1.lerp(ips[0],0.5)
        else:return p2
    elif len(ips) > 1:raise ValueError


# IVE BEEN DOING THIS ALL WRONG...
# THE GOAL IS TO PARTITION THE XY PLANE.. NOT GENERATE BBOXES AND USE POLYGON MATH..
#   DO NOT MAKE A POLYGON OF ALL ROADS IN THE FOOTPRINT...
#   COMPUTE FOOTPRINTS OF THE ROADS AND INTERSECTIONS AND MAKE VERTICES FOR THEM


# generate topology where the first intersection is on the boundary
#   and it is the only resulting intersection with the boundary
def peninsula(fp,ixs,rds):
    ixcnt = 1
    exit = ixs[0][1]


    new = placeroad(fp,ixs,rds,exit,exit.lerp(vec3(0,0,0).com(fp),0.25))
    fi = (ixcnt,new,[0])
    ixs.append(fi)
    ixs[0][2].append(ixcnt)
    rds[(0,ixcnt)] = {'s':0,'e':ixcnt,'type':'generic','lvl':0}
    ixcnt += 1


    tipx = ixcnt-1
    for x in range(2):
        tip = ixs[tipx]
        new = tip[1].cp().trn(vec3(10,0,0))
        new = placeroad(fp,ixs,rds,tip[1],new)

        ni = (ixcnt,new,[tipx])
        ixs.append(ni)
        tip[2].append(ixcnt)
        rds[(tipx,ixcnt)] = {'s':tipx,'e':ixcnt,'type':'generic','lvl':0}

        tipx = ixcnt

        ixcnt += 1


        ax = dtl.plot_axes(100)
        ax = dtl.plot_polygon(fp,ax,lw = 3,col = 'b')
        for y in range(1,len(ixs)):
            ax = dtl.plot_edges((ixs[y-1][1],ixs[y][1]),ax,lw = 3)
            ax = dtl.plot_point(ixs[y-1][1],ax,col = 'r')
            ax = dtl.plot_point(ixs[y][1],ax,col = 'r')
        plt.show()


# based on topology, generate boundaries for roads and intersections
def bboxes(fp,ixs,rds):

    # for each intersection, create a seam for each road
    for i in ixs:
        ring = i[2]
        



    #ibbs = [i[1].pring(3,8) for i in ixs]
    ibbs = [i[1].sq(6,6) for i in ixs]
    #ibbs = []
    rbbs = {}
    for rk in rds:
        r = rds[rk]
        rs,re = ixs[r['s']][1],ixs[r['e']][1]
        rt = rs.tov(re).nrm()
        rn = vec3(0,0,1).crs(rt).uscl(2)
        rs = rs.cp().trn(rt.cp().uscl(-2))
        re = re.cp().trn(rt.cp().uscl( 2))

        #rs = rs.cp().trn(rt.cp().uscl( 1))
        #re = re.cp().trn(rt.cp().uscl(-1))
        #rs = rs.cp()
        #re = re.cp()

        rbb = [
            rs.cp().trn(rn),rs.cp().trn(rn.flp()),
            re.cp().trn(rn),re.cp().trn(rn.flp())]
        #rstn = rt.cp().rot(quat(1,0,0,0).av(numpy.pi/8.0,vec3(0,0,1)))
        #retn = rt.cp().rot(quat(1,0,0,0).av(numpy.pi/8.0,vec3(0,0,1)))
        #l = rbb[0].spline(rbb[3],rstn.cp(),retn.cpf(),3)
        #r = rbb[1].spline(rbb[2],rstn.cp(),retn.cpf(),3)
        rbb = pym.ebdxy(rbb,ibbs[r['s']])
        rbb = pym.ebdxy(rbb,ibbs[r['e']])

        #rbbs.append(rbb)
        rbbs[rk] = rbb

        #ax = dtl.plot_axes(100)
        #ax = dtl.plot_edges(l,ax,lw = 2,col = 'b')
        #ax = dtl.plot_edges(r,ax,lw = 2,col = 'g')
        #plt.show()

    return ibbs,rbbs



    






#def bboxes(fp,ixs,rds):
def _____peninsula(fp,ixs,rds):
    #ibbs = [i[1].pring(5,8) for i in ixs]
    #ibbs = [i[1].sq(2,2) for i in ixs]
    #ibbs = []
    #rbbs = []
    for rk in rds:
        r = rds[rk]
        rs,re = ixs[r['s']][1],ixs[r['e']][1]
        rt = rs.tov(re).nrm()
        rn = vec3(0,0,1).crs(rt).uscl(2)
        rs = rs.cp().trn(rt.cp().uscl(-1))
        re = re.cp().trn(rt.cp().uscl( 1))

        rbb = [
            rs.cp().trn(rn),rs.cp().trn(rn.flp()),
            re.cp().trn(rn),re.cp().trn(rn.flp())]
        rstn = rt.cp().rot(quat(1,0,0,0).av(numpy.pi/8.0,vec3(0,0,1)))
        retn = rt.cpf().rot(quat(1,0,0,0).av(numpy.pi/8.0,vec3(0,0,1)))
        l = rbb[0].spline(rbb[3],rstn.cp(),retn.cpf(),3)
        r = rbb[1].spline(rbb[2],rstn.cp(),retn.cpf(),3)


        #ax = dtl.plot_axes(100)
        #ax = dtl.plot_edges(l,ax,lw = 2,col = 'b')
        #ax = dtl.plot_edges(r,ax,lw = 2,col = 'g')
        #plt.show()


        rds[(sv,ixcnt)] = {'s':sv,'e':ixcnt,'type':ser['type'],'lvl':ser['lvl'],'bb':rbb}
        rds[(ixcnt,ev)] = {'s':ixcnt,'e':ev,'type':ser['type'],'lvl':ser['lvl']}
        del rds[k]

        sv = ixcnt
        ixcnt += 1

    for nix in nixs:
        # connect to exactly two intersections... pick dir and extend toward boundary
        ni = ixs[nix]
        nip = ni[1]
        ip1 = ixs[ni[2][0]][1]
        ip2 = ixs[ni[2][1]][1]
        t1 = ip1.tov(nip).nrm()
        t2 = nip.tov(ip2).nrm()
        t1ct2 = t1.crs(t2)
        if   t1ct2.z > 0:d = 0
        elif t1ct2.z < 0:d = 1
        else:d = random.randint(0,1)
        if   d == 0:nm = vec3(0,0,1).crs(t1)
        elif d == 1:nm = vec3(0,0,1).crs(t1.cp().flp())

        aip = nip.cp().trn(nm.uscl(30))
        aip = placeroad(fp,ixs,rds,nip,aip)

        ixs.append((ixcnt,aip,[nix]))
        ni[2].append(ixcnt)
        rds[(nix,ixcnt)] = {'s':nix,'e':ixcnt,'type':'generic','lvl':0}

        ixcnt += 1

        aip = nip.cp().trn(nm.flp())
        aip = placeroad(fp,ixs,rds,nip,aip)

        ixs.append((ixcnt,aip,[nix]))
        ni[2].append(ixcnt)
        rds[(nix,ixcnt)] = {'s':nix,'e':ixcnt,'type':'generic','lvl':0}

        ixcnt += 1

    print('peninsula generation')
    return
    pdb.set_trace()

def road_graph(fp,exits):
    ixcnt = 0
    ixs = []
    rds = {}

    if len(exits) == 0:
        mc = vec3(0,0,0).com(fp)
        ixs.append((ixcnt,mc,[]))
        ixcnt += 1
        contained(fp,ixs,rds)
    else:
        for e in exits:
            ixs.append((ixcnt,exits[0],[]))
            ixcnt += 1
        if len(exits) == 1:
            peninsula(fp,ixs,rds)
        else:net(fp,ixs,rds)

    plot(fp,ixs,rds)
    plt.show()

    return ixs,rds




    mc = vec3(0,0,0).com(fp)

    ixcnt = 1
    ixs = [(0,mc,[])]
    rds = {}
    for x in range(10):
        i = ixs[-1]
        rdir = growdir(ixs,i)
        rlen = 10 + 2*x
        nip = i[1].cp().trn(rdir.uscl(rlen))
        ixs.append((ixcnt,nip,[i[0]]))
        i[2].append(ixcnt)
        r = {'s':i[0],'e':ixcnt,
            'type':'generic','lvl':0}
        rds[(i[0],ixcnt)] = r
        #rds[(ixcnt,i[0])] = r
        ixcnt += 1



    return ixs,rds


