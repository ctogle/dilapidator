from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat
import dilap.geometry.polymath as pym

import dilap.topology.planargraph as pgr

import dilap.core.plotting as dtl
import matplotlib.pyplot as plt

import pdb



class partitiongraph(pgr.planargraph):

    # break a vertex based on intersection with a line (adds a vertex)
    def bv(self,vx,p,d,**kws):
        b = self.vs[vx][1]['b']
        p = pym.ptob(b[0],p)
        sp1 = p.cp().trn(d.cp().uscl( 10000))
        sp2 = p.cp().trn(d.cp().uscl(-10000))
        nfps = pym.bsegsxy(b[0],sp1,sp2)
        new = []
        l = nfps.pop(0)
        for r in nfps:new.append(self.sv(vx,l,r,**kws))
        return new

    # split a vertex using two new boundaries (adds a vertex)
    def sv(self,vx,l,r,connect_split = True,connect_ring = True):
        sv = self.vs[vx]
        if pym.binbxy(r,l):sv[1]['b'][1].append(r)
        sv[1]['b'][0] = l
        sv[1]['p'] = vec3(0,0,0).com(l)
        vinfo = {'b':[r,[]],'p':vec3(0,0,0).com(r),'l':0}
        for k in sv[1]:
            if not k in vinfo:
                vinfo[k] = sv[1][k]
        new = self.av(**vinfo)
        if pym.binbxy(l,r):self.vs[new][1]['b'][1].append(l)
        ekws = {'p':0.25}
        if connect_ring:
            for ringvx in self.orings[vx]:
                self.ae(new,ringvx,**ekws)
        if connect_split:self.ae(new,vx,**ekws)
        return new

    # insert another graph into a vertex 
    def ig(self,vx,l,r,connect_split = True,connect_ring = True):
        raise NotImplemented

    # verify the geometric viability of each edge of a vertex
    # when an exit resides geometrically between two vertices, add an edge
    def vve(self,vx):
        v = self.vs[vx]
        if v is None:return
        vb = v[1]['b']
        rem = []
        for adj in self.orings[vx]:
            ov = self.vs[adj]
            if ov is None:
                print('vverem')
                rem.append(adj)
                continue
            ob = ov[1]['b']
            aws = pym.badjbxy(vb[0],ob[0])
            if not aws:self.re(adj,vx)
        for rm in rem:self.orings[vx].remove(rm)

    # perform edge verification for all vertices
    def vves(self):
        for vx in range(self.vcnt):self.vve(vx)

    # produce a planar graph such that boundaries of vertices become edges
    # NOTE: loops of this new graph correspond to vertices of self
    def bgraph(self,epsilon = 0.1):
        bg = pgr.planargraph()
        for vx in range(self.vcnt):
            v = self.vs[vx]
            if v is None:continue
            eb = v[1]['b'][0]

            '''#
            # check if intersections happen in between!
            fnd = True
            while fnd:
                fnd = False
                for ox in range(self.vcnt):
                    o = self.vs[ox]
                    if o is None or ox == vx:continue
                    ob = o[1]['b'][0]
                    adjs = pym.badjbxy(eb,ob,0.1)
                    acnt = len(adjs)
                    if acnt == 0:continue
                    else:
                        if acnt > 1:print('multiadjacency!!')
                        for adj in adjs:
                            vbx,obx = adj
                            vb1,vb2 = eb[vbx-1],eb[vbx]
                            ob1,ob2 = ob[obx-1],ob[obx]
                            ips = pym.sintsxyp(vb1,vb2,ob1,ob2,ieb = 0,skew = 0,ie = 0)
                            if ips is None:continue
                            ip1,ip2 = ips
                            if not ip1 in eb or not ip2 in eb:
                                if ip1.onsxy(vb1,vb2,0):
                                    eb.insert(vbx,ip1);fnd = True
                                if ip2.onsxy(vb1,vb2,0):
                                    eb.insert(vbx,ip2);fnd = True
                        if fnd:break
            '''#

            for ebx in range(len(eb)):
                e1,e2 = eb[ebx-1],eb[ebx]
                v1,v2 = bg.fp(e1,epsilon),bg.fp(e2,epsilon)
                if not v2 in bg.rings[v1]:bg.ae(v1,v2)

        return bg

    def plotxy(self,ax = None,scale = 25):
        if ax is None:ax = dtl.plot_axes_xy(scale)
        ax = pgr.planargraph.plotxy(self,ax)
        ekeys = []
        for vx in range(self.vcnt):
            v = self.vs[vx]
            if v is None:continue
            vb = v[1]['b']
            vrtls = '-'
            if 't' in v[1]:vrtls = '--' if 'developed' in v[1]['t'] else '-'
            ax = dtl.plot_polygon_xy(pym.contract(vb[0],1),
                ax,lw = 2,ls = vrtls,col = 'b')
            for ib in vb[1]:
                ax = dtl.plot_polygon_xy(pym.contract(ib,-1),ax,lw = 2,col = 'r')
            vbc = vec3(0,0,0).com(vb[0])
            rs = str(vx)+','+str(self.orings[vx])#+',\n'+str(v[1]['type'])
            ax = dtl.plot_point_xy_annotate(vbc,ax,rs)
            for ve in self.orings[vx]:
                ov = self.vs[ve]
                if ov is None:continue

                if not (vx,ve) in ekeys and not (ve,vx) in ekeys:
                    ekeys.append((vx,ve))
                    aws = pym.badjbxy(vb[0],ov[1]['b'][0],0)
                    for awx in range(len(aws)):
                        vwx,owx = aws[awx]
                        ex,ekws = self.es[self.elook[ekeys[-1]]]
                        door = vb[0][vwx-1].lerp(vb[0][vwx],ekws['p'])
                        ax = dtl.plot_point_xy(door,ax,col = 'r')

                ovbc = vec3(0,0,0).com(self.vs[ve][1]['b'][0])
                ax = dtl.plot_edges_xy((vbc,ovbc),ax,lw = 3,col = 'c')

            #for ep in v[1]['exits']:
            #    ax = dtl.plot_point_xy(ep,ax,mk = 's',col = 'r')

        return ax



    

