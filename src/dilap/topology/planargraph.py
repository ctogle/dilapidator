from dilap.geometry.vec3 import vec3

import dilap.topology.wiregraph as dwg

import dilap.core.plotting as dtl
import matplotlib.pyplot as plt

import math,numpy,pdb



###############################################################################
###############################################################################

# planar wire graph class (topological + xy-projected geometry)
class planargraph(dwg.wiregraph):

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

    # find a vertex within epsilon of p or create one
    def fp(self,p,epsilon,**vkws):
        for j in range(self.vcnt):
            if self.vs[j][1]['p'].d(p) < epsilon:
                return j
        return self.av(p = p.cp(),**vkws)

    ###################################

    # return a full polygon representing the partitioning of 
    # the xy plane affected by the edges of this graph
    def polygon(self,r = 2,d = 'cw'):

        # TEMPORARY HACK -> FIX THIS DAMNIT!!
        import dilap.geometry.polymath as pym

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
        for sx in range(len(seams)):
            if pym.bnrm(seams[sx]).z < 0:seams[sx].reverse()
        #py = (tuple(seams.pop(eseam)),tuple(tuple(s) for s in seams))
        py = [seams.pop(eseam),seams]
        return py

    ###################################

    # plot the vertices and edges of the graph
    def plotxy(self,ax = None,l = 10,s = 1.0):
        if ax is None:ax = dtl.plot_axes_xy(l)
        for j in range(self.vcnt):
            i = self.vs[j]
            if i is None:continue
            ip = i[1]['p']
            ax = dtl.plot_point_xy_annotate(ip,ax,str(j))
            ax = dtl.plot_point_xy(ip,ax,col = 'r')
        for k in self.rings:
            vr = self.rings[k]
            for ov in vr:
                if vr[ov] is None:continue
                re = (self.vs[k][1]['p'].cp(),self.vs[ov][1]['p'].cp())
                rn = vec3(0,0,1).crs(re[0].tov(re[1])).nrm().uscl(s)
                re = (re[0].trn(rn),re[1].trn(rn))
                ax = dtl.plot_edges_xy(re,ax,lw = 2,col = 'g')
        return ax

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




 
