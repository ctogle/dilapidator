from dilap.geometry.vec3 import vec3
import dilap.geometry.tools as gtl

import dilap.topology.wiregraph as dwg

import dilap.core.plotting as dtl
import matplotlib.pyplot as plt

import numpy,pdb



###############################################################################

# planar wire graph class (topological + xy-projected geometry)
class planargraph(dwg.wiregraph):

    # split an edge into two edges
    # make the third vertex by lerping between the other two
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

        etn1 = e1.cp().nrm()
        etn2 = e2.cp().nrm()
        para = gtl.isnear(abs(etn1.dot(etn2)),1)
        if para:sa = numpy.pi
        else:sa = e1.sang(e2,vec3(0,0,1))

        #print('easa1',sa,para,etn1.dot(etn2))
        #print('alkdfj',up,vp,wp)

        if sa <= 0:sa = numpy.pi+sa

        '''#
        ax = dtl.plot_axes_xy(500)
        ax = dtl.plot_edges_xy((vp,wp),ax,lw = 2,col = 'b')
        ax = dtl.plot_edges_xy((up,vp),ax,lw = 2,col = 'g')
        ax = dtl.plot_point_xy_annotate(up,ax,'u '+str(gtl.deg(sa)))
        ax = dtl.plot_point_xy(up,ax)
        ax = dtl.plot_point_xy_annotate(vp,ax,'v ')
        ax = dtl.plot_point_xy(vp,ax)
        ax = dtl.plot_point_xy_annotate(wp,ax,'w ')
        ax = dtl.plot_point_xy(wp,ax)
        plt.show()
        '''#

        #print('easa2',sa)

        return sa

    # given an edge, find the next edge taking the first
    # clockwise turn available
    def cw(self,u,v):
        uor = self.orings[u]
        vor = self.orings[v]
        uori = vor.index(u)
        ror = vor[uori+1:]+vor[:uori]
        if ror:tip = ror[0]
        else:tip = u
        return tip

    # given an edge, find the next edge taking the first
    # counterclockwise turn available
    def ccw(self,u,v):
        uor = self.orings[u]
        vor = self.orings[v]
        uori = vor.index(u)
        ror = vor[uori+1:]+vor[:uori]
        if ror:tip = ror[-1]
        else:tip = u
        return tip

    # find a vertex within epsilon of p or create one
    def fp(self,p,epsilon,**vkws):
        for j in range(self.vcnt):
            if self.vs[j][1]['p'].d(p) < epsilon:
                return j
        return self.av(p = p.cp(),**vkws)

    ###################################

    # plot the vertices and edges of the graph
    def plotxy(self,ax = None,l = 10,s = 1.0,number = True):
        if ax is None:ax = dtl.plot_axes_xy(l)
        for j in range(self.vcnt):
            i = self.vs[j]
            if i is None:continue
            ip = i[1]['p']
            if number:ax = dtl.plot_point_xy_annotate(ip,ax,str(j))
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




 
