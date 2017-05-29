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
        e1 = vp.tov(up)
        e2 = vp.tov(wp)
        etn1 = e1.cp().nrm()
        etn2 = e2.cp().nrm()
        para  = gtl.isnear(etn1.dot(etn2), 1)
        apara = gtl.isnear(etn1.dot(etn2),-1)
        if apara:sa = numpy.pi
        elif para:sa = 0
        else:sa = e1.sang(e2,vec3(0,0,1))
        if sa < 0:sa = 2*numpy.pi+sa
        return sa

    # compute the edge ordering of a vertex v relative to an edge u,v
    def eo(self,u,v):
        uvas,avor = [],[]
        vor = self.orings[v]
        for w in vor:
            if w == u:uvwa = 0
            else:uvwa = self.ea(u,v,w)
            if not uvas or uvwa > uvas[-1]:
                uvas.append(uvwa);avor.append(w)
            else:
                fnd = False
                for uvax in range(len(uvas)):
                    if uvwa < uvas[uvax]:
                        fnd = True
                        break
                if fnd:
                    uvas.insert(uvax,uvwa)
                    avor.insert(uvax,w)
                else:
                    uvas.append(uvwa)
                    avor.append(w)
        uori = uvas.index(min(uvas))
        ror = avor[uori+1:]+avor[:uori]
        return ror

    # given an edge, find the next edge taking the first
    # clockwise turn available
    def cw(self,u,v):
        vor = self.orings[v][:]
        if len(vor) == 2:vor.remove(u)
        if len(vor) == 1:return vor[0]
        ror = self.eo(u,v)
        if ror:tip = ror[0]
        else:tip = u
        return tip

    # given an edge, find the next edge taking the first
    # counterclockwise turn available
    def ccw(self,u,v):
        vor = self.orings[v][:]
        if len(vor) == 2:vor.remove(u)
        if len(vor) == 1:return vor[0]
        ror = self.eo(u,v)
        if ror:tip = ror[-1]
        else:tip = u
        return tip

    # find a vertex within epsilon of p or create one
    # ie : 0 - split edges if sufficiently close to an edge
    def fp(self,p,epsilon,ie = 1,**vkws):
        for j in range(self.vcnt):
            if self.vs[j][1]['p'].d(p) < epsilon:
                return j
        nv = self.av(p = p.cp(),**vkws)
        if ie:
            for u,v in self.elook:
                up,vp = self.vs[u][1]['p'],self.vs[v][1]['p']
                ped = p.dexy(up,vp)
                if ped > -1 and ped < epsilon:
                    self.se(u,v,nv)
                    break
        return nv

    # find an edge between two vertices or create one 
    def fe(self,u,v,**ekws):
        if   (u,v) in self.elook:return self.elook[(u,v)]
        elif (v,u) in self.elook:return self.elook[(v,u)]
        return self.ae(u,v,**ekws)

    ###################################

    def smooth(self,i = 1,w = 1.0):
        for s in range(i):
            dels = []
            for j in range(self.vcnt):
                v = self.vs[j]
                if v is None:continue
                ops = [self.vs[o][1]['p'] for o in self.rings[j]]
                dels.append(v[1]['p'].tov(vec3(0,0,0).com(ops)).uscl(w*v[1]['w']))
            for j in range(self.vcnt):
                v = self.vs[j]
                if v is None:continue
                v[1]['p'].trn(dels[j])
        return self

    def smooth_sticks(self,i = 1,w = 1.0,d = 5.0):
        for s in range(i):
            dels = []
            for j in range(self.vcnt):
                v = self.vs[j]
                if v is None:continue
                ops = []
                for o in self.rings[j]:
                    vod = v[1]['p'].d(self.vs[o][1]['p'])
                    stickp = v[1]['p'].lerp(self.vs[o][1]['p'],1.0-d/vod)
                    ops.append(stickp)
                dels.append(v[1]['p'].tov(vec3(0,0,0).com(ops)).uscl(w*v[1]['w']))
            for j in range(self.vcnt):
                v = self.vs[j]
                if v is None:continue
                v[1]['p'].trn(dels[j])
        return self

    ###################################

    # plot the vertices and edges of the graph
    def plotxy(self,ax = None,l = 10,s = 1.0,number = True,**kws):
        if ax is None:ax = dtl.plot_axes_xy(l,**kws)
        for j in range(self.vcnt):
            i = self.vs[j]
            if i is None:continue
            ip = i[1]['p']
            if number:
                jstr = str(j)+str([v for v in self.rings[j]])
                jstrccw = str([self.ccw(j,v) for v in self.rings[j]])
                jstrcw = str([self.cw(j,v) for v in self.rings[j]])
                jstr += '\n'+jstrccw+'\n'+jstrcw
                ax = dtl.plot_point_xy_annotate(ip,ax,jstr)
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
    def plot(self,ax = None,l = 10,**kws):
        if ax is None:ax = dtl.plot_axes(l)
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




 
