import dilap.core.lsystem as lsy
import dilap.core.context as cx
import dilap.geometry.tools as gtl
from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat
from dilap.geometry.pointset import pointset
import dilap.geometry.polymath as pym
import dilap.modeling.model as dmo
import dilap.topology.trimesh as dtm
import dilap.topology.tree as dtr
import dilap.topology.vert as dvt
import dilap.topology.edge as deg
import dilap.topology.loop as dlp
import dilap.topology.face as dfc

#import dilap.topology.tools.triangulate as dtg
import dilap.geometry.triangulate as dtg

import dilap.core.plotting as dtl
import matplotlib.pyplot as plt

import numpy,math
import pdb



###############################################################################
###
###############################################################################

class building(cx.context):

    def __init__(self,*args,**kwargs):
        self._def('name','buildingcontext',**kwargs)
        cx.context.__init__(self,*args,**kwargs)

    def asurf(self,poly,fm = 'generic',
            p = None,q = None,s = None,
            m = None,gm = None,rv = False):
        if gm is None:
            if m is None:
                m = dmo.model()
                sgv = self.amodel(p,q,s,m,None)
            gm = m.agfxmesh()
        eb,ibs = poly
        hmin,ref,smo = 1,False,False
        eb,ibs = dtg.split_nondelauney_edges(eb,ibs)
        if ref:hmin,eb,ibs = dtg.split_nondelauney_edges_chew1(eb,ibs)
        tris,bnds = dtg.triangulate(eb,ibs,hmin,ref,smo)
        #tris,bnds = [],[]
        for tri in tris:
            p1,p2,p3 = tri
            n = gtl.nrm(p1,p2,p3)
            v1 = gm.avert(*m.avert(p1,n))
            v2 = gm.avert(*m.avert(p2,n))
            v3 = gm.avert(*m.avert(p3,n))
            if rv:f1 = gm.aface(v1,v3,v2,fm) 
            else:f1  = gm.aface(v1,v2,v3,fm) 
        return m

    def awall(self,p1,p2,wh,hs = None,p = None,q = None,s = None,m = None,gm = None):
        eb = [p1,p2,p2.cp().ztrn(wh),p1.cp().ztrn(wh)]
        ibs = []
        hpts = []
        if not hs is None:
            for h in hs:
                # which side of eb,where along it, how wide, how tall, dist from side
                dx,dp,dw,dh,dz = h
                ddp = dw/(2.0*p1.d(p2))
                d1,d2 = p1.lerp(p2,dp-ddp),p1.lerp(p2,dp+ddp)
                h = [d2.cp().ztrn(dz),d2.cp().ztrn(dz+dh),
                    d1.cp().ztrn(dz+dh),d1.cp().ztrn(dz)]
                hpts.append(h)
                if dz == 0:
                    for hp in h:eb.insert(1,hp)
                else:ibs.append(h)
        wpy = (tuple(eb),tuple(tuple(x) for x in ibs))
        m = self.asurf(wpy,m = m,gm = gm)
        return m,hpts

    def aportal(self,oloop,p1,p2,ww,p = None,q = None,s = None,m = None,gm = None):
        pn = vec3(0,0,1).crs(p1.tov(p2).nrm()).uscl(ww)
        iloop = [p.cp().trn(pn) for p in oloop]
        for lx in range(len(oloop)):
            spy = ((oloop[lx-1],iloop[lx-1],iloop[lx],oloop[lx]),())
            m = self.asurf(spy,m = m,gm = gm,rv = True if lx > 0 else False)
        return m

    def aroom(self,**kws):
        m = dmo.model()
        sgv = self.amodel(None,None,None,m,None)
        gm = m.agfxmesh()
        fbnd,cbnd = kws['floor'],kws['ceiling']
        m = self.asurf(fbnd,m = m,gm = gm)
        m = self.asurf(cbnd,m = m,gm = gm,rv = True)
        for x in range(len(fbnd[0])):
            w2,w1 = fbnd[0][x-1],fbnd[0][x]
            hs,wh,ww = kws['wholes'][x-1],kws['wheight'],kws['wwidth']
            wt = kws['wtypes'][x-1]
            m,hpts = self.awall(w1,w2,wh,hs,m = m,gm = gm)
            kws['wholepoints'][x-1].extend(hpts)
            if wt == 'i':
                for hps in hpts:
                    self.aportal(hps,w1,w2,ww,m = m,gm = gm)
        return m

    def contract(self,ps,ds):
        fbnd = []
        for x in range(len(ps)):
            p1,p2,p3 = ps[x-2],ps[x-1],ps[x]
            pn = p2.cp().trn((p2.tov(p1).nrm().mid(p2.tov(p3).nrm())).nrm())
            fbnd.append(p2.lerp(pn,math.sqrt(2)*ds))
        fbnd.append(fbnd.pop(0))
        return fbnd

    def edgehole(self,r1,r2):
        bnd = r1['bound']
        oeb = r2['bound']
        edist,pair = None,None
        for oex in range(len(oeb)):
            oe1,oe2 = oeb[oex-1],oeb[oex]
            otn = oe1.tov(oe2).nrm()
            for ex in range(len(bnd)):
                e1,e2 = bnd[ex-1],bnd[ex]
                etn = e1.tov(e2).nrm()
                eds = oe1.d(e2)+oe2.d(e1)
                if edist is None or eds < edist:
                    edist = eds
                    pair = (ex,oex)
        if not pair is None:
            ex,oex = pair
            r1['wholes'][ex-1].append((0,0.5,1.5,3,0))
            r2['wholes'][oex-1].append((0,0.5,1.5,3,0))
            r1['wtypes'][ex-1] = 'i'
            r2['wtypes'][oex-1] = 'i'

    def roofhole(self,r1,r2):
        f1 = r1['floor']
        fh = vec3(0,0,0).com(f1[0]).sq(4,4)
        r1['floor'] = (f1[0],(fh,))
        c2 = r2['ceiling']
        ch = vec3(0,0,0).com(c2[0]).sq(4,4)
        r2['ceiling'] = (c2[0],(ch,))

    def connect_rooms(self,r1,r2):
        if   r1['level'] == r2['level']:self.edgehole(r1,r2)
        elif r1['level']  > r2['level']:self.roofhole(r1,r2)
        elif r1['level']  < r2['level']:self.roofhole(r2,r1)

    def armkws(self,bnd,es,lvl,ww,wh,sh,ch):
        fbnd = self.contract(bnd,ww)
        floor   = (tuple([p.cp() for p in fbnd]),())
        ceiling = (tuple([p.cp().ztrn(wh) for p in fbnd]),())
        wcnt = len(floor[0])
        rmkws = {
            'bound':bnd,'wheight':wh,'wwidth':ww,'level':lvl,
            'floor':floor,'ceiling':ceiling,'skirt':sh,'crown':ch,
            'wholes':[[] for x in range(wcnt)],
            'wholepoints':[[] for x in range(wcnt)],
            'wtypes':['e' for x in range(wcnt)],
                }

        if es is None:
            rmkws['wholes'][0].append((0,0.5,1.5,3,0))
        else:
            for rmx in es:
                self.connect_rooms(rmkws,self.rmkws[rmx])

        rmkwsx = len(self.rmkws)
        self.rmkws.append(rmkws)
        return rmkwsx

    def boundlevel(self,lvl):
        m = dmo.model()
        sgv = self.amodel(None,None,None,m,None)
        gm = m.agfxmesh()
        bnd,mbnd = [],[]
        for rmkws in self.rmkws:
            if not rmkws['level'] == lvl:continue
            for wx in range(len(rmkws['wtypes'])):
                wt,wh,ww = rmkws['wtypes'][wx-1],rmkws['wheight'],rmkws['wwidth']
                wp1,wp2 = rmkws['bound'][wx-1],rmkws['bound'][wx]
                sh,ch = rmkws['skirt'],rmkws['crown']
                hs = rmkws['wholes'][wx-1]
                if wt == 'e':
                    if not hs:hs.append((0,0.5,1.5,2,1))
                    bnd.append((wp1,wp2))
                    wtn = wp1.tov(wp2).crs(vec3(0,0,1)).nrm().uscl(ww)
                    wp1 = wp1.cp().trn(wtn)
                    wp2 = wp2.cp().trn(wtn)
                    mbnd.append((wp1,wp2))
                    m,hpts = self.awall(wp1,wp2,wh,hs,m = m,gm = gm)
                    rmkws['wholepoints'][wx-1].extend(hpts)
                    for hps in hpts:self.aportal(hps,wp1,wp2,2*ww,m = m,gm = gm)
                    sp1,sp2 = wp1.cp().ztrn(-sh),wp2.cp().ztrn(-sh)
                    cp1,cp2 = wp1.cp().ztrn(wh),wp2.cp().ztrn(wh)
                    m,hpts = self.awall(sp1,sp2,sh,None,m = m,gm = gm)
                    m,hpts = self.awall(cp1,cp2,ch,None,m = m,gm = gm)
        tbnd = [(bnd[0],mbnd[0])]
        while not len(tbnd) == len(bnd):
            lst = tbnd[-1][0][1]
            for mbe,be in zip(mbnd,bnd):
                if (be,mbe) in tbnd:continue
                if lst.isnear(be[0]):
                    tbnd.append((be,mbe))
                    lst = tbnd[-1][0][1]
        self.rims.append([])
        for x in range(len(tbnd)):
            p1,p2 = tbnd[x-1][1]
            p3,p4 = tbnd[x][1]
            if p2.isnear(p3):continue
            cpy = (
                (p2.cp().ztrn(-sh),p3.cp().ztrn(-sh),
                p3.cp().ztrn(wh+ch),p2.cp().ztrn(wh+ch)),())
            m = self.asurf(cpy,m = m,gm = gm)
            self.rims[-1].append(cpy[0][3].cp())
            self.rims[-1].append(cpy[0][2].cp())

    def boundroof(self):
        m = dmo.model()
        sgv = self.amodel(None,None,None,m,None)
        gm = m.agfxmesh()

        #ax = dtl.plot_axes()
        for lx in range(self.floors):
            lvl = self.rims[lx]
            #ax = dtl.plot_polygon(lvl,ax,lw = 2+2*lx)
            if lx == self.floors-1:
                tlvl = lvl
                print('make doors if the level below is not fully eclipsed')
            else:tlvl = pym.ebdxy(lvl,self.rims[lx+1])
            rbnd = (tuple(tlvl),())
            m = self.asurf(rbnd,m = m,gm = gm)
            for x in range(len(tlvl)):
                p1,p2 = tlvl[x-1],tlvl[x]
                self.awall(p1,p2,2.0,None,m = m,gm = gm)
            itlvl = self.contract(tlvl,0.125)
            for x in range(len(itlvl)):
                p1,p2 = itlvl[x-1],itlvl[x]
                self.awall(p2,p1,2.0,None,m = m,gm = gm)

            # must zoffset tlvl and itlvl?
            #m = self.asurf((tlvl,(itlvl,)),m = m,gm = gm)

        #plt.show()

    def rgraph(self):
        print('make a wiregraph of rm bounds, rm edges, and rm levels')
        rmvs = []

        #( boundary, edges, level, wallwidth, wallheight, skirtlength, crownheight )

        self.floors = 1

        rmv = vec3(0,0,0).sq(12,8),None,0,0.125,4.0,0.5,1.0
        rmvs.append(rmv)

        rmv = vec3(10,0,0).sq(8,8),[0],0,0.125,4.0,0.5,1.0
        rmvs.append(rmv)

        rmv = vec3(10,10,0).sq(8,12),[1],0,0.125,4.0,0.5,1.0
        rmvs.append(rmv)

        rmv = vec3(0,10,0).sq(12,12),[0,2],0,0.125,4.0,0.5,1.0
        rmvs.append(rmv)

        #rmv = vec3(0,10,5.5).sq(12,12),[3],1,0.125,4.0,0.5,1.0
        #rmvs.append(rmv)

        for rmv in rmvs:self.rtopo.append(rmv)

    # do something which fills the scenegraph
    def generate(self,worn = 0):
        print('generate with worn',worn)

        #m = dmo.model()
        #sgv = self.amodel(vec3(0,0,1),None,None,m,None)
        #gm = m.atricube('generic')

        self.rmkws,self.rtopo,self.rims = [],[],[]
        self.rgraph()
        for rtp in self.rtopo:self.armkws(*rtp)
        for lx in range(self.floors):self.boundlevel(lx)
        self.boundroof()
        for rmkws in self.rmkws:m = self.aroom(**rmkws)
        return self





