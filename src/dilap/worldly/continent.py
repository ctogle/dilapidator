from dilap.core import *
from dilap.core.model import model
import dilap.worldly.terrain as ter
import dilap.worldly.infrastructure as inf
import dilap.worldly.partitiongraph as ptg
from dilap.geometry import *
import dilap.geometry.polymath as pym


class continent(scenegraph):


    def generate(self):
        for v,depth in self.a.enum():
            if 'world' in v.style:pass
            elif 'ocean' in v.style:
                self.ocean(v,self.a)
            elif 'land' in v.style:
                pass
                #self.natural(w,v,self.a)
            elif 'bounded' in v.style:
                pass
                #self.developed(w,v,self.a)
            elif 'infrastructure' in v.style:
                self.infrastructure(v,self.t,self.a)
        for v,depth in self.a.enum():
            if 'land' == v.style or 'bounded' == v.style:
                print('generate terrain for style: \'%s\'' % v.style)
                self.terrain(v,self.t,self.r,self.a)
        return self


    def __init__(self,b = None,e = None,**kws):
        scenegraph.__init__(self)
        if b is None:b = vec3(0,0,0).pring(1000,8)
        if e is None:
            xpj = vec3(1,0,0).prjps(b)
            e = (xpj[1]-xpj[0])/1000.0
        self.b = b
        self.e = e
        self.t = ter.terrain(b,e,**kws)
        self.r = inf.roadmap(self.t,e,**kws)
        self.a = self.partition(self.t,self.r,e)
        self.generate()


    @staticmethod
    def partition(terr,roads,e):
        '''partition the space spanned by a topography
         - the boundary to the first set of topo loops is the ocean
         - the first set of topo loops are the landmasses
         - partition each landmass based upon infrastructure bounds
        '''
        areas = ptg.partition(
            loop = terr.natural.root.loop,holes = [],style = 'world')
        blgfps = roads.footprints
        rpy = pym.pgtopy(roads.roads,5)

        #ax = dtl.plot_axes(200)
        #ax = dtl.plot_polygon(rpy[0],ax,lw = 3,col = 'b')
        #ax = dtl.plot_points(rpy[0],ax)
        #plt.show()

        rpy = None
        for lmtv in terr.natural.below(terr.natural.root):
            rootpvs,lmpv = areas.split(areas.root,lmtv.loop,None,style = 'land')
            if not rpy is None:
                lmpvs,rdpv = areas.split(lmpv,rpy[0],style = 'infrastructure')
                for fp in blgfps[0]:
                    fnd = False
                    for lmpv in lmpvs:
                        if pym.binbxy(fp[0],lmpv.loop):
                            #blgpvs,blgpv = areas.split(lmpv,fp[0],style = 'structure')
                            #found = True
                            #break
                            pass
                    if fnd:break
                if len(rpy[1]) > 0:
                    irdpvs,irdpv = areas.split(rdpv,rpy[1][0],style = 'bounded')
                if len(rpy[1]) > 1:
                    for irpy in rpy[1][1:]:
                        irdpvs,irdpv = areas.split(irdpvs[0],irpy,style = 'bounded')
                rtv = terr.natural.avert(
                    terr.natural.locloop(rpy[0],1)[0],loop = rpy[0])
                for dpy in rpy[1]:
                    rtv = terr.natural.avert(
                        terr.natural.locloop(dpy,1)[0],loop = dpy)

        print('... the continent ...')
        return areas


    def ocean(self,v,pg,depth = 2,bleed = 5):
        print('ocean vertex',v[0])
        ### create an unrefined flat surface for the bottom of the ocean
        ###  and an unrefined flat surface with contracted holes for the 
        ###  surface of the ocean
        m = model()
        sgv = self.avert(None,None,None,models = [m],parent = self.root)
        tm = m.agfxmesh()
        gb = v[1]['b']
        #lhs = [pym.contract([p.cp().ztrn(depth) for p in ib],bleed,5) for ib in gb[1]]
        lhs = [pym.contract([p.cp() for p in ib],bleed,5) for ib in gb[1]]
        gb = [[p.cp() for p in gb[0]],lhs]
        wb = [[p.cp().ztrn(depth) for p in gb[0]],
            [[p.cp().ztrn(depth) for p in lh] for lh in lhs]]
        ngvs = m.asurf(gb,tm,fm = 'concrete1',ref = False,hmin = 100,zfunc = None)
        ngvs = m.asurf(wb,tm,fm = 'concrete1',ref = False,hmin = 100,zfunc = None)


    def infrastructure(self,v,t,pg):
        print('infrastructure vertex',v.style,len(v.loop),len(v.holes))
        m = model()
        sgv = self.avert(None,None,None,models = [m],parent = self.root)
        tm = m.agfxmesh()

        rh = 0.0
        vb = [b.cp().ztrn(rh) for b in v.loop]
        vibs = [[b.cp().ztrn(rh) for b in h] for h in v.holes]
        vb = self.vstitch(vb,v,pg)
        vibs = [self.vstitch(vib,v,pg) for vib in v.holes]
        v.loop = vb
        v.holes = vibs

        ngvs = m.asurf((vb,vibs),tm,
            fm = 'concrete1',ref = True,hmin = 100,zfunc = t,
            uvstacked = None,autoconnect = True)
        lockf = lambda p : p.onpxy((v.loop,v.holes)) 
        m.subdiv(tm,False,True,lockf)

        return m

        rh = 0.0
        vb = [b.cp().ztrn(rh) for b in v[1]['b'][0]]
        vibs = [[b.cp().ztrn(rh) for b in v[1]['b'][1][x]] 
            for x in range(len(v[1]['b'][1]))]
        vb = self.vstitch(vb,v,pg)
        vibs = [self.vstitch(vib,v,pg) for vib in v[1]['b'][1]]
        v[1]['b'] = vb,vibs

        ngvs = m.asurf((vb,vibs),tm,
            fm = 'concrete1',ref = True,hmin = 100,zfunc = v[1]['tmesh'],
            uvstacked = None,autoconnect = True)
        lockf = lambda p : p.onpxy(v[1]['b']) 
        m.subdiv(tm,False,True,lockf)
        #m.subdiv(tm,False,True,lockf)
        #m.subdiv(tm,False,True,lockf)
        #m.uvs(tm)
        print('generated infrastructure')
        return m


    def terrain(self,v,t,r,a):
        m = model()
        sgv = self.avert(None,None,None,models = [m],parent = self.root)
        tm = m.agfxmesh()
        #vstitch isnt being used!!!
        # need to enforce that refinement is no problem for stitching!!!
        #  need to stitch from holes to interior regions!!!

        #vb = [b.cp() for b in v[1]['b'][0]]
        vb = [b.cp() for b in v.loop]
        vibs = [[b.cp() for b in h] for h in v.holes]
        #vibs = []
        vb = self.stitchareas(vb,v,a)
        vibs = [self.stitchareas(h,v,a) for h in v.holes]
        #v[1]['b'] = vb,vibs
        v.loop = vb
        v.holes = vibs

        print('generating terrain')
        ngvs = m.asurf((vb,vibs),tm,
            fm = 'grass2',ref = True,hmin = 16,zfunc = t,
            rv = pym.bnrm(vb).z < 0,uvstacked = None,autoconnect = True)
        lockf = lambda p : p.onpxy((vb,vibs)) 
        m.subdiv(tm,False,True,lockf)
        m.uvs(tm)
        print('generated terrain')
        return m


    @staticmethod
    def stitchareas(vb,v,a,l = 10):
        #vb = [b.cp() for b in v[1]['b'][0]]
        if not pym.bccw(vb):vb.reverse()
        fnd = True
        while fnd:
            fnd = False
            for o,depth in a.enum():
                if o.ix == v.ix:continue
                #if ox == v.ix:continue
                ob = o.loop
                adjs = pym.badjbxy(vb,ob,0.1)
                for adj in adjs:
                    vbx,obx = adj
                    vb1,vb2,ob1,ob2 = vb[vbx-1],vb[vbx],ob[obx-1],ob[obx]
                    ips = pym.sintsxyp(vb1,vb2,ob1,ob2,ieb = 0,skew = 0,ie = 0)
                    if ips is None:continue
                    ip1,ip2 = ips
                    if not ip1 in vb or not ip2 in vb:
                        if ip1.onsxy(vb1,vb2,0):vb.insert(vbx,ip1);fnd = True
                        if ip2.onsxy(vb1,vb2,0):vb.insert(vbx,ip2);fnd = True
                if fnd:break
        vbx = 0
        while vbx < len(vb):
            vb1,vb2 = vb[vbx-1],vb[vbx]
            if vb1.d(vb2) < l:vbx += 1
            else:vb.insert(vbx,vb1.mid(vb2))
        return vb
