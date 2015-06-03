import dilap.core.base as db
import dilap.core.tools as dpr
import dilap.core.sgraph as dsg
import dilap.core.context as dgc
import dilap.core.mesh as dms
import dilap.io.io as dio
import dilap.primitive.cube as dcu
import dilap.primitive.terrain as dt

import dilap.core.tmesh as dtm

import dp_vector as dpv
import dp_ray as dr
import dp_bbox as dbb

import random,pdb

class landscape(dgc.context):

    def __init__(self,*args,**kwargs):
        dgc.context.__init__(self,*args,**kwargs)
        self._def('boundary',dpr.point_ring(250,6),**kwargs)
        self._def('controls',[],**kwargs)
        self._def('sealevel',-0.5,**kwargs)

    def search_offset_random(self,p1,p2):
        return (2.0*random.random()-1.0)*dpv.distance(p1,p2)/8.0

    def search_offset(self,p,p1,p2,w,thresh1 = 20,thresh2 = 100):
        ccnt = len(self.controls)
        if ccnt < 1:return self.search_offset_random(p1,p2)
        controlx = dpv.find_closest_xy(p,self.controls,ccnt,5.0)
        controld = dpv.distance_xy(p,self.controls[controlx])
        if controld < thresh1:
            z = self.controls[controlx].z-p.z
            w.scale_u(0.0)
            #self.controls.pop(controlx)
        elif controld < thresh2:
            rz = self.search_offset_random(p1,p2)
            pz = self.controls[controlx].z-p.z
            z = pz+(rz-pz)*(controld/thresh2)**2
        else:z = self.search_offset_random(p1,p2)
        return z

    # if theres a coord near p in ps, return its index
    # else append and return the new index
    def search(self,ps,ps0,p1,p2,w):
        p = dpv.midpoint(p1,p2)
        pfnd = p.nearest(ps0)
        if not pfnd == -1:return pfnd
        z = self.search_offset(p,p1,p2,w)
        pcnt = len(ps)
        ps.append((p,z,w))
        ps0.append(p)
        return pcnt

    def split(self,ps,ws,ts):
        level0 = []
        level = []
        after = []
        mpt = dpv.midpoint
        pcnt = len(ps)
        for t in ts:
            if ws[t[0]] == ws[t[1]]:w1 = ws[t[0]].copy()
            else:w1 = dpv.one()
            if ws[t[1]] == ws[t[2]]:w2 = ws[t[1]].copy()
            else:w2 = dpv.one()
            if ws[t[2]] == ws[t[0]]:w3 = ws[t[2]].copy()
            else:w3 = dpv.one()
            m1 = self.search(level,level0,ps[t[0]],ps[t[1]],w1)+pcnt
            m2 = self.search(level,level0,ps[t[1]],ps[t[2]],w2)+pcnt
            m3 = self.search(level,level0,ps[t[2]],ps[t[0]],w3)+pcnt
            after.append((t[0],m1,m3))
            after.append((m1,t[1],m2))
            after.append((m3,m2,t[2]))
            after.append((m1,m2,m3))
        for l in level:
            l[0].z += l[1]
            ps.append(l[0])
            ws.append(l[2])
        return after

    def tmodels(self):
        tris = []
        for x in range(6):
            c1 = 6
            c2 = x
            c3 = x+1
            if c3 == 6:c3 = 0
            tris.append((c1,c2,c3))
        pts = [b.copy() for b in self.boundary]
        wts = [dpv.zero() for p in pts]
        pts.append(dpv.center_of_mass(pts))
        wts.append(dpv.one())
        bnd = dpr.inflate([b.copy() for b in pts],-10.0)

        for p in pts:p.translate_z((2.0*random.random()-1.0)*100)
        for x in range(5):tris = self.split(pts,wts,tris)

        m = dtm.meshme(pts,None,None,wts,[],tris)

        '''#
        pfaces = [(m.vs[x].p,m.vs[y].p,m.vs[z].p) for x,y,z in m.fs]
        mbb = dcu.cube().scale_u(100).scale_z(10)._aaabbb()
        mps = m.gpdata()
        hitfaces = dbb.intersect_tri_filter(mbb,m.fs,mps)
        hbnd = m.order_loop(m.cut_hole(hitfaces))
        patch = m.advfrontmethod(hbnd)
        flatpatch = m.flatten(patch,dpv.vector(0,0,50),dpv.zhat.copy())
        '''#

        #vbnd = [v for v in m.vs if not dpv.inside(v.p,bnd)]
        #for vb in vbnd:
        #    vb.w.x = 0.0
        #    vb.w.y = 0.0
        return m

    def generate(self,other = None,worn = 0):
        random.seed(0)
        m = self.tmodels()
        
        '''#
        pfaces = [(m.vs[x].p,m.vs[y].p,m.vs[z].p) for x,y,z in m.fs]
        mbb = dcu.cube().scale_u(100).scale_z(10)._aaabbb()
        mps = m.gpdata()
        hitfaces = dbb.intersect_tri_filter(mbb,m.fs,mps)
        #hbnd = m.order_loop(m.cut_hole(hitfaces))
        hbnd = m.cut_hole(hitfaces)
        hbps = [m.vs[hx].p for hx in hbnd]

        v = dpv.vector(0,0,50)
        bps = dpv.translate_coords(dpr.dice_edges(dpr.corners(80,60),3),v)
        #bps.reverse()
        bns = [dpv.zhat.copy() for x in bps]
        bus = [dpv.zero2d() for x in bps]
        bws = [dpv.one() for x in bps]

        bvs = m.newvdata(bps,bns,bus,bws)
        for bvx in range(len(bvs)):m.vs[bvx].w.scale_u(0.0)
        bvsloop = [(bvs[x-1],bvs[x]) for x in range(len(bvs))]
        for e in bvsloop:m.vs[e[0]].edge(e[1])

        #patch = m.advfrontmethod(hbnd)
        #patch = m.advfrontmethod(bvs)
        #patch = m.afm_bridge(hbnd,bvs)

        #patch = m.delaunay_bridge(hbnd,bvs)
        #patch = m.delaunaymethod(hbnd)
        patch = m.delaunaymethod(bvs)

        flatpatch = m.flatten(patch,dpv.vector(0,0,50),dpv.zhat.copy())
        '''#

        m.smooths(50,0.1,method = 'uniform')
                                                                            
        # add terrain models to scenegraph
        self._models_to_graph(m.pelt())
        # add water models to scenegraph
        water = dcu.cube().scale_x(100).scale_y(100).scale_z(20)
        water.translate_z(-19.5).translate_z(self.sealevel)
        #wnode = self._node_wrap(water)
        #self._nodes_to_graph(wnode)
        return self


