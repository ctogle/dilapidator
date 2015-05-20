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

    # if theres a coord near p in ps, return its index
    # else append and return the new index
    def search(self,ps,p1,p2):
        p = dpv.midpoint(p1,p2)
        pcnt = len(ps)
        for pdx in range(pcnt):
            if ps[pdx][0].near(p):
                return pdx
        z = (2.0*random.random()-1.0)*dpv.distance(p1,p2)/8.0
        ps.append((p,z))
        return pcnt

    def split(self,ps,ts):
        level = []
        after = []
        mpt = dpv.midpoint
        pcnt = len(ps)
        for t in ts:
            m1 = self.search(level,ps[t[0]],ps[t[1]])+pcnt
            m2 = self.search(level,ps[t[1]],ps[t[2]])+pcnt
            m3 = self.search(level,ps[t[2]],ps[t[0]])+pcnt
            after.append((t[0],m1,m3))
            after.append((m1,t[1],m2))
            after.append((m3,m2,t[2]))
            after.append((m1,m2,m3))
        for l in level:
            l[0].z += l[1]
            ps.append(l[0])
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
        pts.append(dpv.center_of_mass(pts))
        bnd = dpr.inflate([b.copy() for b in pts],-10.0)
        for p in pts:p.translate_z((2.0*random.random()-1.0)*100)
        for x in range(5):tris = self.split(pts,tris)

        m = dtm.meshme(pts,None,None,None,[],tris)

        pfaces = [(m.vs[x].p,m.vs[y].p,m.vs[z].p) for x,y,z in m.fs]
        mbb = dcu.cube().scale_u(100).scale_z(10)._aaabbb()
        #hitfaces,hitcasts = dr.ray_grid(dpv.nzhat,mbb,pfaces,1.0)

        mps = m.gpdata()
        hitfaces = dbb.intersect_tri_filter(mbb,m.fs,mps)

        hbnd = m.cut_hole(hitfaces)
        for hb in hbnd:m.vs[hb].w.scale_u(0.0)
        m.advfrontmethod(hbnd)
        vbnd = [v for v in m.vs if not dpv.inside(v.p,bnd)]
        for vb in vbnd:
            vb.w.x = 0.0
            vb.w.y = 0.0
        m.smooths(10,0.1)
        return [m.pelt()]

    def generate(self,other = None,worn = 0):
        random.seed(0)
        tmds = self.tmodels()
        # add terrain models to scenegraph
        self._models_to_graph(*tmds)
        # add water models to scenegraph
        water = dcu.cube().scale_x(100).scale_y(100).scale_z(20)
        water.translate_z(-19.5).translate_z(self.sealevel)
        #wnode = self._node_wrap(water)
        #self._nodes_to_graph(wnode)


