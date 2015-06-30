import dilap.core.base as db
import dilap.core.tools as dpr
import dilap.core.sgraph as dsg
import dilap.core.context as dgc
import dilap.core.mesh as dms
import dilap.io.io as dio
import dilap.primitive.cube as dcu
import dilap.primitive.terrain as dt

import dilap.core.tmesh as dtm
import dilap.core.mesh.piecewisecomplex as pwc
import dilap.core.mesh.tools as dtl

import dp_vector as dpv
import dp_ray as dr
import dp_bbox as dbb

import random,pdb
import matplotlib.pyplot as plt

class landscape(dgc.context):

    def __init__(self,*args,**kwargs):
        dgc.context.__init__(self,*args,**kwargs)
        self._def('controls',[],**kwargs)
        self._def('holes',[],**kwargs)
        self._def('regions',[],**kwargs)
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
            #z = pz+(rz-pz)*(controld/thresh2)**2
            z = pz+(rz-pz)*(controld/thresh2)**3
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

    def _cover(self,radius = 100,tri_edgelength = 10,mod_edgelength = 250):
        convexcover = dpr.pts_to_convex_xy(self.regions)
        dpr.inflate(convexcover,radius)
        pts,tris = dpr.triangle_cover(convexcover,mod_edgelength)
        wts = [dpv.one() for x in pts]
        self.center = dpv.center_of_mass(convexcover)
        return pts,wts,tris,convexcover

    def tmodels(self,splits = 3):
        pts,wts,tris,convex = self._cover()
        maxh = 200.0
        dtbs = [dpv.distance_to_border_xy(p,convex) for p in pts]
        mdtb = max(dtbs)
        deltas = [(0.0 + maxh*(dtbs[x]/mdtb)**3) 
            if dpv.inside(pts[x],convex) else -20.0 for x in range(len(dtbs))]
        for x in range(len(pts)):pts[x].translate_z(deltas[x])
        for x in range(splits):tris = self.split(pts,wts,tris)
        m = dtm.meshme(pts,None,None,wts,[],tris)
        vbnd = [x for x in range(len(m.vs)) if len(m.vs[x].vring) < 6]
        #dpr.plot_points([m.vs[x].p for x in vbnd])
        for vbx in vbnd:
            m.vs[vbx].w.scale_u(0.0)
            m.vs[vbx].p.z = -20.0
        return m

    def generate(self,other = None,worn = 0):
        random.seed(0)
        m = self.tmodels()
        m.smooths(100,0.1,method = 'uniform')

        flatholes = dpr.flatten(self.holes)
        hbb = dbb.bb_from_ps(flatholes)
        hbb._consume_x(dpv.vector2d(hbb.x.x-5.0,hbb.x.y+5.0))
        hbb._consume_y(dpv.vector2d(hbb.y.x-5.0,hbb.y.y+5.0))
        hbb._consume_z(dpv.vector2d(-1000.0,1000.0))
        mps = m.gpdata()
        mfs = m.intersect_aaabbb(hbb)

        ax = dtl.plot_axes_xy()
        cutrng = []
        for hdx in range(len(self.holes)):
            hole = self.holes[hdx]

            dtl.plot_polygon_xy(hole,ax)
        plt.show()

        '''#
        cutrng = []
        for mfx in range(len(mfs)):
            mfps = m.gfpdat(mfs[mfx])
            for hdx in range(len(self.holes)):
                hole = self.holes[hdx]
                isect = 1 - dpv.separating_axis(mfps,hole)
                if isect:
                    if not mfx in cutrng:
                        cutrng.append(mfs[mfx])
        '''#

        hbnd = m.cut_hole(cutrng)
        for h in hbnd:m.vs[h].w.scale_u(0.0)

        #bns = [dpv.zhat.copy() for x in flatholes]
        #bus = [dpv.zero2d() for x in flatholes]
        #bws = [dpv.one() for x in flatholes]
        #fvs = m.newvdata(flatholes,bns,bus,bws)

        holecover = pwc.model_plc(points = [m.vs[h].p for h in hbnd])
        #dpr.plot_points([m.vs[h].p for h in hbnd],edges = False)

        #for hx in hbnd:
        #    m.vs[hx].p.z = 100
        #    for hxx in m.vs[hx].vring:m.vs[hxx].p.z = 100

        #loops = m.looprank(hbnd)
        #dvs = loops[0][0] + fvs
        #dvs = hbnd[:]

        #patch = m.delaunaymethod(dvs)
        #flatpatch = m.flatten(patch,dpv.vector(0,0,100),dpv.zhat.copy())

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

        self.landmd = m.pelt()
        self.landmd._consume(holecover)
        self.landbb = self.landmd._aaabbb()
        self._models_to_graph(self.landmd)
        return self


