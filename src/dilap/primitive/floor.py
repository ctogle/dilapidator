import dilap.core.model as dmo
import dilap.core.tools as dpr
import dilap.core.vector as dpv

import pdb

class floor(dmo.model):

    def __init__(self,l,w,*args,**kwargs):
        dmo.model.__init__(self,*args,**kwargs)
        self.l = l
        self.w = w
        self._def('h',0.5,**kwargs)
        self._def('gap',None,**kwargs)
        self._def('m','generic',**kwargs)
        self._geo()

    def _geo(self):
        if self.gap:self._geo_gap()
        else:self._geo_gapless()

    def _geo_gap(self):
        l,w,h,g,m = self.l,self.w,self.h,self.gap,self.m
        gp,gl,gw = g
        iloop = dpr.corners(gl,gw,gp)
        oloop = dpr.corners(l,w)
        iloop.append(iloop[0])
        oloop.append(oloop[0])
        iloopb = [c.copy() for c in iloop]
        oloopb = [c.copy() for c in oloop]
        [c.translate_z(-h) for c in iloopb]
        [c.translate_z(-h) for c in oloopb]

        nfs = self._bridge(iloop,oloop,m = m)
        self._project_uv_xy(nfs)
        nfs = self._bridge(oloopb,iloopb,m = m)
        self._project_uv_xy(nfs)
        nfs = self._bridge(iloopb,iloop,m = m)
        self._project_uv_flat(nfs)
        nfs = self._bridge(oloop,oloopb,m = m)
        self._project_uv_flat(nfs)

    def _geo_gapless(self):
        l,w,h,m = self.l,self.w,self.h,self.m
        corners = dpr.corners(l,w)
        us = dpr.polygon(4)
        dpv.translate_coords(us,dpv.one().scale_u(0.5))
        dpv.scale_coords_x(us,l)
        dpv.scale_coords_y(us,w)
        us = [v.xy2d() for v in us]
        self._quad(*corners,us = us,m = m)
        bcorners = [c.copy().translate_z(-h) for c in corners]
        bcorners.reverse()
        self._quad(*bcorners,us = us,m = m)

        bcorners.reverse()
        bcorners.append(bcorners[0])
        corners.append(corners[0])
        nfs = self._bridge(corners,bcorners,m = m)
        self._project_uv_flat(nfs)


