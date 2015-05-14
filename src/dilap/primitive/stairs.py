import dilap.core.model as dmo
import dilap.core.tools as dpr

import dp_vector as dpv

import pdb,math

class stairs(dmo.model):

    def __init__(self,*args,**kwargs):
        dmo.model.__init__(self,*args,**kwargs)
        self._def('steps',8,**kwargs)
        self._def('l',10,**kwargs)
        self._def('w',4,**kwargs)
        self._def('h',8,**kwargs)
        self._geo()

    def _profile(self,stepheight,steplength):
        line = []
        p = dpv.zero()
        for sx in range(self.steps):
            line.append(p.copy())
            p.translate_z(stepheight)
            line.append(p.copy())
            p.translate_y(steplength)
        line.append(p.copy())
        return line

    def _geo_from_profile(self,line,l,w,h,steps,stepheight,steplength):
        topleft = [pt.copy().translate_x(-w/2.0) for pt in line] 
        topright = [pt.copy().translate_x(w/2.0) for pt in line] 
        bottom = dpr.point_line(dpv.zero(),dpv.vector(0,l,h),steps)
        for bdx in range(steps):
            bottom.insert(2*bdx+1,bottom[2*bdx+1].copy())
        dpv.translate_coords_z(bottom[1:],-stepheight)
        bottomleft = [pt.copy().translate_x(-w/2.0) for pt in bottom] 
        bottomright = [pt.copy().translate_x(w/2.0) for pt in bottom] 
        nfs = []
        nfs.extend(self._bridge(topleft,topright))
        nfs.extend(self._bridge(bottomleft,topleft))
        nfs.extend(self._bridge(topright,bottomright))
        nfs.extend(self._bridge(bottomright,bottomleft))
        nfs.extend(self._quad(
            topleft[-1],topright[-1],
            bottomright[-1],bottomleft[-1]))
        #self._project_uv_flat(nfs)

    def _geo(self):
        l,w,h = float(self.l),float(self.w),float(self.h)
        steps = self.steps
        stepheight = h/steps
        steplength = l/steps
        line = self._profile(stepheight,steplength)
        self._geo_from_profile(line,l,w,h,steps,stepheight,steplength)


