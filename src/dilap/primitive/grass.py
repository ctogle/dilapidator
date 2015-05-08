import dilap.core.model as dmo
import dilap.core.tools as dpr

import dp_vector as dpv

import pdb,numpy,random

clump = 1
class grass_clump(dmo.model):

    def __init__(self,*args,**kwargs):
        dmo.model.__init__(self,*args,**kwargs)
        self._def('pieces',20,**kwargs)
        self._def('radius',1,**kwargs)
        self._def('width',0.75,**kwargs)
        self._def('height',0.5,**kwargs)
        self._geo()

    # v is zero, take a genuinely random sample
    # v is nonzero, take a safe sample self.width away from v
    def _point(self,v):
        vmag = v.magnitude()
        if vmag > 0.000000001:
            r = self.width
        else:
            mr = self.radius - self.width
            r = random.random()*mr
        t = random.random()*2*numpy.pi
        x = r*numpy.cos(t)
        y = r*numpy.sin(t)
        p = dpv.vector(x,y,0)
        p.translate(v)
        return p

    def _geo(self):
        global clump
        if clump == 1:
           clump += 1
           self._consume(grass_clump1(radius = self.radius))
        elif clump == 2:
           clump = 1
           self._consume(grass_clump1(radius = self.radius))
           #self._consume(grass_clump2(radius = self.radius))

class grass_clump1(grass_clump):

    def _geo(self):
        for pdx in range(self.pieces):
            zero = dpv.zero()
            p1 = self._point(zero)
            p2 = self._point(p1)
            ps = [p1,p2]
            t = random.random()*2*numpy.pi
            dpv.rotate_z_coords_about(ps,dpv.center_of_mass(ps),t)
            p3 = p2.copy().translate_z(self.height)
            p4 = p1.copy().translate_z(self.height)
            nfs = self._quad(p4,p1,p2,p3,m = 'grass4')
            self._scale_uv_v(nfs,0.75)
            self._translate_uv_v(nfs,0.05)

class grass_clump2(grass_clump):

    def _geo(self):
        n = 8
        bottom = dpr.point_ring(self.radius,n)
        top = dpv.center_of_mass(bottom).translate_z(self.height)
        bottom.reverse()
        self._tripie(bottom)
        bottom.append(bottom[0])
        bottom.reverse()
        self._trifan(top,bottom)


