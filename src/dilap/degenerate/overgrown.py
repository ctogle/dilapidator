import dilap.core.dilapidor as dd
import dilap.core.model as dmo
import dilap.primitive.vine as dv

import dp_vector as dpv
import dp_ray as dr

import pdb,numpy,random

class ivy(dd.dilapidor):

    def __init__(self,*args,**kwargs):
        dd.dilapidor.__init__(self,*args,**kwargs)
        self._def('z_max',10,**kwargs)
        self.withers.append('ivy')

    # desire a model in world space of ALL nodes/children of the context
    def ivy(self,model,years):
        growth = dmo.model()

        mfaces = model._face_positions(model.faces)
        mbb = model._aaabbb()

        dx = 1
        xax = numpy.arange(mbb.x.x,mbb.x.y,dx)
        yax = numpy.arange(mbb.y.x,mbb.y.y,dx)
        z = mbb.z.y + 1
        pts = [(x,y) for x in xax for y in yax]
        raygrid = [dr.ray(dpv.vector(x,y,z),dpv.nzhat) for x,y in pts]

        hitfaces = []
        hitcasts = []
        if mfaces:
            for zray in raygrid:
                if random.random() < 0.9:continue
                hf,hc = dr.intersect_hits_closest(zray,mfaces)
                #if not hf == -1 and not hf in hitfaces:
                if not hf == -1:
                    hitfaces.append(hf)
                    hitcasts.append(hc)
        
        for hdx in range(len(hitfaces)):
            hf = hitfaces[hdx]
            hc = hitcasts[hdx]

            v0,v1,v2 = mfaces[hf]
            v0,v1,v2 = v0.copy(),v1.copy(),v2.copy()
            u,v = hc.y,hc.z
            bcc = v0.scale_u(1-u-v)+v1.scale_u(u)+v2.scale_u(v)
            ivy = dv.vine().grow(years).translate(bcc)

            growth._consume(ivy)

        return growth


