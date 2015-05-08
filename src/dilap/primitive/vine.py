import dilap.core.model as dmo
import dilap.core.tools as dpr

import dilap.primitive.cube as dcu

import dp_vector as dpv
import dp_quaternion as dq

import pdb,random,numpy

class vine(dmo.model):

    def __init__(self,*args,**kwargs):
        dmo.model.__init__(self,*args,**kwargs)

    # growth begins at position seed
    # the models geometry should determine the vine path
    # the number of growth steps taken depends on years
    def grow(self,seed,pfaces,nfaces,years):
        gface,seedposition = seed
        gpt = seedposition.copy()
        pactive = pfaces[gface]
        tactive = dpr.tangent(*pactive)
        nactive = dpr.normal(*pactive)
        #nactive = nfaces[seedface]
        ctrlpts = [gpt]
        ctrlpts.append(gpt.copy().translate_z(0.25))
        pi2 = numpy.pi/2.0
        while years:
            years -= 1
            r = dpr.clamp(random.random(),0.5,1.0)
            t = random.random()*numpy.pi-pi2
            q = dq.quaternion(t,*nactive)
            #dgpt = tactive.copy().rotate(q).scale_u(r)
            dgpt = q.rotate_vector(tactive.copy()).scale_u(r)
            gpt = ctrlpts[-1].copy().translate(dgpt)
            ctrlpts.append(gpt)
            tactive = dpv.v1_v2(ctrlpts[-2],ctrlpts[-1]).normalize()
            # consider if the active face should change

        splined = []
        splined.append(ctrlpts[0])
        for x in range(3,len(ctrlpts)):
            v1,v2,v3,v4 = ctrlpts[x-3],ctrlpts[x-2],ctrlpts[x-1],ctrlpts[x]
            splined.extend(dpv.vector_spline(v1,v2,v3,v4,5))
        #splined.append(ctrlpts[-1])

        loop = dpr.point_ring(0.05,8)
        loop.append(loop[0].copy())
        dpv.translate_coords(loop,splined[0])

        self._extrude(loop,splined)
        #self._extrude(loop,ctrlpts)
        #for cp in ctrlpts:
        for cp in splined:
            leaf = dcu.cube().translate_z(0.5)
            leaf.scale_u(0.1).translate(cp)
            #self._consume(leaf)

        #    #pactive = pfaces[gface]
        #    #tactive = dpr.tangent(*pactive)
        #    #nactive = dpr.normal(*pactive)

        return self


