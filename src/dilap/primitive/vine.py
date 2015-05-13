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
        ctrlpts = [gpt.copy()]
        ctrlpts.append(ctrlpts[-1].copy().translate_z(0.25))
        pi2 = numpy.pi/2.0
        while years:
            years -= 1

            r = 1.0
            t = -pi2/3.0 if years % 2 == 0 else pi2/3.0

            q = dq.q_from_av(t,nactive)
            #dgpt is where to go from the last point
            # it must be aware of structures to creep on
            dgpt = tactive.copy().rotate(q).scale_u(r)

            ctrlpts.append(ctrlpts[-1].copy().translate(dgpt))
            tactive = dpv.v1_v2(ctrlpts[-2],ctrlpts[-1]).normalize()
            # consider if the active face should change

        ctrlpts.append(ctrlpts[-1].copy().translate_z(-0.25))

        #splined = []
        #splined.append(ctrlpts[0])
        #for x in range(3,len(ctrlpts)):
        #    v1,v2,v3,v4 = ctrlpts[x-3],ctrlpts[x-2],ctrlpts[x-1],ctrlpts[x]
        #    splined.extend(dpv.vector_spline(v1,v2,v3,v4,5))
        ##splined.append(ctrlpts[-1])
        #splined = ctrlpts[:]

        loop = dpr.point_ring(0.1,8)
        ctrl = dpv.zero()

        loop.append(loop[0].copy())
        nfs = self._extrude(loop,ctrlpts,dpv.zero())
        return self


