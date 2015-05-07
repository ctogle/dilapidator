import dilap.core.dilapidor as dd
import dilap.core.model as dmo
import dilap.primitive.cylinder as dcyl
import dilap.primitive.cone as dco
import dilap.primitive.vine as dv
import dilap.primitive.grass as dg

import dp_vector as dpv
import dp_ray as dr

import pdb,numpy,random

class ivy(dd.dilapidor):

    def __init__(self,*args,**kwargs):
        dd.dilapidor.__init__(self,*args,**kwargs)
        self._def('z_max',10,**kwargs)
        self.withers.append('init')
        self.withers.append('trees')
        self.withers.append('ivy')
        self.withers.append('grass')

    # perform some calcs that are useful to more than one wither
    def init(self,model,years):
        pfaces = model._face_positions(model.faces)
        nfaces = model._face_normals(model.faces)
        mbb = model._aaabbb()
        # pfaces must not be empty...
        hitfaces,hitcasts = dr.ray_grid(dpv.nzhat,mbb,pfaces,1.0)
        for hdx in range(len(hitfaces)):
            hf = hitfaces[hdx]
            hc = hitcasts[hdx]
            v0,v1,v2 = pfaces[hf]
            u,v = hc.y,hc.z
            bcc = dpv.barymetric_to_world(u,v,v0,v1,v2)
            hitcasts[hdx] = bcc
        self.hitdata = {
            'hitfaces':hitfaces,
            'hitcasts':hitcasts,
            'pfaces':pfaces,
            'nfaces':nfaces,
                }

    # desire a model in world space of ALL nodes/children of the context
    def trees(self,model,years):
        growth = dmo.model()
        hitfaces = self.hitdata['hitfaces']
        hitcasts = self.hitdata['hitcasts']
        pfaces = self.hitdata['pfaces']
        #nfaces = self.hitdata['nfaces']
        for hdx in range(len(hitfaces)):
            if random.random() < 0.99:continue
            hf = hitfaces[hdx]
            hc = hitcasts[hdx]
            dtb = dpv.distance_to_border_xy(hc,pfaces[hf])
            rad = dtb
            if rad < 1:continue
            tree = dcyl.cylinder().translate_z(0.5).scale_z(10)
            growth._consume(tree.translate(hc))
        return growth

    # desire a model in world space of ALL nodes/children of the context
    def ivy(self,model,years):
        growth = dmo.model()
        hitfaces = self.hitdata['hitfaces']
        hitcasts = self.hitdata['hitcasts']
        pfaces = self.hitdata['pfaces']
        #nfaces = self.hitdata['nfaces']
        seeds = []
        for hdx in range(len(hitfaces)):
            if random.random() < 0.999:continue
            hf = hitfaces[hdx]
            hc = hitcasts[hdx]
            seeds.append(hc)

        for sd in seeds:
            vine = dv.vine()
            vine.grow(sd,model,years)
            #fake = dco.cone().translate_z(0.5).scale_z(2)
            #growth._consume(fake.translate(sd))
            growth._consume(vine)

        return growth

    # desire a model in world space of ALL nodes/children of the context
    def grass(self,model,years):
        growth = dmo.model()
        hitfaces = self.hitdata['hitfaces']
        hitcasts = self.hitdata['hitcasts']
        pfaces = self.hitdata['pfaces']
        #nfaces = self.hitdata['nfaces']
        for hdx in range(len(hitfaces)):
            if random.random() < 0.9:continue
            hf = hitfaces[hdx]
            hc = hitcasts[hdx]
            dtb = dpv.distance_to_border_xy(hc,pfaces[hf])
            rad = dtb
            if rad < 0.4:continue
            clump = dg.grass_clump(width = 0.75,radius = rad)
            growth._consume(clump.translate(hc))
        return growth


