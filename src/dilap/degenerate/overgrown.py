import dilap.core.dilapidor as dd
import dilap.core.model as dmo
import dilap.primitive.vine as dv

import dp_vector as dpv

import pdb

class ivy(dd.dilapidor):

    def __init__(self,*args,**kwargs):
        self._def('z_max',10,**kwargs)

    def grow(self,model,years):
        growth = dmo.model()
        for f in model.faces:
            ns = [model.ncoords[fx] for fx in f]
            ncom = dpv.center_of_mass(ns).normalize()
            pitched = dpv.distance(ncom,dpv.zhat) > 0.5
            if not pitched:
                ps = [model.pcoords[fx] for fx in f]
                pcom = dpv.center_of_mass(ps)
                vpts = [dpv.midpoint(pcom,v) for v in ps]
                for p in vpts:
                    ivy = dv.vine().grow(years).translate(p)
                    growth._consume(ivy)
                print('grow ivy for',years,'years at',pcom)
        model._consume(growth)
        return model

    def wither_node(self,node,years):
        for tf in node.tform.children:self.wither_node(tf.owner,years)
        for m in node.models:self.grow(m,years)

    def wither(self,context,years):
        for n in context.sgraph.nodes:
            self.wither_node(n,years)


