import dilap.core.context as cx
import dilap.modeling.model as dmo

from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat

import dilap.worldly.treeskin as ltr
import dilap.worldly.building as blg

import dilap.core.plotting as dtl
import matplotlib.pyplot as plt

import pdb

###############################################################################
###
###############################################################################

class world(cx.context):

    def __init__(self,*args,**kwargs):
        self._def('name','buildingcontext',**kwargs)
        cx.context.__init__(self,*args,**kwargs)

    def generate(self,worn = 0):
        print('generate world with worn',worn)

        #m = dmo.model()
        #sgv = self.amodel(vec3(0,0,10),None,None,m,None)
        #gm = m.atricube('generic')

        bfa = blg.blgfactory()

        fp = vec3(0,0,0).sq(48,64)
        ex = [vec3(0,-32,0),vec3(24,-16,0)]

        blgcx = bfa.new(fp,ex)
        blgcx.generate(worn)
        self.achild(blgcx)

        '''#
        bp,bq,bs = vec3(-20,-20,0),quat(1,0,0,0),vec3(1,1,1)
        blgcx = bfa.new(fp)
        blgcx.generate(worn,bp,bq,bs)
        self.achild(blgcx)
        '''#

        '''#
        p,d = vec3(0,0,0),vec3(0,0,1)
        ax = dtl.plot_axes()
        cx = ltr.tree(p,d,ax = ax)
        cx.generate(worn)
        self.achild(cx)
        '''#

        return self





