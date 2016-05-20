import dilap.core.context as cx
import dilap.modeling.model as dmo

from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat

import dilap.worldly.treeskin as ltr
import dilap.worldly.building as blg

import dilap.core.plotting as dtl
import matplotlib.pyplot as plt

import random,pdb

###############################################################################
###
###############################################################################

def ushape(l,w,z):
    fp = [vec3(0,0,z),
        vec3(l/6,0,z),vec3(l/6,-w/3,z),vec3(l/3,-w/2,z),vec3(l/2,-w/2,z),
        #vec3(l/2,2*w/3,z),vec3(-l/2,2*w/3,z),vec3(-l/2,-w/2,z),vec3(-l/3,-w/2,z),
        vec3(l/2,w/2,z),vec3(-l/2,w/2,z),vec3(-l/2,-w/2,z),vec3(-l/3,-w/2,z),
        vec3(-l/6,-w/3,z),vec3(-l/6,0,z)]
    stack = [
        (0,(vec3(-5*l/6,w/6,z),vec3(5*l/6,w/6,z))),
        (1,(vec3(-l/4,-5*w/3,z),vec3(-l/4,5*w/3,z))),
        (2,(vec3(l/4,-5*w/3,z),vec3(l/4,5*w/3,z))),
        (1,(vec3(-5*l/6,0,z),vec3(5*l/6,0,z))),
        (3,(vec3(-5*l/6,0,z),vec3(5*l/6,0,z))),
        (0,(vec3(-l/24,-5*w/3,z),vec3(-l/24,5*w/3,z))),
        (6,(vec3(l/24,-5*w/3,z),vec3(l/24,5*w/3,z))),
        (0,(vec3(-l/4,-5*w/3,z),vec3(-l/4,5*w/3,z))),
        (7,(vec3(l/4,-5*w/3,z),vec3(l/4,5*w/3,z))),
            ]
    exs = [vec3(0,0,0),vec3(-l/2,w/6,0)]
    rtstack = [(1,4),(3,5),(0,8),(8,6),(6,7),(7,9)]
    #atstack = [(1,4),(3,5),(0,8),(8,6),(6,7),(7,9)]
    atstack = []
    shafts = [6]
    return fp,stack,exs,rtstack,atstack,shafts

class world(cx.context):

    def __init__(self,*args,**kwargs):
        self._def('name','buildingcontext',**kwargs)
        self._def('bfa',blg.blgfactory(),**kwargs)
        cx.context.__init__(self,*args,**kwargs)

    def generate(self,worn = 0):
        print('generate world with worn',worn)

        #m = dmo.model()
        #sgv = self.amodel(vec3(0,0,10),None,None,m,None)
        #gm = m.atricube('generic')

        p,q,s = vec3(0,0,0),quat(1,0,0,0),vec3(1,1,1)
        for m in range(1):

            lvls = random.randint(1,10)
            l = random.uniform(80.0,160.0)
            w = random.uniform(40.0,80.0)
            z = random.uniform(0,10)

            #lvls = 3
            #l,w,z = 120.0,60.0,0
            #shaft = 6

            print('lzod',l,w,z)

            fh = 7

            fp,stack,exs,rtstack,atstack,shafts = ushape(l,w,z)

            exs = []

            fps = [[p.cp().ztrn(fh*x) for p in fp] for x in range(lvls)]
            stacks = [[
                (s[0],(s[1][0].cp().ztrn(fh*x),s[1][1].cp().ztrn(fh*x)))
                    for s in stack] for x in range(lvls)]
            rtstacks = [[s[:] for s in rtstack] for x in range(lvls)]
            atstacks = [[s[:] for s in atstack] for x in range(lvls)]
            alshafts = [shafts for x in range(lvls)]

            blgcx = self.bfa.new(exs,lvls,fps,stacks,
                rtstacks,atstacks,alshafts,p,q,s)
            blgcx.generate(worn)
            self.achild(blgcx)

            p,q,s = p.cp().trn(vec3(160,0,0)),q.cp(),s.cp()

        # must solve issue of floor to floor topological connections...

        return self





