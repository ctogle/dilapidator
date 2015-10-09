from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat

import dilap.modeling.model as dmo

import dilap.core.tools as dpr

import pdb










__doc__ = '''dilapidator\'s implementation of a model factory'''
# dilapidators implementation of a model factory
# NOTE: most methods are meant to be overloaded by subclasses
# NOTE: this base class just makes a cube instead
class factory:

    modelclass = dmo.model

    def __str__(self):return 'cube factory:'

    def __init__(self,*args,**kwargs):
        pass

    # return a new independent model based on the 
    # current parameters of the factory
    # what does a new model require to be complete?
    #   sufficient data to build each set of trimeshes
    #   trimeshes can be made directly
    #   trimeshes can be made from other trimeshes
    #   trimeshes can be made from polygonmeshes
    #
    def new(self,*args,**kwargs):
        model = self.modelclass(*args,**kwargs)
        mpoly = model.apolymesh()

        v1 = mpoly.avert(vec3(-1,-1,-1))
        v2 = mpoly.avert(vec3( 1,-1,-1))
        v3 = mpoly.avert(vec3( 1, 1,-1))
        v4 = mpoly.avert(vec3(-1, 1,-1))
        v5 = mpoly.avert(vec3(-1,-1, 1))
        v6 = mpoly.avert(vec3( 1,-1, 1))
        v7 = mpoly.avert(vec3( 1, 1, 1))
        v8 = mpoly.avert(vec3(-1, 1, 1))
        bot = [v1,v2,v3,v4]
        top = [v5,v6,v7,v8]
        #l1 = mpoly.aloop(bot)
        #l2 = mpoly.aloop(top)
        #f1 = mpoly.aface(l1,())
        #f2 = mpoly.aface(l2,())

        return model

        mesh = mpoly
        ax = dtl.plot_axes()
        for f in mesh.faces:
            ps = mesh.gvps(f)
            ax = dtl.plot_polygon(ps,ax)
        plt.show()






 



