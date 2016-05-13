import dilap.core.base as db
from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat

import dilap.modeling.model as dmo

import pdb



__doc__ = '''dilapidator\'s implementation of a model factory'''
# dilapidators implementation of a simple factory
class factory(db.base):

    def __str__(self):return 'cube factory:'
    def __init__(self,**kws):
        self._def('bclass',dmo.model,**kws)
    # return a new instance of self.bclass
    def new(self,*ags,**kws):
        n = self.bclass(*ags,**kws)
        return n

 



