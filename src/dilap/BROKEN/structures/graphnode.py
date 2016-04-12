import dilap.core.base as db
import dilap.core.tools as dpr
import dilap.core.vector as dpv
import dilap.core.graphnode as gnd

import dilap.mesh.tools as dtl

import matplotlib.pyplot as plt
import pdb



class node(gnd.node):

    def __init__(self,p,**kwargs):
        gnd.node.__init__(self,p,**kwargs)
        self._def('height',10.0,**kwargs)
        self._def('gap',1.0,**kwargs)










