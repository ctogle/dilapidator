from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat
import dilap.core.base as db
import dilap.core.qtgui as qg
import dilap.modeling.scenegraph as dsg

import random,pdb


###############################################################################
### context corresponds to a procedurally generated world
### it includes a scenegraph
### it read/writes .dfest files
###
### it uses randomness to generate the scenegraph/models
###############################################################################


class context(db.base):

    def display(self,**kws):
        if qg.figure is None:qg.init_figure()
        kws['window_title'] = self.name+' display window'
        p = qg.displaycontext(self,**kws)
        if not p is None:p.join()


    def __init__(self,*args,**kwargs):
        self.sgraph = dsg.scenegraph()
        self._def('name','context',**kwargs)
        self._def('children',[],**kwargs)
        self._def('dilapidors',[],**kwargs)
        self._def('params',(),**kwargs)


    def transform(self,t,q,s):
        '''transform all nodes in the sgraph by t,q,s'''
        self.sgraph.root.scl(s).rot(q).trn(t)


    def amodel(self,p = None,q = None,s = None,m = None,par = None):
        if p is None:p = vec3(0,0,0)
        if q is None:q = quat(1,0,0,0)
        if s is None:s = vec3(1,1,1)
        return self.sgraph.avert(p,q,s,models = [m],parent = par)


    def achild(self,o):
        self.children.append(o)
