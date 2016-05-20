import dilap.core.base as db
import dilap.core.qtgui as qg

from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat

import dilap.modeling.scenegraph as dsg
import dilap.modeling.model as dmo

import dilap.io.io as dio

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
        #self.sgraph = dsg.sgraph()
        self.sgraph = dsg.scenegraph()
        self._def('name','context',**kwargs)
        self._def('children',[],**kwargs)
        self._def('iotype','obj',**kwargs)
        self._def('dilapidors',[],**kwargs)
        self._def('params',(),**kwargs)

    # transform all nodes in the sgraph by t,q,s
    def transform(self,t,q,s):
        self.sgraph.root.scl(s).rot(q).trn(t)
        #for n in self.sgraph.verts:
        #    n.scl(s).rot(q).trn(t)

    def amodel(self,p = None,q = None,s = None,m = None,par = None):
        if p is None:p = vec3(0,0,0)
        if q is None:q = quat(1,0,0,0)
        if s is None:s = vec3(1,1,1)
        return self.sgraph.avert(p,q,s,models = [m],parent = par)

    def achild(self,o):
        self.children.append(o)

    def graph(self,io = None,**kws):
        if not io is None:self.iotype = io
        iotype = self.iotype
        if type(iotype) is type(''):
            iotype = dio.iotypes[iotype]
        for ch in self.children:ch.graph(io,**kws)
        self.sgraph.graph(iotype,**kws)
        return self

    # do something which fills the scenegraph
    def generate(self,worn = 0):
        print('fill the scenegraph with models')
        print('generate with worn',worn)
        return self





