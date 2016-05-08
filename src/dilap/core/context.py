import dilap.core.base as db
import dilap.core.qtgui as qg

from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat

#import dilap.core.vector as dpv
#import dilap.core.sgraph as dsg
import dilap.modeling.scenegraph as dsg
#import dilap.core.model as dmo
import dilap.modeling.model as dmo

import dilap.io.io as dio

#import dilap.primitive.cone as dco

import random

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
        if q is None:q = quat(0,0,0,1)
        if s is None:s = vec3(1,1,1)
        return self.sgraph.avert(p,q,s,models = [m],parent = par)

    def graph(self,io = None,**kws):
        if not io is None:self.iotype = io
        iotype = self.iotype
        if type(iotype) is type(''):
            iotype = dio.iotypes[iotype]
        self.sgraph.graph(iotype,**kws)
        return self



    # do something which fills the scenegraph
    def generate(self,worn = 0):
        print('generate with worn',worn)
        #dcone = dco.cone()
        #self._models_to_graph(dcone)
        return self

    def passtime(self,years):
        self._nodes_to_world()

        # must acquire structural data about the sgraph!!

        proxy = self._nodal_proxy()
        for d in self.dilapidors:
            d.wither(self,years,proxy)
        self._nodes_to_local()



    #### THIS CLASS MUST BE HEAVILY WORKED ON


    # for one or many nodes, return a node consuming them
    #def _node_consume(self,*nodes):
    #    #consumenode = dsg.node(children = list(nodes),consumption = True)
    #    consumenode = dsg.scenevert(children = list(nodes),consumption = True)
    #    return consumenode

    # for one or many models, return a node containing them
    # can specify node properties such as consumption via kwargs
    #def _node_wrap(self,*models,**kwargs):
    #    #nodewrap = dsg.node(models = list(models),**kwargs)
    #    nodewrap = dsg.scenevert(models = list(models),**kwargs)
    #    return nodewrap

    # for one or many models, add them to the scenegraph
    def _models_to_graph(self,*models,**kwargs):
        #self.sgraph.verts.append(self._node_wrap(*models,**kwargs))
        #self.sgraph.avert(self._node_wrap(*models,**kwargs))
        self.sgraph.avert(0,tform,*models,**kwargs)

    # for one or many nodes, add them to the scenegraph
    def _nodes_to_graph(self,*nodes,**kwargs):
        for n in nodes:self.sgraph.nodes.append(n)

    # move all nodes to world space
    def _nodes_to_world(self):
        for n in self.sgraph.nodes:n._to_space('world')

    # move all nodes to local space
    def _nodes_to_local(self):
        for n in self.sgraph.nodes:n._to_space('local')

    # combine this context with another in place
    def _consume(self,other):
        self._nodes_to_graph(*other.sgraph.nodes)

    def _nodal_proxy_single(self,proxy,node):
        for c in node.tform.children:
            self._nodal_proxy_single(proxy,c.owner)
        for m in node.models:proxy._consume_preserving(m)

    def _nodal_proxy(self):
        proxy = dmo.model()
        for n in self.sgraph.nodes:
            self._nodal_proxy_single(proxy,n)
        return proxy






