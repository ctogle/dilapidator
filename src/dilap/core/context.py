import dilap.core.base as db
import dilap.core.sgraph as dsg
import dilap.io.io as dio
import dilap.primitive.cube as dcu
import dilap.primitive.cone as dco
import dilap.primitive.cylinder as dcyl
import dilap.primitive.wall as dw
import dilap.primitive.floor as df

import dp_vector as dpv

import random

###############################################################################
### context corresponds to a procedurally generated world
### it includes a scenegraph
### it read/writes .dfest files
###
### it uses randomness to generate the scenegraph/models
###############################################################################

class context(db.base):

    def __init__(self,*args,**kwargs):
        self.sgraph = dsg.sgraph()
        self._def('iotype','obj',**kwargs)

    # for one or many nodes, return a node consuming them
    def _node_consume(self,*nodes):
        consumenode = dsg.node(children = list(nodes),consumption = True)
        return consumenode

    # for one or many models, return a node containing them
    # can specify node properties such as consumption via kwargs
    def _node_wrap(self,*models,**kwargs):
        nodewrap = dsg.node(models = list(models),**kwargs)
        return nodewrap

    # for one or many models, add them to the scenegraph
    def _models_to_graph(self,*models,**kwargs):
        self.sgraph.nodes.append(self._node_wrap(*models,**kwargs))

    # for one or many nodes, add them to the scenegraph
    def _nodes_to_graph(self,*nodes,**kwargs):
        for n in nodes:self.sgraph.nodes.append(n)

    # combine this context with another in place
    def _consume(self,other):
        self._nodes_to_graph(*other.sgraph.nodes)

    # do something which fills the scenegraph
    def generate(self,worn = 0):
        print('generate with worn',worn)
        
        dcone = dco.cone()
        self._models_to_graph(dcone)
        return self

    def graph(self):
        iotype = self.iotype
        if type(iotype) is type(''):
            iotype = dio.iotypes[iotype]
        self.sgraph.graph(iotype)
        return self


