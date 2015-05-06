import dilap.core.base as db
import dilap.core.model as dmo

###############################################################################
### dilapidor modifies a context to apply some time-dep effect
###
### this should include things like ivy overgrowing sturctures
### rain wearing away stone or land
### decay and collapse of non-concrete structures
### sinking and deforming of structures and their supports
###
### it should effectively simulate the passage of time on a context
###############################################################################

class dilapidor(db.base):

    def __init__(self,*args,**kwargs):
        self._def('pernode',False,**kwargs)
        self._def('withers',[],**kwargs)

    # desire a model in world space of ALL nodes/children of the context
    def wither_model(self,model,years):
        change = dmo.model()
        for w in self.withers:
            change._consume(self.__getattribute__(w)(model,years))
        return change

    # recursively wither nodes and their children instead of using the proxy
    def wither_node(self,node,years):
        for tf in node.tform.children:self.wither_node(tf.owner,years)
        for m in node.models:self.wither_model(m,years)

    # strength of effect should somehow be proportional to years
    # effects are applied to the nodes in the scenegraph of context
    def wither(self,context,years,proxy = None):
        if self.pernode:
            for n in context.sgraph.nodes:
                self.wither_node(n,years)

        withering = self.wither_model(proxy,years)
        context._nodes_to_graph(context._node_wrap(withering))
        print('dilapidor withered',years,'years on context',context)


