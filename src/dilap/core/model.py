import dilap.core.base as db

###############################################################################
### model is the basic unit of geometry for dilap
###
### it corresponds to a model space object
### it contains vertex data to build a model
### it can inherit transforms from a scenegraph
###
### it can write itself as a blender model
### it can write itself as a fbx model
###############################################################################

unused_model_id = 0
class model(base):

    def _dpid(self):
        global unused_model_id
        self.dpid = unused_model_id
        unused_model_id += 1

    def __init__(self,*args,**kwargs):
        self._dpid()


