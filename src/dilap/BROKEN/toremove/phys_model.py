import dilap.core.base as db
import dilap.core.tools as dpr
import dilap.core.model as dmo

import dp_vector as dpv
import dp_bbox as dbb


### THIS IS CURRENTLY UNUSED

#
# should be necessary at least because model normal coords are not always flat
# should be completely worthwhile because of the unreal engine
#

###############################################################################
### phys_model is the basic unit of structural geometry for dilap
###
### one phys_model corresponds to exactly one model object
### it contains colliders to accompany the model
### the colliders are convex polyhedrons, which should be 
### outputtable to the unreal engine via fbx...
###
### it can inherit transforms from a scenegraph
### its colliders inherit from it
###
### it can write itself as a blender model
### it can write itself as a obj model
### it can write itself as a fbx model
###############################################################################

unused_phys_model_id = 0
class phys_model(dmo.model):

    def _dpid(self):
        global unused_phys_model_id
        self.dpid = unused_phys_model_id
        unused_phys_model_id += 1

    def __init__(self,*args,**kwargs):
        self._dpid()



        # geometric data
        self._def('colliders',[],**kwargs)
        # non geometric data
        self._def('reps',{},**kwargs)
        self._def('filename','physmodel.mesh',**kwargs)

    # return 3d bounding box for this model
    def _aaabbb(self):
        raise NotImplemented
        #xproj = dpv.project_coords(self.pcoords,dpv.xhat)
        #yproj = dpv.project_coords(self.pcoords,dpv.yhat)
        #zproj = dpv.project_coords(self.pcoords,dpv.zhat)
        #bb = dbb.bbox(xproj,yproj,zproj)
        return bb


