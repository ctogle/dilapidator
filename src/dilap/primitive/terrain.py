import dilap.core.model as dmo
import dilap.core.tools as dpr

import dp_vector as dpv

import pdb

###############################################################################
### terrain represents a piece of deterministic terrain
### terrain is a triangle when xy-projected
###############################################################################

class terrain(dmo.model):

    def __init__(self,data,*args,**kwargs):
        dmo.model.__init__(self,*args,**kwargs)
        self._geo(data)

    # from data, generate necessary model data
    # data is a list of lists of terrain_points
    # each sublist represents one face in the mesh
    def _geo(self,data):
        for fdx,fdat in enumerate(data):
            if not fdat:continue
            fposs = [f.position.copy() for f in fdat]
            #ns = [f.calculate_smooth_normal() for f in fdat]
            m = 'generic'
            #nfs = self._triangle(*fposs,ns = ns,m = m)
            nfs = self._triangle(*fposs,m = m)
            self._project_uv_xy(nfs)


