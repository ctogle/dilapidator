import dilap.core.tools as dpr

import dilap.geometry.tools as gtl
from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat
from dilap.geometry.pointset import pointset

from dilap.topology.trimesh import trimesh
from dilap.topology.polygonmesh import polygonmesh

import dilap.modeling.model as dmo

import pdb










__doc__ = '''terrain model factory class'''
# dilapidators implementation of a general terrain factory
class factory(dmo.model):

    def __init__(self,*args,**kwargs):
        dmo.model.__init__(self,*args,**kwargs)
        self.filename = 'terrain.model.mesh'










