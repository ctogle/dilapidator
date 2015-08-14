import dilap.core.base as db
import dilap.core.vector as dpv
import dilap.core.tools as dpr
import dilap.core.sgraph as dsg
import dilap.core.context as dgc
import dilap.primitive.cube as dcu

import random as rm

import pdb

class roof(dgc.context):

    def __init__(self,floorplan,*args,**kwargs):
        dgc.context.__init__(self,*args,**kwargs)
        self.fplan = floorplan

    # generate and return roof pieces
    def generate(self,worn = 0):
        roofpieces = []
        rps,eps,ips,sps = self.fplan.allplans
        for rp in rps:

            # make a list of outter and inner points
            # walk around the outter loop making tris or quads 
            # to fill in the roof

            x,y,l,w = rp[1]['x'],rp[1]['y'],rp[1]['l'],rp[1]['w']
            z = self.fplan.bldg._roof_height()+0.5
            l += 1.0
            w += 1.0
            rpiece = dcu.cube().scale_x(l).scale_y(w)
            rpiece.translate(dpv.vector(x,y,z))

            roofpieces.append(rpiece)
        self._nodes_to_graph(self._node_wrap(*roofpieces))


