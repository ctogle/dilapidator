import dilap.core.base as db
import dilap.core.tools as dpr
import dilap.core.sgraph as dsg
import dilap.core.context as dgc
import dilap.io.io as dio
import dilap.primitive.cube as dcu
import dilap.primitive.cone as dco
import dilap.primitive.cylinder as dcyl
import dilap.primitive.wall as dw
import dilap.primitive.floor as df
import dilap.primitive.road as dr
import dilap.generate.floorplan as dfp

import dp_vector as dpv

import random

class house(dgc.context):

    def __init__(self,l,w,*args,**kwargs):
        dgc.context.__init__(self,*args,**kwargs)
        self.l = l
        self.w = w
        self._def('stories',3,**kwargs)
        self._def('fheights',[0.25 for x in range(self.stories)],**kwargs)
        self._def('cheights',[0.25 for x in range(self.stories)],**kwargs)
        self._def('wheights',[4.0 for x in range(self.stories)],**kwargs)
        self._def('stoop',(5,3,2),**kwargs)
        self.fplan = dfp.floorplan(self)

    def _terrain_points(self):
        tpts = self.fplan._terrain_points()
        return tpts

    def _story_height(self,level):
        z = 0.0
        for x in range(level):
            z += self.wheights[level]
        return z

    def _roof_height(self):
        z = self._story_height(self.stories-1)
        z += self.wheights[-1]
        return z

    def generate_stoop(self,worn = 0):
        s = self.stoop
        cub = dcu.cube()
        cub.scale_x(s[0]).scale_y(s[1]).scale_z(s[2])
        cub.translate_y(-(10+1.5)).translate_z(-s[2]/2.0)
        stoopnode = self._node_wrap(cub)
        stoopnode = self._node_consume(self._node_wrap(cub))
        self._nodes_to_graph(stoopnode)

    def generate_story(self,level,worn = 0):
        offset = dpv.zero().translate_z(self._story_height(level))
        storynodes = self.fplan.sgraph.nodes
        for n in storynodes:n.translate(offset)
        self._nodes_to_graph(*storynodes)

    def generate(self,worn = 0):
        self.fplan.plan()
        for x in range(self.stories):
            self.fplan.plan_specific(x)
            self.fplan.sgraph.nodes = []
            # need to get bboxes for the floorplan?
            self.fplan.generate(self.wheights[x],worn)
            if x == 0:self.generate_stoop(worn)
            self.generate_story(x,worn)
        for s in self.fplan.generate_shafts(worn):self._consume(s)
        self._consume(self.fplan.generate_roof(worn))
        return self


