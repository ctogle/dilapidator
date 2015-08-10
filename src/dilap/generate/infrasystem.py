import dilap.core.base as db
import dilap.core.context as dgc
import dilap.core.tools as dpr
import dilap.core.tmesh as dtm
import dilap.core.lsystem as dls
import dilap.core.mesh.tools as dtl

import dilap.primitive.cube as dcu
import dilap.primitive.road as dr

import dp_vector as dpv
import dp_quaternion as dpq

import pdb,numpy,math
import matplotlib.pyplot as plt

class infrastructure(dgc.context):

    def __init__(self,igraph,*args,**kwargs):
        self.igraph = igraph
        dgc.context.__init__(self,*args,**kwargs)

    def _terrain_points(self):
        self.tpts = []
        [self.tpts.extend(r.tpts) for r in self.roads]
        [self.tpts.extend(i.tpts) for i in self.intersections]
        return self.tpts

    def _hole_points(self):
        self.hpts = []
        [self.hpts.extend(r.hpts) for r in self.roads]
        [self.hpts.extend(i.hpts) for i in self.intersections]
        return self.hpts

    def generate(self,seed,boundary,worn = 0):

        #rgraph = infragraph(seed,boundary)

        ax = self.igraph.plot()
        plt.show()


        self.intersections = []
        self.roads = []
        rcs,ris = dr.graph_roads(self.igraph)
        #rcs,ris = dr.circle(dr.highway)
        
        self.intersections.extend([i.translate(seed) for i in ris])
        self.roads.extend([r.translate(seed) for r in rcs])

        #p = dpv.zero()
        #d = dpv.xhat.copy()
        #rmeshmodel = rmesh(-1)._realize(p,d).model
        #rmeshmodel.translate(seed)
        #self._nodes_to_graph(self._node_wrap(rmeshmodel))

        [self._nodes_to_graph(self._node_wrap(i)) for i in self.intersections]
        [self._nodes_to_graph(self._node_wrap(r)) for r in self.roads]
        return self



