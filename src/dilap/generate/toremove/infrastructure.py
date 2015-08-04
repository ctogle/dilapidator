import dilap.core.base as db
import dilap.core.context as dgc
import dilap.core.tools as dpr
import dilap.core.tmesh as dtm
import dilap.core.lsystem as dls
import dilap.core.mesh.tools as dtl

import dilap.primitive.cube as dcu
import dilap.primitive.road as dr

import dilap.generate.region as drg

import dp_vector as dpv
import dp_quaternion as dpq

import pdb,numpy,math
import matplotlib.pyplot as plt

# a region should represent a part of a continent, 
# connected to other regions by roads only!
#
# roads may only intersect the boundary of a region by covering
# it perfectly, or by intersecting it at some angle
# 
# it has a possibly concave boundary 
# a list of plugs, which are intersections on the boundary
#   which are always connected to the regions road graph
# a list of edges, which are pairs of plugs between which a road
#   should be drawn which perfectly overlaps the boundary
#
class region(dgc.context):

    def plot_xy(self,ax = None):
        dtl.plot_polygon_xy(self.bound,ax,True)
        return ax 

    def __init__(self,*args,**kwargs):
        self._def('bound',None,**kwargs)
        self._def('plugs',[],**kwargs)
        self._def('edges',[],**kwargs)
        dgc.context.__init__(self,*args,**kwargs)

        #self._define()

    def _define(self):

        pdb.set_trace()

# an infragraph lays out roads within a region
class infragraph(db.base):

    def __init__(self,sealevel,**kwargs):
        self.sealevel = sealevel

# a regiongraph lays out regions and how theyre connected
class regiongraph(db.base):

    def __init__(self,sealevel,**kwargs):
        self.sealevel = sealevel

        rbnds = self._define_region_boundaries()
        self._define_region_interiors(rbnds)
        
    def _define_region_boundaries(self):
        rpts1 = dpr.polygon(8)
        rpts2 = dpr.polygon(8)

        r = 250.0
        scale = r/(0.5/math.tan(dpr.rad(22.5)))
        dpv.scale_coords(rpts1,dpv.one().scale_u(scale))
        dpv.scale_coords(rpts2,dpv.one().scale_u(scale))
        dpv.translate_coords_x(rpts2,2.0*r)

        rbnds = [rpts1,rpts2]
        return rbnds

    def _define_region_interiors(self,rbnds):
        rpts1,rpts2 = rbnds
        rsd1 = dpv.center_of_mass(rpts1)
        rsd2 = dpv.center_of_mass(rpts2)

        rg1 = region(bound = rpts1)
        rg2 = region(bound = rpts2)

        ax = dtl.plot_axes_xy()
        rg1.plot_xy(ax)
        rg2.plot_xy(ax)
        plt.show()

        pdb.set_trace()

class infrastructure(dgc.context):

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

        rgraph = infragraph(seed,boundary)

        self.intersections = []
        self.roads = []
        rcs,ris = dr.circle(dr.highway)
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



