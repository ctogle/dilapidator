import dilap.core.vector as dpv
import dilap.core.quaternion as dpq
import dilap.core.tools as dpr
import dilap.core.context as dgc
import dilap.core.profiler as dprf
import dilap.mesh.piecewisecomplex as pwc
import dilap.mesh.tools as dtl
import dilap.infrastructure.graphregion as grg
import dilap.infrastructure.infragraph as ifg
import dilap.generate.landscape as dls
import dilap.generate.city as dcy
import dilap.foliage.plant as dflg

import dilap.generate.lot as dlt
import dilap.primitive.cylinder as dcyl

import matplotlib.pyplot as plt
import pdb

# an area is a region bounded by roads, 
# whose extent is defined by self.boundary
class area(dgc.context):

    def __init__(self,*args,**kwargs):
        dgc.context.__init__(self,*args,**kwargs)
        self._def('boundary',None,**kwargs)

    # populate the boundary polygon with holes for buildings and such
    def populate(self):
        def clear(b):
            b = tuple(dpr.inflate([bp.copy() for bp in b],2))
            clr = True
            '''#
            if dpr.concaves_intersect(bnd,b):return False
            if not dpr.inconcave(b[0],bnd):return False
            for h in holes:
                if dpr.concaves_intersect(h,b):return False
                if dpr.inconcave(b[0],h):return False
            for p in prints:
                if dpr.concaves_intersect(p,b):return False
                if dpr.inconcave(b[0],p):return False
            '''#
            if dpr.concaves_intersect(bnd,b):clr = False
            if not dpr.inconcave_xy(b[0],bnd):clr = False
            for h in holes:
                if dpr.concaves_intersect(h,b):clr = False
                if dpr.inconcave_xy(b[0],h):clr = False
            for p in prints:
                if dpr.concaves_intersect(p,b):clr = False
                if dpr.inconcave_xy(b[0],p):clr = False

            #if not clr:
            if False:
                ax = dtl.plot_axes_xy()
                ax = dtl.plot_polygon_xy(list(bnd),ax,lw = 1)
                ax = dtl.plot_polygon_xy(list(b),ax,lw = 5)
                for h in holes:
                    ax = dtl.plot_polygon_xy(list(h),ax,lw = 2)
                for p in prints:
                    ax = dtl.plot_polygon_xy(list(p),ax,lw = 1)
                plt.show()

            #return True
            return clr
        bnd = self.boundary[0]
        holes = self.boundary[1]

        '''#
        prints = []
        for segx in range(1,len(bnd)):
            p1,p2 = bnd[segx-1],bnd[segx]
            tn = dpv.v1_v2(p1,p2).normalize()
            nn = tn.copy().rotate_z(-dpr.PI/2.0)
            bl,bw,bh = 50,50,20
            segp = tn.copy().scale_u(bl/2.0)+nn.copy().scale_u(bw)
            p,r = p1 + segp,dpr.angle_from_xaxis_xy(tn)
            lfootprint = tuple(dpr.square(bl,bw,p,r))
            if clear(lfootprint):prints.append(lfootprint)
        '''#
        prints = [tuple(dpr.square(50,50,dpv.center_of_mass(list(bnd))))]

        self.lfootprints = prints
        #self.boundary = (bnd,holes+tuple(prints))

    def place_lots(self):
        bfootprints = []
        for prt in self.lfootprints:
            #dlot = dlt.lot(boundary = prt).generate(7)
            #for s in dlot.structures:
            #    bfootprints.append(tuple(b.copy() for b in s.boundary))
            #self._consume(dlot)
            pass
        bnd,holes = self.boundary
        self.boundary = (bnd,holes+tuple(bfootprints))
        
    def generate(self):
        self.populate()
        self.place_lots()
        tplc = pwc.piecewise_linear_complex(refine = True,smooth = False)
        tplc.add_polygons(self.boundary)

        ax = dtl.plot_axes()
        tplc.plot(ax)
        plt.show()

        dprf.profile_function(tplc.triangulate)
        tpelt = tplc.pelt()
        tnode = self._node_wrap(tpelt)
        self._nodes_to_graph(tnode)

        plantmodel = dflg.plant().model()
        plantnode = self._node_wrap(plantmodel)
        self._nodes_to_graph(plantnode)

        return self



        

