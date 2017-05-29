import dilap.core.base as db
import dilap.core.lsystem as lsy
import dilap.geometry.tools as gtl
from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat
import dilap.geometry.polymath as pym

import dilap.topology.tree as dtr
import dilap.topology.planargraph as pgr

import dilap.worldly.polygen as pyg

import dilap.core.plotting as dtl
import matplotlib.pyplot as plt

import math,numpy,random,pdb



class roadmap(db.base):

    @staticmethod
    def layroads(t,e):
        '''generate a planargraph of a network of roads'''

        #return pgr.planargraph()

        drho,i = 20,3
        p,d = vec3(0,0,0),vec3(0,1,0)
        #axiom = 'X'
        #rules = dict([('X','F{[X}{]X}X'),('F','FF')])
        #params = dict(dazimuthal = gtl.rad(90),drho = drho)
        axiom = 'X'
        rules = dict([('X','F]{{X}[X}[F{[FX}]X'),('F','FF')])
        params = dict(dazimuthal = gtl.rad(90),drho = drho)
        pg = pgr.planargraph()
        for piece in lsy.lgen(p,d,axiom,rules,i,**params):
            if isinstance(piece,tuple):
                p1,p2 = piece
                v1,v2 = pg.fp(p1,drho/2),pg.fp(p2,drho/2)
                rdkws = {'style':'urban'}
                e12 = pg.fe(v1,v2,**rdkws)
            elif isinstance(piece,vec3):pass

        #for v in pg.vs:
        #    if v is None:continue
        #    p = v[1]['p']
        #    p.ztrn(t(p.x,p.y))

        #pyscale = vec3(0,1,0).prjps(py)
        #tipscale = vec3(0,1,0).prjps(tip.loop)
        #pyscale = pyscale[1]-pyscale[0]
        #tipscale = tipscale[1]-tipscale[0]
        #scale = tipscale/pyscale
        #com = vec3(0,0,0).com(tip.loop).tov(vec3(0,0,0).com(py)).uscl(-1)
        #for p in py:p.scl(vec3(scale,scale,0)).trn(com)

        print('ROADS')
        ax = dtl.plot_axes_xy(200)
        ax = pg.plotxy(ax)
        plt.show()
        #pdb.set_trace()
        return pg

    @staticmethod
    def layfootprints(r,t,e):
        eseam,seams,loops = pym.pgtopy(r,5,e,findeseam = True)
        ebnd = seams[eseam]
        eloop = loops[eseam]
        k = 0
        easement = None
        for j in range(len(ebnd)):
            road = r.es[r.elook[(eloop[k-1],eloop[k])]]
            edge = (ebnd[j-1],ebnd[j])
            tang = edge[0].tov(edge[1])
            norm = tang.crs(vec3(0,0,1)).nrm()
            p1 = edge[0].cp().trn(tang.cp().uscl(0.0))
            p2 = edge[0].cp().trn(tang.cp().uscl(1.0))
            p3 = p2.cp().trn(norm.uscl(10))
            p4 = p1.cp().trn(norm)
            newfp = [p1,p2,p3,p4]
            newfp = pym.contract(newfp,2)
            newfp = pym.ebdxy(newfp,ebnd)[0]
            newfp = (newfp,[])
            preasement = easement
            if easement is None:easement = newfp
            else:easement = pym.epuxy(easement,newfp)[0]
            if not eloop[k] == eloop[k-2]:k += 1

            if False and j > 85:
                print('walk',j,len(eloop),len(ebnd))
                ax = dtl.plot_axes_xy(300)
                ax = dtl.plot_polygon_xy(ebnd,ax,lw = 3,col = 'b')
                ax = dtl.plot_points_xy(ebnd,ax,number = True)
                ax = dtl.plot_edges_xy(edge,ax,lw = 7,col = 'r')
                ax = dtl.plot_polygon_full_xy(easement,ax,lw = 5,col = 'g')
                if preasement:
                    ax = dtl.plot_polygon_full_xy(preasement,ax,lw = 3,col = 'm')
                plt.show()

        #easement = pym.bsuxy(easement,e)
        fps = ([easement],[[] for l in range(len(loops)-1)])
        # generate a set of footprints connected to road
        # footprint must not overlap interiors with ebnd
        ax = dtl.plot_axes_xy(300)
        ax = dtl.plot_polygon_xy(ebnd,ax,lw = 3,col = 'b')
        for fp in fps[0]:
            ax = dtl.plot_polygon_full_xy(fp,ax,lw = 3,col = 'g')
        plt.show()

        return fps

    def __init__(self,t,e,**kws):
        self.terrain = t
        self.e = e
        self.roads = self.layroads(t,e)
        self.bounds = pym.pgtopy(self.roads,5)
        self.footprints = self.layfootprints(self.roads,t,e)
        self.layers = {
            'roads':self.roads,
            'footprints':self.footprints,
                }

    def __call__(self,x,y):
        '''identify the height of terrain if x,y overlap any infrastructure layers'''
        pdb.set_trace()


