import dilap.core.base as db
import dilap.core.model as dmo
import dilap.core.tools as dpr

import dilap.primitive.cube as dcu

import dp_vector as dpv
import dp_quaternion as dpq

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import pdb,numpy,random

class segment(dmo.model):

    def _calculate(self):
        t = cv.v1_v2(self.start,self.end).normalize()
        l = cv.distance(self.start,self.end)
        n = t.copy().xy().rotate_z(dpr.rad(90)).normalize()
        sn = n.copy().rotate_z(-1*self.a1).normalize()
        en = n.copy().rotate_z(   self.a2).normalize()
        self.t = t
        self.l = l
        self.n = n
        self.sn = sn
        self.en = en

    def __init__(self,start,end,**kwargs):
        dmo.model.__init__(self,*args,**kwargs)
        self.start = start
        self.end = end
        self._def('a1',0.0,**kwargs)
        self._def('a2',0.0,**kwargs)
        self._def('lanecount',2,**kwargs)

        self._calculate()
        self._geo()

    def _geo(self):
        lc,lw,m = self.lanecount,self.lw,'generic'
        lanes = [x-(lc-1)/2.0 for x in range(lc)]
        for lan in lanes:
            lanw = lan*lw
            self._build_lane(lanw,lw,project = True,m = m)

    def _geo_lane(self):
        pass
    def _build_lane(self,l,lw,lh = 0,rh = 0,
            project = False,m = 'gridmat'):
        lt1 = self.sn.copy().scale_u(l)
        lt2 = self.en.copy().scale_u(l)
        tr = self.n.copy().scale_u(lw/2.0)
        v1 = self.start.copy().translate(lt1)
        v2 = v1.copy()
        v3 = self.end.copy().translate(lt2)
        v4 = v3.copy()
        v1.translate(tr).translate_z(lh)
        v2.translate(tr.flip()).translate_z(rh)
        v3.translate(tr).translate_z(rh)
        v4.translate(tr.flip()).translate_z(lh)
        mbp.rotate_pair([v1,v2],-1*self.a1)
        mbp.rotate_pair([v3,v4],self.a2)
        nfs = self._quad(v1,v2,v3,v4,m = m)
        if project:self._project_uv_xy(nfs)
        return nfs

class lane(db.base):

    def __init__(self,*args,**kwargs):
        self._def('road',None,**kwargs)
        self._def('w',5.0,**kwargs)
        # alignment designates the number of lanes between
        # this lane and center of the road to which it belongs
        self._def('alignment',0,**kwargs)
        # direction can be -1,1 or 0
        self._def('direction',1,**kwargs)

        # sequence is a an encoding of the topology
        self._def('sequence','s',**kwargs)

    def _tolane(self,ptx):
        nsoffs = [self.road.lanes[lx].w for lx in self.lnstocenter]
        offdst = self.direction*(self.w/2.0)+sum(nsoffs)
        offset = self.road.normals[ptx].copy().scale_u(offdst)
        pt = self.road.vertices[ptx].copy().translate(offset)
        return pt

    def _fromlane(self,ptx):
        return pt
        
    def _lanestocenter(self):
        lanestocenter = []
        for lx in range(len(self.road.lanes)):
            l = self.road.lanes[lx]
            if not l.direction == self.direction:continue
            if not l.alignment < self.alignment:continue
            lanestocenter.append(lx)
        return lanestocenter

    def _geo_pair(self,x):
        inframe = self._tolane(x)
        norm = self.road.normals[x]
        offset = norm.copy().scale_u(self.w/2.0)
        pt1 = inframe.copy().translate(offset)
        pt2 = inframe.copy().translate(offset.flip())
        return pt1,pt2

    # return a model containing this lane
    def _geo(self):
        self.lnstocenter = self._lanestocenter()
        #m = dtm.meshme(pts,None,None,None,[],tris)
        m = dmo.model()
        
        lpt1,lpt2 = self._geo_pair(0)
        m._consume(dcu.cube().translate(lpt1).translate_z(0.5))
        m._consume(dcu.cube().translate(lpt2).translate_z(0.5))

        for pdx in range(1,len(self.road.vertices)):
            pt1,pt2 = self._geo_pair(pdx)

            m._quad(pt1,pt2,lpt2,lpt1)

            lpt1,lpt2 = pt1,pt2

            m._consume(dcu.cube().translate(pt1).translate_z(0.5))
            m._consume(dcu.cube().translate(pt2).translate_z(0.5))
        return m
        
# a road is a topological structure where lanes (vertices) are connected
# by sharing edges and merging/splitting

class road(dmo.model):

    def calculate_tips(self,controls):
        eleng = 3.0
        start_tip = self.start.copy().translate(
                self.tail.copy().scale_u(eleng))
        end_tip = self.end.copy().translate(
            self.tip.copy().flip().scale_u(eleng))
        self.controls.insert(1,end_tip)
        for cont in controls:self.controls.insert(1,cont)
        self.controls.insert(1,start_tip)

    def calculate_count(self,v1,v2):
        ds = dpv.distance(v1,v2)
        seglen = 3
        self.lastsegcnt = int(ds/seglen)
        return self.lastsegcnt

    def calculate_tangents(self):
        verts = self.vertices
        tangs = []
        for sgdx in range(1,len(verts)):            
            p1,p2 = verts[sgdx-1],verts[sgdx]
            tangs.append(dpv.v1_v2(p1,p2).normalize())
        tangs.append(self.tip.copy())
        self.tangents = tangs

    def calculate_normals(self):
        thats = self.tangents
        rnms = [that.cross(dpv.zhat).normalize() for that in thats]
        self.normals = rnms

    def calculate_vertices(self):
        verts = [self.start.copy()]
        for dx in range(len(self.controls)-3):
            if len(verts) > 1:verts.pop(-1)
            v1 = self.controls[dx]
            v2 = self.controls[dx+1]
            v3 = self.controls[dx+2]
            v4 = self.controls[dx+3]
            scnt = self.calculate_count(v2,v3)
            verts.extend(dpv.vector_spline(v1,v2,v3,v4,scnt))
        verts.append(self.end.copy())
        self.vertices = verts
        self.calculate_tangents()
        self.calculate_normals()

        #self.set_angles(self.tangents)
        #self.set_corners()
        #self.total_length = self.set_arc_length()
        #self.set_safe_vertices()

    def calculate(self,controls = []):
        self.tipnormal = self.tip.copy().flip().xy()
        self.tipnormal.rotate_z(dpr.rad(-90)).normalize()
        self.tailnormal = self.tail.copy().xy()
        self.tailnormal.rotate_z(dpr.rad(90)).normalize()

        self.controls = [self.start.copy(),self.end.copy()]
        self.calculate_tips(controls)
        self.calculate_vertices()

        self.plot_vertices()
        #self.xy_bbox()                                

    def plot_vertices(self):
        xs = [v.x for v in self.vertices]
        ys = [v.y for v in self.vertices]
        zs = [v.z for v in self.vertices]
        fig = plt.figure()
        ax = fig.add_subplot(111,projection='3d')
        ax.plot(xs,ys,zs = zs,marker = 'o')
        plt.show()

    def _terrain_points(self):
        return [v.copy() for v in self.vertices]

    def _lotspace(self,bbs):

        rdtangs = self.tangents
        rdnorms = self.normals
        #rdwidth = self.total_width
        rdvts = self.vertices
        segcnt = len(rdvts) - 1

        rdist = -1.0*sum([l.w for l in self._rightside()])
        which = random.randint(0,segcnt)
        rvt,rnm = rdvts[which],rdnorms[which].copy()
        lotp = rvt.copy().translate(rnm.scale_u(rdist))
        lotzrot = dpv.angle_from_xaxis_xy(rdtangs[which])

        #####
        #lotzrot = 0.0 
        #####

        l,w,p,q = 50,50,lotp,dpq.q_from_av(lotzrot,dpv.zhat)
        return l,w,p,q

    def __init__(self,start,end,tip,tail,**kwargs):
        dmo.model.__init__(self,**kwargs)
        self.start = start
        self.end = end
        self.tip = tip
        self.tail = tail
        self.calculate()
        self._lanes()
        self._geo()

    def _leftside(self):
        return self._side(-1)

    def _rightside(self):
        return self._side(1)

    def _side(self,side):
        lanes = []
        for l in self.lanes:
            if l.direction == side:
                lanes.append(l)
        return lanes

    # create a collection of lane objects representing how this road
    # should evolve (merges, shoulders, gutters, sidewalks, etc...)
    def _lanes(self):
        same = {'road':self}
        self.lanes = [
            lane(direction = -1,alignment = 0,**same),
            lane(direction =  1,alignment = 0,**same),
            lane(direction =  1,alignment = 1,**same)]

    # use the lane objects to build a model
    def _geo(self):
        m = dmo.model()
        for l in self.lanes:m._consume(l._geo())
        self._consume(m)


