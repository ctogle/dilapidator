import dilap.core.base as db
import dilap.core.tools as dpr
import dilap.core.sgraph as dsg
import dilap.core.context as dgc
import dilap.io.io as dio
import dilap.primitive.cube as dcu
import dilap.primitive.terrain as dt

import dp_vector as dpv

import random,pdb

class terrain_point:

    def __init__(self,position):
        self.position = position
        self.weights = dpv.one().scale_u(0.05)
        self.neighbors = []
        self.owners = []
        self.owner_count = 0
        self.boundary = False
        self.hole_boundary = False
        self.is_corner = False
    
    def set_neighbor_count(self):
        if self.boundary:
            if self.owner_count == 1:
                self.neighbor_count = 4
            elif self.is_corner == True:
                self.neighbor_count = 6
            else: self.neighbor_count = 6
        else: self.neighbor_count = 6

    def calculate_smooth_normal(self):
        spos = self.position
        nposs = self.neighbor_positions()
        nnorms = []
        #for ndx in range(self.neighbor_count):
        for ndx in range(len(self.neighbors)):
            n1 = self.neighbors[ndx-1].position
            n2 = self.neighbors[ndx].position

            #if near_xy(n1,n2):pdb.set_trace()

            vn1 = dpv.v1_v2(spos,n1)
            vn2 = dpv.v1_v2(spos,n2)
            nnorm = vn1.cross(vn2)
            if nnorm.z < 0.0: nnorm.flip().normalize()
            nnorms.append(nnorm)

        ncom = dpv.center_of_mass(nnorms)
        ncom.normalize()
        return ncom

    def neighbor_positions(self):
        return [n.position for n in self.neighbors]

    def reneighbors(self,pts,dthresh = 50):
        renew = []
        dt2 = dthresh**2
        spos = self.position

        target_ncnt = self.neighbor_count
        for npt in self.neighbors:
            tpos = npt.position
            ndist = dpv.distance_xy(tpos,spos)
            if ndist < dthresh:
                renew.append(npt)

        rdx = len(renew)
        if rdx < target_ncnt:
            for tp in pts:
                tpos = tp.position
                dx2 = (tpos.x - spos.x)**2
                if dx2 > dt2: continue
                dy2 = (tpos.y - spos.y)**2
                if dy2 > dt2: continue
                if tp in renew: continue
                if tp is self: continue
                if dx2 + dy2 < dt2:
                    renew.append(tp)
                    rdx += 1
                    if rdx == target_ncnt: break
        #else: print 'already reneighbored'
        self.neighbors = renew

class ttri(object):

    def _face_data(self,depth = None,max_depth = None):
        data = []
        if not depth is None:
            if depth == max_depth:
                data.append(self.local_points)
                return data
            else:depth += 1
        if self.children:
            [data.extend(ch._face_data(depth,max_depth)) 
                for ch in self.children]
        else:data.append(self.local_points)
        return data

    def _geo_data(self,lod = False):
        depth = None
        max_depth = None
        if lod:
            depth = 0
            max_depth = self.splits - 1
        data = self._face_data(depth,max_depth)
        return data

class landscape(dgc.context):

    def __init__(self,*args,**kwargs):
        dgc.context.__init__(self,*args,**kwargs)
        self._def('controls',dpr.point_ring(50,6),**kwargs)
        self._def('sealevel',-0.5,**kwargs)

    # should port algorithms from make_places for this...
    # 
    # self.controls represents fixed, in-mesh verts, forming a convex loop
    # landscape will generate a set of models with terrain that spans
    # the xy-projection of the polygon defined by self.controls
    #
    # first dice the loop so that the average distance between points
    # is at some threshold
    #
    # use the advancing front method to fill the diced loop with triangles
    # 
    # create a collection terrain_points from this process, 
    # ideally preserving the groups of 3...
    def generate_terrain_points(self,other = None,worn = 0):
        tptstack = []

        tpts = self.controls[:]
        tcom = dpv.center_of_mass(tpts)
        for x in range(len(tpts)):
            c1 = terrain_point(tcom.copy())
            c2 = terrain_point(tpts[x-1].copy())
            c3 = terrain_point(tpts[x].copy())
            tptstack.extend([c1,c2,c3])

        return tptstack

    # from a set of points that should appear in the mesh
    # create a list of lists of points corresponding to faces
    # for now assuming simplest stack of points organization
    def generate_data_from_points(self,tpts):
        tdata = [[]]
        for x in range(len(tpts))[::3]:
            face = [tpts[x],tpts[x+1],tpts[x+2]]
            tdata[-1].append(face)
        return tdata

    # from a list of terrain data lists, generate terrain models
    def generate_terrain_models(self,tdata):
        tmodels = []
        for tp in tdata:
            tmodels.append(dt.terrain(tp))
        return tmodels

    # ITS TIME TO ADD SOME SIMPLE TERRAIN, SO THAT IVY CAN BE ADDED
    def generate(self,other = None,worn = 0):
        tpts = self.generate_terrain_points(other,worn)
        terrain_data = self.generate_data_from_points(tpts)
        tmds = self.generate_terrain_models(terrain_data)

        dcube = dcu.cube().scale_x(100).scale_y(100)
        dcube.translate_z(-0.5).scale_z(20)
        dcube.translate_z(self.sealevel)
        dnode = self._node_wrap(dcube)
        tmods = self._node_wrap(*tmds)
        self._nodes_to_graph(dnode)
        self._nodes_to_graph(tmods)


