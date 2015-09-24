import dilap.core.tools as dpr

import dilap.geometry.tools as gtl
from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat
from dilap.geometry.pointset import pointset

from dilap.topology.trimesh import trimesh
from dilap.topology.polygonmesh import polygonmesh

import pdb










__doc__ = '''dilapidator\'s implementation of an all purpose model'''
# dilapidators implementation of an all purpose model
# a model contains geometric information and 
#   topological information which references it
# a model generates application content:
#   a set of trimeshes with position,normal,uv coordinates
#   a set of convex disjoint trimeshes composing a collision hull
#   a set of lod trimeshes
# it contains other objects which can operate on its trimeshes
class model:

    def __str__(self):return 'model:'

    '''#
    # return 3d bounding box for this model
    def aaabbb(self):
        xproj = dpv.project_coords(self.pcoords,dpv.xhat)
        yproj = dpv.project_coords(self.pcoords,dpv.yhat)
        zproj = dpv.project_coords(self.pcoords,dpv.zhat)
        bb = dbb.bbox(xproj,yproj,zproj)
        return bb
    '''#

    def __init__(self,*args,**kwargs):
        self.gfxmeshes = []
        self.colmeshes = []
        self.lodmeshes = []
        self.polymeshes = []
        self.pset = pointset()
        self.nset = pointset()
        self.uset = pointset()

        #self._def('reps',{},**kwargs)
        #self._def('filename','model.mesh',**kwargs)
        self.filename = 'model.mesh'

    # create new gfx trimesh
    def agfxmesh(self):
        ngm = trimesh()
        self.gfxmeshes.append(ngm)
        return ngm

    # create new col trimesh
    def acolmesh(self):
        ncm = trimesh()
        self.colmeshes.append(ncm)
        return ncm

    # create new lod trimesh
    def alodmesh(self):
        nlm = trimesh()
        self.lodmeshes.append(nlm)
        return nlm

    # create new polygonmesh
    def apolymesh(self):
        npm = polygonmesh()
        self.polymeshes.append(npm)
        return npm

    # given a polygonmesh, add a gfx trimesh to the 
    # model representing the polygonmesh
    def gfx(self,pmesh):
        ngm = self.agfxmesh()

        print('GEN TRIMESH FROM POLYGONMESH PLEASE!!!')
        print('GEN TRIMESH FROM POLYGONMESH PLEASE!!!')
        print('GEN TRIMESH FROM POLYGONMESH PLEASE!!!')

        return ngm

    # create new vertex tuple not present in any 
    # current mesh of the model (new p,n, and/or u)
    # NOTE: a polygonmesh might only need points, so 
    #   it should probably not require this function
    def avert(self,p = None,n = None,u = None):
        if p is None:p = vec3(0,0,0)
        if n is None:n = vec3(0,0,1)
        if u is None:u = vec3(0,0,0)
        px = self.pset.ap(p)
        nx = self.nset.ap(n)
        ux = self.uset.ap(u)
        return px,nx,ux

    # given the indices of some vertices, 
    # return their position vectors
    def gvps(self,mesh,vxs):
        pxs = (mesh.verts[vx][0] for vx in vxs)
        ps = self.pset.gps(pxs)
        return ps

    #
    # i really want two fundamental concepts abstractly added
    #   splitting operators which add vertices to the mesh
    #       and maintain topological correctness
    #   smoothing operators which converge to limit surfaces
    #
    # perform a sqrt(3) subdivision on a trimesh
    #   topological splitting of the face
    #   also geometric smooothing afterwards
    def subdiv(self,mesh):
        oldvs = list(mesh.ve_rings.keys())
        oldes = list(mesh.ef_rings.keys())
        for f in list(mesh.fs_mats.keys()):
            v,w,x = f
            mp = gtl.com(self.gvps(mesh,f))
            px,nx,ux = self.avert(mp)
            u = mesh.avert(px,nx,ux)
            mesh.sface(u,v,w,x)
        for u,v in oldes:mesh.fedge(u,v)
        dels = []
        for v in oldvs:
            p = self.pset.ps[v[0]]
            vns = mesh.mask(0,v,None,None)
            pns = self.pset.gps((v[0] for v in vns))
            alpha = mesh.alphan(len(pns))
            sdel = p.tov(gtl.com(pns)).uscl(alpha)
            dels.append((p,sdel))
        for sd in dels:sd[0].trn(sd[1])

    # translate the position pointset of the model
    def trn(self,v):
        self.pset.trn(v)
        return self

    # rotate the position and normal pointsets of the model
    def rot(self,q):
        self.pset.rot(q)
        self.nset.rot(q)
        return self

    # scale the position pointset of the model
    def scl(self,s):
        self.pset.mul(s)




 




