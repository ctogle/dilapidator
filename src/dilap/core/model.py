import dilap.core.base as db

import dp_vector as dpv

###############################################################################
### model is the basic unit of geometry for dilap
###
### it corresponds to a model space object
### it contains vertex data to build a model
### it can inherit transforms from a scenegraph
###
### it can write itself as a blender model
### it can write itself as a fbx model
###############################################################################

unused_model_id = 0
class model(base):

    def _dpid(self):
        global unused_model_id
        self.dpid = unused_model_id
        unused_model_id += 1

    def __init__(self,*args,**kwargs):
        self._dpid()
        self._def('pcoords',[],**kwargs)
        self._def('ncoords',[],**kwargs)
        self._def('ucoords',[],**kwargs)
        self._def('faces',[],**kwargs)

    #######################################################
    #methods for modifying the models geometry data
    #######################################################

    # given three points, add new triangle face
    def _triangle(self,v1,v2,v3,ns = None,us = None,m = None,pm = None):
        nfstart = len(self.faces)
        nps = [v1.copy(),v2.copy(),v3.copy()]
        if ns is None:
            n = normal(*nps)
            nns = [n,n,n]
        else: nns = ns
        if us is None:
            nus = [cv.vector2d(0,1),cv.vector2d(0,0),cv.vector2d(1,0)]
        else: nus = us
        self._add_vdata(nps,nns,nus)

        foffset = len(self.pcoords) - len(nps)
        nfs = mpu.offset_faces([[0,1,2]],foffset)

        m = self._lookup_mat(m)
        pm = self._lookup_pmat(pm)

        nfms = [m]
        nfpms = [pm]

        self._add_fdata(nfs,nfms,nfpms)
        nfend = len(self.faces)
        return range(nfstart,nfend)

    #######################################################

    #######################################################
    #methods for transforming the model in world space
    #######################################################

    def translate_x(self,dx):
        dpv.translate_coords_x(self.pcoords,dx)
        return self

    def translate_y(self,dy):
        dpv.translate_coords_y(self.pcoords,dy)
        return self

    def translate_z(self,dz):
        dpv.translate_coords_z(self.pcoords,dz)
        return self

    def translate_u(self,u):
        trn = dpv.vector(u,u,u)
        dpv.translate_coords(self.pcoords,trn)
        return self

    def translate(self,v):
        dpv.translate_coords(self.pcoords,v)
        return self

    #######################################################

    def scale_x(self, sx):
        cv.scale_coords_x(self.coords, sx)
        if self._scale_uvs_: self.scale_uvs(cv.vector(sx,0,0))
        self.modified = True
        return self

    def scale_y(self,sy):
        dpv.scale_coords_y(self.pcoords,sy)
        #if self._scale_uvs_: self.scale_uvs(cv.vector(0,sy,0))
        return self

    def scale_z(self,sz):
        dpv.scale_coords_z(self.pcoords,sz)
        #if self._scale_uvs_: self.scale_uvs(cv.vector(0,0,sz))
        return self

    def scale_u(self,u):
        scl = dpv.vector(u,u,u)
        dpv.scale_coords(self.pcoords,scl)
        #if self._scale_uvs_: self.scale_uvs(vect)
        return self

    def scale(self,v):
        dpv.scale_coords(self.pcoords,v)
        #if self._scale_uvs_: self.scale_uvs(vect)
        return self

    #######################################################


