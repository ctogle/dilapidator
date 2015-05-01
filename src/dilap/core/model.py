import dilap.core.base as db
import dilap.primitive.tools as dpr

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
class model(db.base):

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
        self._def('face_mats',[],**kwargs)
        self._def('mats',['generic'],**kwargs)

    #######################################################
    #methods for modifying the models material data
    #######################################################

    # provide the index of a material
    def _lookup_mat(self,m):
        if m is None:m = 0
        else:
            if m in self.mats:m = self.mats.index(m)
            else:
                self.mats.append(m)
                m = len(self.mats) - 1
        return m

    #######################################################

    #######################################################
    #methods for modifying the models geometry data
    #######################################################

    # add vertex data given coords,normals,uvs
    def _add_vdata(self,ps,ns,us):
        self.pcoords.extend(ps)
        self.ncoords.extend(ns)
        self.ucoords.extend(us)

    # add face data given face indices,materials
    def _add_fdata(self,fs,fms):
        self.faces.extend(fs)
        self.face_mats.extend(fms)

    # given 3 verts(vs) and the passed in normals(ns)
    # return a list of certainly acceptable normals
    def _def_normals(self,vs,ns):
        if ns is None:
            if len(vs) == 3:
                n = dpr.normal(*vs)
                nns = [n,n,n]
            elif len(vs) == 4:
                n = dpr.normal(*vs[:-1])
                nns = [n,n,n,n]
            else:
                print '_def_normals requires 3 or 4 vertices only'
                raise ValueError
        else:nns = ns
        return nns

    # given 3 verts(vs) and the passed in uvs(us)
    # return a list of certainly acceptable uvs
    def _def_uvs(self,vs,us):
        if us is None:
            if len(vs) == 3:
                nus = [dpv.vector2d(0,1),dpv.vector2d(0,0),dpv.vector2d(1,0)]
            elif len(vs) == 4:
                nus = [dpv.vector2d(0,1),dpv.vector2d(0,0),
                      dpv.vector2d(1,0),dpv.vector2d(1,1)]
            else:
                print '_def_uvs requires 3 or 4 vertices only'
                raise ValueError
        else:nus = us
        return nus

    # given three points, add new triangle face
    def _triangle(self,v1,v2,v3,ns = None,us = None,m = None):
        nfstart = len(self.faces)
        nps = [v1.copy(),v2.copy(),v3.copy()]
        nns = self._def_normals(nps,ns)
        nus = self._def_uvs(nps,us)
        self._add_vdata(nps,nns,nus)
        foffset = len(self.pcoords) - len(nps)
        nfs = [[foffset,foffset+1,foffset+2]]
        m = self._lookup_mat(m)
        nfms = [m]
        self._add_fdata(nfs,nfms)
        nfend = len(self.faces)
        return range(nfstart,nfend)

    # given four points, add two new triangle faces
    def _quad(self,v1,v2,v3,v4,ns = None,us = None,m = None):
        nfstart = len(self.faces)
        vs = [v1,v2,v3,v4]
        nns = self._def_normals(vs,ns)
        nus = self._def_uvs(vs,us)
        us1 = [nus[0],nus[1],nus[2]]
        us2 = [nus[0],nus[2],nus[3]]
        ns1 = [nns[0],nns[1],nns[2]]
        ns2 = [nns[0],nns[2],nns[3]]

        self._triangle(v1,v2,v3,ns = ns1,us = us1,m = m)
        self._triangle(v1,v3,v4,ns = ns2,us = us2,m = m)
        nfend = len(self.faces)
        return range(nfstart,nfend)

    #######################################################

    #######################################################
    #methods for transforming the model in world space
    #######################################################

    def center(self):
        com = dpv.center_of_mass(self.pcoords)
        self.translate(com.flip())
        return self

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


