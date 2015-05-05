import dilap.core.base as db
import dilap.core.tools as dpr

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
        self._def('reps',{},**kwargs)
        self._def('filename','modelfile.mesh',**kwargs)

    # consume another model in place, adding all its data
    def _consume(self,other):
        ofmats = other.face_mats
        ofacnt = len(ofmats)
        for dx in range(len(other.mats)):
            omat = other.mats[dx]
            if not omat in self.mats:
                self.mats.append(omat)
                mdx = len(self.mats) - 1
                for fdx in range(ofacnt):
                    if ofmats[fdx] == dx:
                        ofmats[fdx] = mdx
            else:
                mdx = self.mats.index(omat)
                if not mdx == dx:
                    for fdx in range(ofacnt):
                        if ofmats[fdx] == dx:
                            ofmats[fdx] = mdx

        other_offset = len(self.pcoords)
        dpr.offset_faces(other.faces,other_offset)
        self.pcoords.extend(other.pcoords)
        self.ncoords.extend(other.ncoords)
        self.ucoords.extend(other.ucoords)
        self.faces.extend(other.faces)
        self.face_mats.extend(other.face_mats)
        self.reps = {}
        return self

    # return geometry data organized as dict of materials
    def _face_dict(self):
        mcnt = len(self.mats)
        fcnt = len(self.face_mats)
        fa = {}
        for mdx in range(mcnt):
            ma = self.mats[mdx]
            fa[ma] = []
            for fmdx in range(fcnt):
                if self.face_mats[fmdx] == mdx:
                    fa[ma].append(self.faces[fmdx])
            if not fa[ma]:del fa[ma]
        return fa

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

    # assign material m to range of faces rng
    def _assign_material(self,m,rng = None):
        m = self._lookup_mat(m)
        if rng is None:rng = range(len(self.faces))
        for dx in rng:
            self.face_mats[dx] = m

    #######################################################

    #######################################################
    #methods for modifying uv coordinate data
    #######################################################

    # for range of faces rng, project uvs xy
    def _project_uv_xy(self,rng):
        for nf in rng:
            face = self.faces[nf]
            for fdx in face:
                p = self.pcoords[fdx]
                nu = p.copy().xy2d()
                self.ucoords[fdx] = nu

    # for range of faces rng, project uvs yz
    def _project_uv_yz(self,rng):
        for nf in rng:
            face = self.faces[nf]
            for fdx in face:
                p = self.pcoords[fdx]
                nu = p.copy().yz2d()
                self.ucoords[fdx] = nu

    # for range of faces rng, project uvs xz
    def _project_uv_xz(self,rng):
        for nf in rng:
            face = self.faces[nf]
            for fdx in face:
                p = self.pcoords[fdx]
                nu = p.copy().xz2d()
                self.ucoords[fdx] = nu

    # for range of faces rng, project uvs flat
    def _project_uv_flat(self,rng):
        for nf in rng:
            face = self.faces[nf]
            for fdx in face:
                p = self.pcoords[fdx]
                n = self.ncoords[fdx]
                if dpv.near(n,dpv.nxhat) or dpv.near(n,dpv.xhat):
                    nu = p.copy().yz2d()
                elif dpv.near(n,dpv.nyhat) or dpv.near(n,dpv.yhat):
                    nu = p.copy().xz2d()
                elif dpv.near(n,dpv.nzhat) or dpv.near(n,dpv.zhat):
                    nu = p.copy().xy2d()
                else:continue
                self.ucoords[fdx] = nu

    # for range of faces rng, 
    # scale uvs along u coordinate by du
    def _scale_uv_u(self,rng,du):
        for nf in rng:
            face = self.faces[nf]
            for fdx in face:
                u = self.ucoords[fdx]
                u.scale_x(du)

    # for range of faces rng, 
    # scale uvs along v coordinate by dv
    def _scale_uv_v(self,rng,dv):
        for nf in rng:
            face = self.faces[nf]
            for fdx in face:
                u = self.ucoords[fdx]
                u.scale_y(dv)

    # given a tform for world space, scale uvs to local space
    def _uvs_to_local(self,uv_ttf):
        sx = uv_ttf.scl.x
        sy = uv_ttf.scl.y
        for uvc in self.ucoords:
            uvc.x *= sx
            uvc.y *= sy

    # given a tform for world space, scale uvs to world space
    def _uvs_to_world(self,uv_ttf):
        sx = uv_ttf.scl.x
        sy = uv_ttf.scl.y
        for uvc in self.ucoords:
            uvc.x /= sx
            uvc.y /= sy

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
                print('_def_normals requires 3 or 4 vertices only')
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
                print('_def_uvs requires 3 or 4 vertices only')
                raise ValueError
        else:nus = us
        return nus

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

    # given a point apex and a list of points blade, add fan of tris
    def _trifan(self,apex,blade,ns = None,us = None,m = None):
        nfstart = len(self.faces)
        tcnt = len(blade) - 1
        for trdx in range(tcnt):
            c2dx = trdx
            c3dx = trdx+1
            c1 = apex.copy()
            c2 = blade[c2dx].copy()
            c3 = blade[c3dx].copy()
            self._triangle(c1,c2,c3,ns,us,m)
        nfend = len(self.faces)
        return range(nfstart,nfend)

    # given one loop, fill with a fan of triangles
    def _tripie(self,loop,ns = None,us = None,m = None):
        nfstart = len(self.faces)
        lcom = dpv.center_of_mass(loop)
        tcnt = len(loop)
        for trdx in range(tcnt):
            c2dx = trdx
            c3dx = trdx+1
            if c3dx == tcnt: c3dx = 0
            c1 = lcom.copy()
            c2 = loop[c2dx].copy()
            c3 = loop[c3dx].copy()
            self._triangle(c1,c2,c3,ns,us,m)
        nfend = len(self.faces)
        return range(nfstart,nfend)

    # raise ValueError if l1 and l2 differ in length
    def _check_loop_equality(self,l1,l2):
        if not len(l1) == len(l2):
            print('_bridge loops must have equal length')
            raise ValueError

    # given two loops of equal length, bridge with quads
    def _bridge(self,loop1,loop2,ns = None,us = None,m = None):
        self._check_loop_equality(loop1,loop2)
        nfstart = len(self.faces)
        lcnt = len(loop1)
        for ldx in range(1,lcnt):
            v1 = loop1[ldx-1]
            v2 = loop2[ldx-1]
            v3 = loop2[ldx]
            v4 = loop1[ldx]
            self._quad(v1,v2,v3,v4,ns,us,m)
        nfend = len(self.faces)
        return range(nfstart,nfend)

    # given two loops of equal length bridge with a spline extrusion
    def _bridge_spline(self,loop1,loop2,n = 3,
            n1 = None,n2 = None,ns = None,us = None,m = None):
        self._check_loop_equality(loop1,loop2)
        nfstart = len(self.faces)

        if n1 is None:n1 = normal(*loop1[:3])
        if n2 is None:n2 = normal(*loop2[:3]).flip()

        curves = []
        lcnt = len(loop1)
        for x in range(lcnt):
            v2 = loop1[x].copy()
            v3 = loop2[x].copy()
            v1 = v2.copy().translate(n1)
            v4 = v3.copy().translate(n2)
            curve = dpv.spline(v1,v2,v3,v4,n)
            curves.append(curve)

        ccnt = len(curves)
        for y in range(1,ccnt):
            lp2 = curves[y-1]
            lp1 = curves[y]
            self._bridge(lp1,lp2,ns,us,m)
        
        nfend = len(self.faces)
        return range(nfstart,nfend)

    # given a loop, incrementally fill with loops 
    # and then seal with tripie
    def _bridge_patch(self,loop,n = 10,m = None):
        nfstart = len(self.faces)
        def move(oloop):
            iloop = [l.copy() for l in oloop]
            [l.translate(ray) for l,ray in zip(iloop,rays)]
            self._bridge(iloop,oloop,m = m)
            return iloop

        com = dpv.center_of_mass(loop)
        loop.append(loop[0])
        rays = [dpv.v1_v2(l,com).scale_u(1.0/n) for l in loop]
        for x in range(n):loop = move(loop)
        self._tripie(loop,m = m)
        nfend = len(self.faces)
        return range(nfstart,nfend)

    # given a line of points make n faces between angles a1 and a2
    def _revolve_z(self,loop,a1,a2,n,ns = None,us = None,m = None):
        nfstart = len(self.faces)
        rotstep = (a2-a1)/float(n)
        for step in range(n):
            ta1 = a1+step*rotstep
            ta2 = a1+(step+1)*rotstep
            loop1 = [p.copy().rotate_z(ta1) for p in loop]
            loop2 = [p.copy().rotate_z(ta2) for p in loop]
            self._bridge(loop1,loop2,ns = ns,us = us,m = m)
        nfend = len(self.faces)
        return range(nfstart,nfend)

    # given a curve of points make faces to extrude loop along the curve
    def _extrude(self,loop,curve,ns = None,us = None,m = None):
        nfstart = len(self.faces)
        tailloop = [l.copy() for l in loop]
        n = len(curve)
        for step in range(1,n):
            tn = dpv.v1_v2(curve[step-1],curve[step])
            tiploop = [l.copy().translate(tn) for l in tailloop]
            # i want to rotate tiploop according to tangents - need quat
            # rotation...
            self._bridge(tailloop,tiploop,ns = ns,us = us,m = m)
            tailloop = [l.copy() for l in tiploop]
        nfend = len(self.faces)
        return range(nfstart,nfend)

    # for range of faces rng, flip each face and its normals
    def _flip_faces(self,rng):
        for nf in rng:
            face = self.faces[nf]
            face.reverse()
            for fdx in face:
                self.ncoords[fdx].flip()

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

    def scale_x(self,sx):
        dpv.scale_coords_x(self.pcoords,sx)
        #if self._scale_uvs_: self.scale_uvs(dpv.vector(sx,0,0))
        return self

    def scale_y(self,sy):
        dpv.scale_coords_y(self.pcoords,sy)
        #if self._scale_uvs_: self.scale_uvs(dpv.vector(0,sy,0))
        return self

    def scale_z(self,sz):
        dpv.scale_coords_z(self.pcoords,sz)
        #if self._scale_uvs_: self.scale_uvs(dpv.vector(0,0,sz))
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

    def rotate(self,q):
        print('its time to try the quaternion')
        pdb.set_trace()

        dpv.scale_coords_x(self.pcoords,sx)
        #if self._scale_uvs_: self.scale_uvs(dpv.vector(sx,0,0))
        return self

    def rotate_x(self,rx):
        dpv.rotate_x_coords(self.pcoords,rx)
        return self

    def rotate_y(self,ry):
        dpv.rotate_y_coords(self.pcoords,ry)
        return self

    def rotate_z(self,rz):
        dpv.rotate_z_coords(self.pcoords,rz)
        return self

    #######################################################


