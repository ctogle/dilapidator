import dilap.core.base as db
import dilap.core.tools as dpr

import dilap.core.model as dmo
#import dilap.primitive.terrain as dt

import dp_vector as dpv
import dp_quaternion as dpq
import dp_bbox as dbb
import dp_ray as dr

import numpy,random,pdb

# meshdata is a collection of vectors
# one meshdata object is consumed by meshes
class meshdata(db.base):

    # add coordinate data, to be consumed by vertices
    def _add_data(self,ps = None,ns = None,us = None):
        if not ps is None:
            self.pcoords.extend(ps)
            self.pcnt = len(self.pcoords)
        if not ns is None:
            self.ncoords.extend(ns)
            self.ncnt = len(self.ncoords)
        if not us is None:
            self.ucoords.extend(us)
            self.ucnt = len(self.ucoords)

    # provide the index of a material, given its name
    def _lookup_mat(self,m):
        if m is None:m = 0
        else:
            if m in self.materials:
                m = self.materials.index(m)
            else:
                self.materials.append(m)
                m = len(self.materials) - 1
        return m

    def __init__(self,*args,**kwargs):
        # data consumed by topology
        self._def('pcoords',[],**kwargs)  # list of position coord vectors
        self._def('ncoords',[],**kwargs)  # list of normal coord vectors
        self._def('ucoords',[],**kwargs)  # list of uv coord vectors
        self.pcnt = len(self.pcoords)
        self.ncnt = len(self.ncoords)
        self.ucnt = len(self.ucoords)

        self._def('meshes',[],**kwargs)   # list of mesh objects
        self._def('vertices',[],**kwargs) # list of vertex tuples
        self.mcnt = len(self.meshes)
        self.vcnt = len(self.vertices)

        self._def('materials',['generic'],**kwargs)

    # add a new vertex and return its index
    # while incrementing the total count of vertices
    def _add_vertex(self,v):
        newv = self.vcnt
        self.vertices.append(v)
        self.vcnt += 1
        return newv

    # return a sharp vertex, whose data is not shared
    def _coord_sharp(self,p,n,u):
        self.pcoords.append(p.copy())
        self.ncoords.append(n.copy())
        self.ucoords.append(u.copy())
        sharp = (self.pcnt,self.ncnt,self.ucnt)
        self.pcnt += 1
        self.ncnt += 1
        self.ucnt += 1
        return self._add_vertex(sharp)

    # look in coords for a vector near to c
    # return its index, or None if none is found
    # coord must also be present in sg
    def _coord_search(self,sg,coords,c):
        for cdx in range(len(coords)):
            if cdx in sg:
                e = coords[cdx]
                if e.near(c):return cdx

    # return a smooth vertex, whose data may be shared
    def _coord_smooth(self,sg,p,n,u):
        psg = [self.vertices[x][0] for x in sg]
        pindex = self._coord_search(psg,self.pcoords,p)
        if pindex is None:return self._coord_sharp(p,n,u)
        else:
            nsg = [self.vertices[x][1] for x in sg 
                if self.vertices[x][0] == pindex]
            nindex = self._coord_search(nsg,self.ncoords,n)
            if nindex is None:
                self.ncoords.append(n.copy())
                self.ucoords.append(u.copy())
                self.ncnt += 1
                self.ucnt += 1
                return self._add_vertex((pindex,self.ncnt-1,self.ucnt-1))
            else:
                self.ucoords.append(u.copy())
                self.ucnt += 1
                return self._add_vertex((pindex,nindex,self.ucnt-1))

    # given a smoothing group,position,normal,uv vector return vertex tuple
    # this must somehow know how to deal with data duplication
    # to achieve proper topological info and shared coords
    # return an index within self.vertices
    def _vertex(self,sg,p,n,u,smooth = True):
        if smooth:v = self._coord_smooth(sg,p,n,u)
        else:v = self._coord_sharp(p,n,u)
        sg.append(v)
        return v

    # create and return new mesh object
    def _mesh(self):
        new = mesh(self)
        self.meshes.append(new)
        return new

    # return kwargs to create models from meshes [this will go away...]
    def _modeldata(self):
        mdatas = []
        for m in self.meshes:
            mdatas.append(m._modeldata())
        return mdatas

# a mesh is a topological object, referencing meshdata indices
unused_mesh_id = 0
class mesh(db.base):

    def _dpid(self):
        global unused_mesh_id
        self.dpid = unused_mesh_id
        unused_mesh_id += 1

    #######################################################
    #methods for modifying the models material data
    #######################################################

    # assign material m to range of faces rng
    def _assign_material(self,m,rng = None):
        m = self.meshdata._lookup_mat(m)
        if rng is None:rng = range(len(self.faces))
        for dx in rng:self.face_mats[dx] = m

    #######################################################

    #######################################################
    #methods for modifying uv coordinate data
    #######################################################

    # for range of faces rng, project uvs flat
    # if a vertex is used more than once, it will get
    #   the projected uv from the last face which uses 
    #   it and which has a flat normal near the x,y,z axis
    def _project_uv_flat(self,rng = None):
        if rng is None:rng = range(len(self.faces))
        pcs = self.meshdata.pcoords
        ucs = self.meshdata.ucoords
        #vcs = self.vertices
        vcs = self.meshdata.vertices
        for nf in rng:
            face = self.faces[nf]
            #v1,v2,v3 = face
            #p1,p2,p3 = pcs[vcs[v1][0]],pcs[vcs[v2][0]],pcs[vcs[v3][0]]
            ps = [pcs[vcs[vx][0]] for vx in face]
            #n = dpr.normal(p1,p2,p3)
            n = dpr.normal(ps[0],ps[1],ps[2])
            if   dpv.near(n,dpv.nxhat) or dpv.near(n,dpv.xhat):natt = 'yz2d'
            elif dpv.near(n,dpv.nyhat) or dpv.near(n,dpv.yhat):natt = 'xz2d'
            elif dpv.near(n,dpv.nzhat) or dpv.near(n,dpv.zhat):natt = 'xy2d'
            else:continue
            for v in face:
                px,nx,ux = vcs[v]
                ucs[ux] = pcs[px].__getattribute__(natt)()
                print('projected uv coord',ucs[ux])
        return self

    #######################################################

    def __init__(self,meshdata,*args,**kwargs):
        self._dpid()
        self.meshdata = meshdata

        # smoothing groups are smoothing groups:
        # lists of verts which can share data with meshdata
        self._def('groups',[[]],**kwargs)
        self._def('group',0,**kwargs)

        # indices in the vertex list on meshdata
        self._def('vertices',[],**kwargs)
        # indices of vertices in self.vertices
        self._def('faces',[],**kwargs)
        # 1-1 with faces, specifying a material by index
        self._def('face_mats',[],**kwargs)
        self.vcnt = len(self.vertices)
        self.fcnt = len(self.faces)

        # topological data
        self._def('tverts',[],**kwargs)
        self._def('twghts',[],**kwargs)
        self._def('tfaces',[],**kwargs)
        self._def('vrings',[],**kwargs)
        self._def('frings',[],**kwargs)
        self.tvcnt = len(self.tverts)
        self.tfcnt = len(self.tfaces)

    # given a vertex, 
    # create/store topovertex, 
    # return topovertex index
    def _topovertex(self,v):
        vx = self.meshdata.vertices[v]
        for tdx in range(self.tvcnt):
            tv = self.tverts[tdx]
            vrep = self.meshdata.vertices[tv[0]]
            if vrep[0] == vx[0] and vrep[1] == vx[1]:
                self.tverts[tdx].append(v)
                return tdx
        self.tverts.append([v])
        self.twghts.append(dpv.one().scale_u(0.1))
        self.vrings.append([])
        self.frings.append([])
        newtv = self.tvcnt
        self.tvcnt += 1
        return newtv

    # create topology of an edge between v1,v2
    def _topoedge(self,v1,v2):
        if not v2 in self.vrings[v1]:
            self.vrings[v1].append(v2)
            self.vrings[v2].append(v1)

    # compute topology of this tface, return its index
    def _topoface(self,vs,f):
        if len(f) == 3:return self._topoface_triangle(vs,f)
        if len(f) == 4:return self._topoface_quad(vs,f)

    # compute topology of this triangular tface, return its index
    def _topoface_triangle(self,vs,f):
        t1,t2,t3 = (vs[f[0]],vs[f[1]],vs[f[2]])
        self.tfaces.append((t1,t2,t3))
        self._topoedge(t1,t2)
        self._topoedge(t1,t3)
        self._topoedge(t2,t3)
        newtf = self.tfcnt
        if not newtf in self.frings[t1]:self.frings[t1].append(newtf)
        if not newtf in self.frings[t2]:self.frings[t2].append(newtf)
        if not newtf in self.frings[t3]:self.frings[t3].append(newtf)
        self.tfcnt += 1
        return newtf

    # compute topology of this triangular tface, return its index
    def _topoface_quad(self,vs,f):
        t1,t2,t3,t4 = (vs[f[0]],vs[f[1]],vs[f[2]],vs[f[3]])
        self.tfaces.append((t1,t2,t3,t4))
        self._topoedge(t1,t2)
        self._topoedge(t2,t3)
        self._topoedge(t3,t4)
        self._topoedge(t4,t1)
        newtf = self.tfcnt
        if not newtf in self.frings[t1]:self.frings[t1].append(newtf)
        if not newtf in self.frings[t2]:self.frings[t2].append(newtf)
        if not newtf in self.frings[t3]:self.frings[t3].append(newtf)
        if not newtf in self.frings[t4]:self.frings[t4].append(newtf)
        self.tfcnt += 1
        return newtf

    # erase topology of an edge between v1,v2
    # create a new point between them, and two new edges
    def _topoedge_bisect(self,v1,v2):
        mdvs = self.meshdata.vertices
        mdps = self.meshdata.pcoords
        vcs = self.vertices
        vts = self.tverts
        if v2 in self.vrings[v1]:
            print('bisection in ring')
            self.vrings[v1].remove(v2)
            self.vrings[v2].remove(v1)
        if True:#THIS LINE IS TO BE REMOVE....
            p1 = mdps[mdvs[vcs[vts[v1][0]]][0]]
            p2 = mdps[mdvs[vcs[vts[v2][0]]][0]]
            #deltam = dpv.vector(0,0,random.random()*10)
            deltam = dpv.zero()
            m = dpv.midpoint(p1,p2).translate(deltam)
        else:
            sharedv = None
            for vr1 in self.vrings[v1]:
                if vr1 in self.vrings[v2]:
                    sharedv = vr1
                    break
            if not sharedv:pdb.set_trace()
            m = mdps[mdvs[vcs[vts[sharedv][0]]][0]]
        return m

    # swap all quads for pairs of triangles
    def _triangulate(self):
        mdvs = self.meshdata.vertices
        mdps = self.meshdata.pcoords
        vcs = self.vertices
        vts = self.tverts
        toremove = []
        for tfx in range(self.tfcnt):
            tf = self.tfaces[tfx]
            if len(tf) == 3:continue
            toremove.append(tfx)
            ps = [mdps[mdvs[vcs[vts[tv][0]]][0]] for tv in tf]
            p1,p2,p3,p4 = ps
            self._triangle(p1,p2,p3)
            self._triangle(p1,p3,p4)
        self._rem_faces(toremove)
        return self

    # dice each triangle in a topologically intelligent way
    # for each tface, replace with 4 new tfaces
    # the vertices just need to be from the same smoothing group
    #   to inherit data correctly
    def _decimate(self):
        mdvs = self.meshdata.vertices
        mdps = self.meshdata.pcoords
        vcs = self.vertices
        vts = self.tverts
        toremove = []
        for tfx in range(self.tfcnt):
            tf = self.tfaces[tfx]
            toremove.append(tfx)
            vs = [vts[tv][0] for tv in tf]
            ps = [mdps[mdvs[vcs[vt]][0]] for vt in vs]
            if len(ps) == 3:
                v1,v2,v3 = tf
                p1,p2,p3 = ps
                m1 = self._topoedge_bisect(v1,v2)
                m2 = self._topoedge_bisect(v2,v3)
                m3 = self._topoedge_bisect(v3,v1)
                self._triangle(p1,m1,m3)
                self._triangle(p2,m2,m1)
                self._triangle(p3,m3,m2)
                self._triangle(m1,m2,m3)
            elif len(ps) == 4:
                v1,v2,v3,v4 = tf
                p1,p2,p3,p4 = ps
                m1 = self._topoedge_bisect(v1,v2)
                m2 = self._topoedge_bisect(v2,v3)
                m3 = self._topoedge_bisect(v3,v4)
                m4 = self._topoedge_bisect(v4,v1)
                cp = dpv.center_of_mass(ps)
                self._quad(p1,m1,cp,m4)
                self._quad(m1,p2,m2,cp)
                self._quad(cp,m2,p3,m3)
                self._quad(m4,cp,m3,p4)
        self._rem_faces(toremove)
        return self

    def _decimated(self,dcnt = 2):
        for x in range(dcnt):self._decimate()
        return self

    # use topolgy to incrementally smooth the mesh
    def _smooth(self):
        mdvs = self.meshdata.vertices
        pcs = self.meshdata.pcoords
        vs = self.vertices
        ts = self.tverts
        deltas = []
        tomove = []
        for tv in range(self.tvcnt):
            v = mdvs[vs[ts[tv][0]]]
            vring = self.vrings[tv]
            twght = self.twghts[tv]
            neighbors = [mdvs[vs[ts[vr][0]]] for vr in vring]
            ncom = dpv.center_of_mass([pcs[n[0]] for n in neighbors])
            #delta = dpv.v1_v2(pcs[v[0]],ncom).scale_u(0.02)
            delta = dpv.v1_v2(pcs[v[0]],ncom).scale(twght)
            deltas.append(delta)
            tomove.append(pcs[v[0]])
        for tdx in range(self.tvcnt):
            tomove[tdx].translate(deltas[tdx])
        return self

    # use topolgy to incrementally smooth the mesh iteratively
    def _smooths(self,scnt):
        for x in range(scnt):self._smooth()
        return self

    # recalculate normals for each vertex based on topology
    # faces added as sharp should appear flat
    # faces added as smooth should appear smooth
    def _calculate_normals(self):
        pcs = self.meshdata.pcoords
        ncs = self.meshdata.ncoords
        vcs = self.meshdata.vertices
        for vdx in range(self.tvcnt):
            v = vcs[self.vertices[self.tverts[vdx][0]]]
            fring = self.frings[vdx]
            if not fring:continue
            fnorms = []
            for vf in fring:
                aface = self.faces[vf]
                p1 = pcs[vcs[aface[0]][0]]
                p2 = pcs[vcs[aface[1]][0]]
                p3 = pcs[vcs[aface[2]][0]]
                fnorm = dpr.normal(p1,p2,p3)
                skip = False
                for fn in fnorms:
                    if fnorm.near(fn):
                        skip = True
                        break
                if not skip:fnorms.append(fnorm)
            ncs[v[1]] = dpv.center_of_mass(fnorms).normalize()
        return self

    def _add_vertex(self,v):
        newv = self.vcnt
        self.vertices.append(v)
        self.vcnt += 1
        return newv

    def _add_face(self,f,fm):
        newf = self.fcnt
        self.faces.append(f)
        self.face_mats.append(fm)
        self.fcnt += 1
        return newf

    # this messes up frings
    def _rem_face(self,fdx):
        self.tfaces.pop(fdx)
        self.faces.pop(fdx)
        self.face_mats.pop(fdx)
        self.tfcnt -= 1
        self.fcnt -= 1
        for fr in self.frings:
            torem = []
            for fdx in range(len(fr)):
                fv = fr[fdx]
                if fdx == fv:torem.append(fv)
                elif fdx < fv:fr[fdx] -= 1
                else:pass
            for trm in torem:
                fr.remove(trm)

    def _rem_faces(self,rem):
        rem.reverse()
        for trm in rem:
            self._rem_face(trm)

    # update topology using new face data
    # vs is a list of indices to elements of meshdata.vertices
    # fs is a list of faces, with indices in the space of vs
    #
    # for each v, find its topological entry, or create one
    # translate fs to that space?
    def _topology(self,vs,fs):
        vs = [self._topovertex(v) for v in vs]
        for f in fs:tface = self._topoface(vs,f)
        #return self
        return vs

    # given vs, a list of indices of vertices in meshdata.vertices
    # given fs, a list of tuples of indices to elements of self.vertices
    def _add_data(self,vs,fs,fms):
        vis = []
        fis = []
        for v in vs:vis.append(self._add_vertex(v))
        for f in fs:fis.append(tuple(vis[fpt] for fpt in f))
        for fdx in range(len(fis)):
            f,fm = fis[fdx],fms[fdx]
            self._add_face(f,fm)
        tvs = self._topology(vs,fs)
        return tvs
        #return self

    # given three points, add new triangle face
    def _add_triangle(self,p1,p2,p3,n1,n2,n3,u1,u2,u3,w1,w2,w3,m):
        v1 = self.meshdata._vertex(self.groups[self.group],p1,n1,u1)
        v2 = self.meshdata._vertex(self.groups[self.group],p2,n2,u2)
        v3 = self.meshdata._vertex(self.groups[self.group],p3,n3,u3)
        vs = [v1,v2,v3]
        fs = [(0,1,2)]
        fm = self.meshdata._lookup_mat(m)
        fms = [fm]
        tvs = self._add_data(vs,fs,fms)
        self.twghts[tvs[0]] = w1
        self.twghts[tvs[1]] = w2
        self.twghts[tvs[2]] = w3

    # given three points, add new triangle face
    def _add_quad(self,p1,p2,p3,p4,n1,n2,n3,n4,u1,u2,u3,u4,w1,w2,w3,w4,m):
        v1 = self.meshdata._vertex(self.groups[self.group],p1,n1,u1)
        v2 = self.meshdata._vertex(self.groups[self.group],p2,n2,u2)
        v3 = self.meshdata._vertex(self.groups[self.group],p3,n3,u3)
        v4 = self.meshdata._vertex(self.groups[self.group],p4,n4,u4)
        vs = [v1,v2,v3,v4]
        fs = [(0,1,2,3)]
        fm = self.meshdata._lookup_mat(m)
        fms = [fm]
        tvs = self._add_data(vs,fs,fms)
        self.twghts[tvs[0]] = w1
        self.twghts[tvs[1]] = w2
        self.twghts[tvs[2]] = w3
        self.twghts[tvs[3]] = w4

    # given three points, add new triangle face
    def _triangle(self,p1,p2,p3,ns = None,us = None,ws = None,m = None):
        if ns is None:
            n = dpr.normal(p1,p2,p3)
            ns = (n.copy(),n.copy(),n)
        if us is None:
            us = (dpv.vector2d(0,1),dpv.vector2d(0,0),dpv.vector2d(1,0))
        if ws is None:
            ws = (dpv.one(),dpv.one(),dpv.one())
        n1,n2,n3 = ns
        u1,u2,u3 = us
        w1,w2,w3 = ws
        if m is None:m = 'generic'
        self._add_triangle(p1,p2,p3,n1,n2,n3,u1,u2,u3,w1,w2,w3,m)
        return self

    # given four points, add two new triangle faces
    def _quad(self,p1,p2,p3,p4,ns = None,us = None,ws = None,m = None):
        if ns is None:
            #n = dpr.normal(p1,p2,p3)
            n = dpv.zhat.copy()
            ns = (n.copy(),n.copy(),n.copy(),n)
        if us is None:
            #us = (dpv.vector2d(0,1),dpv.vector2d(0,0),
            #      dpv.vector2d(1,0),dpv.vector2d(1,1))
            us = (dpv.vector2d(0,0),dpv.vector2d(0,0),
                  dpv.vector2d(0,0),dpv.vector2d(0,0))
        if ws is None:
            n = dpv.one()
            ns = (n.copy(),n.copy(),n.copy(),n)
        n1,n2,n3,n4 = ns
        u1,u2,u3,u4 = us
        w1,w2,w3,w4 = ws
        if m is None:m = 'generic'
        self._add_quad(p1,p2,p3,p4,n1,n2,n3,n4,u1,u2,u3,u4,w1,w2,w3,w4,m)
        return self

    # return geometry data organized as dict of materials
    def _face_dict(self):
        mats = self.meshdata.materials
        mcnt = len(mats)
        fcnt = len(self.face_mats)
        fa = {}
        for mdx in range(mcnt):
            ma = mats[mdx]
            fa[ma] = []
            for fmdx in range(fcnt):
                if self.face_mats[fmdx] == mdx:
                    fa[ma].append(self.faces[fmdx])
            if not fa[ma]:del fa[ma]
        return fa

    def _modeldata(self):
        ps,ns,us,fs,fms = [],[],[],[],[]
        pcs = self.meshdata.pcoords
        ncs = self.meshdata.ncoords
        ucs = self.meshdata.ucoords
        vcs = self.meshdata.vertices
        vs = self.vertices

        fnum = 0
        for fdx in range(self.fcnt):
            f = self.faces[fdx]
            fvn = len(f)
            fps = [pcs[vcs[vs[fv]][0]] for fv in f]
            fns = [ncs[vcs[vs[fv]][1]] for fv in f]
            fus = [ucs[vcs[vs[fv]][2]] for fv in f]
            ps.extend(fps)
            ns.extend(fns)
            us.extend(fus)
            fs.append([fnum+x for x in range(fvn)])
            fms.append(self.face_mats[fdx])
            fnum += fvn
        margs = {
            'pcoords':ps,
            'ncoords':ns,
            'ucoords':us,
            'faces':fs,
            'face_mats':fms,
                }
        return margs



#####
class terrain(dmo.model):
    def __init__(self,*args,**kwargs):
        dmo.model.__init__(self,*args,**kwargs)
        self._def('mesh',None,**kwargs)

    # return geometry data organized as dict of materials
    def _face_dict(self):
        pxy = self.mesh._modeldata()
        self.pcoords = pxy['pcoords']
        self.ncoords = pxy['ncoords']
        self.ucoords = pxy['ucoords']
        self.faces = pxy['faces']
        self.face_mats = pxy['face_mats']
        return dmo.model._face_dict(self)
#####


# if two points share a position, they can share a normal
# if two points positions differ, they may not share normals/uvs

# if two points share a normal, they can share a uv
# if two points normals differ, they may not share uv

# if two points share position, normal and uv, they are the same vertex
#   two uses of the same vertex will reference the same vertex tuple
#   they form a smooth corner, where the normal should be recalculated
#     once the topology is no longer changing
# if two points share position but not normal (and thus not uv)
#   they form a sharp corner, where edge may not appear smooth
# if two points share 

# the reuse of vertex tuples is what constitutes topological identity
# two vertices which share a face form an edge

def test():
    import dilap.construct as dlc

    mdata = meshdata()
    p = [dpv.vector(0,0,0),dpv.vector(1,0,0),dpv.vector(1,1,0),dpv.vector(0,1,0),
         dpv.vector(0,0,1),dpv.vector(1,0,1),dpv.vector(1,1,1),dpv.vector(0,1,1)]
    dpv.scale_coords(p,dpv.one().scale_u(10))

    m = mdata._mesh()
    m._quad(p[3],p[2],p[1],p[0])._quad(p[4],p[5],p[6],p[7])
    m._quad(p[0],p[1],p[5],p[4])._quad(p[2],p[3],p[7],p[6])
    m._quad(p[1],p[2],p[6],p[5])._quad(p[3],p[0],p[4],p[7])
    m._decimate()._decimate()._smooths(100)._triangulate()._calculate_normals()._project_uv_flat()

    tmodels = [terrain(mesh = msh) for msh in mdata.meshes]
    dlc.build(*tmodels)

    pdb.set_trace()



