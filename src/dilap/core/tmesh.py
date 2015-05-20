import dilap.core.base as db
import dilap.core.tools as dpr

import dilap.core.model as dmo
import dilap.primitive.cube as dcu
import dilap.primitive.cylinder as dcy

import dp_vector as dpv
import dp_quaternion as dpq

import numpy,random,pdb

class vertex:

    # given the index of a vert with which this vert forms an edge
    # append this index to self.vring if its not already there
    def edge(self,other):
        if not other in self.vring:
            self.vring.append(other)
    # given a face (list of indices), the index of this vertex
    # and the index of the new face within the mesh
    # properly connect edges and append facedex to self.frings
    def face(self,fs,selfdex,facedex):
        flen = len(fs)
        if not facedex in self.fring:
            self.fring.append(facedex)
        for fdx in range(flen):
            if fdx == selfdex:
                if fdx == 0:e1dx = flen-1
                else:e1dx = fdx-1
                if fdx == flen-1:e2dx = 0
                else:e2dx = fdx+1
                self.edge(fs[e1dx])
                self.edge(fs[e2dx])
                return
    def __init__(self,p,n,u,w,m = None):
        self.p = p # position vector
        self.n = n # normal vector
        self.u = u # uv vector2d
        self.w = w # movement weight vector
        self.vring = [] # list of indices in the vertex list of m
        self.fring = [] # list of indices in the face list of m
        self.mesh = m
        if not m is None:m.vdat(self)

    def smooth(self,selfdex):
        ns = [self.mesh.vs[x].p for x in self.vring]
        ns.append(self.mesh.vs[selfdex].p)
        ncom = dpv.center_of_mass(ns)
        return dpv.v1_v2(self.p,ncom).scale(self.w)


class mesh:
    
    def vdat(self,v):
        self.vs.append(v)
    def vdata(self,vs):
        for v in vs:
            self.vs.append(v)
    def gpdata(self,rng = None):
        if rng is None:rng = range(len(self.vs))
        gps = [self.vs[x].p for x in rng]
        return gps
    def edata(self,es):
        for e in es:
            self.es.append(e)
    def fdata(self,fs):
        for f in fs:
            self.fs.append(f)
    def rmfdat(self,fdx):
        for vdx in self.fs[fdx]:
            self.vs[vdx].fring.remove(fdx)
        self.fs[fdx] = None
        # these indices can be recycled if i keep a stack of them!!
    def rmfdata(self,fdxs):
        for fdx in fdxs:
            self.rmfdat(fdx)
    def __init__(self):
        self.vs = [] # list of vertex objects
        self.es = [] # list of vertex indices forming edges
        self.fs = [] # list of vertex indices forming faces
        
    def cut_hole(self,rng):
        bound = []
        for rx in rng:
            f = self.fs[rx]
            if f is None:continue
            for vx in f:
                v = self.vs[vx]
                for vfr in v.fring:
                    if not vfr in rng:
                        if not vx in bound:
                            bound.append(vx)
            self.rmfdat(rx)
        return bound

    # eloop is a list of topologically connected (edges) points
    def advfrontmethod(self,eloop):
        #self.advfront(eloop)
        #while len(eloop) > 2:self.advfront(eloop)
        center = dpv.center_of_mass([self.vs[x] for x in eloop])
        for x in range(200):self.advfront(eloop,center)

    def advfront(self,eloop,attract):
        vs = self.vs
        eacnt = 0
        whiche = (0,0,0,10.0)
        for el in eloop:
            ev = vs[el]
            ns = [vx for vx in ev.vring if vx in eloop]
            if len(ns) != 2:
                eloop.pop(eacnt)
                return
            e1 = dpv.v1_v2(ev.p,vs[ns[0]].p).normalize()
            e2 = dpv.v1_v2(ev.p,vs[ns[1]].p).normalize()
            ea = numpy.arccos(e1.dot(e2))
            if ea < whiche[3]:
                whiche = (eacnt,ns[0],ns[1],ea)
            eacnt += 1
        self.advfrontpoint(eloop,whiche,attract)

    def advfrontpoint(self,eloop,which,attract):
        dpvd = dpv.distance

        vs = self.vs
        we = eloop[which[0]]
        w1 = which[1]
        w2 = which[2]
        e1d = dpvd(vs[we].p,vs[w1].p)
        e2d = dpvd(vs[we].p,vs[w2].p)
        e3d = dpvd(vs[w1].p,vs[w2].p)
        avd = (e1d+e2d)/2.0

        #if e3d > avd:
        if False:

            #mpt = dpv.midpoint(vs[w1].p,vs[w2].p)
            #mdist = avd - dpv.distance(vs[we].p,mpt)
            #mdelt = dpv.v1_v2(vs[we].p,mpt).normalize().scale_u(mdist)
            #mpt.translate(mdelt)

            eg1 = dpv.v1_v2(vs[we].p,vs[w1].p).normalize()
            eg2 = dpv.v1_v2(vs[we].p,vs[w2].p).normalize()
            wwc = dpv.cross(eg2,eg1)
            qa = numpy.arcsin(wwc.magnitude())/2.0
            qv = wwc.normalize()
            if dpvd(qv,vs[we].n) > dpvd(qv.copy().flip(),vs[we].n):
                qv.flip()
            q = dpq.q_from_av(qa,qv)
            mpt = dpv.v1_v2(vs[we].p,vs[w2].p).rotate(q).translate(vs[we].p)

            #if dpv.distance(mpt,attract) > dpv.distance(vs[we].p,attract):
            #    mpt.translate(dpv.v1_v2(vs[we].p,mpt).flip().scale_u(2.0))

            newv = vertex(mpt,dpv.zhat.copy(),dpv.zero2d(),dpv.zero(),self)

            w3 = len(self.vs)-1
            fs = [(w1,we,w3),(we,w2,w3)]
            eloop.append(w3)
        else:
            fs = [(w1,we,w2)]

        self.fdata(fs)
        for fx in range(len(fs)):
            f = fs[fx]
            [self.vs[f[x]].face(f,x,fx) for x in range(len(f))]
        eloop.pop(which[0])

    def smooths(self,scnt,sscl):
        for s in range(scnt):
            self.smooth(sscl)
        return self
    def smooth(self,sscl):
        deltas = []
        vs = self.vs
        for vdx in range(len(vs)):
            deltas.append(vs[vdx].smooth(vdx))
        for vdx in range(len(vs)):
            vs[vdx].p.translate(deltas[vdx].scale_u(sscl))
        return self

    def pelt(self):
        s = dmo.model()
        for f in self.fs:
            if f is None:continue
            v1 = self.vs[f[0]]
            v2 = self.vs[f[1]]
            v3 = self.vs[f[2]]
            s._triangle(v1.p,v2.p,v3.p,[v1.n,v2.n,v3.n])
        return s
    def skeleton(self):
        s = dmo.model()
        for v in self.vs:
            m = dcu.cube().scale(dpv.one().scale_u(0.1)).translate(v.p)
            s._consume(m)
        for e in self.es:
            cq = dpq.q_from_uu(dpv.zhat,
                dpv.v1_v2(self.vs[e[0]].p,self.vs[e[1]].p))
            c = dcy.cylinder().translate_z(0.5)
            c.scale_x(0.05).scale_y(0.05).rotate(cq)
            c.translate(self.vs[e[0]].p)
            s._consume(c)
        return s

# given points, and edges/faces on them
# create a mesh with the proper vertex topology
# CANNOT detect preexisting topology
# successive calls will create topologically unconnected mesh
def meshme(ps,ns = None,us = None,ws = None,es = [],fs = []):
    if ns is None:ns = [dpv.zhat.copy() for p in ps]
    if us is None:us = [dpv.zero2d() for p in ps]
    if ws is None:ws = [dpv.one() for p in ps]
    nvs = []
    m = mesh()
    for vdx in range(len(ps)):
        newv = vertex(ps[vdx],ns[vdx],us[vdx],ws[vdx],m)
        nvs.append(newv)
    for e in es:
        nvs[e[0]].edge(e[1])
        nvs[e[1]].edge(e[0])
    m.edata(es)
    m.fdata(fs)
    for fx in range(len(fs)):
        f = fs[fx]
        [nvs[f[x]].face(f,x,fx) for x in range(len(f))]
    return m                   

    


