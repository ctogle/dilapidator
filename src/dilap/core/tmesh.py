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

    def smooth_uniform(self,selfdex):
        ns = [self.mesh.vs[x].p for x in self.vring]
        ns.append(self.mesh.vs[selfdex].p)
        ncom = dpv.center_of_mass(ns)
        return dpv.v1_v2(self.p,ncom).scale(self.w)

    def smooth_radial(self,selfdex):
        ns = [self.mesh.vs[x].p for x in self.vring]
        ns.append(self.mesh.vs[selfdex].p)
        ncom = dpv.center_of_mass(ns)
        # weight the location of ncom by the distance of
        # each 1-ring vertex to the self.p
        # 
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
    # but the points are not in a coherent order
    # return a new list with the same contents ordered as they
    # are topologically ordered
    # sdx will be the starting vertex of the new loop
    def order_loop(self,eloop,sdx = 0):
        vs = self.vs
        ellen = len(eloop)
        nllen = 1
        nloop = [eloop[sdx]]
        while nllen < ellen:
            last = nloop[-1]
            ns = [vr for vr in vs[last].vring 
                if vr in eloop and not vr in nloop]
            nxt = ns[0]
            nloop.append(nxt)
            nllen += 1
        return nloop

    def checkfront(self,eloop):
        for x in range(len(eloop)):
            el = x-1
            en = x+1 if x < len(eloop)-1 else 0
            if not eloop[el][1] == eloop[x][0]:
                return x
            if not eloop[en][0] == eloop[x][1]:
                return x

    # eloop is a list of topologically connected (edges) points
    def advfrontmethod(self,eloop):
        fcnt = len(self.fs)
        eloop = [(eloop[x-1],eloop[x]) for x in range(len(eloop))]
        #for x in range(500):
        #    if len(eloop) < 3:break
        #    self.advfront(eloop)
        while len(eloop) > 2:self.advfront(eloop)
        frng = range(fcnt,len(self.fs))
        return frng

    def advfront(self,eloop):
        vs = self.vs
        eacnt = 0
        whch = (0,10.0)
        for ex in range(len(eloop)):
            el1,el2,el3 = eloop[ex-1][0],eloop[ex][0],eloop[ex][1]
            v1,v2,v3 = vs[el1],vs[el2],vs[el3]
            e1 = dpv.v1_v2(v2.p,v1.p).normalize()
            e2 = dpv.v1_v2(v2.p,v3.p).normalize()
            ea = dpv.signed_angle_between(e2,e1,v2.n)
            if ea < 0.01:ea = numpy.pi
            if ea <= whch[1]:whch = (eacnt,ea)
            eacnt += 1
        self.advfrontpoint(eloop,whch)

    def advfrontpoint(self,eloop,point):
        ex,ea = point
        el1,el2,el3 = eloop[ex-1][0],eloop[ex][0],eloop[ex][1]
        v1,v2,v3 = self.vs[el1],self.vs[el2],self.vs[el3]
        e1 = dpv.v1_v2(v2.p,v1.p)
        e2 = dpv.v1_v2(v2.p,v3.p)
        e1len,e2len = e1.magnitude(),e2.magnitude()
        avd = (e1len+e2len)/2.0
        e1.normalize()
        e2.normalize()
        if ea < numpy.pi/2.0:
            fs = [(el1,el2,el3)]
            eloop.pop(ex)
            eloop.pop(ex-1)
            eloop.insert((0 if ex == 0 else ex-1),(el1,el3))
        elif ea < 2.0*numpy.pi/3.0:
            mpt = dpv.midpoint(v1.p,v3.p)
            mdist = avd - dpv.distance(v2.p,mpt)
            mdelt = dpv.v1_v2(v2.p,mpt).normalize().scale_u(mdist)
            mpt.translate(mdelt)
            newv = vertex(mpt,dpv.zhat.copy(),dpv.zero2d(),dpv.one(),self)
            el4 = len(self.vs)-1
            fs = [(el1,el2,el4),(el2,el3,el4)]
            eloop.pop(ex)
            eloop.pop(ex-1)
            eloop.insert((0 if ex == 0 else ex-1),(el4,el3))
            eloop.insert((0 if ex == 0 else ex-1),(el1,el4))
        else:
            mdist = min((e1len,e2len))
            mdelt = dpv.v1_v2(v2.p,v1.p).cross(v2.n)
            mdelt.normalize().scale_u(mdist)
            mpt = dpv.midpoint(v2.p,v1.p).translate(mdelt)
            newv = vertex(mpt,dpv.zhat.copy(),dpv.zero2d(),dpv.one(),self)
            el4 = len(self.vs)-1
            fs = [(el1,el2,el4)]
            eloop.pop(ex-1)
            eloop.insert((0 if ex == 0 else ex-1),(el4,el2))
            eloop.insert((0 if ex == 0 else ex-1),(el1,el4))
        
        check = self.checkfront(eloop)
        if not check is None:pdb.set_trace()

        self.fdata(fs)
        for fx in range(len(fs)):
            f = fs[fx]
            [self.vs[f[x]].face(f,x,fx) for x in range(len(f))]

    def smooths(self,scnt,sscl,rng = None,method = 'uniform'):
        for s in range(scnt):self.smooth(sscl,rng,method)
        return self
    def smooth(self,sscl,rng = None,method = 'uniform'):
        vs = self.vs
        if rng is None:rng = range(len(vs))
        deltas = []
        for vdx in rng:
            deltas.append(vs[vdx].__getattribute__('smooth_'+method)(vdx))
            #deltas.append(vs[vdx].smooth(vdx))
        for vdx in range(len(rng)):
            vs[rng[vdx]].p.translate(deltas[vdx].scale_u(sscl))
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

    


