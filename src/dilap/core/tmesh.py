import dilap.core.base as db
import dilap.core.tools as dpr

import dilap.core.model as dmo
import dilap.primitive.cube as dcu
import dilap.primitive.cylinder as dcy

import dp_vector as dpv
import dp_quaternion as dpq
import dp_ray as dr

import pdb,numpy,random,sys,traceback
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D

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

    def plot(self,ax):
        ax.plot([self.p.x],[self.p.y],zs = [self.p.z],marker = 'o')
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
        if not self.vring:return dpv.zero()
        ns = [self.mesh.vs[x].p for x in self.vring]
        ds = [dpv.distance(n,self.p) for n in ns]
        ws = [dpr.clamp(d/max(ds),0.5,1.0) for d in ds]
        ncom = dpv.center_of_mass_weighted(ns,ws)
        return dpv.v1_v2(self.p,ncom).scale(self.w)

# vs is a list of vertex objects, to be plotted one by one
def plot_vertices(vs):
    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    for v in vs:v.plot(ax)
    plt.show()

# vs is a list of vertex objects referenced by the tuples
# representing edges, found in es
def plot_edges(vs,es):
    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    for e in es:
        v1,v2 = e
        vs[v1].plot(ax)
        vs[v2].plot(ax)
        ax.plot([vs[v1].p.x,vs[v2].p.x],
                [vs[v1].p.y,vs[v2].p.y],
            zs =[vs[v1].p.z,vs[v2].p.z])
    plt.show()
    
# e1 and e2 are edges (tuples of ints)
# return True if they are the same, irrespective of direction
def sameedge(e1,e2):
    if e1[0] == e2[0] and e1[1] == e2[1]:return True
    if e1[0] == e2[1] and e1[1] == e2[0]:return True
    return False

class mesh:
    
    def newvdat(self,p,n,u,w):
        new = vertex(p,n,u,w,self)
    def newvdata(self,ps,ns,us,ws):
        vst = len(self.vs)
        for vdx in range(len(ps)):
            newv = vertex(ps[vdx],ns[vdx],us[vdx],ws[vdx],self)
        vrng = [x for x in range(vst,len(self.vs))]
        return vrng
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
        fdxs = sorted(fdxs)
        fdxs.reverse()
        for fdx in fdxs:
            self.rmfdat(fdx)
    def gfdata(self,rng = None):
        if rng is None:rng = range(len(self.fs))
        gfs = [self.fs[x] for x in rng]
        return gfs
    def gfpdat(self,fx):
        gfps = [self.vs[x].p for x in self.fs[fx]]
        return gfps
    def gfpdata(self,rng = None):
        if rng is None:rng = range(len(self.fs))
        gfps = [[self.vs[y].p for y in self.fs[x]] for x in rng]
        return gfps
    def gfpdata_aaabbb(self,bb,rng = None):
        if rng is None:rng = range(len(self.fs))
        gfps = []
        for x in rng:
            ptri = [self.vs[y].p for y in self.fs[x]]
            if bb.intersect_tri(ptri):gfps.append(ptri)
        return gfps
    def intersect_aaabbb(self,bb,rng = None):
        if rng is None:rng = range(len(self.fs))
        gfs = []
        for x in rng:
            ptri = [self.vs[y].p for y in self.fs[x]]
            if bb.intersect_tri(ptri):gfs.append(x)
        return gfs

    def plot(self,rng = None):
        if rng is None:rng = range(len(self.vs))
        vs = [self.vs[x] for x in rng]
        plot_vertices(vs)
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

    def delaunay_bridge(self,eloop1,eloop2):
        fcnt = len(self.fs)
        eloop1 = self.edgeloop(eloop1)
        eloop2 = self.edgeloop(eloop2)
        eloop = self.edgeloop_bridge(eloop1,eloop2)

        #vrng = [x[0] for x in eloop]
        vrng = [e[0] for e in eloop1]
        patch = self.delaunaymethod(vrng)

        frng = range(fcnt,len(self.fs))
        return frng

    # order the indices found in rng based on the positions of the points
    # of the vertices which they correspond to
    def delaunay_lexicographic_range(self,rng):
        ps = [self.vs[x].p for x in rng]
        lrng = [rng[x] for x in dpv.lexicographic_order_xy(ps)]
        return lrng

    # rng is a set of vertex indices which are relevant to the patching
    # rng is assumed to be coplanar??
    def delaunaymethod(self,rng):
        fcnt = len(self.fs)
        ps = [self.vs[x].p for x in rng]
        bb = (dpv.project_coords(ps,dpv.xhat),dpv.project_coords(ps,dpv.yhat))
        bb[0].x -= 100
        bb[0].y += 100
        bb[1].x -= 100
        bb[1].y += 100

        ps.extend([
            dpv.vector(bb[0].x,bb[1].x,0),dpv.vector(bb[0].y,bb[1].x,0),
            dpv.vector(bb[0].y,bb[1].y,0),dpv.vector(bb[0].x,bb[1].y,0)])
        
        rgcnt = len(rng)
        cover = [(0+rgcnt,1+rgcnt,2+rgcnt),(0+rgcnt,2+rgcnt,3+rgcnt)]
        plane = (dpv.zero(),dpv.zhat.copy())

        print('prerange',rng)
        rng = self.delaunay_lexicographic_range(rng)
        print('postrange',rng)

        for x in rng:
            self.delaunaystep(ps,cover,x,plane)

            xs = [v.x for v in ps]
            ys = [v.y for v in ps]
            zs = [v.z for v in ps]
            fig = plt.figure()
            #ax = fig.add_subplot(111,projection = '3d')
            ax = fig.add_subplot(111)
            #[ax.plot([x],[y],[z],marker = 'x',zdir = 'z') for x,y,z in zip(xs,ys,zs)]
            #[ax.plot([x],[y],[0.0],marker = 'x',zdir = 'z') for x,y,z in zip(xs,ys,zs)]
            [ax.plot([x],[y],marker = 'x') for x,y,z in zip(xs,ys,zs)]
            for c in cover:
                ptri = [ps[cx] for cx in c]
                ptri.append(ptri[0])
                xs = [v.x for v in ptri]
                ys = [v.y for v in ptri]
                #zs = [v.z for v in ptri]
                zs = [0.0 for v in ptri]
                #ax.plot(xs,ys,zs = zs,marker = 'o',zdir = 'z')
                ax.plot(xs,ys,marker = 'o')
            plt.show()

        #pdb.set_trace()

        extras = []
        for cvx in range(len(cover)):
            for c in cover[cvx]:
                if not c < rgcnt:
                    extras.append(cvx)
                    break
        extras = sorted(extras)
        extras.reverse()
        for x in extras:cover.pop(x)
        cover = [tuple([rng[c] for c in cover[cx]]) for cx in range(len(cover))]

        self.fdata(cover)
        for cx in range(len(cover)):
            f = cover[cx]
            [self.vs[f[x]].face(f,x,cx) for x in range(len(f))]
        frng = range(fcnt,len(self.fs))
        return frng

    # ps is a list of points, 
    # where cover specifies triangles by referencing points in ps
    # cover is the current list of triangles composing the triangulation
    # pt is the next candidate point (its index) for this delaunay step
    def delaunaystep(self,ps,cover,ptx,plane):
        control = self.vs[ptx].p
        dtris = []
        for cvx in range(len(cover)):
            ptri = [ps[c] for c in cover[cvx]]
            circump,circumr = dpr.circumscribe_tri(
                    ptri[0],ptri[1],ptri[2],plane)
            #try:
            #    circump,circumr = dpr.circumscribe_tri(
            #            ptri[0],ptri[1],ptri[2],plane)
            #except:
            #    traceback.print_exc(file=sys.stdout)
            #    print('delaunay step failure')
            #    return
            if dpr.inside_circle(control,circump,circumr,plane):
                dtris.append(cvx)
        dtris = sorted(dtris)
        dtris.reverse()

        alledges = []
        for dtri in dtris:
            alledges.append((cover[dtri][0],cover[dtri][1]))
            alledges.append((cover[dtri][1],cover[dtri][2]))
            alledges.append((cover[dtri][2],cover[dtri][0]))
            cover.pop(dtri)

        eloop = []
        for e1 in alledges:
            same = False
            for e2 in alledges:
                if e2 is e1:continue
                same = sameedge(e1,e2)
                if same:break
            if not same:eloop.append(e1)

        bnd = [eloop[0][0],eloop[0][1]]
        nxt = bnd[-1]
        while len(bnd) < len(eloop):
            for e in eloop:
                if   nxt == e[0] and not e[1] in bnd:
                    nxt = e[1]
                    break
                elif nxt == e[1] and not e[0] in bnd:
                    nxt = e[0]
                    break
            bnd.append(nxt)
        bnd = list(set(bnd))

        xbnd = bnd[:]
        xbnd.append(xbnd[0])

        fig = plt.figure()
        ax = fig.add_subplot(111)
        xs = [ps[b].x for b in xbnd]
        ys = [ps[b].y for b in xbnd]
        ax.plot(xs,ys,marker = 'o')
        plt.show()
        
        tcnt = len(bnd)
        ptdx = ps.index(control)
        for trdx in range(tcnt):
            c2 = trdx
            c3 = trdx+1
            if c3 == tcnt:c3 = 0
            cover.append((ptdx,bnd[c2],bnd[c3]))

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
            try:nxt = ns[0]
            except IndexError:raise ValueError
            nloop.append(nxt)
            nllen += 1
        return nloop

    def afm_bridge(self,eloop1,eloop2):
        eloop1 = self.edgeloop(eloop1)
        eloop2 = self.edgeloop(eloop2)
        
        #need to somehow plot this to check it..
        eloop = self.edgeloop_bridge(eloop1,eloop2)

        fcnt = len(self.fs)
        for x in range(500):
            if len(eloop) < 3:break
            self.advfront(eloop)
        #while len(eloop) > 2:self.advfront(eloop)
        frng = range(fcnt,len(self.fs))
        return frng

    def checkfront(self,eloop):
        for x in range(len(eloop)):
            el = x-1
            en = x+1 if x < len(eloop)-1 else 0
            if not eloop[el][1] == eloop[x][0]:return x
            if not eloop[en][0] == eloop[x][1]:return x

    # given an unordered list of vertex indices constituting a 
    # vertex ring, produce ordered edges between them
    def edgeloop(self,loop):
        oloop = self.order_loop(loop)
        looprng = range(len(oloop))
        eloop = [(oloop[x-1],oloop[x]) for x in looprng]
        return eloop

    # given two proper edgeloops, sever intelligently and
    # bridge the gap such that the result is a single edgeloop
    def edgeloop_bridge(self,eloop1,eloop2):
        vxs1 = [eloop1[x-1][1] for x in range(len(eloop1))]
        vxs2 = [eloop2[x-1][1] for x in range(len(eloop2))]
        vs1 = [self.vs[x] for x in vxs1]
        vs2 = [self.vs[x] for x in vxs2]
        vps2 = [v.p for v in vs2]
        minx1 = 0
        minx2 = dpv.find_closest(vs1[minx1].p,vps2,len(vps2),0.0)

        pbrg = dpr.point_line(vs1[minx1].p,vs2[minx2].p,5)
        preb = len(self.vs)
        for pb in pbrg[1:-1]:
            newv = vertex(pb,dpv.zhat.copy(),dpv.zero2d(),dpv.one(),self)
        newvrng = list(range(preb,len(self.vs)))
        newvrng.reverse()
        brg1 = [(newvrng[x-1],newvrng[x]) for x in range(1,len(newvrng))]
        brg2 = [(y,x) for x,y in brg1[::-1]]

        #THIS IS WHERE THE REAL PROBLEM IS!!
        brg1.append((brg1[-1][1],eloop1[minx1][0]))
        brg1.insert(0,(eloop2[minx2][0],brg1[0][0]))
        brg2.append((brg2[-1][1],eloop2[minx2][0]))
        brg2.insert(0,(eloop1[minx1][0],brg2[0][0]))

        bridged = brg1+eloop1[:minx1]+eloop1[minx1:]+\
                    brg2+eloop2[minx2:]+eloop2[:minx2]

        check = self.checkfront(bridged)
        #if not check is None:pdb.set_trace()

        #plot_edges(self.vs,bridged)

        #pdb.set_trace()

        return bridged

    # eloop is a list of topologically connected (edges) points
    def advfrontmethod(self,eloop):
        fcnt = len(self.fs)
        eloop = self.edgeloop(eloop)
        for x in range(500):
            if len(eloop) < 3:break
            self.advfront(eloop)
        #while len(eloop) > 2:self.advfront(eloop)
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

            #ea = dpv.signed_angle_between(e2,e1,v2.n)

            ea = numpy.arccos(dpv.dot(e1,e2))
            vn = e1.cross(e2)
            if vn.dot(v2.n) < 0.0:ea *= -1.0

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
            # shouldnt this always default to e1 and not just use whichever is shorter...
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

    # flatten the points on the interior of the loop to 
    # the plane defined by planepoint and planenormal
    # interfaces is a list of faces whose verts are moved
    # return the indices of the vertices that were moved
    def flatten(self,interfaces,planepoint,planenormal):
        projected = []
        for iface in interfaces:
            for vx in self.fs[iface]:
                if not vx in projected:
                    projected.append(vx)
                    v = self.vs[vx]
                    vp = v.p.project_plane(planepoint,planenormal)
                    v.p.translate(dpv.v1_v2(v.p,vp))
                    v.w.scale_u(0.0)
        return projected

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





###############################################################################
###############################################################################
###############################################################################

def afmtest():
    ps = [
        dpv.vector(0,0,0),dpv.vector(5,0,0),dpv.vector(5,5,0),
        dpv.vector(.5,5,0),dpv.vector(.5,10,0),
        dpv.vector(10,10,0),dpv.vector(10,-10,0),
        dpv.vector(0,-10,0),
        dpv.vector(-10,-10,0),dpv.vector(-10,10,0),
        dpv.vector(-.5,10,0),dpv.vector(-.5,5,0),
        dpv.vector(-5,5,0),dpv.vector(-5,0,0),
            ]
    ps1 = [
        dpv.vector(0,0,0),dpv.vector(5,0,0),dpv.vector(5,5,0),
        dpv.vector(.5,5,0),dpv.vector(.5,10,0),
        dpv.vector(10,10,0),dpv.vector(10,-10,0),
        dpv.vector(0,-10,0),
            ]
    ps2 = [
        dpv.vector(0,-10,0),
        dpv.vector(-10,-10,0),dpv.vector(-10,10,0),
        dpv.vector(-.5,10,0),dpv.vector(-.5,5,0),
        dpv.vector(-5,5,0),dpv.vector(-5,0,0),dpv.vector(0,0,0),
            ]

    hole = [
        dpv.vector(-5,0,0),dpv.vector(5,0,0),
        dpv.vector(5,5,0),dpv.vector(-5,5,0),
            ]

    #ps = dpr.dice_edges(ps,1)
    #ps1 = dpr.dice_edges(ps1,1)
    #ps2 = dpr.dice_edges(ps2,1)
    ps.reverse()
    ps1.reverse()
    ps2.reverse()
    prng = [x for x in range(len(ps))]
    prng1 = [x for x in range(len(ps1))]
    prng2 = [x for x in range(len(ps2))]

    def edgeloop(vrng):
        return [(vrng[x-1],vrng[x]) for x in vrng]

    es = edgeloop(prng)
    es1 = edgeloop(prng1)
    es2 = edgeloop(prng2)
    #m = meshme(ps,es = es)
    m = meshme(ps1,es = es1)

    arng1 = range(len(ps1))
    #arng2 = range(len(ps2))

    plot_edges(m.vs,es1)
    #plot_edges(m.vs,es2)
    #pdb.set_trace()

    patch = m.delaunaymethod(arng1)
    ptchr = []
    for p in patch:
        ptri = [m.vs[c].p for c in m.fs[p]]
        if dpv.inside(dpv.center_of_mass(ptri),hole):ptchr.append(p)
    m.rmfdata(ptchr)

    #patch = m.advfrontmethod(range(len(m.vs)))
    #patch = m.advfrontmethod(arng1)
    #patch = m.advfrontmethod(arng2)
    return m.pelt()
    


