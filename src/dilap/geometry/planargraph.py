from .tools import isnear, epsilon
from .vec3 import vec3
from .polymath import sintsxyp, binbxy, bintselfxy, bccw
from dilap.topology import wiregraph
from dilap.core.plotting import plot_axes_xy, plot_polygon_xy, plot_point_xy, plt
import numpy
import pdb


# planar wire graph class (topological + xy-projected geometry)
class planargraph(wiregraph):

    @staticmethod
    def orderkey(x, y):
        return (x, y) if x < y else ((y, x) if y < x else None)

    def check_disconnected(self, x, y):
        key = self.orderkey(x, y)
        if key in self.elook:
            return self.re(x, y)

    def check_connected(self, x, y):
        key = self.orderkey(x, y)
        if key and not key in self.elook:
            return self.fe(x, y)

    def addintersection(self, u, v, p, q, i):
        self.check_disconnected(u, v)
        self.check_disconnected(u, i)
        self.check_disconnected(i, v)
        self.check_disconnected(p, q)
        self.check_disconnected(p, i)
        self.check_disconnected(i, q)
        self.check_connected(u, i)
        self.check_connected(i, v)
        self.check_connected(p, i)
        self.check_connected(i, q)
        return self

    def addintersections(self, e=0.001):
        found = True
        while found:
            found = False
            edges = self.edgekeys()
            for u, v in edges:
                for p, q in edges:
                    if u in (p, q) or v in (p, q):
                        continue
                    up, vp, qp, pp = self[u], self[v], self[q], self[p]
                    ip = sintsxyp(up, vp, qp, pp, e=e)
                    if isinstance(ip, tuple):
                        raise ValueError
                    elif isinstance(ip, vec3):
                        nv = self.fp(ip, e)
                        self.addintersection(u, v, q, p, nv)
                        found = True
                        break

    def interiorpoints(self):
        loops = [[self[p] for p in l] for l in self.uloops('ccw')]
        marked = []
        for j, v in self:
            for l in loops:
                if self[j].inbxy(l):
                    marked.append(j)
                    break
        return marked

    @staticmethod
    def zipb(py, e=0.01):
        c = len(py)
        ss = [(py[j - 1], py[j]) for j in range(c)]
        pg = planargraph.segstopg(ss, e)
        pg.addintersections(e)
        for j in pg.interiorpoints():
            pg.rv(j)
        uloops = pg.uloops('ccw')
        if len(uloops) > 1:
            loopls = [len(l) for l in uloops]
            uloops.pop(loopls.index(max(loopls)))
        ploops = [[pg[v] for v in l] for l in uloops]
        '''
        import dilap.core.plotting as plt
        from dilap.geometry.polymath import longestb
        ax = plt.plot_axes_xy(60, (-60, 20))
        for j, l in enumerate(ploops):
            plt.plot_polygon_xy(l, ax, lw=2, col='b')
            plt.plot_point_xy_annotate(vec3(0, 0, 0).com(l), ax, 'polygon: %d' % j)
        best = longestb(ploops)
        plt.plot_polygon_xy(best, ax, lw=8, col='g')
        plt.plt.show()
        '''
        return ploops

    @staticmethod
    def aggb(py, r, e=0.01, a=numpy.pi/3.0):
        g = planargraph.pytopg(py)
        found = True
        while found:
            found = False
            edgekeys = g.edgekeys()
            if not edgekeys:
                return
            edges, lengths = [], []
            for u, v in edgekeys:
                edges.append((u, v))
                lengths.append(g[u].d(g[v]))
            shortest = lengths.index(min(lengths))
            #print('shortest', r, lengths[shortest], [('%.2f' % l) for l in lengths])
            #print('shortest', r, lengths[shortest], 'pcnt', len(lengths))
            if lengths[shortest] < r:
                g.de(*edges[shortest])
                found = True
            if not found:
                for j, v in g:
                    ves = [(v, o) for o in g.rings[j]]
                    assert(len(ves) == 2)
                    x, y, z = g[ves[0][1]], g[j], g[ves[1][1]]
                    t = y.tov(x).angxy(y.tov(z))
                    if t < a:
                        g.miterv(j, r, e)
                        found = True
                        break
            if found and len(list(g)) < 3:
                #print('polygon collapsed')
                return
        py = g.boundary()
        return py

    def miterv(self, v, r, e):
        x, y = tuple(self.rings[v].keys())
        p, q = self[v].tov(self[x]), self[v].tov(self[y])
        a = numpy.arccos(p.cp().nrm().dot(q.cp().nrm()))
        d = r / numpy.tan(a)
        d = min(d, p.mag(), q.mag())
        dx = self[v].cp().trn(p.cp().nrm().uscl(d))
        dy = self[v].cp().trn(q.cp().nrm().uscl(d))
        nx = self.fp(dx, e)
        ny = self.fp(dy, e)
        ne = self.fe(nx, ny)
        self.rv(v)

    def boundary(self):
        loops = [[self[p] for p in l] for l in self.uloops('ccw')]
        # find loop such that none of its points are inside any other loop
        foundloops = []
        for j, u in enumerate(loops):
            found = False
            for p in u:
                for i, v in enumerate(loops):
                    if i == j:
                        continue
                    if p.inbxy(v):
                        found = True
                        break
            if not found:
                foundloops.append(u)
        loopls = [len(l) for l in foundloops]
        u = foundloops[loopls.index(max(loopls))]
        return u

    @staticmethod
    def merge_pgs(u, v, e=0.1):
        ues, uevs = u.pgtosegs()
        ves, vevs = v.pgtosegs()
        ips = []
        for (a, b), ue in zip(ues, uevs):
            for (c, d), ve in zip(ves, vevs):
                ip = sintsxyp(a, b, c, d)
                if isinstance(ip, tuple):
                    raise ValueError
                elif ip:
                    ips.append((ip, ue, ve))
        new = planargraph.segstopg(ues + ves)
        for ip, (p, q), (r, s) in ips:
            ipv = new.fp(ip, e)
            up, uq = u[p], u[q]
            vr, vs = v[r], v[s]
            nup, nuq = new.fp(up, e), new.fp(uq, e)
            nvr, nvs = new.fp(vr, e), new.fp(vs, e)
            if not (ipv == new.vcnt - 1):
                new.addintersection(nup, nuq, nvr, nvs, ipv)
        return new

    @classmethod
    def segstopg(cls, segs, epsilon=0.1):
        pg = cls()
        for e1, e2 in segs:
            v1, v2 = pg.fp(e1, epsilon, ie=False), pg.fp(e2, epsilon, ie=False)
            if v1 == v2:
                #print('seg is smaller than epsilon')
                pass
            elif not v2 in pg.rings[v1]:
                pg.fe(v1, v2)
            elif not v1 in pg.rings[v2]:
                pg.fe(v2, v1)
        return pg

    @classmethod
    def pytopg(cls, py, epsilon=0.1):
        segs = [(py[j - 1], py[j]) for j in range(len(py))]
        return planargraph.segstopg(segs, epsilon)

    def pgtosegs(self):
        es = set()
        for j, x in self:
            for o in self.rings[j]:
                if j < o:
                    es.add((j, o))
        segs = [(self[x], self[y]) for x, y in es]
        es = list(es)
        return segs, es

    def aggregate(self, e=1.0, a=0.98):
        found = True
        while found:
            found = False
            for j, v in self:
                ring = self.rings[j]
                os = set(ring.keys())
                if len(os) > 0 and len(os) < 3:
                    vp = self[j]
                    ops = [self[o] for o in os]
                    es = [vp.tov(op) for op in ops]
                    ods = [t.mag() for t in es]
                    if any(o < e for o in ods):
                        if len(os) == 1:
                            self.rv(j)
                            found = True
                            break
                        elif len(os) == 2:
                            e1, e2 = es
                            alpha = abs(e1.ang(e2))
                            if alpha > a * numpy.pi:
                                self.rv(j)
                                self.ae(*os)
                                found = True
                                break

    @staticmethod
    def aggregate_polygon(polygon, e, a):
        pg = planargraph()
        u = pg.av(p=polygon[-1])
        for j in range(len(polygon)):
            v = pg.av(p=polygon[j])
            uv = pg.ae(u, v)
            u = v

        xs, ys, zs = zip(*polygon)

        ax = plot_axes_xy((max(xs) - min(xs)), (sum(xs) / len(xs), sum(ys) / len(ys)))
        plot_polygon_xy(polygon, ax, col='b', lw=3)
        pg.plotxy(ax)

        plt.show()

        pg.aggregate(e, a)
        py = pg.polygon(0)[0]

        ax = plot_axes_xy((max(xs) - min(xs)), (sum(xs) / len(xs), sum(ys) / len(ys)))
        plot_polygon_xy(polygon, ax, col='b', lw=4)
        plot_polygon_xy(py, ax, col='y', lw=2)
        pg.plotxy(ax)

        plt.show()

        return py

    def fitbxy(self, b, w=1.0):
        ps = [v[1]['p'] for v in self.vs if v is not None]
        bpjx = vec3(1, 0, 0).prjps(b)
        bpjy = vec3(0, 1, 0).prjps(b)
        pjx = vec3(1, 0, 0).prjps(ps)
        pjy = vec3(0, 1, 0).prjps(ps)
        cx = sum(pjx)-sum(bpjx)
        cy = sum(pjy)-sum(bpjy)
        recenter = vec3(cx,cy,0).uscl(-0.5)
        for v in self.vs:
            if v is not None:
                v[1]['p'].trn(recenter)
        d1x = (pjx[1]-pjx[0])+(pjy[1]-pjy[0])
        d2x = (bpjx[1]-bpjx[0])+(bpjy[1]-bpjy[0])
        scale = w * d2x / d1x
        vec3(scale,scale,0).sclps(ps)

    def project(self, v):
        ps = [v[1]['p'] for v in self]
        minv = 100000000
        maxv = -10000000
        for p in ps:
            p = p.dot(v)
            if p < minv:
                minv = p
            if p > maxv:
                maxv = p
        return minv, maxv

    def intersectb(self, b, oobs=[], plot=False):
        if plot:
            ax = plot_axes_xy(700, aspect='equal')
            plot_polygon_xy(b, ax, col='b')
            for oob in oobs:
                plot_polygon_xy(oob, ax, col='r', lw=1)
            self.plotxy(ax, number=False)
            plt.show()

        marked = []
        for vx, v in enumerate(self):
            if not v[1]['p'].inbxy(b):
                marked.append(vx)
        for oob in oobs:
            for vx, v in enumerate(self):
                if v[1]['p'].inbxy(oob):
                    marked.append(vx)
        for vx, v in enumerate(self):
            if self.orphan(vx):
                marked.append(vx)
        for vx in set(marked):
            self.rv(vx)

        if plot:
            ax = plot_axes_xy(700, aspect='equal')
            plot_polygon_xy(b, ax, col='b')
            for oob in oobs:
                plot_polygon_xy(oob, ax, col='r', lw=1)
            self.plotxy(ax, number=False)
            plt.show()

        return self

    def polygon(self,r,epsilon = 0.1,z = vec3(0,0,1),findeseam = False):
        loops, seams, eseam = self.uloops('ccw'), [], None

        '''
        for u in loops:
            ax = plot_axes_xy(50)
            for v in loops:
                v = [self[x] for x in v]
                ax = plot_polygon_xy(v, ax, lw=2, col='b')
            u = [self[x] for x in u]
            ax = plot_polygon_xy(u, ax, lw=4, col='g')
            plt.show()
        '''

        #print('planargraph loop count: %d' % len(loops))

        for lp in loops:
            seam = []
            seams.append(seam)
            for lpx in range(len(lp)):
                lp0, lp1, lp2 = lp[lpx-2], lp[lpx-1], lp[lpx]
                lpp0, lpp1, lpp2 = self[lp0], self[lp1], self[lp2]
                lptn1, lptn2 = lpp0.tov(lpp1).nrm(), lpp1.tov(lpp2).nrm()
                lpnm1, lpnm2 = z.crs(lptn1).uscl(r), z.crs(lptn2).uscl(r)
                if lp0 == lp2:
                    stemoffset = lptn1.cp().uscl(r)
                    seam.append(lpp1.cp().trn(lpnm1+stemoffset))
                    seam.append(lpp1.cp().trn(lpnm2+stemoffset))
                else:
                    if lpnm1.isnear(lpnm2):
                        cnm = lpnm1
                    else:
                        s1 = lpp0.cp().trn(lpnm1).trn(lptn1.cp().uscl(-1000))
                        s2 = lpp1.cp().trn(lpnm1).trn(lptn1.cp().uscl( 1000))
                        s3 = lpp1.cp().trn(lpnm2).trn(lptn2.cp().uscl(-1000))
                        s4 = lpp2.cp().trn(lpnm2).trn(lptn2.cp().uscl( 1000))
                        ip = sintsxyp(s1, s2, s3, s4, col=0)
                        cnm = lpp1.tov(ip)
                    seam.append(lpp1.cp().trn(cnm))

            #if bintselfxy(seam):
            #    print('self intersecting polygon!!')

            if not bccw(seam):
                seam.reverse()
            if eseam is None:
                print('eseam was none; asserted 0')
                eseam = 0
            else:
                try:
                    if binbxy(seams[eseam], seam):
                        eseam = len(seams) - 1
                        print('eseam rechosen as latest loop;', eseam)
                except:

                    print('crap', len(loops))
                    ax = plot_axes_xy(100, (50, 100))
                    ax = self.plotxy(ax,number = True, col='g')
                    plt.show()

                    for u in loops:
                        ax = plot_axes_xy(100, (50, 100))
                        #for v in loops:
                        #    v = [self.vs[x][1]['p'] for x in v]
                        #    ax = plot_polygon_xy(v, ax, lw=2, col='b')
                        u = [self.vs[x][1]['p'] for x in u]
                        ax = plot_polygon_xy(u, ax, lw=4, col='g')
                        plt.show()

                    quit()
                #many = (sum(1 if x == y else 0 for x, y in zip(seam, seams[eseam])))
                #print(many, many / len(seam), many / len(seams[eseam]))
                #print(seam[0] == seams[eseam][0])
                #print(seam[-1], seam[0], seam[1])
                #print(seams[eseam][-1], seams[eseam][0], seams[eseam][1])

        if eseam is None:
            print('no eseam!')
            self.plotxy(l=700)
            plt.show()
            return None

        if findeseam:
            return eseam, seams, loops
        else:
            py = [seams.pop(eseam),seams]
            return py

    def trn(self, t):
        for x in range(self.vcnt):
            v = self.vs[x]
            if not v is None:
                v[1]['p'].trn(t)
        return self

    # split an edge into two edges
    # make the third vertex by lerping between the other two
    def be(self,u,v,f = 0.5,**kws):
        if not 'p' in kws:
            up,vp = self.vs[u][1]['p'],self.vs[v][1]['p']
            kws['p'] = up.lerp(vp,f)
        nv = self.av(**kws)
        nr1,nr2 = self.se(u,v,nv)
        return nv,nr1,nr2

    # compute the counterclockwise angle between two edges
    def ea(self,u,v,w):
        up = self.vs[u][1]['p']
        vp = self.vs[v][1]['p']
        wp = self.vs[w][1]['p']
        e1 = vp.tov(up)
        e2 = vp.tov(wp)
        etn1 = e1.cp().nrm()
        etn2 = e2.cp().nrm()
        #para  = gtl.isnear(etn1.dot(etn2), 1, gtl.epsilon)
        #apara = gtl.isnear(etn1.dot(etn2),-1, gtl.epsilon)
        para  = isnear(etn1.dot(etn2), 1, epsilon)
        apara = isnear(etn1.dot(etn2),-1, epsilon)
        if apara:sa = numpy.pi
        elif para:sa = 0
        else:sa = e1.sang(e2,vec3(0,0,1))
        if sa < 0:sa = 2*numpy.pi+sa
        return sa

    # compute the edge ordering of a vertex v relative to an edge u,v
    def eo(self,u,v):
        uvas,avor = [],[]
        vor = self.orings[v]
        for w in vor:
            if w == u:uvwa = 0
            else:uvwa = self.ea(u,v,w)
            if not uvas or uvwa > uvas[-1]:
                uvas.append(uvwa);avor.append(w)
            else:
                fnd = False
                for uvax in range(len(uvas)):
                    if uvwa < uvas[uvax]:
                        fnd = True
                        break
                if fnd:
                    uvas.insert(uvax,uvwa)
                    avor.insert(uvax,w)
                else:
                    uvas.append(uvwa)
                    avor.append(w)
        uori = uvas.index(min(uvas))
        ror = avor[uori+1:]+avor[:uori]
        return ror

    # given an edge, find the next edge taking the first
    # clockwise turn available
    def cw(self,u,v):
        vor = self.orings[v][:]
        if len(vor) == 2:vor.remove(u)
        if len(vor) == 1:return vor[0]
        ror = self.eo(u,v)
        if ror:tip = ror[0]
        else:tip = u
        return tip

    # given an edge, find the next edge taking the first
    # counterclockwise turn available
    def ccw(self,u,v):
        vor = self.orings[v][:]
        if len(vor) == 2:vor.remove(u)
        if len(vor) == 1:return vor[0]
        ror = self.eo(u,v)
        if ror:tip = ror[-1]
        else:tip = u
        return tip

    # find a vertex within epsilon of p or create one
    # ie : 0 - split edges if sufficiently close to an edge
    def fp(self,p,epsilon,ie = 1,**vkws):
        for j, v in self:
            if self[j].d(p) < epsilon:
                nv = j
                break
        else:
            nv = self.av(p = p.cp(),**vkws)
        if ie and nv == self.vcnt - 1:
            found = True
            while found:
                found = False
                for u, v in self.elook:
                    up,vp = self[u], self[v]
                    ped = p.dexy(up, vp)
                    if ped > -1 and ped < epsilon:
                        self.se(u, v, nv)
                        found = True

                        #print('se', ped)
                        #ax = plot_axes_xy(5, (3, 0))
                        #self.plotxy(ax, lw=3, col='b', s=0)
                        #plt.show()

                        break
        return nv

    # find an edge between two vertices or create one 
    def fe(self,u,v,**ekws):
        if u == v:
            print('attempted fe with u == v??')
        if   (u,v) in self.elook:return self.elook[(u,v)]
        elif (v,u) in self.elook:return self.elook[(v,u)]
        return self.ae(u,v,**ekws)

    def smooth(self,i = 1,w = 1.0):
        for s in range(i):
            dels = []
            for j in range(self.vcnt):
                v = self.vs[j]
                if v is None:continue
                ops = [self.vs[o][1]['p'] for o in self.rings[j]]
                dels.append(v[1]['p'].tov(vec3(0,0,0).com(ops)).uscl(w*v[1]['w']))
            for j in range(self.vcnt):
                v = self.vs[j]
                if v is None:continue
                v[1]['p'].trn(dels[j])
        return self

    def smooth_sticks(self,i = 1,w = 1.0,d = 5.0):
        for s in range(i):
            dels = []
            for j in range(self.vcnt):
                v = self.vs[j]
                if v is None:continue
                ops = []
                for o in self.rings[j]:
                    vod = v[1]['p'].d(self.vs[o][1]['p'])
                    stickp = v[1]['p'].lerp(self.vs[o][1]['p'],1.0-d/vod)
                    ops.append(stickp)
                dels.append(v[1]['p'].tov(vec3(0,0,0).com(ops)).uscl(w*v[1]['w']))
            for j in range(self.vcnt):
                v = self.vs[j]
                if v is None:continue
                v[1]['p'].trn(dels[j])
        return self

    # plot the vertices and edges of the graph
    def plotxy(self,ax = None,l = 10,s = 1.0,number = True, col='g', marker=None, **kws):
        import dilap.core.plotting as dtl
        if ax is None:ax = dtl.plot_axes_xy(l,**kws)
        for j in range(self.vcnt):
            i = self.vs[j]
            if i is None:
                continue
            ip = i[1]['p']
            #if number and (len(tuple(self.rings[j].keys())) > 2):
            if number:
                jstr = str(j)+str([v for v in self.rings[j]])
                #jstrccw = str([self.ccw(j,v) for v in self.rings[j]])
                #jstrcw = str([self.cw(j,v) for v in self.rings[j]])
                #jstr += '\n'+jstrccw+'\n'+jstrcw
                ax = dtl.plot_point_xy_annotate(ip,ax,jstr)
                ax = dtl.plot_point_xy(ip,ax,col = 'r')
            if marker:
                ax = dtl.plot_point_xy(ip, ax, col='r', mk=marker)
        for k in self.rings:
            vr = self.rings[k]
            for ov in vr:
                if vr[ov] is None:
                    continue
                re = (self.vs[k][1]['p'].cp(),self.vs[ov][1]['p'].cp())
                rn = vec3(0,0,1).crs(re[0].tov(re[1])).nrm().uscl(s)
                re = (re[0].trn(rn),re[1].trn(rn))
                ax = dtl.plot_edges_xy(re, ax, lw=4, col=col)
        return ax

    # plot the vertices and edges of the graph
    def plot(self,ax = None,l = 10,**kws):
        import dilap.core.plotting as dtl
        if ax is None:ax = dtl.plot_axes(l)
        for j in range(self.vcnt):
            i = self.vs[j]
            if i is None:continue
            ip = i[1]['p']
            ax = dtl.plot_point(ip,ax,col = 'r')
        for k in self.rings:
            vr = self.rings[k]
            for ov in vr:
                if vr[ov] is None:continue
                re = (self.vs[k][1]['p'].cp(),self.vs[ov][1]['p'].cp())
                rn = vec3(0,0,1).crs(re[0].tov(re[1])).nrm()
                re = (re[0].trn(rn),re[1].trn(rn))
                ax = dtl.plot_edges(re,ax,lw = 2,col = 'g')
        return ax
