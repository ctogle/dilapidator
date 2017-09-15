from dilap.core import *
from dilap.core.lsystem import lgrammer
from dilap.topology import tree
import dilap.geometry.tools as gtl
import functools
import pdb


class treegrammer(lgrammer):

    def edge(ls):
        st = ls.tip.p.cp()

        tw = quat(0,0,0,0).av(ls.tip.t, ls.tip.ld)
        d = ls.tip.d.cp().rot(tw)
        ls.tip.ld = d.cp()

        #ed = st.cp().trn(d.uscl(ls.drho))

        # return two points: st & ed
        # ed is on the ball ls.drho away from st
        # such that D is minimized
        # where D is a function representing proximity to other verts
        #os = list(ls.enum(False))
        #ds = [o.p.d(ed) for o,d in os]
        #ts = [o.p.tov(ed) for o,d in os]
        #net = functools.reduce(lambda x,y : x.trn(y), ts).nrm()
        #print('SPOOOKY', len(os), ds, net)

        #ed = st.cp().trn(net.uscl(ls.drho))
        xy = ls.root.p.tov(st).cpxy().ztrn(5).nrm()
        dd = d.nrm()

        ed = st.cp().trn((xy + dd).nrm().uscl(ls.drho))

        ls.tip.p.trn(st.tov(ed))
        return st,ed

#treegrammer.dic['F'] = treegrammer.edge


def treeskin(m, gm, lsys, n=8, r=1.0):


    def branch(a, b, c, d, n):
        q12,q23 = quat(1,0,0,0),quat(1,0,0,0)
        if a:
            q12.uu(vec3(0,0,1),a.tov(b).cp().nrm())
        q23.uu(vec3(0,0,1),b.tov(c).cp().nrm())
        r1 = radius(d)
        r2 = radius(d+1)
        pts1 = b.trnps(q12.rotps(vec3(0,0,0).pring(r1,n)))
        pts2 = c.trnps(q23.rotps(vec3(0,0,0).pring(r2,n)))
        uvs = (vec3(0,0,0),vec3(0,1,0),vec3(1,1,0),vec3(1,0,0))
        for x in range(n):
            sq = ((pts1[x-1],pts1[x],pts2[x],pts2[x-1]),())
            m.asurf(sq, gm, uvstacked=(u.cp() for u in uvs))


    def frond(a, b, d, n):
        q12 = quat(1,0,0,0)
        q12.uu(vec3(0,0,1),a.tov(b).cp().nrm())
        pts = b.trnps(q12.rotps(vec3(0,0,0).pring(10,n)))
        for x in range(len(pts)):
            m.asurf(((pts[x-1],pts[x],b),()), gm)
        qt = quat(0,0,0,0).av(pi/2, a.tov(b).crs(b.tov(pts[0])).nrm())
        pts = b.trnps(qt.rotps(q12.rotps(vec3(0,0,0).pring(10,n))))
        for x in range(len(pts)):
            m.asurf(((pts[x-1],pts[x],b),()), gm)


    def walk(v, d=0, n=8):
        d += 1
        a = lsys.above(v)
        if v.term:
            frond(a.p, v.p, d, n)
        bel = lsys.below(v)
        if bel:
            for b in bel:
                branch(a.p, v.p, b.p, d, n)
                walk(b, d=d, n=n)


    totaldepth = lsys.depth()
    radius = lambda d : (totaldepth - d) * r
    branch(lsys.rootp, lsys.rootp, lsys.root.p, 0, 8)
    walk(lsys.root)
    m.normals(gm)


def oak(i, p, d):
    #i,p,d = 4,vec3(0,0,0),vec3(0,0,1)
    axiom,rules = 'X',dict([('X','{^(FQX}{_)FQX}'),('F','AF'),('A','BF'),('B','F')])
    params = dict(
        dpolar=gtl.rad(20), 
        dazimuthal=gtl.rad(15), 
        dtwist=gtl.rad(30), 
        drho=10)
    lsys = lsystem(i, p, d, axiom, rules, grammer=treegrammer, **params)()
    return lsys


def beach(i, p, d):
    #i,p,d = 4,vec3(0,0,0),vec3(0,0,1)
    i = 8
    axiom = 'FX'
    rules = dict([
        ('X','^F{[(_F^)X})]A'),
        ('A','F^F_{)F(X}'), 
            ])
    params = dict(
        dpolar=gtl.rad(20), 
        dazimuthal=gtl.rad(20), 
        dtwist=gtl.rad(20), 
        drho=10)
    lsys = lsystem(i, p, d, axiom, rules, grammer=treegrammer, **params)()
    return lsys
