from dilap.core import *
from dilap.topology import tree
import dilap.geometry.tools as gtl
import pdb


def treeskin(m, gm, lsys, n=8, r=0.5):


    def branch(a, b, c, d, n):
        q12,q23 = quat(1,0,0,0),quat(1,0,0,0)
        if a:
            q12.uu(vec3(0,0,1),a.p.tov(b.p).cp().nrm())
        q23.uu(vec3(0,0,1),b.p.tov(c.p).cp().nrm())
        r1 = radius(d)
        r2 = radius(d+1)
        pts1 = b.p.trnps(q12.rotps(vec3(0,0,0).pring(r1,n)))
        pts2 = c.p.trnps(q23.rotps(vec3(0,0,0).pring(r2,n)))
        uvs = (vec3(0,0,0),vec3(0,1,0),vec3(1,1,0),vec3(1,0,0))
        for x in range(n):
            sq = ((pts1[x-1],pts1[x],pts2[x],pts2[x-1]),())
            m.asurf(sq, gm, uvstacked=(u.cp() for u in uvs))


    def frond(a, b, d, n):
        q12 = quat(1,0,0,0)
        q12.uu(vec3(0,0,1),a.p.tov(b.p).cp().nrm())
        pts = b.p.trnps(q12.rotps(vec3(0,0,0).pring(10,n)))
        for x in range(len(pts)):
            m.asurf(((pts[x-1],pts[x],b.p),()), gm)
        qt = quat(0,0,0,0).av(pi/2, a.p.tov(b.p).crs(b.p.tov(pts[0])).nrm())
        pts = b.p.trnps(qt.rotps(q12.rotps(vec3(0,0,0).pring(10,n))))
        for x in range(len(pts)):
            m.asurf(((pts[x-1],pts[x],b.p),()), gm)


    def walk(v, d=0, n=8):
        d += 1
        a = lsys.above(v)
        if v.term:
            frond(a, v, d, n)
        bel = lsys.below(v)
        if bel:
            for b in bel:
                branch(a, v, b, d, n)
                walk(b, d=d, n=n)


    totaldepth = lsys.depth()
    radius = lambda d : (totaldepth - d) * r
    walk(lsys.root)
    m.normals(gm)


def oak():
    i,p,d = 4,vec3(0,0,0),vec3(0,0,1)
    axiom,rules = 'X',dict([('X','{^(FQX}{_)FQX}'),('F','AF'),('A','BF'),('B','F')])
    params = dict(
        dpolar=gtl.rad(20), 
        dazimuthal=gtl.rad(15), 
        dtwist=gtl.rad(30), 
        drho=10)
    lsys = lsystem(i, p, d, axiom, rules, **params)()
    return lsys


def beach():
    i,p,d = 4,vec3(0,0,0),vec3(0,0,1)
    axiom,rules = 'X',dict([('X','F{^^(F_QX}{_)FQX}'),('F','AF'),('A','BF'),('B','F')])
    params = dict(
        dpolar=gtl.rad(20), 
        dazimuthal=gtl.rad(15), 
        dtwist=gtl.rad(30), 
        drho=10)
    lsys = lsystem(i, p, d, axiom, rules, **params)()
    return lsys
