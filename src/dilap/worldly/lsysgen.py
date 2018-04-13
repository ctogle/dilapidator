from dilap.core import lsystem
from dilap.geometry import vec3, planargraph
from dilap.geometry.tools import rad, deg

def lsys_to_pg(lsys):
    pg = planargraph()
    for piece in lsys:
        if isinstance(piece, tuple):
            p1, p2 = piece
            v1 = pg.fp(p1, 10, w=0.5)
            v2 = pg.fp(p2, 10, w=0.5)
            e12 = pg.fe(v1, v2)
        elif isinstance(piece, vec3):
            pass
        else:
            raise ValueError('unknown piece type: "%s"' % type(piece))
    return pg

def grass_lsys(i=3, o=None, d=None):
    o = vec3(0, 0, 0) if o is None else o
    d = vec3(0, 1, 0) if d is None else d
    axiom = 'X' 
    rules = dict([('X','F]{{X}[X}[F{[FX}]X'),('F','FF')])
    params = dict(dazimuthal=rad(25), drho=20)
    return lsystem(i, o, d, axiom, rules, **params)
    
def tight_sprawl(i=6, o=None, d=None):
    o = vec3(0, 0, 0) if o is None else o
    d = vec3(0, 1, 0) if d is None else d
    axiom = 'X'
    rules = dict([('X','{[[X}{]X}F]X'),('F','FA'),('A','F')])
    params = dict(dazimuthal=rad(25.7), drho=20)
    return lsystem(i, o, d, axiom, rules, **params)
    
def loose_sprawl(i=3, o=None, d=None):
    o = vec3(0, 0, 0) if o is None else o
    d = vec3(0, 1, 0) if d is None else d
    axiom = 'X'
    rules = dict([('X','F]{{X}[X}[F{[FX}]X'),('F','FF')])
    params = dict(dazimuthal=rad(25), drho=20)
    return lsystem(i, o, d, axiom, rules, **params)



def waterway(b):
    i = 3
    p,d = vec3(0,0,0),vec3(0,1,0)
    axiom = 'X'
    rules = dict([('X','F]{{X}[X}[F{[FX}]X'),('F','FF')])
    params = dict(dazimuthal=gtl.rad(25), drho=32)
    lsys = lsystem(i, p, d, axiom, rules, **params)
    pg = planargraph()
    for piece in lsys:
        if isinstance(piece, tuple):
            p1, p2 = piece
            v1 = pg.fp(p1, 16, w=0.5)
            v2 = pg.fp(p2, 16, w=0.5)
            e12 = pg.fe(v1, v2)
        elif isinstance(piece, vec3):
            pass
    pg.smooth_sticks()
    pg.fitbxy(b, 1.1)
    return pg


