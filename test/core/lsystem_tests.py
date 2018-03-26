import dilap.core.plotting as dtl
import dilap.core.lsystem as lsy
#import dilap.worldly.treeskin as ltr
from dilap.geometry.vec3 import vec3
import dilap.geometry.tools as gtl
import dilap.geometry.polymath as pym
import dilap.geometry.planargraph as dpg

import matplotlib.pyplot as plt
import itertools as it
import unittest,numpy,math

import pdb

######################################################
### python3 -m unittest discover -v ./ "*tests.py" ###
######################################################

def realize(i, p, d, axiom, rules, **params):
    ax = dtl.plot_axes(10)
    for piece in lsy(i, p, d, axiom, rules, **params):
        if isinstance(piece, tuple):
            ax = dtl.plot_edges(piece, ax, lw=2, col='k')
        elif isinstance(piece, vec3):
            ax = dtl.plot_point(piece, ax)
    plt.show()

class lsystem(unittest.TestCase):

    def test_pythagorastree(self):
        i = 5
        p,d = vec3(0,0,0),vec3(1,0,0)
        axiom = 'Q'
        rules = dict([('F','FF'),('Q','F{[Q}]Q')])
        params = dict(dazimuthal = numpy.pi/12)
        realize(i,p,d,axiom,rules,**params)

    def test_dragoncurve(self):
        i = 9
        p,d = vec3(0,0,0),vec3(0,1,0)
        axiom = 'FX'
        rules = dict([('X','X[YF['),('Y',']FX]Y')])
        params = dict(dazimuthal = gtl.rad(90))
        realize(i,p,d,axiom,rules,**params)

    def test_axialtree(self):
        i = 5
        p,d = vec3(0,0,0),vec3(0,1,0)
        axiom = 'X'
        rules = dict([('X','F{[X}{]X}FX'),('F','FF')])
        params = dict(dazimuthal = gtl.rad(25.7))
        realize(i,p,d,axiom,rules,**params)

    def test_plant(self):
        i = 3
        p,d = vec3(0,0,0),vec3(0,1,0)
        axiom = 'X'
        rules = dict([('X','F]{{X}[X}[F{[FX}]X'),('F','FF')])
        params = dict(dazimuthal = gtl.rad(25))
        realize(i,p,d,axiom,rules,**params)

    def test_grass(self):
        i = 5
        p,d = vec3(0,0,0),vec3(0,1,0)
        axiom = 'F'
        rules = dict([('F','F{[F}F{]F}F')])
        params = dict(dazimuthal=gtl.rad(25.7))
        realize(i,p,d,axiom,rules,**params)

    def aatest_trees(self):
        i = 5
        p, d = vec3(0,0,0), vec3(0,1,0)
        axiom = 'A'
        rules = dict((
            ('F','F~[~FA]F'),
            ('A','[--////FQA]-<F[>>\\\\FQA]>+F[-\>>QA\[X]][</++QA/[Y]]Q'),
            ('X','F+Q'),
            ('Y','F-Q')))
        params = dict(dazimuthal=gtl.rad(25.7))
        realize(i,p,d,axiom,rules,**params)
        #loadouts.append(('Q',[
        #    ('Q','<FF[)}+FQ][){Q]')]))

class axialmtn(unittest.TestCase):

    keepers = {
        'lake':(
                (6, vec3(0,0,0), vec3(0,1,0),
                    'X', dict([('X', '{[[X}{]X}F]X'), ('F', 'FA'), ('A', 'F')])
                ),
                dict(dazimuthal=gtl.rad(25.7), drho=20)),
            }

    def pg(self,*ags,**kws):
        pg = dpg()
        for piece in lsy(*ags,**kws):
            if isinstance(piece,tuple):
                p1,p2 = piece
                v1,v2 = pg.fp(p1,10),pg.fp(p2,10)
                e12 = pg.fe(v1,v2)
            elif isinstance(piece,vec3):
                pass

        py = pg.polygon(1)[0]
        py = pym.bisectb(py)
        py = pym.smoothxy(py,0.5,2)

        ax = dtl.plot_axes_xy(50)
        ax = dtl.plot_polygon_xy(py,ax,lw = 3)
        plt.show()

    def test_lake(self):
        lakeags,params = self.keepers['lake']
        self.pg(*lakeags,**params)

    def aaatest(self):
        #i = 5
        i = 6
        p,d = vec3(0,0,0),vec3(0,1,0)
        axiom = 'X'
        rules = dict([('X','{[[X}{]X}F]X'),('F','FA'),('A','F')])
        #rules = dict([('X','F{[X}{]X}FX'),('F','FF')])
        #axiom = 'F'
        #rules = dict([('F','F{[F}F{]F}F')])
        params = dict(dazimuthal = gtl.rad(25.7),drho = 20)

        pg = dpg()
        for piece in lsy(i, p, d, axiom, rules, **params):
            if isinstance(piece,tuple):
                p1,p2 = piece
                v1,v2 = pg.fp(p1,10),pg.fp(p2,10)
                e12 = pg.fe(v1,v2)
            elif isinstance(piece,vec3):pass

        py = pym.pgtopy(pg,1)[0]
        py = pym.bisectb(py)
        py = pym.smoothxy(py,0.5,2)
        #py = pym.aggregate(py,2)

        ax = dtl.plot_axes_xy(20)
        ax = dtl.plot_polygon_xy(py,ax,lw = 3)
        plt.show()

if __name__ == '__main__':
    unittest.main()
