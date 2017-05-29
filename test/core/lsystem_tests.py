import dilap.core.plotting as dtl
import dilap.core.lsystem as lsy
#import dilap.worldly.treeskin as ltr
from dilap.geometry.vec3 import vec3
import dilap.geometry.tools as gtl
import dilap.geometry.polymath as pym
import dilap.topology.planargraph as dpg

import matplotlib.pyplot as plt
import itertools as it
import unittest,numpy,math

import pdb

#python3 -m unittest discover -v ./ "*tests.py"

###############################################################################

def realize(i,p,d,axiom,rules,**params):
    ax = dtl.plot_axes(10)
    for piece in lsy.lgen(p,d,axiom,rules,i,**params):
        if isinstance(piece,tuple):
            ax = dtl.plot_edges(piece,ax,lw = 2,col = 'k')
        elif isinstance(piece,vec3):
            ax = dtl.plot_point(piece,ax)
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
        params = dict(dazimuthal = gtl.rad(25.7))
        realize(i,p,d,axiom,rules,**params)


class axialmtn(unittest.TestCase):

    keepers = {
        'lake':((vec3(0,0,0),vec3(0,1,0),'X',
                dict([('X','{[[X}{]X}F]X'),('F','FA'),('A','F')]),6),
                dict(dazimuthal = gtl.rad(25.7),drho = 20)),
            }

    def pg(self,*ags,**kws):
        pg = dpg.planargraph()
        for piece in lsy.lgen(*ags,**kws):
            if isinstance(piece,tuple):
                p1,p2 = piece
                v1,v2 = pg.fp(p1,10),pg.fp(p2,10)
                e12 = pg.fe(v1,v2)
            elif isinstance(piece,vec3):pass
        py = pym.pgtopy(pg,1)[0]
        py = pym.bisectb(py)
        py = pym.smoothxy(py,0.5,2)

        #ax = pg.plotxy(l = 20)
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
        
        pg = dpg.planargraph()
        for piece in lsy.lgen(p,d,axiom,rules,i,**params):
            if isinstance(piece,tuple):
                p1,p2 = piece
                v1,v2 = pg.fp(p1,10),pg.fp(p2,10)
                e12 = pg.fe(v1,v2)
            elif isinstance(piece,vec3):pass

        py = pym.pgtopy(pg,1)[0]
        py = pym.bisectb(py)
        py = pym.smoothxy(py,0.5,2)
        #py = pym.aggregate(py,2)

        #ax = pg.plotxy(l = 20)
        ax = dtl.plot_axes_xy(20)
        ax = dtl.plot_polygon_xy(py,ax,lw = 3)
        plt.show()

if __name__ == '__main__':
    unittest.main()

'''#
l,w,n = 25,25,3
ax = dtl.plot_axes(50)
pd = it.product(range(n),range(n),range(1))
pstack = [vec3(*p).scl(vec3(l,w,1)) for p in pd]
dstack = [vec3(0,0,1) for x in range(len(pstack))]

class ___KEEPTOGETEXAMPLEStest_lsystem(unittest.TestCase):

    def atest_ltree_loadout(self):
        p,d = pstack.pop(0),dstack.pop(0)
        l = ltr.ltree(0)._realize(p,d,ax)
        p,d = pstack.pop(0),dstack.pop(0)
        l = ltr.ltree(1)._realize(p,d,ax)

    def test_ltree(self):
        ''#
        kws = {
            'axiom':'F',
            #'axiom':'G',
            'rules':[
                ('F','F=![+++++F][-------F]-![++++F][------F]-![+++F][-----F]-!F')],
                #('G','<FF[)}+\\FG][)-{G]')],
            'iterations':1,'seed':0,'polar':gtl.rad(10),'azimuthal':gtl.rad(55),
                }
        ''#
        ldots = [
                {
            'axiom':'F',
            'rules':[
                #('F','F=![+++++F][-------F]-![++++F][------F]-![+++F][-----F]-!F')],
                ('F','F[+F]F[-F]F')],
            'iterations':3,'seed':0,'polar':gtl.rad(10),'azimuthal':gtl.rad(25),
                },
                {
            'axiom':'G',
            'rules':[
                ('G','<FF[)}+\\FG][)-{G]')],
            'iterations':5,'seed':0,'polar':gtl.rad(10),'azimuthal':gtl.rad(55),
                },
                {
            'axiom':'G',
            'rules':[
                ('G','FG[<}GF][))-{FG][F]')],
            'iterations':3,'seed':0,'polar':gtl.rad(20),'azimuthal':gtl.rad(75),
                },
                {
            'axiom':'Q',
            'rules':[
                ('Q','FQ[<}Q][))-{FQ][F]')],
            'iterations':3,'seed':0,'polar':gtl.rad(20),'azimuthal':gtl.rad(75),
                },
                {
            'axiom':'A',
            'rules':[
                ('F','F~[~FA]F'),
                ('A','[--////FQA]-<F[>>\\\\FQA]>+F[-\>>QA\[X]][</++QA/[Y]]Q'),
                ('X','F+Q'),('Y','F-Q')],
            'iterations':3,'seed':0,'polar':gtl.rad(20),'azimuthal':gtl.rad(75),
                }
                ]
        for x in range(len(ldots)):
            p,d = pstack.pop(0),dstack.pop(0)
            kws = ldots[x]
            l = ltr.ltree(0,**kws)._realize(p,d,ax)

        #p,d = pstack.pop(0),dstack.pop(0)
        #l = ltr.ltree(0,**kws)._realize(p,d,ax)

    def atest_treeskin(self):
        p,d = pstack.pop(0),dstack.pop(0)
        l = ltr.treeskin(p,d,ax = ax)
'''#
