import dilap.core.plotting as dtl
import dilap.core.lsystem as lsy
import dilap.topology.worldly.treeskin as ltr
from dilap.geometry.vec3 import vec3
import dilap.geometry.tools as gtl

import matplotlib.pyplot as plt
import itertools as it
import unittest,numpy,math

import pdb

#python3 -m unittest discover -v ./ "*tests.py"

###############################################################################

l,w,n = 25,25,3
ax = dtl.plot_axes(50)
pd = it.product(range(n),range(n),range(1))
pstack = [vec3(*p).scl(vec3(l,w,1)) for p in pd]
dstack = [vec3(0,0,1) for x in range(len(pstack))]

class test_lsystem(unittest.TestCase):

    def tearDownClass():plt.show()
    #def setUpClass():

    def atest_pythagoras_tree(self):
        p,d = pstack.pop(0),dstack.pop(0)
        l = lsy.pythagoras_tree()._realize(p,d,ax)

    def atest_dragon_curve(self):
        p,d = pstack.pop(0),dstack.pop(0)
        l = lsy.dragon_curve()._realize(p,d,ax)

    def atest_axial_tree(self):
        p,d = pstack.pop(0),dstack.pop(0)
        l = lsy.axial_tree()._realize(p,d,ax)

    def atest_plant(self):
        p,d = pstack.pop(0),dstack.pop(0)
        l = lsy.plant()._realize(p,d,ax)

    def atest_grass(self):
        kws = {
            'axiom':'F',
            'rules':[('F','F[+F]F[-F]F')],
            'iterations':5,'seed':0,
            'angle':gtl.rad(25.7),
                }
        p,d = pstack.pop(0),dstack.pop(0)
        l = lsy.lsystem(**kws)._realize(p,d,ax)

    def atest_ltree_loadout(self):
        p,d = pstack.pop(0),dstack.pop(0)
        l = ltr.ltree(0)._realize(p,d,ax)

    def test_ltree(self):
        kws = {
            'axiom':'G',
            'rules':[
                #('F','F=![+++++F][-------F]-![++++F][------F]-![+++F][-----F]-!F'),
                #('F','F{{F'),
                ('G','FF[)}FG][){G]')],
            'iterations':1,'seed':0,
            'polar':gtl.rad(12),
            'azimuthal':gtl.rad(75),
                }
        for x in range(5):
            kws['iterations'] += 1
            p,d = pstack.pop(0),dstack.pop(0)
            l = ltr.ltree(0,**kws)._realize(p,d,ax)
        #p,d = pstack.pop(0),dstack.pop(0)
        #l = ltr.ltree(0,**kws)._realize(p,d,ax)

    def atest_treeskin(self):
        p,d = pstack.pop(0),dstack.pop(0)
        l = ltr.treeskin(p,d,ax = ax)

###############################################################################

if __name__ == '__main__':
    unittest.main()

###############################################################################





