import dilap.core.base as db

from dilap.geometry.vec3 import vec3
import dilap.geometry.tools as gtl
import dilap.core.lsystem as lsy
#import dilap.topology.worldly.treeskin as ltr
#import dilap.topology.worldly.building as blg
import dilap.worldly.world as dwo

import dilap.core.plotting as dtl

import matplotlib.pyplot as plt
import itertools as it
import unittest,numpy,math

import pdb

#python3 -m unittest discover -v ./ "*tests.py"

###############################################################################

l,w,n = 25,25,3
#ax = dtl.plot_axes(50)
pd = it.product(range(n),range(n),range(1))
pstack = [vec3(*p).scl(vec3(l,w,1)) for p in pd]
dstack = [vec3(0,0,1) for x in range(len(pstack))]

class test_context(unittest.TestCase):

    #def tearDownClass():plt.show()
    #def setUpClass():

    def atest_grass(self):
        kws = {
            'axiom':'F',
            'rules':[('F','F[+F]F[-F]F')],
            'iterations':3,'seed':0,
            'angle':gtl.rad(25.7),
                }
        p,d = pstack.pop(0),dstack.pop(0)
        cx = ltr.tree(p,d,ax = ax)
        cx.display()
        #plt.show()

    def atest_building(self):
        cx = blg.building().generate()

    def test_world(self):
        #cx = dwo.world().generate()
        wfa = dwo.worldfactory()
        wcx = wfa.new()
        #db.profile_function(wfa.new,*(),**{})

###############################################################################

if __name__ == '__main__':
    unittest.main()

###############################################################################





