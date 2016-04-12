from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat

import dilap.geometry.triangulate as dtg

#import dilap.mesh.tools as dtl
import dilap.core.plotting as dtl
import matplotlib.pyplot as plt

import unittest,numpy,pdb

#python3 -m unittest discover -v ./ "*tests.py"

class test_triangulate(unittest.TestCase):

    def setUp(self):pass

    def test_triangulate(self):
        eb = (vec3(-2,-2,0),vec3(2,-2,0),vec3(2,2,0),vec3(-2,2,0))
        ibs = ()
        hmin,ref,smo = 1,False,False

        tris,bnds = dtg.triangulate(eb,ibs,hmin,ref,smo)

        pdb.set_trace()

if __name__ == '__main__':
    unittest.main()



