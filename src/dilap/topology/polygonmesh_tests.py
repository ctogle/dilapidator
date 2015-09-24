import dilap.core.tools as dpr
from dilap.geometry.vec3 import vec3

import dilap.topology.polygonmesh as dmsh

import dilap.mesh.tools as dtl
import matplotlib.pyplot as plt

import unittest,numpy,math

import pdb

#python3 -m unittest discover -v ./ "*tests.py"

class test_polygonmesh(unittest.TestCase):

    def setUp(self):
        self.mesh = dmsh.polygonmesh()

    def test_init(self):
        self.assertEqual(self.mesh.vcnt(),1)
        self.assertEqual(self.mesh.ecnt(),1)
        self.assertEqual(self.mesh.lcnt(),1)
        self.assertEqual(self.mesh.fcnt(),1)

    def test_avert(self):
        self.mesh.avert(self.mesh.verts[0].p.cp())
        self.assertEqual(self.mesh.vcnt(),2)
        self.assertEqual(self.mesh.ecnt(),1)
        self.assertEqual(self.mesh.lcnt(),1)
        self.assertEqual(self.mesh.fcnt(),1)

    def test_sedge(self):
        v2,e1,e2 = self.mesh.sedge(self.mesh.edges[0],vec3(2,0,0))
        v3,e2,e3 = self.mesh.sedge(e2,vec3(0,2,0))
        self.assertEqual(self.mesh.vcnt(),3)
        self.assertEqual(self.mesh.ecnt(),3)
        self.assertEqual(self.mesh.lcnt(),1)
        self.assertEqual(self.mesh.fcnt(),1)

    def test_mask(self):
        tr = dmsh.polygonmesh()
        print('mask woot for polygonmesh!')

if __name__ == '__main__':
    unittest.main()









 


