from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat

import dilap.geometry.tools as gtl

#import dilap.mesh.tools as dtl
import dilap.core.plotting as dtl
import matplotlib.pyplot as plt

import unittest,numpy

#python3 -m unittest discover -v ./ "*tests.py"

class test_tools(unittest.TestCase):

    def setUp(self):
        self.pslist = [
            vec3(-1,-1,0),vec3( 1,-1,0),
            vec3( 1, 1,0),vec3(-1, 1,0)]
        self.pstupl = (
            vec3(-1,-1,0),vec3( 1,-1,0),
            vec3( 1, 1,0),vec3(-1, 1,0))

    def test_isnear(self):
        self.assertEqual(gtl.isnear(0.0001,0),1)
        self.assertEqual(gtl.isnear(0.1,0),0)

    def test_near(self):
        self.assertEqual(gtl.near(0.0001,0),0)
        self.assertEqual(gtl.isnear(gtl.near(0.1,0),0.1),1)

    def test_rad(self):
        self.assertEqual(gtl.isnear(gtl.rad(180),gtl.PI),1)

    def test_deg(self):
        self.assertEqual(gtl.isnear(gtl.deg(gtl.PI),180),1)

    def test_clamp(self):
        self.assertEqual(gtl.clamp(180,0,90),90)
        self.assertEqual(gtl.clamp(-90,0,90), 0)
        self.assertEqual(gtl.clamp( 45,0,90),45)

    def test_wrap(self):
        self.assertEqual(gtl.wrap(120,0,90),30)
        self.assertEqual(gtl.wrap(-20,0,90),70)
        self.assertEqual(gtl.wrap( 45,0,90),45)

    def test_inrng(self):
        self.assertEqual(gtl.inrng(   0,0,1),0)
        self.assertEqual(gtl.inrng(   1,0,1),0)
        self.assertEqual(gtl.inrng( 0.1,0,1),1)
        self.assertEqual(gtl.inrng( 0.9,0,1),1)
        self.assertEqual(gtl.inrng(-0.1,0,1),0)
        self.assertEqual(gtl.inrng( 1.1,0,1),0)

    def test_adist(self):
        deg10 = gtl.rad(10)
        self.assertEqual(gtl.isnear(gtl.adist(deg10*2,deg10*6),deg10*4),1)
        self.assertEqual(gtl.isnear(gtl.adist(deg10*6,deg10*2),deg10*4),1)
        self.assertEqual(gtl.isnear(gtl.adist(deg10*6,deg10*22),deg10*16),1)
        self.assertEqual(gtl.isnear(gtl.adist(deg10*6,deg10*32),deg10*10),1)

    def test_orient2d(self):
        p1,p2,p3,p4 = vec3(1,1,0),vec3(0,1,0),vec3(0,0,0),vec3(0,2,0)
        self.assertTrue(gtl.orient2d(p1,p2,p3) > 0)
        self.assertTrue(gtl.orient2d(p2,p1,p3) < 0)
        self.assertTrue(gtl.orient2d(p2,p4,p3) == 0)

if __name__ == '__main__':
    unittest.main()



