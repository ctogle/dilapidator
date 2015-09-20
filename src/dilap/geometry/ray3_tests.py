from dilap.geometry.vec3 import vec3
from dilap.geometry.ray3 import ray3

import dilap.core.tools as dpr

import dilap.mesh.tools as dtl
import matplotlib.pyplot as plt

import unittest,numpy,math

#python3 -m unittest discover -v ./ "*tests.py"

class test_ray3(unittest.TestCase):

    def test_cp(self):
        r1 = ray3(vec3(1,1,2),vec3(0,1,0))
        r2 = ray3(vec3(1,1,2),vec3(0,1,0))
        self.assertTrue(r1 is r1)
        self.assertFalse(r1 is r1.cp())
        self.assertTrue(r1 == r1.cp())
        self.assertFalse(r1 is r2)
        self.assertFalse(r1 is r2.cp())
        self.assertTrue(r1 == r2.cp())

    def test_isnear(self):
        r1 = ray3(vec3(1,1,2),vec3(0,1,0))
        r2 = ray3(vec3(1,1,2),vec3(0,5,0))
        r3 = ray3(vec3(1,1,2),vec3(1,1,0))
        r4 = ray3(vec3(2,1,2),vec3(0,1,0))
        self.assertTrue(r1.isnear(r1))
        self.assertTrue(r1.isnear(r1.cp()))
        self.assertTrue(r1.isnear(r2))
        self.assertFalse(r1.isnear(r3))
        self.assertFalse(r1.isnear(r4))

    def test_hittri(self):
        r1 = ray3(vec3(1,1,2),vec3(0,1,0))
        tri1 = (vec3(0,2,0),vec3(3,2,1),vec3(3,3,4))
        tri2 = (vec3(0,2,0),vec3(3,2,1),vec3(2,3,4))
        tri3 = (vec3(0,2,0),vec3(3,2,0),vec3(1,3,4))
        tri4 = (vec3(0,2,0),vec3(3,2,1),vec3(1,3,4))
        self.assertFalse(r1.hittri(*tri1))
        self.assertTrue(r1.hittri(*tri2))
        self.assertTrue(r1.hittri(*tri3))
        self.assertTrue(r1.hittri(*tri4))

    def test_hittris(self):
        r1 = ray3(vec3(1,1,2),vec3(0,1,0))
        tri1 = (vec3(0,2,0),vec3(3,2,1),vec3(3,3,4))
        tri2 = (vec3(0,2,0),vec3(3,2,1),vec3(2,3,4))
        tri3 = (vec3(0,2,0),vec3(3,2,0),vec3(1,3,4))
        tri4 = (vec3(0,2,0),vec3(3,2,1),vec3(1,3,4))
        self.assertEqual(r1.hittris([tri1,tri2,tri3,tri4]),[1,2,3])

    def test_hittri_close(self):
        r1 = ray3(vec3(1,1,2),vec3(0,1,0))
        tri1 = (vec3(0,3,0),vec3(3,2,1),vec3(3,3,4))
        tri2 = (vec3(0,2,0),vec3(3,4,1),vec3(2,4,4))
        tri3 = (vec3(0,2,0),vec3(3,2,0),vec3(1,1,4))
        tri4 = (vec3(0,2,0),vec3(3,3,1),vec3(1,3,4))
        self.assertEqual(r1.hittri_close([tri1,tri2,tri3,tri4]),2)

    def test_hitpln(self):
        r1 = ray3(vec3(1,1,2),vec3(0,1,0))
        p01,p02 = vec3(-2,6,12),vec3(3,-2,0)
        pn1,pn2,pn3 = vec3(0,1,0),vec3(0,-1,0),vec3(1,0,0)
        self.assertTrue(r1.hitpln(p01,pn1))
        self.assertTrue(r1.hitpln(p01,pn2))
        self.assertFalse(r1.hitpln(p01,pn3))
        self.assertFalse(r1.hitpln(p02,pn1))
        self.assertFalse(r1.hitpln(p02,pn2))
        self.assertFalse(r1.hitpln(p02,pn3))

if __name__ == '__main__':
    unittest.main()









 

