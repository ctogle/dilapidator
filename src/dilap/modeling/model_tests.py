import dilap.core.tools as dpr
from dilap.geometry.vec3 import vec3
from dilap.geometry.pointset import pointset

import dilap.topology.trimesh as tmsh

import dilap.mesh.tools as dtl
import matplotlib.pyplot as plt

import unittest,numpy,math

import pdb

#python3 -m unittest discover -v ./ "*tests.py"

class test_trimesh(unittest.TestCase):

    def avert(self,p):
        px = self.pset.ap(p)
        nx = self.nset.ap(vec3(0,0,1))
        ux = self.uset.ap(vec3(0,0,0))
        return self.mesh.avert(px,nx,ux)

    def quad(self):
        self.v1 = self.avert(vec3(-1,-1,-1))
        self.v2 = self.avert(vec3( 1,-1,-1))
        self.v3 = self.avert(vec3( 1, 1,-1))
        self.v4 = self.avert(vec3(-1, 1,-1))
        self.f1 = self.mesh.aface(self.v1,self.v2,self.v3) 
        self.f2 = self.mesh.aface(self.v1,self.v3,self.v4) 

    def cube(self):
        self.v1  = self.avert(vec3(-1,-1,-1))
        self.v2  = self.avert(vec3( 1,-1,-1))
        self.v3  = self.avert(vec3( 1, 1,-1))
        self.v4  = self.avert(vec3(-1, 1,-1))
        self.v5  = self.avert(vec3(-1,-1, 1))
        self.v6  = self.avert(vec3( 1,-1, 1))
        self.v7  = self.avert(vec3( 1, 1, 1))
        self.v8  = self.avert(vec3(-1, 1, 1))
        self.f1  = self.mesh.aface(self.v1,self.v3,self.v2) 
        self.f2  = self.mesh.aface(self.v1,self.v4,self.v3) 
        self.f3  = self.mesh.aface(self.v5,self.v6,self.v7) 
        self.f4  = self.mesh.aface(self.v5,self.v7,self.v8) 
        self.f5  = self.mesh.aface(self.v1,self.v2,self.v5) 
        self.f6  = self.mesh.aface(self.v1,self.v5,self.v4) 
        self.f7  = self.mesh.aface(self.v2,self.v3,self.v6) 
        self.f8  = self.mesh.aface(self.v2,self.v6,self.v5) 
        self.f9  = self.mesh.aface(self.v3,self.v4,self.v7) 
        self.f10 = self.mesh.aface(self.v3,self.v7,self.v6) 
        self.f11 = self.mesh.aface(self.v4,self.v1,self.v8) 
        self.f12 = self.mesh.aface(self.v4,self.v8,self.v7) 

    def assert_counts(self,vcnt,ecnt,fcnt):
        self.assertEqual(self.mesh.vcnt(),vcnt)
        self.assertEqual(self.mesh.ecnt(),ecnt)
        self.assertEqual(self.mesh.fcnt(),fcnt)

    def setUp(self):
        self.mesh = tmsh.trimesh()
        self.pset = pointset()
        self.nset = pointset()
        self.uset = pointset()

    def test_init(self):
        self.assert_counts(0,0,0)

    def test_avert(self):
        self.v1 = self.avert(vec3(-1,-1,-1))
        self.v2 = self.avert(vec3( 1,-1,-1))
        self.v3 = self.avert(vec3( 1, 1,-1))
        self.v4 = self.avert(vec3(-1, 1,-1))
        self.assert_counts(4,0,0)

    #def test_rvert(self):
    #def test_aedge(self):
    #def test_redge(self):

    def test_aface(self):
        self.v1 = self.avert(vec3(-1,-1,-1))
        self.v2 = self.avert(vec3( 1,-1,-1))
        self.v3 = self.avert(vec3( 1, 1,-1))
        self.v4 = self.avert(vec3(-1, 1,-1))
        fx = self.mesh.aface(self.v1,self.v2,self.v3)
        self.assert_counts(4,3,1)

    #def test_rface(self):
    #def test_sface(self):

    def test_mask_v1(self):
        self.quad()
        v10m = self.mesh.mask(0,self.mesh.verts[self.v1],None,None)
        self.assertTrue(self.mesh.verts[self.v2] in v10m)
        self.assertTrue(self.mesh.verts[self.v3] in v10m)
        self.assertTrue(self.mesh.verts[self.v4] in v10m)
        v11m = self.mesh.mask(1,self.mesh.verts[self.v1],None,None)
        self.assertTrue((0,1) in v11m)
        self.assertTrue((3,0) in v11m)
        self.assertTrue((0,2) in v11m)
        self.assertTrue((2,0) in v11m)
        v12m = self.mesh.mask(2,self.mesh.verts[self.v1],None,None)
        self.assertTrue((0,1,2) in v12m)
        self.assertTrue((0,2,3) in v12m)
        self.assertEqual(len(v10m),3)
        self.assertEqual(len(v11m),4)
        self.assertEqual(len(v12m),2)

    def test_mask_v4(self):
        self.quad()
        v40m = self.mesh.mask(0,self.mesh.verts[self.v4],None,None)
        self.assertTrue(self.mesh.verts[self.v1] in v40m)
        self.assertFalse(self.mesh.verts[self.v2] in v40m)
        self.assertTrue(self.mesh.verts[self.v3] in v40m)
        v41m = self.mesh.mask(1,self.mesh.verts[self.v4],None,None)
        self.assertTrue((3,0) in v41m)
        self.assertTrue((2,3) in v41m)
        v42m = self.mesh.mask(2,self.mesh.verts[self.v4],None,None)
        self.assertTrue((0,2,3) in v42m)
        self.assertEqual(len(v40m),2)
        self.assertEqual(len(v41m),2)
        self.assertEqual(len(v42m),1)

    #def test_vonb(self):
    #def test_eonb(self):

    def test_quad(self):
        self.quad()
        self.assert_counts(4,6,2)

    def test_cube(self):
        self.cube()
        self.assert_counts(8,36,12)

if __name__ == '__main__':
    unittest.main()








 

