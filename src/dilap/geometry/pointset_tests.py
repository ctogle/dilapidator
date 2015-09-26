from dilap.geometry.quat import quat
from dilap.geometry.vec3 import vec3
from dilap.geometry.pointset import pointset

import dilap.core.tools as dpr

import dilap.mesh.tools as dtl
import matplotlib.pyplot as plt

import unittest,numpy,math

import pdb





#python3 -m unittest discover -v ./ "*tests.py"

class test_pointset(unittest.TestCase):

    def setUp(self):
        self.pset = pointset()
        self.pset.ap(vec3(0,0,0))
        self.pset.ap(vec3(1,0,0))
        self.pset.ap(vec3(1,1,0))
        self.pset.ap(vec3(1,1,1))
        self.pset.ap(vec3(0,1,1))
        self.pset.ap(vec3(0,0,1))
        self.pset.ap(vec3(1,0,1))
        self.pset.ap(vec3(0,1,0))

    def test_gpscp(self):
        gps = self.pset.gpscp((1,3,4))
        self.assertEqual(gps[0],vec3(1,0,0))
        self.assertEqual(gps[0],self.pset.ps[1])
        self.assertFalse(gps[0] is self.pset.ps[1])
        self.assertEqual(gps[1],vec3(1,1,1))
        self.assertEqual(gps[1],self.pset.ps[3])
        self.assertFalse(gps[1] is self.pset.ps[3])
        self.assertEqual(gps[2],vec3(0,1,1))
        self.assertEqual(gps[2],self.pset.ps[4])
        self.assertFalse(gps[2] is self.pset.ps[4])

    def test_gps(self):
        gps = self.pset.gps((1,3,4))
        self.assertEqual(gps[0],vec3(1,0,0))
        self.assertEqual(gps[0],self.pset.ps[1])
        self.assertTrue(gps[0] is self.pset.ps[1])
        self.assertEqual(gps[1],vec3(1,1,1))
        self.assertEqual(gps[1],self.pset.ps[3])
        self.assertTrue(gps[1] is self.pset.ps[3])
        self.assertEqual(gps[2],vec3(0,1,1))
        self.assertEqual(gps[2],self.pset.ps[4])
        self.assertTrue(gps[2] is self.pset.ps[4])

    def test_ap(self):
        pset = pointset()
        np = vec3(0,4,0)
        x = pset.ap(np)
        self.assertEqual(pset.ps[x],np)
        self.assertTrue(pset.ps[x] is np)
        np = vec3(2,-1,3.2)
        x = pset.ap(np)
        self.assertEqual(pset.ps[x],np)
        self.assertTrue(pset.ps[x] is np)

    def test_aps(self):
        np1 = vec3(0,4,0)
        np2 = vec3(2,-1,3.2)
        nps = [np1,np2]
        nxs = self.pset.aps(nps)
        for x in range(len(nxs)):
            nx = nxs[x]
            self.assertEqual(self.pset.ps[nx],nps[x])
            self.assertTrue(self.pset.ps[nx] is nps[x])

    def test_np(self):
        np = vec3(1,1,0)
        x = self.pset.np(np)
        self.assertEqual(x,2)
        self.assertFalse(self.pset.ps[x] is np)
        self.assertEqual(self.pset.ps[x],np)
        np = vec3(0,1,2)
        x = self.pset.np(np)
        self.assertEqual(x,8)
        self.assertTrue(self.pset.ps[x] is np)
        self.assertEqual(self.pset.ps[x],np)

    def test_nps(self):
        np1 = vec3(1,1,0)
        np2 = vec3(0,1,2)
        x,y = self.pset.nps([np1,np2])
        self.assertEqual(x,2)
        self.assertEqual(y,8)
        self.assertFalse(self.pset.ps[x] is np1)
        self.assertEqual(self.pset.ps[x],np1)
        self.assertTrue(self.pset.ps[y] is np2)
        self.assertEqual(self.pset.ps[y],np2)

    def test_fp(self):
        fp = vec3(1,1,0)
        x = self.pset.fp(fp)
        self.assertEqual(x,2)
        self.assertFalse(self.pset.ps[x] is fp)
        self.assertEqual(self.pset.ps[x],fp)
        fp = vec3(2,1,0)
        x = self.pset.fp(fp)
        self.assertEqual(x,-1)
        self.assertFalse(fp in self.pset.ps)

    def test_fps(self):
        fp1 = vec3(1,1,0)
        fp2 = vec3(2,1,0)
        x,y = self.pset.fps([fp1,fp2])
        self.assertEqual(x,2)
        self.assertEqual(y,-1)
        self.assertFalse(self.pset.ps[x] is fp1)
        self.assertEqual(self.pset.ps[x],fp1)
        self.assertFalse(fp2 in self.pset.ps)

    def test_disjoint(self):
        pset = pointset()
        self.assertTrue(pset.disjoint(self.pset))
        pset.ap(vec3(3,0,0))
        self.assertTrue(pset.disjoint(self.pset))
        pset.ap(vec3(0,0,0))
        self.assertFalse(pset.disjoint(self.pset))

    def test_trn(self):
        self.pset.trn(vec3(1,0,0))
        self.assertFalse(vec3(0,0,0) in self.pset.ps)
        self.assertTrue(vec3(1,0,0) in self.pset.ps)
        self.assertTrue(vec3(2,1,1) in self.pset.ps)
        self.assertFalse(vec3(1,2,0) in self.pset.ps)

    def test_rot(self):
        q = quat(0,0,0,0).av(dpr.PI2,vec3(0,0,1))
        self.pset.rot(q)
        self.assertTrue(vec3( 0,0,0) in self.pset.ps)
        self.assertTrue(vec3(-1,0,0) in self.pset.ps)
        self.assertTrue(vec3(-1,1,0) in self.pset.ps)
        self.assertTrue(vec3( 0,1,0) in self.pset.ps)
        self.assertTrue(vec3( 0,0,1) in self.pset.ps)
        self.assertTrue(vec3(-1,0,1) in self.pset.ps)
        self.assertTrue(vec3(-1,1,1) in self.pset.ps)
        self.assertTrue(vec3( 0,1,1) in self.pset.ps)
        self.assertFalse(vec3(1,0,0) in self.pset.ps)
        self.assertFalse(vec3(1,1,0) in self.pset.ps)
        self.assertFalse(vec3(1,0,1) in self.pset.ps)
        self.assertFalse(vec3(1,1,1) in self.pset.ps)

    def test_scl(self):
        sc = vec3(2,3,-1)
        self.pset.scl(sc)
        self.assertTrue(vec3(0,0, 0) in self.pset.ps)
        self.assertTrue(vec3(2,0, 0) in self.pset.ps)
        self.assertTrue(vec3(2,3, 0) in self.pset.ps)
        self.assertTrue(vec3(0,3, 0) in self.pset.ps)
        self.assertTrue(vec3(0,0,-1) in self.pset.ps)
        self.assertTrue(vec3(2,0,-1) in self.pset.ps)
        self.assertTrue(vec3(2,3,-1) in self.pset.ps)
        self.assertTrue(vec3(0,3,-1) in self.pset.ps)
        self.assertFalse(vec3(1,0,0) in self.pset.ps)
        self.assertFalse(vec3(1,1,0) in self.pset.ps)
        self.assertFalse(vec3(1,0,1) in self.pset.ps)
        self.assertFalse(vec3(1,1,1) in self.pset.ps)

    def test_uscl(self):
        self.pset.uscl(-2.3)
        self.assertTrue(vec3(   0,   0,   0) in self.pset.ps)
        self.assertTrue(vec3(-2.3,   0,   0) in self.pset.ps)
        self.assertTrue(vec3(-2.3,-2.3,   0) in self.pset.ps)
        self.assertTrue(vec3(   0,-2.3,   0) in self.pset.ps)
        self.assertTrue(vec3(   0,   0,-2.3) in self.pset.ps)
        self.assertTrue(vec3(-2.3,   0,-2.3) in self.pset.ps)
        self.assertTrue(vec3(-2.3,-2.3,-2.3) in self.pset.ps)
        self.assertTrue(vec3(   0,-2.3,-2.3) in self.pset.ps)
        self.assertFalse(vec3(1,0,0) in self.pset.ps)
        self.assertFalse(vec3(1,1,0) in self.pset.ps)
        self.assertFalse(vec3(1,0,1) in self.pset.ps)
        self.assertFalse(vec3(1,1,1) in self.pset.ps)

if __name__ == '__main__':
    unittest.main()









 


