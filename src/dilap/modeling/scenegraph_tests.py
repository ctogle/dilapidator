import dilap.core.tools as dpr

from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat
from dilap.geometry.tform import tform

import dilap.modeling.scenegraph as dsg

import dilap.mesh.tools as dtl
import matplotlib.pyplot as plt

import unittest,numpy,math

import pdb

#python3 -m unittest discover -v ./ "*tests.py"

class test_scenegraph(unittest.TestCase):

    def test_init(self):
        tr = dsg.scenegraph()
        self.assertEqual(tr.above(tr.root),None)
        self.assertEqual(tr.below(tr.root),[])
        self.assertEqual(tr.vertcount,1)
        self.assertEqual(tr.edgecount,0)

    #def test_mask(self):
    #def test_aroot(self):

    def test_avert(self):
        p1,r1,s1 = vec3(0,0,0),quat(0,0,0,1),vec3(1,1,1)
        p2,r2,s2 = vec3(2,0,1),quat(0,0,0,1).av(dpr.PI2,vec3(0,0,1)),vec3(1,1,1)
        p3,r3,s3 = vec3(0,3,-1),quat(0,0,0,1).av(dpr.PI4,vec3(0,0,1)),vec3(1,1,1)
        tr = dsg.scenegraph()
        one = tr.avert(p1,r1,s1,tr.root)
        two = tr.avert(p2,r2,s2,tr.root)
        three = tr.avert(p3,r3,s3,two)
        self.assertEqual(tr.below(tr.root),[one,two])
        self.assertEqual(tr.above(one),tr.root)
        self.assertEqual(tr.above(two),tr.root)
        self.assertEqual(tr.below(one),[])
        self.assertEqual(tr.below(two),[three])
        self.assertEqual(tr.above(three),two)

    #def test_aedge(self):
    #def test_above(self):
    #def test_below(self):

    def test_graph(self):
        p1,r1,s1 = vec3(0,0,0),quat(0,0,0,1).av(0,vec3(0,0,1)),vec3(1,1,1)
        p2,r2,s2 = vec3(2,0,1),quat(0,0,0,1).av(dpr.PI2,vec3(0,0,1)),vec3(1,1,1)
        p3,r3,s3 = vec3(0,3,-1),quat(0,0,0,1).av(dpr.PI4,vec3(0,0,1)),vec3(1,1,1)
        tr = dsg.scenegraph()
        one = tr.avert(p1,r1,s1,tr.root)
        two = tr.avert(p2,r2,s2,tr.root)
        thr = tr.avert(p3,r3,s3,two)

        tr.graph()

if __name__ == '__main__':
    unittest.main()









 


