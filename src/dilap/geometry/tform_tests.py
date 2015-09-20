from dilap.geometry.vec3 import vec3 
from dilap.geometry.quat import quat
from dilap.geometry.tform import tform

import dilap.core.tools as dpr

import dilap.mesh.tools as dtl
import matplotlib.pyplot as plt

import unittest,numpy,math

import pdb

#python3 -m unittest discover -v ./ "*tests.py"

class test_tform(unittest.TestCase):

    def test_cp(self):
        p1,s1 = vec3(1,1,0),vec3(1,1,2)
        r1 = quat(0,1,0,0).av(dpr.PI4,vec3(0,0,1))
        tf1 = tform(p1,r1,s1)
        p2,s2 = vec3(1,1,0),vec3(1,1,2)
        r2 = quat(0,1,0,0).av(dpr.PI4,vec3(0,0,1))
        tf2 = tform(p2,r2,s2)
        self.assertTrue(tf1 is tf1)
        self.assertFalse(tf1 is tf1.cp())
        self.assertTrue(tf1 == tf1.cp())
        self.assertFalse(tf1 is tf2)
        self.assertFalse(tf1 is tf2.cp())
        self.assertTrue(tf1 == tf2.cp())

    def test_true(self):
        p1,s1 = vec3(1,1,0),vec3(1,1,2)
        r1 = quat(0,1,0,0).av(dpr.PI4,vec3(0,0,1))
        tf1 = tform(p1,r1,s1)
        p2,s2 = vec3(1,1,0),vec3(1,1,2)
        r2 = quat(0,1,0,0).av(dpr.PI4,vec3(0,0,1))
        tf2 = tform(p2,r2,s2)
        p3,s3 = vec3(1,1+math.sqrt(2)/2.0,0),vec3(1,1,4)
        r3 = quat(0,1,0,0).av(dpr.PI2,vec3(0,0,1))
        tf3 = tform(p3,r3,s3)
        tf4 = tf1.true(tf2)
        self.assertEqual(tf3,tf4)

if __name__ == '__main__':
    unittest.main()









 

