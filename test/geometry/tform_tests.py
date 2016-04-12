from dilap.geometry.vec3 import vec3 
from dilap.geometry.quat import quat
from dilap.geometry.tform import tform

import dilap.geometry.tools as dpr

import matplotlib.pyplot as plt

import unittest,numpy,math

import pdb

#python3 -m unittest discover -v ./ "*tests.py"

class test_tform(unittest.TestCase):

    def setUp(self):
        a1 = dpr.PI2
        v1,v2,v3 = vec3(1,0,0),vec3(0,1,0),vec3(0,0,1)
        q0 = quat(0,0,0,0).av( 0,v3)
        q1 = quat(0,0,0,0).av(a1,v1)
        q2 = quat(0,0,0,0).av(a1,v2)
        q3 = quat(0,0,0,0).av(a1,v3)
        self.tf1 = tform(vec3(1,1,0),q3.cp(),vec3(1,2,1))
        self.tf2 = tform(vec3(1,1,0),q3.cp(),vec3(1,2,1))
        self.tf3 = tform(vec3(0,1,0),q1,vec3(1,1,1))
        a2 = dpr.PI
        q4 = quat(0,0,0,0).av(a2,v3)
        self.tf4 = tform(vec3(0,2,0),q4,vec3(1,4,1))

    def test_cp(self):
        self.assertTrue(self.tf1 is self.tf1)
        self.assertFalse(self.tf1 is self.tf1.cp())
        self.assertTrue(self.tf1 == self.tf1.cp())
        self.assertFalse(self.tf1 is self.tf2)
        self.assertTrue(self.tf1 == self.tf2)

    def test_true(self):
        tf4 = self.tf1.true(self.tf2)
        self.assertEqual(self.tf4,tf4)

if __name__ == '__main__':
    unittest.main()









 

