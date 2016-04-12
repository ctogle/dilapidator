from dilap.geometry.vec3 import vec3
from dilap.geometry.curve import curve

import matplotlib.pyplot as plt

import unittest,numpy,math

import pdb

#python3 -m unittest discover -v ./ "*tests.py"

class test_curve(unittest.TestCase):

    def setUp(self):
        self.c1 = curve(vec3(-1,-1,-1),vec3(1,1,1))
        self.c1.segs = 5
        self.c1.calc()

    def test_cp(self):
        c2 = self.c1.cp()
        self.assertFalse(self.c1 is c2)
        self.assertEqual(self.c1,c2)

    def test_cpxy(self):
        c2 = self.c1.cpxy()
        self.assertFalse(self.c1 is c2)
        self.assertFalse(self.c1 == c2)

    def test_cln(self):
        self.assertTrue(len(self.c1.pts)>0)
        self.assertTrue(len(self.c1.tns)>0)
        self.assertTrue(len(self.c1.nms)>0)
        self.c1.cln()
        self.assertTrue(len(self.c1.pts)==0)
        self.assertTrue(len(self.c1.tns)==0)
        self.assertTrue(len(self.c1.nms)==0)

    #def test_calcone(self):

    def test_calc(self):
        self.assertTrue(len(self.c1.pts)==self.c1.segs-1)
        self.assertTrue(len(self.c1.tns)==self.c1.segs-1)
        self.assertTrue(len(self.c1.nms)==self.c1.segs-1)
        print(self.c1)

if __name__ == '__main__':
    unittest.main()









 


