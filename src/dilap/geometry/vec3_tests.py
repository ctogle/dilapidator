from dilap.geometry.vec3 import vec3 

import dilap.core.tools as dpr
import dilap.core.quaternion as dpq

import dilap.mesh.tools as dtl
import matplotlib.pyplot as plt

import unittest,numpy,math

#python3 -m unittest discover -v ./ "*tests.py"

class test_vec3(unittest.TestCase):

    def test_cp(self):
        v1,v2 = vec3(1,1,0),vec3(1,2,0)
        c1,c2 = vec3(1,1,0),vec3(1,2,0)
        self.assertTrue(v1 is v1)
        self.assertFalse(v1 is v1.cp())
        self.assertTrue(v1 == v1.cp())
        self.assertFalse(v1 is v2)
        self.assertFalse(v1 is v2.cp())
        self.assertFalse(v1 == v2.cp())

    #def test_cpxy(self):

    def test_d(self):
        v1,v2,v3,v4 = vec3(1,1,0),vec3(1,2,0),vec3(1,2,1),vec3(1,1,1)
        self.assertEqual(dpr.isnear(v1.d(v1),0.0000000),1)
        self.assertEqual(dpr.isnear(v1.d(v2),1.0000000),1)
        self.assertEqual(dpr.isnear(v1.d(v2), v2.d(v1)),1)
        self.assertEqual(dpr.isnear(v1.d(v2),1.0100000),0)
        self.assertEqual(dpr.isnear(v1.d(v3),1.4142356),1)
        self.assertEqual(dpr.isnear(v1.d(v3),1.0100000),0)
        self.assertEqual(dpr.isnear(v1.d(v4),1.0000000),1)
        self.assertEqual(dpr.isnear(v1.d(v4),1.0100000),0)

    #def test_dxy(self):

    def test_ang(self):
        v1,v2,v3,v4 = vec3(1,1,0),vec3(-1,1,0),vec3(-1,0,0),vec3(1,1,1)
        self.assertTrue(dpr.isnear(v1.ang(v1),0))
        self.assertTrue(dpr.isnear(v1.ang(v2),dpr.PI2))
        self.assertTrue(dpr.isnear(v2.ang(v1),dpr.PI2))
        self.assertTrue(dpr.isnear(v1.ang(v3),3*dpr.PI4))
        self.assertTrue(dpr.isnear(v3.ang(v1),3*dpr.PI4))
        self.assertTrue(dpr.isnear(v1.ang(v4),numpy.arctan(1.0/math.sqrt(2))))
        self.assertTrue(dpr.isnear(v4.ang(v1),numpy.arctan(1.0/math.sqrt(2))))

    #def test_angxy(self):

    def test_dot(self):
        v1,v2,v3,v4 = vec3(1,1,0),vec3(-1,1,0),vec3(-1,0,0),vec3(1,1,1)
        self.assertEqual(dpr.isnear(v1.dot(v1), 2),1)
        self.assertEqual(dpr.isnear(v1.dot(v2), 0),1)
        self.assertEqual(dpr.isnear(v1.dot(v3),-1),1)
        self.assertEqual(dpr.isnear(v1.dot(v4), 2),1)
        self.assertEqual(dpr.isnear(v2.dot(v2), 2),1)
        self.assertEqual(dpr.isnear(v2.dot(v3), 1),1)
        self.assertEqual(dpr.isnear(v2.dot(v4), 0),1)
        self.assertEqual(dpr.isnear(v3.dot(v3), 1),1)
        self.assertEqual(dpr.isnear(v3.dot(v4),-1),1)
        self.assertEqual(dpr.isnear(v4.dot(v4), 3),1)

    def test_crs(self):
        v1,v2,v3,v4 = vec3(1,1,0),vec3(-1,1,0),vec3(-1,0,0),vec3(1,1,1)
        self.assertEqual(v1.crs(v1),vec3(0,0,0))
        self.assertEqual(v2.crs(v2),vec3(0,0,0))
        self.assertEqual(v3.crs(v3),vec3(0,0,0))
        self.assertEqual(v4.crs(v4),vec3(0,0,0))
        self.assertEqual(v1.crs(v2),vec3(0,0,2))
        self.assertEqual(v1.crs(v3),vec3(0,0,1))
        self.assertEqual(v1.crs(v4),vec3(1,-1,0))
        self.assertEqual(v2.crs(v3),vec3(0,0,1))
        self.assertEqual(v2.crs(v4),vec3(1,1,-2))
        self.assertEqual(v3.crs(v4),vec3(0,1,-1))
        self.assertEqual(v1.crs(v2),v2.crs(v1).flp())
        self.assertEqual(v1.crs(v3),v3.crs(v1).flp())
        self.assertEqual(v1.crs(v4),v4.crs(v1).flp())

    #def test_prj(self):

    def test_inneighborhood(self):
        v1,v2,v3,v4 = vec3(1,1,0),vec3(1,2,0),vec3(1,2,1),vec3(1,1,1)
        self.assertEqual(v1.inneighborhood(v2,1.00),0)
        self.assertEqual(v1.inneighborhood(v2,1.01),1)
        self.assertEqual(v1.inneighborhood(v3,1.00),0)
        self.assertEqual(v1.inneighborhood(v3,1.01),0)
        self.assertEqual(v1.inneighborhood(v4,1.00),0)
        self.assertEqual(v1.inneighborhood(v4,1.01),1)

    def test_isnear(self):
        v1,v2,v3,v4 = vec3(1,1,0),vec3(1,1,0.1),vec3(1,1,1),vec3(1.000001,1,1)
        self.assertEqual(v1.isnear(v1),1)
        self.assertEqual(v3.isnear(v3),1)
        self.assertEqual(v1.isnear(v2),0)
        self.assertEqual(v2.isnear(v1),0)
        self.assertEqual(v1.isnear(v3),0)
        self.assertEqual(v2.isnear(v3),0)
        self.assertEqual(v2.isnear(v4),0)
        self.assertEqual(v3.isnear(v4),1)

    #def test_isnearxy(self):

    def test_mag2(self):
        v1,v2,v3 = vec3(1,0,0),vec3(1,1,1),vec3(2,5,11)
        self.assertEqual(dpr.isnear(v1.mag2(),1),1)
        self.assertEqual(dpr.isnear(v2.mag2(),3),1)
        self.assertEqual(dpr.isnear(v3.mag2(),150),1)

    def test_mag(self):
        v1,v2,v3 = vec3(1,0,0),vec3(1,1,1),vec3(2,5,11)
        self.assertEqual(dpr.isnear(v1.mag(),1),1)
        self.assertEqual(dpr.isnear(v2.mag(),math.sqrt(3)),1)
        self.assertEqual(dpr.isnear(v3.mag(),math.sqrt(150)),1)

    def test_nrm(self):
        v1,v2,v3 = vec3(1,0,0),vec3(1,2,5),vec3(10,20,50)
        self.assertTrue(v1.nrm() == v1)
        self.assertTrue(v1.nrm() is v1)
        self.assertTrue(v2.nrm() == v3.nrm())
        self.assertFalse(v2.nrm() is v3.nrm())
        self.assertFalse(v2.nrm() == v1.nrm())

    #def test_trn(self):
    #def test_xtrn(self):
    #def test_ytrn(self):
    #def test_ztrn(self):

    def test_scl(self):
        v1,v2 = vec3(-1,2,5),vec3(-12,24,60)
        self.assertTrue(v1.scl(12) == v2)
        self.assertTrue(v1.scl(12) is v1)
        self.assertFalse(v1.scl(12) is v2)
        self.assertTrue(v1.scl(12) == v1)

    #def test_xscl(self):
    #def test_yscl(self):
    #def test_zscl(self):

    def test_rot(self):
        v1,v2 = vec3(0,2,0),vec3(-1,0,1)

    def test_xrot(self):
        v1,v2 = vec3(0,2,0),vec3(-1,0,1)

    def test_yrot(self):
        v1,v2 = vec3(0,2,0),vec3(-1,0,1)

    def test_zrot(self):
        v1,v2 = vec3(0,2,0),vec3(-1,0,1)

    def test_flp(self):
        v1,v2 = vec3(-1,-2,-5),vec3(1,2,5)
        self.assertTrue(v1.flp() == v1)
        self.assertFalse(v1.cp() == v1.flp())
        self.assertTrue(v1.flp() is v1)
        self.assertFalse(v1.flp() is v2)
        self.assertTrue(v1.flp() == v2)

    #def test_tov(self):
    #def test_mid(self):

    def test_lterp(self):
        v1,v2 = vec3(0,2,0),vec3(-1,0,1)
        v3,v4 = vec3(-0.75,0.5,0.75),vec3(0,2,0)
        self.assertEqual(v1.lterp(v2,0.75),v3)
        self.assertFalse(v1.lterp(v2,0) is v1)
        self.assertEqual(v1.lterp(v2,0),v1)
        self.assertFalse(v1.lterp(v2,1) is v2)
        self.assertEqual(v1.lterp(v2,1),v2)

if __name__ == '__main__':
    unittest.main()









 
