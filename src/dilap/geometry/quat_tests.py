from dilap.geometry.quat import quat
from dilap.geometry.vec3 import vec3

import dilap.core.tools as dpr

import dilap.mesh.tools as dtl
import matplotlib.pyplot as plt

import unittest,numpy,math

import pdb

#python3 -m unittest discover -v ./ "*tests.py"

class test_quat(unittest.TestCase):

    def test_av(self):
        a = 3*dpr.PI4
        u1,u2,u3 = vec3(1,0,0),vec3(0,-1,0),vec3(0,0,1)
        q1,q2 = quat(0,0,0,0).av(a,u1),quat(0,0,0,0).av(a,u2)
        q3,q4 = quat(0,0,0,0).av(-a,u3),quat(0,0,0,0).av(-a,u2)
        self.assertTrue(q1.w > 0.1)
        self.assertTrue(q1.x > 0.1)
        self.assertTrue(dpr.isnear(q1.y,0))
        self.assertTrue(dpr.isnear(q1.z,0))
        self.assertTrue(q2.w > 0.1)
        self.assertTrue(dpr.isnear(q2.x,0))
        self.assertTrue(q2.y < -0.1)
        self.assertTrue(dpr.isnear(q2.z,0))
        self.assertTrue(q3.w > 0.1)
        self.assertTrue(dpr.isnear(q3.x,0))
        self.assertTrue(dpr.isnear(q3.y,0))
        self.assertTrue(q3.z < -0.1)
        self.assertFalse(q2 == q4.cp().flp())
        self.assertTrue(q2 == q4.cnj())

    def test_uu(self):
        u1,u2,u3 = vec3(1,0,0),vec3(0,-1,0),vec3(0,0,1)
        q1,q2 = quat(0,0,0,0).uu(u1,u2),quat(0,0,0,0).uu(u1,u3)
        q3,q4 = quat(0,0,0,0).uu(u2,u3),quat(0,0,0,0).uu(u3,u2)
        self.assertTrue(q1.w >  0.1)
        self.assertTrue(dpr.isnear(q1.x,0))
        self.assertTrue(dpr.isnear(q1.y,0))
        self.assertTrue(q1.z < -0.1)
        self.assertTrue(q2.w >  0.1)
        self.assertTrue(dpr.isnear(q2.x,0))
        self.assertTrue(q2.y < -0.1)
        self.assertTrue(dpr.isnear(q2.z,0))
        self.assertTrue(q3 == q4.cnj())

    def test_cp(self):
        q1 = quat(1,2,3,4)
        self.assertTrue(q1 is q1)
        self.assertFalse(q1 is q1.cp())
        self.assertTrue(q1 == q1.cp())

    def test_isnear(self):
        q1,q2 = quat(1,1,1,0),quat(1,1,1,0.1)
        q3,q4 = quat(1,1,1,1),quat(1,1.000001,1,1)
        self.assertEqual(q1.isnear(q1),1)
        self.assertEqual(q3.isnear(q3),1)
        self.assertEqual(q1.isnear(q2),0)
        self.assertEqual(q2.isnear(q1),0)
        self.assertEqual(q1.isnear(q3),0)
        self.assertEqual(q2.isnear(q3),0)
        self.assertEqual(q2.isnear(q4),0)
        self.assertEqual(q3.isnear(q4),1)

    def test_mag2(self):
        q1,q2,q3 = quat(1,0,0,0),quat(1,1,1,0),quat(0,2,5,11)
        self.assertEqual(dpr.isnear(q1.mag2(),1),1)
        self.assertEqual(dpr.isnear(q2.mag2(),3),1)
        self.assertEqual(dpr.isnear(q3.mag2(),150),1)

    def test_mag(self):
        q1,q2,q3 = quat(1,0,0,0),quat(1,1,1,0),quat(0,2,5,11)
        self.assertEqual(dpr.isnear(q1.mag(),1),1)
        self.assertEqual(dpr.isnear(q2.mag(),math.sqrt(3)),1)
        self.assertEqual(dpr.isnear(q3.mag(),math.sqrt(150)),1)

    def test_nrm(self):
        q1,q2,q3 = quat(1,0,0,0),quat(1,1,1,0),quat(0,2,5,11)
        self.assertEqual(dpr.isnear(q1.cp().nrm().mag(),1),1)
        self.assertEqual(dpr.isnear(q2.cp().nrm().mag(),1),1)
        self.assertEqual(dpr.isnear(q3.cp().nrm().mag(),1),1)
        self.assertTrue(q1.cp().nrm().mag() == q1.mag())
        self.assertTrue(q1.nrm() is q1)
        self.assertFalse(q2.cp().nrm().mag() == q2.mag())
        self.assertTrue(q2.nrm() is q2)
        self.assertFalse(q3.cp().nrm().mag() == q3.mag())
        self.assertTrue(q3.nrm() is q3)

    def test_flp(self):
        q1,q2 = quat(1,0,0,0),quat(1,1,1,0)
        q3,q4 = quat(0,2,5,11),quat(-1,1,1,0)
        self.assertFalse(q1.cp().flp() == q1)
        self.assertFalse(q2.cp().flp() == q2)
        self.assertTrue(q3.cp().flp() == q3)
        self.assertFalse(q4.cp().flp() == q4)
        self.assertTrue(q2.cp().flp() == q4)
        self.assertTrue(q1.flp() is q1)
        self.assertTrue(q2.flp() is q2)
        self.assertTrue(q3.flp() is q3)
        self.assertTrue(q4.flp() is q4)

    def test_uscl(self):
        q1,q2 = quat(1,0,0,0),quat(1,1,1,0)
        q3,q4 = quat(0,2,5,11),quat(0,1,2.5,5.5)
        self.assertTrue(q1.cp().uscl(1) == q1)
        self.assertFalse(q1.cp().uscl(3) == q1)
        self.assertTrue(q2.cp().uscl(1) == q2)
        self.assertFalse(q2.cp().uscl(3) == q2)
        self.assertTrue(q3.cp().uscl(0.5) == q4)
        self.assertTrue(q1.uscl(1) is q1)

    def test_cnj(self):
        q1,q2 = quat(1,0,0,0),quat(1,1,1,0)
        q3,q4 = quat(-1,2,5,11),quat(1,-2,-5,-11)
        self.assertTrue(q1.cp().cnj() == q1)
        self.assertTrue(q1.cnj() is q1)
        self.assertFalse(q2.cp().cnj() == q2)
        self.assertFalse(q3.cnj() == q4)
               
    def test_inv(self):
        a1,v1 = dpr.PI4,vec3(0,0,1)
        a2,v2 = dpr.threePI4,vec3(0,0,1)
        q1,q2 = quat(1,0,0,0).av(a1,v1),quat(1,1,1,0).av(a2,v2)
        self.assertEqual(q1.cp().cnj(),q1.inv())
        self.assertEqual(q2.cp().cnj(),q2.inv())
        self.assertFalse(q1.inv() is q1)
        
    def test_add(self):
        q1,q2 = quat(0.5,0.3,-2.2,3),quat(1,1.1,2,-0.5)
        q3 = quat(1.5,1.4,-0.2,2.5)
        self.assertEqual(q1.add(q2),q3)
        self.assertFalse(q1.add(q2) is q1)

    def test_sub(self):
        q1,q2 = quat(0.5,0.3,-2.2,3),quat(1,1.1,2,-0.5)
        q3 = quat(-0.5,-0.8,-4.2,3.5)
        self.assertEqual(q1.sub(q2),q3)
        self.assertFalse(q1.sub(q2) is q1)

    def test_mul(self):
        a1,v1 = dpr.PI4,vec3(0,0,1)
        a2,v2 = dpr.threePI4,vec3(0,0,1)
        q1,q2 = quat(1,0,0,0).av(a1,v1),quat(1,1,1,0).av(a2,v1)
        q3 = quat(0,1,0,0).av(a1+a2,v2)
        self.assertTrue(q1.mul(q2) == q3)
        self.assertFalse(q1.mul(q2) is q1)

    def test_rot(self):
        a1,v1 = dpr.PI4,vec3(0,0,1)
        a2,v2 = dpr.PI2,vec3(0,0,1)
        q1,q2 = quat(1,0,0,0).av(a1,v1),quat(1,1,1,0).av(a1,v1)
        q3 = quat(0,1,0,0).av(a2,v2)
        self.assertTrue(q1.rot(q2) == q3)
        self.assertTrue(q1.rot(q2) is q1)

    def test_dot(self):
        a1,v1 = dpr.PI4,vec3(0,0,1)
        a2,v2 = dpr.PI2,vec3(0,1,0)
        q1,q2 = quat(1,0,0,0).av(a1,v1),quat(1,1,1,0).av(a1,v1)
        q3 = quat(0,1,0,0).av(a2,v2)
        q4 = quat(0,1,0,0).av(0,v1)
        self.assertTrue(dpr.isnear(q1.dot(q2),q1.mag2()))
        self.assertFalse(dpr.isnear(q1.dot(q3),0))
        self.assertTrue(dpr.isnear(q3.dot(q4),q3.w))

    def test_slerp(self):
        a1,v1 = dpr.PI4,vec3(0,0,1)
        a2,v2 = dpr.PI,vec3(0,0,1)
        q1,q2 = quat(1,0,0,0).av(0,v1),quat(1,1,1,0).av(a1,v1)
        q3 = quat(0,1,0,0).av(a2,v2)
        self.assertEqual(q1.slerp(q3,0.25),q2)
        self.assertFalse(q1.slerp(q3,0.25) is q1) 

if __name__ == '__main__':
    unittest.main()









 

