import dilap.geometry.tools as dpr

from dilap.geometry.vec3 import vec3 
from dilap.geometry.quat import quat

import matplotlib.pyplot as plt

import unittest,numpy,math,random

#python3 -m unittest discover -v ./ "*tests.py"

class test_vec3(unittest.TestCase):

    # given a vec3, an op, and a res, verify that the op
    # does not return an independent object, and that the result is correct
    # NOTE: this is for methods which return vec3 objects
    def same(self,op,one,res,*args,**kwargs):
        self.assertTrue(one is one)
        opres = one.__getattribute__(op)(*args,**kwargs)
        self.assertTrue(opres is one)
        self.assertEqual(one,res)

    # given a vec3, an op, and a res, verify that the op
    # does return an independent object, and that the result is correct
    # NOTE: this is for methods which return vec3 objects
    def diff(self,op,one,res,*args,**kwargs):
        self.assertTrue(one is one)
        opres = one.__getattribute__(op)(*args,**kwargs)
        self.assertFalse(opres is one)
        self.assertEqual(opres,res)

    # given a vec3, an op, and a res, 
    # verify that the op does return the correct result
    # verify the op does not modify the input vector
    def comp(self,op,one,res,*args,**kwargs):
        cp = one.cp()
        opres = one.__getattribute__(op)(*args,**kwargs)
        self.assertTrue(dpr.isnear(opres,res))
        self.assertEqual(one,cp)

    def setUp(self):
        self.origin = vec3(0,0,0)
        self.one = vec3(1,1,1)
        self.x = vec3(1,0,0)
        self.y = vec3(0,1,0)
        self.z = vec3(0,0,1)
        self.basis = [self.x,self.y,self.z]
        rd = random.random
        self.r1 = vec3(rd()*10,rd()*10,rd()*10)
        self.r2 = vec3(rd()*10,rd()*10,rd()*10)
        self.r3 = vec3(rd()*10,rd()*10,rd()*10)
        self.r4 = vec3(rd()*10,rd()*10,rd()*10)
        self.rds = [self.r1,self.r2,self.r3,self.r4]
        self.each = [self.origin,self.one]+self.basis+self.rds

    def test_cp(self):
        for e in self.each:self.diff('cp',e,e)

    def test_cpxy(self):
        for e in self.each:self.diff('cpxy',e,vec3(e.x,e.y,0))

    #def test_cpr(self):
    #def test_cpf(self):

    def test_d(self):
        for e in self.each:self.comp('d',e,e.mag(),self.origin)
        for e in self.each:self.comp('d',e,0,e)
        self.comp('d',self.x,math.sqrt(2),self.y)
        self.comp('d',self.y,math.sqrt(2),self.z)
        self.comp('d',self.z,math.sqrt(2),self.x)

    def test_dxy(self):
        v1,v2,v3,v4 = vec3(1,1,0),vec3(1,1,2),vec3(1,2,1),vec3(1,2,4)
        self.assertTrue(dpr.isnear(v1.dxy(v1),0))
        self.assertTrue(dpr.isnear(v1.dxy(v2),0))
        self.assertTrue(dpr.isnear(v1.dxy(v3),1))
        self.assertTrue(dpr.isnear(v1.dxy(v4),1))
        self.assertTrue(dpr.isnear(v2.dxy(v3),1))
        self.assertTrue(dpr.isnear(v2.dxy(v4),1))
        self.assertTrue(dpr.isnear(v3.dxy(v4),0))

    def test_ang(self):
        v1,v2,v3,v4 = vec3(1,1,0),vec3(-1,1,0),vec3(-1,0,0),vec3(1,1,1)
        self.assertTrue(dpr.isnear(v1.ang(v1),0))
        self.assertTrue(dpr.isnear(v1.ang(v2),dpr.PI2))
        self.assertTrue(dpr.isnear(v2.ang(v1),dpr.PI2))
        self.assertTrue(dpr.isnear(v1.ang(v3),3*dpr.PI4))
        self.assertTrue(dpr.isnear(v3.ang(v1),3*dpr.PI4))
        self.assertTrue(dpr.isnear(v1.ang(v4),numpy.arctan(1.0/math.sqrt(2))))
        self.assertTrue(dpr.isnear(v4.ang(v1),numpy.arctan(1.0/math.sqrt(2))))

    def test_sang(self):
        p1,p2,p3,p4 = vec3(1,1,0),vec3(0,1,0),vec3(0,-1,0),vec3(0,0,1)
        pn = vec3(0,0,1)
        self.assertEqual(dpr.isnear(p1.sang(p2,pn), dpr.PI4),1)
        self.assertEqual(dpr.isnear(p2.sang(p1,pn), dpr.PI4),0)
        self.assertEqual(dpr.isnear(p2.sang(p1,pn),-dpr.PI4),1)
        self.assertEqual(dpr.isnear(p2.sang(p3,pn), dpr.PI ),1)
        self.assertEqual(dpr.isnear(p3.sang(p1,pn),dpr.threePI4),1)

    def test_angxy(self):
        v1,v2,v3,v4 = vec3(1,1,2),vec3(-1,1,-1),vec3(-1,0,0),vec3(1,1,1)
        self.assertTrue(dpr.isnear(v1.angxy(v1),0))
        self.assertTrue(dpr.isnear(v1.angxy(v2),dpr.PI2))
        self.assertTrue(dpr.isnear(v2.angxy(v1),dpr.PI2))
        self.assertTrue(dpr.isnear(v1.angxy(v3),3*dpr.PI4))
        self.assertTrue(dpr.isnear(v3.angxy(v1),3*dpr.PI4))
        self.assertTrue(dpr.isnear(v1.angxy(v4),0))
        self.assertTrue(dpr.isnear(v4.angxy(v1),0))

    #def test_sang_xy(self):
    #def test_xang_xy(self):

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

    def test_prj(self):
        p1,pn1 = vec3(-1,1,0),vec3(-1,0,0)
        p2,pn2 = vec3(3,-5,2),vec3(0,1,0)
        v1,v2,v3 = vec3(1,1,0),vec3(2,-10,5),vec3(0,1,-1)
        self.assertEqual(v1.cp().prj(p1,pn1),vec3(-1,1,0))
        self.assertEqual(v2.cp().prj(p1,pn1),vec3(-1,-10,5))
        self.assertEqual(v3.cp().prj(p1,pn1),vec3(-1,1,-1))
        self.assertEqual(v1.cp().prj(p2,pn2),vec3(1,-5,0))
        self.assertEqual(v2.cp().prj(p2,pn2),vec3(2,-5,5))
        self.assertEqual(v3.cp().prj(p2,pn2),vec3(0,-5,-1))
        self.assertTrue(v1.prj(p1,pn1) is v1)

    #def test_bary_xy(self):

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

    def test_isnearxy(self):
        v1,v2,v3,v4 = vec3(1,1,0),vec3(1.1,1,0.1),vec3(1,1,1),vec3(1.000001,1,1)
        self.assertEqual(v1.isnearxy(v1),1)
        self.assertEqual(v3.isnearxy(v3),1)
        self.assertEqual(v1.isnearxy(v2),0)
        self.assertEqual(v2.isnearxy(v1),0)
        self.assertEqual(v1.isnearxy(v3),1)
        self.assertEqual(v2.isnearxy(v3),0)
        self.assertEqual(v2.isnearxy(v4),0)
        self.assertEqual(v3.isnearxy(v4),1)

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

    def test_trn(self):
        v1,v2 = vec3(-1,2,5),vec3(-12,24,60)
        self.assertEqual(v1.cp().trn(v2),vec3(-13,26,65))
        self.assertEqual(v2.cp().trn(v1),vec3(-13,26,65))
        self.assertTrue(v1.trn(v2) is v1)

    def test_xtrn(self):
        v1 = vec3(-1,2,5)
        self.assertEqual(v1.cp().xtrn(2),vec3(1,2,5))
        self.assertEqual(v1.xtrn(2),vec3(1,2,5))
        self.assertEqual(v1.xtrn(-5),vec3(-4,2,5))
        self.assertTrue(v1.xtrn(2) is v1)

    def test_ytrn(self):
        v1 = vec3(-1,2,5)
        self.assertEqual(v1.cp().ytrn(2),vec3(-1,4,5))
        self.assertEqual(v1.ytrn(2),vec3(-1,4,5))
        self.assertEqual(v1.ytrn(-5),vec3(-1,-1,5))
        self.assertTrue(v1.ytrn(2) is v1)

    def test_ztrn(self):
        v1 = vec3(-1,2,5)
        self.assertEqual(v1.cp().ztrn(2),vec3(-1,2,7))
        self.assertEqual(v1.ztrn(2),vec3(-1,2,7))
        self.assertEqual(v1.ztrn(-5),vec3(-1,2,2))
        self.assertTrue(v1.ztrn(2) is v1)

    def test_scl(self):
        v1,v2,v3 = vec3(1,1,0),vec3(2,-10,5),vec3(0,1,-1)
        self.assertEqual(v1*v1,vec3(1,1,0))
        self.assertEqual(v2*v2,vec3(4,100,25))
        self.assertEqual(v3*v3,vec3(0,1,1))
        self.assertEqual(v1*v2,vec3(2,-10,0))
        self.assertEqual(v1*v3,vec3(0,1,0))
        self.assertEqual(v2*v1,vec3(2,-10,0))
        self.assertEqual(v2*v3,vec3(0,-10,-5))
        self.assertEqual(v3*v1,vec3(0,1,0))
        self.assertEqual(v3*v2,vec3(0,-10,-5))
        self.assertTrue(v1.scl(v2) is v1)
        self.assertEqual(v1,vec3(2,-10,0))

    def test_uscl(self):
        v1,v2 = vec3(-1,2,5),vec3(-12,24,60)
        self.assertTrue(v1.uscl(12) == v2)
        self.assertTrue(v1.uscl(12) is v1)
        self.assertFalse(v1.uscl(12) is v2)
        self.assertTrue(v1.uscl(12) == v1)

    def test_xscl(self):
        self.same('xscl',self.one,vec3(4,1,1),4)
        self.same('xscl',self.origin,vec3(0,0,0),4)
        self.same('xscl',self.z,vec3(0,0,1),4)

    def test_yscl(self):
        self.same('yscl',self.one,vec3(1,4,1),4)
        self.same('yscl',self.origin,vec3(0,0,0),4)
        self.same('yscl',self.z,vec3(0,0,1),4)

    def test_zscl(self):
        self.same('zscl',self.one,vec3(1,1,4),4)
        self.same('zscl',self.origin,vec3(0,0,0),4)
        self.same('zscl',self.z,vec3(0,0,4),4)

    def test_rot(self):
        v1,v2 = vec3(0,2,0),vec3(-2,0,0)
        q1 = quat(0,0,0,0).av(dpr.PI2,vec3(0,0,1))
        q2 = quat(0,0,0,0).av(0,vec3(0,0,1))
        self.assertEqual(v1.rot(q1),v2)
        self.assertEqual(v1.cp().rot(q2),v1)

    #def test_fulc(self):

    def test_xrot(self):
        self.same('xrot',self.origin,vec3(0,0,0),dpr.PI2)
        self.same('xrot',self.one,vec3(1,-1,1),dpr.PI2)
        self.same('xrot',self.x,vec3(1,0,0),dpr.PI2)
        self.same('xrot',self.y,vec3(0,0,1),dpr.PI2)
        self.same('xrot',self.z,vec3(0,-1,0),dpr.PI2)

    def test_yrot(self):
        self.same('yrot',self.origin,vec3(0,0,0),dpr.PI2)
        self.same('yrot',self.one,vec3(1,1,-1),dpr.PI2)
        self.same('yrot',self.x,vec3(0,0,-1),dpr.PI2)
        self.same('yrot',self.y,vec3(0,1,0),dpr.PI2)
        self.same('yrot',self.z,vec3(1,0,0),dpr.PI2)

    def test_zrot(self):
        self.same('zrot',self.origin,vec3(0,0,0),dpr.PI2)
        self.same('zrot',self.one,vec3(-1,1,1),dpr.PI2)
        self.same('zrot',self.x,vec3(0,1,0),dpr.PI2)
        self.same('zrot',self.y,vec3(-1,0,0),dpr.PI2)
        self.same('zrot',self.z,vec3(0,0,1),dpr.PI2)

    def test_flp(self):
        v1,v2 = vec3(-1,-2,-5),vec3(1,2,5)
        self.assertTrue(v1.flp() == v1)
        self.assertFalse(v1.cp() == v1.flp())
        self.assertTrue(v1.flp() is v1)
        self.assertFalse(v1.flp() is v2)
        self.assertTrue(v1.flp() == v2)

    def test_tov(self):
        v1,v2 = vec3(1,-2,1),vec3(1,2,5)
        self.assertEqual(v1.tov(v1),vec3(0,0,0))
        self.assertEqual(v1.tov(v2),vec3(0, 4, 4))
        self.assertEqual(v2.tov(v1),vec3(0,-4,-4))
        self.assertEqual(v2.tov(v2),vec3(0,0,0))

    def test_tovxy(self):
        v1,v2 = vec3(1,-2,1),vec3(1,2,5)
        self.assertEqual(v1.tovxy(v1),vec3(0,0,0))
        self.assertEqual(v1.tovxy(v2),vec3(0, 4,0))
        self.assertEqual(v2.tovxy(v1),vec3(0,-4,0))
        self.assertEqual(v2.tovxy(v2),vec3(0,0,0))

    def test_mid(self):
        v1,v2 = vec3(0,2,0),vec3(-1,0,1)
        v3,v4 = vec3(-0.5,1,0.5),vec3(-0.75,0.5,0.75)
        self.assertEqual(v1.mid(v2),v3)
        self.assertEqual(v2.mid(v3),v4)

    def test_lerp(self):
        v1,v2 = vec3(0,2,0),vec3(-1,0,1)
        v3,v4 = vec3(-0.75,0.5,0.75),vec3(0,2,0)
        self.assertEqual(v1.lerp(v2,0.75),v3)
        self.assertFalse(v1.lerp(v2,0) is v1)
        self.assertEqual(v1.lerp(v2,0),v1)
        self.assertFalse(v1.lerp(v2,1) is v2)
        self.assertEqual(v1.lerp(v2,1),v2)

    def test_pline(self):
        pline = self.origin.pline(self.one,2)
        d1 = self.origin.d(pline[0])
        d2 = pline[0].d(pline[1])
        d3 = pline[1].d(self.one)
        self.assertEqual(len(pline),2)
        self.assertTrue(dpr.isnear(d1,d2))
        self.assertTrue(dpr.isnear(d2,d3))
        self.assertTrue(self.origin.mid(self.one),pline[0].mid(pline[1]))

    #def test_pring(self):
    #def test_sq(self):
    #def test_com(self):

if __name__ == '__main__':
    unittest.main()









 
