from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat

import dilap.geometry.tools as gtl

#import dilap.mesh.tools as dtl
import dilap.core.plotting as dtl
import matplotlib.pyplot as plt

import unittest,numpy

#python3 -m unittest discover -v ./ "*tests.py"

class test_tools(unittest.TestCase):

    def setUp(self):
        self.pslist = [
            vec3(-1,-1,0),vec3( 1,-1,0),
            vec3( 1, 1,0),vec3(-1, 1,0)]
        self.pstupl = (
            vec3(-1,-1,0),vec3( 1,-1,0),
            vec3( 1, 1,0),vec3(-1, 1,0))

    def test_isnear(self):
        self.assertEqual(gtl.isnear(0.0001,0),1)
        self.assertEqual(gtl.isnear(0.1,0),0)

    def test_near(self):
        self.assertEqual(gtl.near(0.0001,0),0)
        self.assertEqual(gtl.isnear(gtl.near(0.1,0),0.1),1)

    def test_rad(self):
        self.assertEqual(gtl.isnear(gtl.rad(180),gtl.PI),1)

    def test_deg(self):
        self.assertEqual(gtl.isnear(gtl.deg(gtl.PI),180),1)

    def test_clamp(self):
        self.assertEqual(gtl.clamp(180,0,90),90)
        self.assertEqual(gtl.clamp(-90,0,90), 0)
        self.assertEqual(gtl.clamp( 45,0,90),45)

    def test_wrap(self):
        self.assertEqual(gtl.wrap(120,0,90),30)
        self.assertEqual(gtl.wrap(-20,0,90),70)
        self.assertEqual(gtl.wrap( 45,0,90),45)

    def test_inrng(self):
        self.assertEqual(gtl.inrng(   0,0,1),0)
        self.assertEqual(gtl.inrng(   1,0,1),0)
        self.assertEqual(gtl.inrng( 0.1,0,1),1)
        self.assertEqual(gtl.inrng( 0.9,0,1),1)
        self.assertEqual(gtl.inrng(-0.1,0,1),0)
        self.assertEqual(gtl.inrng( 1.1,0,1),0)

    def test_adist(self):
        deg10 = gtl.rad(10)
        self.assertEqual(gtl.isnear(gtl.adist(deg10*2,deg10*6),deg10*4),1)
        self.assertEqual(gtl.isnear(gtl.adist(deg10*6,deg10*2),deg10*4),1)
        self.assertEqual(gtl.isnear(gtl.adist(deg10*6,deg10*22),deg10*16),1)
        self.assertEqual(gtl.isnear(gtl.adist(deg10*6,deg10*32),deg10*10),1)

    #def test_circumscribe_tri(self):

    '''#
    def test_ang(self):
        p1,p2,p3 = vec3(1, 1,0),vec3(0, 1,0),vec3(0,-1,0)
        meth,nr = gtl.ang,gtl.isnear
        self.assertEqual(nr(meth(p1,p2),gtl.PI4),1)
        self.assertEqual(nr(meth(p2,p1),gtl.PI4),1)
        self.assertEqual(nr(meth(p2,p3),gtl.PI ),1)
        self.assertEqual(nr(meth(p3,p1),gtl.threePI4),1)
    
    def test_sang(self):
        p1,p2,p3,p4 = vec3(1, 1,0),vec3(0, 1,0),vec3(0,-1,0),vec3(0,0,1)
        pn = vec3(0,0,1)
        meth,nr = gtl.sang,gtl.isnear
        self.assertEqual(nr(meth(p1,p2,pn), gtl.PI4),1)
        self.assertEqual(nr(meth(p2,p1,pn), gtl.PI4),0)
        self.assertEqual(nr(meth(p2,p1,pn),-gtl.PI4),1)
        self.assertEqual(nr(meth(p2,p3,pn), gtl.PI ),1)
        self.assertEqual(nr(meth(p3,p1,pn),gtl.threePI4),1)

    def test_ang_xy(self):
        p1,p2 = vec3(1,1,1),vec3(0,1,0)
        meth,nr = gtl.ang_xy,gtl.isnear
        self.assertEqual(nr(meth(p1,p2), gtl.PI4),1)
        self.assertEqual(nr(meth(p2,p1), gtl.PI4),1)
        self.assertEqual(nr(meth(p2,p1),-gtl.PI4),0)

    def test_sang_xy(self):
        p1,p2 = vec3(1,1,1),vec3(0,1,0)
        meth,nr = gtl.sang_xy,gtl.isnear
        self.assertEqual(nr(meth(p1,p2), gtl.PI4),1)
        self.assertEqual(nr(meth(p2,p1), gtl.PI4),0)
        self.assertEqual(nr(meth(p2,p1),-gtl.PI4),1)

    def test_xang_xy(self):
        p1,p2,p3 = vec3(1, 1,1),vec3(0, 1,0),vec3(0,-1,0)
        meth,nr = gtl.xang_xy,gtl.isnear
        self.assertEqual(nr(meth(p1),gtl.PI4),1)
        self.assertEqual(nr(meth(p2),gtl.PI4),0)
        self.assertEqual(nr(meth(p2),gtl.PI2),1)
        self.assertEqual(nr(meth(p3),gtl.threePI2),1)
    
    #def test_bary_xy(self):
    '''#

    def test_onseg_xy(self):
        p1,p2,p3 = vec3(1,1,0),vec3(0,1,0),vec3(2,2,0)
        s1,s2 = vec3(0,0,0),vec3(2,2,0)
        self.assertEqual(gtl.onseg_xy(p1,s1,s2),1)
        self.assertEqual(gtl.onseg_xy(p2,s1,s2),0)
        self.assertEqual(gtl.onseg_xy(p3,s1,s2),1)

    def test_inseg_xy(self):
        p1,p2,p3 = vec3(1,1,0),vec3(0,1,0),vec3(2,2,0)
        s1,s2 = vec3(0,0,0),vec3(2,2,0)
        self.assertEqual(gtl.inseg_xy(p1,s1,s2),1)
        self.assertEqual(gtl.inseg_xy(p2,s1,s2),0)
        self.assertEqual(gtl.inseg_xy(p3,s1,s2),0)

    def test_orient2d(self):
        p1,p2,p3,p4 = vec3(1,1,0),vec3(0,1,0),vec3(0,0,0),vec3(0,2,0)
        self.assertTrue(gtl.orient2d(p1,p2,p3) > 0)
        self.assertTrue(gtl.orient2d(p2,p1,p3) < 0)
        self.assertTrue(gtl.orient2d(p2,p4,p3) == 0)

    #def test_incircle(self):
    #def test_intri(self):
    #def test_onconcave_xy(self):
    #def test_inconcave_xy(self):
    #def test_segs_isect_perp(self):
    #def test_segs_isect_int(self):
    #def test_pline(self):
    #def test_pring(self):
    #def test_rot_poly(self):
    #def test_q_to_xy(self):
    #def test_nrm(self):
    #def test_tng(self):

    def test_poly_nrm(self):
        fp = [vec3(0,0,0),
            vec3(20,0,0),vec3(20,-20,0),vec3(40,-30,0),vec3(60,-30,0),
            vec3(60,40,0),vec3(-60,40,0),vec3(-60,-30,0),vec3(-40,-30,0),
            vec3(-20,-20,0),vec3(-20,0,0)]
        self.assertTrue(gtl.poly_nrm(tuple(fp)).isnear(vec3(0,0,1)))

    #def test_winding(self):

if __name__ == '__main__':
    unittest.main()








'''#

import dilap.core.tools as dpr
import dilap.core.vector as dpv
import dilap.core.quaternion as dpq

import dilap.mesh.tools as dtl
import matplotlib.pyplot as plt

import unittest,numpy,decimal

decimal.getcontext().prec = 5
#python3 -m unittest discover -v ./ "*tests.py"

class test_tools(unittest.TestCase):

    def test_locate_smallest(self):
        vals = [8,1,2,4,10,12,3]
        self.assertEqual(dpr.locate_smallest(vals),1)

    def test_locate_largest(self):
        vals = [8,1,2,4,10,12,3]
        self.assertEqual(dpr.locate_largest(vals),5)

    def test_order_ascending(self):
        vals = [8,1,2,4,10,12,3]
        self.assertEqual(dpr.order_ascending(vals),[1,2,6,3,0,4,5])

    def test_cyclic_permutation(self):
        seq1 = [8,1,2,4,10,12,3]
        seq2 = [2,4,10,12,3,8,1]
        seq3 = [2,8,1,4,10,12,3]
        self.assertEqual(dpr.cyclic_permutation(seq1,seq2),1)
        self.assertEqual(dpr.cyclic_permutation(seq2,seq3),0)

    def test_distance_to_line(self):
        p1 = dpv.vector(1,0,0)
        p2 = dpv.vector(0,0,0)
        p3 = dpv.vector(0,1,0)
        pn = dpv.ny()
        meth = dpr.distance_to_line
        self.assertEqual(dpr.isnear(meth(p3,p2,p1,pn),1.0),1)

    def test_distance_to_border(self):
        p1 = dpv.vector(-2,-2,0)
        p2 = dpv.vector( 2,-2,0)
        p3 = dpv.vector( 2, 2,0)
        p4 = dpv.vector(-2, 2,0)
        tp = dpv.vector( 1, 0,0)
        pn = dpv.vector( 0, 0,1)
        ps = [p1,p2,p3,p4]
        self.assertTrue(dpr.isnear(dpr.distance_to_border(tp,ps),1))
        
    def test_tangents(self):
        p1 = dpv.vector(1,0,0)
        p2 = dpv.vector(0,0,0)
        p3 = dpv.vector(0,1,0)
        self.assertEqual(dpr.tangents([p1,p2,p3]),
            [dpv.vector(-1,0,0),dpv.vector(0,1,0),
                dpv.vector(1,-1,0).normalize()])

    def test_normals(self):
        p1 = dpv.vector(1,0,0)
        p2 = dpv.vector(0,0,0)
        p3 = dpv.vector(0,1,0)
        pn = dpv.z()
        norms = dpr.normals([p1,p2,p3],pn)
        rtres = [
            dpv.vector(0,1,0),dpv.vector(1,0,0),
            dpv.vector(-1,-1,0).normalize()]
        self.assertEqual(norms,rtres)

    def test_revolve_about_edge(self):
        p1 = dpv.vector(2,3,1)
        p2 = dpv.vector(0,0,0)
        p3 = dpv.vector(0,0,5)
        pr = dpr.revolve_about_line(p1,p2,p3,dpr.PI4)
        rr = dpv.vector(2,3,1).rotate_z(dpr.PI4)
        self.assertEqual(pr,rr)

    def test_rotate_coords(self):
        p1 = dpv.vector(2,3,1)
        p2 = dpv.vector(0,0,0)
        p3 = dpv.vector(0,0,5)
        q = dpq.q_from_av(dpr.PI2,dpv.y())
        ps = [p1,p2,p3]
        rps = dpr.rotate_coords(ps,q)
        cps = [p1.rotate(q),p2.rotate(q),p3.rotate(q)]
        self.assertEqual(rps,cps)

    # should be interior only?? or should those guys be on border??
    # should be interior only?? or should those guys be on border??
    # should be interior only?? or should those guys be on border??
    # should be interior only?? or should those guys be on border??
    # should be interior only?? or should those guys be on border??
    def test_intriangle(self):
        t1 = dpv.vector(0,0,0)
        t2 = dpv.vector(1,0,0)
        t3 = dpv.vector(0,1,0)
        self.assertTrue(dpr.intriangle_xy(dpv.vector(0,0,0),t1,t2,t3))
        self.assertTrue(dpr.intriangle_xy(dpv.vector(1,0,0),t1,t2,t3))
        self.assertTrue(dpr.intriangle_xy(dpv.vector(0,1,0),t1,t2,t3))
        self.assertTrue(dpr.intriangle_xy(dpv.vector(0.5,0.5,0),t1,t2,t3))
        self.assertTrue(dpr.intriangle_xy(dpv.vector(0.5,0,0),t1,t2,t3))
        self.assertTrue(dpr.intriangle_xy(dpv.vector(0,0.5,0),t1,t2,t3))
        self.assertTrue(dpr.intriangle_xy(dpv.vector(0.75,0.25,0),t1,t2,t3))
        self.assertFalse(dpr.intriangle_xy(dpv.vector(1,1,0),t1,t2,t3))
        pt = dpv.vector(-22.847272872924805, -22.847267150878906, 0.0)
        t1 = dpv.vector(-21.464466094970703, -28.53553009033203, 0.0)
        t2 = dpv.vector(-17.15900993347168, -24.230073928833008, 0.0)
        t3 = dpv.vector(-24.230077743530273, -17.159006118774414, 0.0)
        self.assertTrue(dpr.intriangle_xy(pt,t1,t2,t3))

    def test_segments_intersect_noncolinear(self):
        p1,p2 = dpv.vector(-10,0,0),dpv.vector(10,0,0)
        p3,p4 = dpv.vector(-5,0,0),dpv.vector(5,0,0)
        p5,p6 = dpv.vector(-5,-2,0),dpv.vector(5,2,0)
        self.assertEqual(dpr.segments_intersect_noncolinear(p1,p2,p3,p4),0)
        self.assertEqual(dpr.segments_intersect_noncolinear(p2,p1,p5,p6),1)

    def test_segments_intersect_at(self):
        p1,p2 = dpv.vector(-10,0,0),dpv.vector(10,0,0)
        p3,p4 = dpv.vector(-5,0,0),dpv.vector(5,0,0)
        isect1 = dtl.segments_intersect_at(p1,p2,p3,p4)
        #isect2 = dtl.segments_intersect_at(p1,p2,p4,p3)
        #isect3 = dtl.segments_intersect_at(p2,p1,p4,p3)
        #isect4 = dtl.segments_intersect_at(p2,p1,p3,p4)
        self.assertTrue(type(isect1) == type(()))

    def test_inconcave_xy(self):
        py = (
            dpv.vector( 0, 0,0),dpv.vector(10,0,0),
            dpv.vector(10,10,0),dpv.vector(0,10,0))
        p1 = dpv.vector(5,5,0)
        p2 = dpv.vector(0,5,0)
        p3 = dpv.vector(0,0,0)
        p4 = dpv.vector(20,5,0)
        self.assertEqual(dpr.inconcave_xy(p1,py),1)
        self.assertEqual(dpr.inconcave_xy(p2,py),0)
        self.assertEqual(dpr.inconcave_xy(p3,py),0)
        self.assertEqual(dpr.inconcave_xy(p4,py),0)

    def test_onconcave_xy(self):
        py = (
            dpv.vector( 0, 0,0),dpv.vector(10,0,0),
            dpv.vector(10,10,0),dpv.vector(0,10,0))
        p1 = dpv.vector(5,5,0)
        p2 = dpv.vector(0,5,0)
        p3 = dpv.vector(0,0,0)
        p4 = dpv.vector(20,5,0)
        self.assertEqual(dpr.onconcave_xy(p1,py),1)
        self.assertEqual(dpr.onconcave_xy(p2,py),1)
        self.assertEqual(dpr.onconcave_xy(p3,py),1)
        self.assertEqual(dpr.onconcave_xy(p4,py),0)
        q = dpv.vector(8.899495124816895, 11.899495124816895, -2.9999990463256836)
        ey = (
            dpv.vector(-15.899495124816895, 11.899495124816895, -2.9999990463256836),
            dpv.vector(-12.79962158203125, 11.899495124816895, -2.9999990463256836),
            dpv.vector(-9.699748039245605, 11.899495124816895, -2.9999990463256836),
            dpv.vector(-6.599874496459961, 11.899495124816895, -2.9999990463256836),
            dpv.vector(-3.5, 11.899495124816895, -2.9999990463256836),
            dpv.vector(-0.40012621879577637, 11.899495124816895, -2.9999990463256836),
            dpv.vector(2.6997475624084473, 11.899495124816895, -2.9999990463256836),
            dpv.vector(5.79962158203125, 11.899495124816895, -2.9999990463256836),
            dpv.vector(8.899495124816895, 11.899495124816895, -2.9999990463256836),
            dpv.vector(8.899495124816895, 8.79962158203125, -2.9999992847442627),
            dpv.vector(8.899495124816895, 5.6997480392456055, -2.999999523162842),
            dpv.vector(8.899495124816895, 2.5998740196228027, -2.999999761581421),
            dpv.vector(8.899495124816895, -0.4999997317790985, -3.0),
            dpv.vector(8.899495124816895, -3.5998735427856445, -3.000000238418579),
            dpv.vector(8.899495124816895, -6.699747085571289, -3.000000476837158),
            dpv.vector(8.899495124816895, -9.79962158203125, -3.0000009536743164),
            dpv.vector(8.899495124816895, -12.899495124816895, -3.0000011920928955),
            dpv.vector(5.919596195220947, -12.899495124816895, -3.0000011920928955),
            dpv.vector(2.939697265625, -12.899495124816895, -3.0000011920928955),
            dpv.vector(-0.040201663970947266, -12.899495124816895, -3.0000011920928955),
            dpv.vector(-3.0201005935668945, -12.899495124816895, -3.0000011920928955),
            dpv.vector(-5.999999523162842, -12.899495124816895, -3.0000011920928955),
            dpv.vector(-8.979898452758789, -12.899495124816895, -3.0000011920928955),
            dpv.vector(-11.959796905517578, -12.899495124816895, -3.0000011920928955),
            dpv.vector(-14.939695358276367, -12.899495124816895, -3.0000011920928955),
            dpv.vector(-17.919593811035156, -12.899495124816895, -3.0000011920928955),
            dpv.vector(-20.899494171142578, -12.899495124816895, -3.0000011920928955),
            dpv.vector(-20.899494171142578, -10.04962158203125, -3.0000009536743164),
            dpv.vector(-20.899494171142578, -7.199747562408447, -3.0000007152557373),
            dpv.vector(-20.899494171142578, -4.349874019622803, -3.000000476837158),
            dpv.vector(-20.899494171142578, -1.499999761581421, -3.000000238418579),
            dpv.vector(-20.899494171142578, 1.3498740196228027, -3.0),
            dpv.vector(-20.899494171142578, 4.1997480392456055, -2.999999523162842),
            dpv.vector(-20.899494171142578, 7.049622058868408, -2.9999992847442627),
            dpv.vector(-20.899494171142578, 9.899495124816895, -2.9999990463256836),
            dpv.vector(-18.399494171142578, 9.899495124816895, -2.9999990463256836),
            dpv.vector(-15.899495124816895, 9.899495124816895, -2.9999990463256836))
        self.assertEqual(dpr.onconcave_xy(q,ey),1)

    def test_concaves_contains(self):
        py1 = (
            dpv.vector( 0, 0,0),dpv.vector(10,0,0),
            dpv.vector(10,10,0),dpv.vector(0,10,0))
        py2 = (
            dpv.vector( 1, 0,0),dpv.vector(11,0,0),
            dpv.vector(11,10,0),dpv.vector(1,10,0))
        py3 = (
            dpv.vector(10, 0,0),dpv.vector(20, 0,0),
            dpv.vector(20,10,0),dpv.vector(10,10,0))
        py4 = (
            dpv.vector(12,0,0),dpv.vector(18,2,0),
            dpv.vector(18,8,0),dpv.vector(12,8,0))
        self.assertEqual(dpr.concaves_contains(py1,py1),1)
        self.assertEqual(dpr.concaves_contains(py1,py3),0)
        self.assertEqual(dpr.concaves_contains(py1,py2),0)
        self.assertEqual(dpr.concaves_contains(py2,py3),0)
        self.assertEqual(dpr.concaves_contains(py3,py4),1)
        self.assertEqual(dpr.concaves_contains(py4,py3),0)
        ey = (
            dpv.vector(-15.899495124816895, 11.899495124816895, -2.9999990463256836),
            dpv.vector(-12.79962158203125, 11.899495124816895, -2.9999990463256836),
            dpv.vector(-9.699748039245605, 11.899495124816895, -2.9999990463256836),
            dpv.vector(-6.599874496459961, 11.899495124816895, -2.9999990463256836),
            dpv.vector(-3.5, 11.899495124816895, -2.9999990463256836),
            dpv.vector(-0.40012621879577637, 11.899495124816895, -2.9999990463256836),
            dpv.vector(2.6997475624084473, 11.899495124816895, -2.9999990463256836),
            dpv.vector(5.79962158203125, 11.899495124816895, -2.9999990463256836),
            dpv.vector(8.899495124816895, 11.899495124816895, -2.9999990463256836),
            dpv.vector(8.899495124816895, 8.79962158203125, -2.9999992847442627),
            dpv.vector(8.899495124816895, 5.6997480392456055, -2.999999523162842),
            dpv.vector(8.899495124816895, 2.5998740196228027, -2.999999761581421),
            dpv.vector(8.899495124816895, -0.4999997317790985, -3.0),
            dpv.vector(8.899495124816895, -3.5998735427856445, -3.000000238418579),
            dpv.vector(8.899495124816895, -6.699747085571289, -3.000000476837158),
            dpv.vector(8.899495124816895, -9.79962158203125, -3.0000009536743164),
            dpv.vector(8.899495124816895, -12.899495124816895, -3.0000011920928955),
            dpv.vector(5.919596195220947, -12.899495124816895, -3.0000011920928955),
            dpv.vector(2.939697265625, -12.899495124816895, -3.0000011920928955),
            dpv.vector(-0.040201663970947266, -12.899495124816895, -3.0000011920928955),
            dpv.vector(-3.0201005935668945, -12.899495124816895, -3.0000011920928955),
            dpv.vector(-5.999999523162842, -12.899495124816895, -3.0000011920928955),
            dpv.vector(-8.979898452758789, -12.899495124816895, -3.0000011920928955),
            dpv.vector(-11.959796905517578, -12.899495124816895, -3.0000011920928955),
            dpv.vector(-14.939695358276367, -12.899495124816895, -3.0000011920928955),
            dpv.vector(-17.919593811035156, -12.899495124816895, -3.0000011920928955),
            dpv.vector(-20.899494171142578, -12.899495124816895, -3.0000011920928955),
            dpv.vector(-20.899494171142578, -10.04962158203125, -3.0000009536743164),
            dpv.vector(-20.899494171142578, -7.199747562408447, -3.0000007152557373),
            dpv.vector(-20.899494171142578, -4.349874019622803, -3.000000476837158),
            dpv.vector(-20.899494171142578, -1.499999761581421, -3.000000238418579),
            dpv.vector(-20.899494171142578, 1.3498740196228027, -3.0),
            dpv.vector(-20.899494171142578, 4.1997480392456055, -2.999999523162842),
            dpv.vector(-20.899494171142578, 7.049622058868408, -2.9999992847442627),
            dpv.vector(-20.899494171142578, 9.899495124816895, -2.9999990463256836),
            dpv.vector(-18.399494171142578, 9.899495124816895, -2.9999990463256836),
            dpv.vector(-15.899495124816895, 9.899495124816895, -2.9999990463256836))
        qy = (
            dpv.vector(8.899495124816895, 8.79962158203125, -2.9999992847442627),
            dpv.vector(8.899495124816895, 11.899495124816895, -2.9999990463256836),
            dpv.vector(5.79962158203125, 11.899495124816895, -2.9999990463256836))
        self.assertEqual(dpr.concaves_contains(ey,qy),1)

if __name__ == '__main__':
    unittest.main()

'''#












