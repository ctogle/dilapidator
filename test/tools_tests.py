import dilap.core.tools as dpr
import dilap.core.vector as dpv
import dilap.core.quaternion as dpq

import dilap.mesh.tools as dtl

import unittest,numpy

#python3 -m unittest discover -v ./ "*tests.py"

class test_tools(unittest.TestCase):

    def test_near(self):
        self.assertEqual(dpr.near(0.00001,0),0)

    def test_isnear(self):
        self.assertEqual(dpr.isnear(0.1,0),0)

    def test_rad(self):
        self.assertEqual(dpr.isnear(dpr.rad(180),dpr.PI),1)

    def test_deg(self):
        self.assertEqual(dpr.isnear(dpr.deg(dpr.PI),180),1)

    def test_clamp(self):
        self.assertEqual(dpr.clamp(180,0,90),90)
        self.assertEqual(dpr.clamp(-90,0,90), 0)

    def test_clamp_periodic(self):
        self.assertEqual(dpr.clamp_periodic(120,0,90),30)
        self.assertEqual(dpr.clamp_periodic(-20,0,90),70)

    def test_inrange(self):
        self.assertEqual(dpr.inrange(  0,0,1),0)
        self.assertEqual(dpr.inrange(0.1,0,1),1)

    def test_adist(self):
        deg10 = dpr.rad(10)
        self.assertEqual(dpr.isnear(dpr.adist(deg10*2,deg10*6),deg10*4),1)
        self.assertEqual(dpr.isnear(dpr.adist(deg10*6,deg10*2),deg10*4),1)
        self.assertEqual(dpr.isnear(dpr.adist(deg10*6,deg10*22),deg10*16),1)
        self.assertEqual(dpr.isnear(dpr.adist(deg10*6,deg10*32),deg10*10),1)

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

    def test_insegment_xy(self):
        p1 = dpv.vector(1,1,0)
        p2 = dpv.vector(0,1,0)
        s1 = dpv.vector(0,0,0)
        s2 = dpv.vector(2,2,0)
        self.assertEqual(dpr.insegment_xy(p1,s1,s2),1)
        self.assertEqual(dpr.insegment_xy(p2,s1,s2),0)

    def test_orient2d(self):
        p1 = dpv.vector(1,1,0)
        p2 = dpv.vector(0,1,0)
        p3 = dpv.vector(0,0,0)
        p4 = dpv.vector(0,2,0)
        self.assertTrue(dpr.orient2d(p1,p2,p3) > 0)
        self.assertTrue(dpr.orient2d(p2,p1,p3) < 0)
        self.assertTrue(dpr.orient2d(p2,p4,p3) == 0)

    def test_angle_between_xy(self):
        p1 = dpv.vector(1,1,1)
        p2 = dpv.vector(0,1,0)
        meth = dpr.angle_between_xy
        self.assertEqual(dpr.isnear(meth(p1,p2), dpr.PI4),1)
        self.assertEqual(dpr.isnear(meth(p2,p1), dpr.PI4),1)
        self.assertEqual(dpr.isnear(meth(p2,p1),-dpr.PI4),0)

    def test_signed_angle_between_xy(self):
        p1 = dpv.vector(1,1,1)
        p2 = dpv.vector(0,1,0)
        meth = dpr.signed_angle_between_xy
        self.assertEqual(dpr.isnear(meth(p1,p2), dpr.PI4),1)
        self.assertEqual(dpr.isnear(meth(p2,p1), dpr.PI4),0)
        self.assertEqual(dpr.isnear(meth(p2,p1),-dpr.PI4),1)

    def test_angle_from_xaxis_xy(self):
        p1 = dpv.vector(1, 1,1)
        p2 = dpv.vector(0, 1,0)
        p3 = dpv.vector(0,-1,0)
        meth = dpr.angle_from_xaxis_xy
        self.assertEqual(dpr.isnear(meth(p1),dpr.PI4),1)
        self.assertEqual(dpr.isnear(meth(p2),dpr.PI4),0)
        self.assertEqual(dpr.isnear(meth(p2),dpr.PI2),1)
        self.assertEqual(dpr.isnear(meth(p3),dpr.threePI2),1)

    def test_angle_between(self):
        p1 = dpv.vector(1, 1,0)
        p2 = dpv.vector(0, 1,0)
        p3 = dpv.vector(0,-1,0)
        meth = dpr.angle_between
        self.assertEqual(dpr.isnear(meth(p1,p2),dpr.PI4),1)
        self.assertEqual(dpr.isnear(meth(p2,p1),dpr.PI4),1)
        self.assertEqual(dpr.isnear(meth(p2,p3),dpr.PI ),1)
        self.assertEqual(dpr.isnear(meth(p3,p1),dpr.threePI4),1)

    def test_signed_angle_between(self):
        p1 = dpv.vector(1, 1,0)
        p2 = dpv.vector(0, 1,0)
        p3 = dpv.vector(0,-1,0)
        pn = dpv.z()
        meth = dpr.signed_angle_between
        self.assertEqual(dpr.isnear(meth(p1,p2,pn), dpr.PI4),1)
        self.assertEqual(dpr.isnear(meth(p2,p1,pn), dpr.PI4),0)
        self.assertEqual(dpr.isnear(meth(p2,p1,pn),-dpr.PI4),1)
        self.assertEqual(dpr.isnear(meth(p2,p3,pn), dpr.PI ),1)
        self.assertEqual(dpr.isnear(meth(p3,p1,pn),dpr.threePI4),1)

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

if __name__ == '__main__':
    unittest.main()





