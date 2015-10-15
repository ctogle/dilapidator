import dilap.core.tools as dpr

from dilap.geometry.vec3 import vec3
from dilap.geometry.pointset import pointset

import dilap.topology.meshes.polygonmesh as dmsh

import dilap.mesh.tools as dtl
import matplotlib.pyplot as plt

import unittest,numpy,math

import pdb

#python3 -m unittest discover -v ./ "*tests.py"

class test_polygonmesh(unittest.TestCase):

    def plotmesh(self):
        mesh = self.mesh
        ax = dtl.plot_axes()
        for f in mesh.faces:                    
            for l in mesh.mask(2,None,None,None,f):
                hes = mesh.mask(6,None,None,l,None)
                p1 = self.pset.ps[hes[0].one[1]]
                p2 = self.pset.ps[hes[0].two[1]]
                ps = [p1,p2]
                for hx in range(1,len(hes)):
                    onep = self.pset.ps[hes[hx].one[1]]
                    twop = self.pset.ps[hes[hx].two[1]]
                    if   ps[-1].isnear(onep):ps.append(twop)
                    elif ps[-1].isnear(twop):ps.append(onep)
                    else:raise ValueError
                ax = dtl.plot_edges(ps,ax)
        plt.show()

    def assert_counts(self,v,e,l,f,fg,s):
        self.assertEqual(self.mesh.vcnt(),v)
        self.assertEqual(self.mesh.ecnt(),e)
        self.assertEqual(self.mesh.lcnt(),l)
        self.assertEqual(self.mesh.fcnt(),f)
        self.assertEqual(self.mesh.fgcnt(),fg)
        self.assertEqual(self.mesh.scnt(),s)

    def cube_points(self):
        self.pset.ap(vec3(0,0,0));self.pset.ap(vec3(2,0,0))
        self.pset.ap(vec3(2,2,0));self.pset.ap(vec3(0,2,0))
        self.pset.ap(vec3(0,0,2));self.pset.ap(vec3(2,0,2))
        self.pset.ap(vec3(2,2,2));self.pset.ap(vec3(0,2,2))

    def setUp(self):
        self.mesh = dmsh.polygonmesh()
        self.pset = pointset()

    #def test_mask_0vnnnnn(self):
    #def test_mask_0nennnn(self):
    #def test_mask_0nnlnnn(self):
    #def test_mask_0nnnfnn(self):
    #def test_mask_0nnnngn(self):
    #def test_mask_0nnnngs(self):

    def test_mask_1vnnnnn(self):
        v1,f1 = self.mesh.mbfv(vgx = 0)
        v2,e1 = self.mesh.mev(v1,vgx = 0)
        v3,e2 = self.mesh.mev(v2,vgx = 0)
        v4,e3 = self.mesh.mev(v3,vgx = 0)
        v5,e4 = self.mesh.mev(v4,vgx = 0)
        v6,e5 = self.mesh.mev(v5,vgx = 0)
        v7,e6 = self.mesh.mev(v6,vgx = 0)
        v8,e7 = self.mesh.mev(v7,vgx = 0)
        e8,f2 = self.mesh.mefl(v8,v1)

        mask1 = self.mesh.mask(1,v1)
        self.assertTrue(e1 in mask1)
        self.assertFalse(e2 in mask1)
        self.assertFalse(e3 in mask1)
        self.assertFalse(e4 in mask1)
        self.assertFalse(e5 in mask1)
        self.assertFalse(e6 in mask1)
        self.assertFalse(e7 in mask1)
        self.assertTrue(e8 in mask1)

        mask2 = self.mesh.mask(1,v2)
        self.assertTrue(e1 in mask2)
        self.assertTrue(e2 in mask2)
        self.assertFalse(e3 in mask2)
        self.assertFalse(e4 in mask2)
        self.assertFalse(e5 in mask2)
        self.assertFalse(e6 in mask2)
        self.assertFalse(e7 in mask2)
        self.assertFalse(e8 in mask2)

        mask8 = self.mesh.mask(1,v8)
        self.assertFalse(e1 in mask8)
        self.assertFalse(e2 in mask8)
        self.assertFalse(e3 in mask8)
        self.assertFalse(e4 in mask8)
        self.assertFalse(e5 in mask8)
        self.assertFalse(e6 in mask8)
        self.assertTrue(e7 in mask8)
        self.assertTrue(e8 in mask8)

    #def test_mask_1nennnn(self):
    #def test_mask_1nnlnnn(self):
    #def test_mask_1nnnfnn(self):
    #def test_mask_1nnnngn(self):
    #def test_mask_1nnnnns(self):

    def test_mask_2vnnnnn(self):
        self.cube_points()

        v1,f1 = self.mesh.mbfv(vgx = 0)
        v2,e1 = self.mesh.mev(v1,vgx = 1)
        v3,e2 = self.mesh.mev(v2,vgx = 2)
        v4,e3 = self.mesh.mev(v3,vgx = 3)
        v5,e4 = self.mesh.mev(v4,vgx = 7)
        v6,e5 = self.mesh.mev(v5,vgx = 6)
        v7,e6 = self.mesh.mev(v6,vgx = 5)
        v8,e7 = self.mesh.mev(v7,vgx = 4)
        e8,f2 = self.mesh.mefl(v8,v1)
        v3hes = self.mesh.mask(6,v3,None,None,f1)
        v6hes = self.mesh.mask(6,v6,None,None,f1)
        oe1,oe2 = v6hes[1],v3hes[1]
        e9,f3 = self.mesh.mefl(v3,v6,oe1,oe2)

        mask1 = self.mesh.mask(2,v1)
        self.assertTrue(self.mesh.loops[0] in mask1)
        self.assertTrue(self.mesh.loops[1] in mask1)
        self.assertFalse(self.mesh.loops[2] in mask1)
        mask2 = self.mesh.mask(2,v2)
        self.assertTrue(self.mesh.loops[0] in mask2)
        self.assertTrue(self.mesh.loops[1] in mask2)
        self.assertFalse(self.mesh.loops[2] in mask2)
        mask3 = self.mesh.mask(2,v3)
        self.assertTrue(self.mesh.loops[0] in mask3)
        self.assertTrue(self.mesh.loops[1] in mask3)
        self.assertTrue(self.mesh.loops[2] in mask3)
        self.assertEqual(self.mesh.lcnt(),3)

    def test_mask_2nennnn(self):
        v1,f1 = self.mesh.mbfv()
        v2,e1 = self.mesh.mev(v1)
        v3,e2 = self.mesh.mev(v2)
        v4,e3 = self.mesh.mev(v3)
        mask1 = self.mesh.mask(2,None,e1)
        mask2 = self.mesh.mask(2,None,e2)
        mask3 = self.mesh.mask(2,None,e3)
        self.assertEqual(mask1,mask2)
        self.assertEqual(mask1,mask3)

    #def test_mask_2nnlnnn(self):
    #def test_mask_2nnnfnn(self):
    #def test_mask_2nnnngn(self):
    #def test_mask_2nnnnns(self):

    #def test_mask_3vnnnnn(self):
    #def test_mask_3nennnn(self):
    #def test_mask_3nnlnnn(self):
    #def test_mask_3nnnfnn(self):
    #def test_mask_3nnnngn(self):
    #def test_mask_3nnnnns(self):

    #def test_mask_4vnnnnn(self):
    #def test_mask_4nennnn(self):
    #def test_mask_4nnlnnn(self):
    #def test_mask_4nnnfnn(self):
    #def test_mask_4nnnngn(self):
    #def test_mask_4nnnnns(self):

    #def test_mask_5vnnnnn(self):
    #def test_mask_5nennnn(self):
    #def test_mask_5nnlnnn(self):
    #def test_mask_5nnnfnn(self):
    #def test_mask_5nnnngn(self):
    #def test_mask_5nnnnns(self):

    #def test_mask_6vnnnnn(self):
    #def test_mask_6nennnn(self):
    #def test_mask_6nnlnnn(self):
    #def test_mask_6nnnfnn(self):
    #def test_mask_6nnnngn(self):
    #def test_mask_6nnnnns(self):

    def test_init(self):
        self.assert_counts(0,0,0,0,0,0)

    def test_avert(self):
        v1,f1 = self.mesh.mbfv()
        v2 = self.mesh.avert()
        self.assert_counts(2,0,0,1,1,1)

    def test_aedge(self):
        v1,f1 = self.mesh.mbfv()
        v2 = self.mesh.avert()
        self.mesh.aedge(v1,v2)
        self.assert_counts(2,1,0,1,1,1)

    #def test_aloop(self):
    #def test_aface(self):
    #def test_afgroup(self):
    #def test_ashell(self):

    def test_mbfv(self):
        self.mesh.mbfv()
        self.assert_counts(1,0,0,1,1,1)

    def test_mev(self):
        v1,f1 = self.mesh.mbfv()
        v2,e1 = self.mesh.mev(v1)
        self.assert_counts(2,1,1,1,1,1)
        v3,e2 = self.mesh.mev(v2)
        self.assert_counts(3,2,1,1,1,1)
        v4,e3 = self.mesh.mev(v3)
        self.assert_counts(4,3,1,1,1,1)

    def test_mefl(self):
        v1,f1 = self.mesh.mbfv()
        v2,e1 = self.mesh.mev(v1)
        v3,e2 = self.mesh.mev(v2)
        v4,e3 = self.mesh.mev(v3)
        e4,f2 = self.mesh.mefl(v4,v1)
        self.assert_counts(4,4,2,2,1,1)

    def test_cube(self):
        self.cube_points()

        v1,f1 = self.mesh.mbfv(vgx = 0)
        v2,e1 = self.mesh.mev(v1,vgx = 1)
        v3,e2 = self.mesh.mev(v2,vgx = 2)
        v4,e3 = self.mesh.mev(v3,vgx = 3)
        v5,e4 = self.mesh.mev(v4,vgx = 7)
        v6,e5 = self.mesh.mev(v5,vgx = 6)
        v7,e6 = self.mesh.mev(v6,vgx = 5)
        v8,e7 = self.mesh.mev(v7,vgx = 4)
        e8,f2 = self.mesh.mefl(v8,v1)

        oe1 = self.mesh.hefrom(v6,f1)
        oe2 = self.mesh.hefrom(v3,f1)
        e9,f3 = self.mesh.mefl(v3,v6,oe1,oe2)
        oe1 = self.mesh.hefrom(v7,f1)
        oe2 = self.mesh.hefrom(v2,f1)
        e10,f4 = self.mesh.mefl(v2,v7,oe1,oe2)
        oe1 = self.mesh.hefrom(v1,f2)
        oe2 = self.mesh.hefrom(v4,f2)
        e11,f5 = self.mesh.mefl(v4,v1,oe1,oe2)
        oe1 = self.mesh.hefrom(v8,f2)
        oe2 = self.mesh.hefrom(v5,f2)
        e12,f6 = self.mesh.mefl(v5,v8,oe1,oe2)

        self.assert_counts(8,12,6,6,1,1)

        self.plotmesh()

if __name__ == '__main__':
    unittest.main()









 


