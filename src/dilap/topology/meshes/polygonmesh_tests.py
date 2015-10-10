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
                hes = mesh.mask(4,None,None,l,None)
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

    def setUp(self):
        self.mesh = dmsh.polygonmesh()
        self.pset = pointset()

    #def test_mask_0nnln(self):

    def test_mask_1vnnn(self):
        v1,f1 = self.mesh.mfv(vgx = 0)
        v2,e1 = self.mesh.mev(v1,vgx = 0)
        v3,e2 = self.mesh.mev(v2,vgx = 0)
        v4,e3 = self.mesh.mev(v3,vgx = 0)
        v5,e4 = self.mesh.mev(v4,vgx = 0)
        v6,e5 = self.mesh.mev(v5,vgx = 0)
        v7,e6 = self.mesh.mev(v6,vgx = 0)
        v8,e7 = self.mesh.mev(v7,vgx = 0)
        e8,f2 = self.mesh.mfe(v8,v1)

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

    def test_mask_2vnnn(self):
        self.pset.ap(vec3(0,0,0));self.pset.ap(vec3(2,0,0))
        self.pset.ap(vec3(2,2,0));self.pset.ap(vec3(0,2,0))
        self.pset.ap(vec3(0,0,2));self.pset.ap(vec3(2,0,2))
        self.pset.ap(vec3(2,2,2));self.pset.ap(vec3(0,2,2))

        v1,f1 = self.mesh.mfv(vgx = 0)
        v2,e1 = self.mesh.mev(v1,vgx = 1)
        v3,e2 = self.mesh.mev(v2,vgx = 2)
        v4,e3 = self.mesh.mev(v3,vgx = 3)
        v5,e4 = self.mesh.mev(v4,vgx = 7)
        v6,e5 = self.mesh.mev(v5,vgx = 6)
        v7,e6 = self.mesh.mev(v6,vgx = 5)
        v8,e7 = self.mesh.mev(v7,vgx = 4)
        e8,f2 = self.mesh.mfe(v8,v1)
        v3hes = self.mesh.mask(4,v3,None,None,f1)
        v6hes = self.mesh.mask(4,v6,None,None,f1)
        oe1,oe2 = v6hes[1],v3hes[1]
        e9,f3 = self.mesh.mfe(v3,v6,oe1,oe2)

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

    def test_mask_2nenn(self):
        v1,f1 = self.mesh.mfv()
        v2,e1 = self.mesh.mev(v1)
        v3,e2 = self.mesh.mev(v2)
        v4,e3 = self.mesh.mev(v3)
        mask1 = self.mesh.mask(2,None,e1)
        mask2 = self.mesh.mask(2,None,e2)
        mask3 = self.mesh.mask(2,None,e3)
        self.assertEqual(mask1,mask2)
        self.assertEqual(mask1,mask3)

    #def test_mask_2nnnf(self):

    #def test_mask_4nnln(self):
    #def test_mask_4vnln(self):

    def test_init(self):
        self.assertEqual(self.mesh.vcnt(),0)
        self.assertEqual(self.mesh.ecnt(),0)
        self.assertEqual(self.mesh.lcnt(),0)
        self.assertEqual(self.mesh.fcnt(),0)

    def test_avert(self):
        v1,f1 = self.mesh.mfv()
        v2 = self.mesh.avert()
        self.assertEqual(self.mesh.vcnt(),2)
        self.assertEqual(self.mesh.ecnt(),0)
        self.assertEqual(self.mesh.lcnt(),0)
        self.assertEqual(self.mesh.fcnt(),1)

    def test_aedge(self):
        v1,f1 = self.mesh.mfv()
        v2 = self.mesh.avert()
        self.mesh.aedge(v1,v2)
        self.assertEqual(self.mesh.vcnt(),2)
        self.assertEqual(self.mesh.ecnt(),1)
        self.assertEqual(self.mesh.lcnt(),0)
        self.assertEqual(self.mesh.fcnt(),1)

    #def test_aloop(self):
    #def test_aface(self):

    def test_mfv(self):
        self.mesh.mfv()
        self.assertEqual(self.mesh.vcnt(),1)
        self.assertEqual(self.mesh.ecnt(),0)
        self.assertEqual(self.mesh.lcnt(),0)
        self.assertEqual(self.mesh.fcnt(),1)

    def test_mev(self):
        v1,f1 = self.mesh.mfv()
        v2,e1 = self.mesh.mev(v1)
        self.assertEqual(self.mesh.vcnt(),2)
        self.assertEqual(self.mesh.ecnt(),1)
        self.assertEqual(self.mesh.lcnt(),1)
        self.assertEqual(self.mesh.fcnt(),1)
        v3,e2 = self.mesh.mev(v2)
        self.assertEqual(self.mesh.vcnt(),3)
        self.assertEqual(self.mesh.ecnt(),2)
        self.assertEqual(self.mesh.lcnt(),1)
        self.assertEqual(self.mesh.fcnt(),1)
        v4,e3 = self.mesh.mev(v3)
        self.assertEqual(self.mesh.vcnt(),4)
        self.assertEqual(self.mesh.ecnt(),3)
        self.assertEqual(self.mesh.lcnt(),1)
        self.assertEqual(self.mesh.fcnt(),1)

    def test_mfe(self):
        v1,f1 = self.mesh.mfv()
        v2,e1 = self.mesh.mev(v1)
        v3,e2 = self.mesh.mev(v2)
        v4,e3 = self.mesh.mev(v3)
        e4,f2 = self.mesh.mfe(v4,v1)
        self.assertEqual(self.mesh.vcnt(),4)
        self.assertEqual(self.mesh.ecnt(),4)
        self.assertEqual(self.mesh.lcnt(),2)
        self.assertEqual(self.mesh.fcnt(),2)

    def test_cube(self):
        def getoe(v,f):
            hes = self.mesh.mask(4,v,None,None,f)
            oe = tuple(x for x in hes if x.one is v)[0]
            return oe

        self.pset.ap(vec3(0,0,0));self.pset.ap(vec3(2,0,0))
        self.pset.ap(vec3(2,2,0));self.pset.ap(vec3(0,2,0))
        self.pset.ap(vec3(0,0,2));self.pset.ap(vec3(2,0,2))
        self.pset.ap(vec3(2,2,2));self.pset.ap(vec3(0,2,2))

        v1,f1 = self.mesh.mfv(vgx = 0)
        v2,e1 = self.mesh.mev(v1,vgx = 1)
        v3,e2 = self.mesh.mev(v2,vgx = 2)
        v4,e3 = self.mesh.mev(v3,vgx = 3)
        v5,e4 = self.mesh.mev(v4,vgx = 7)
        v6,e5 = self.mesh.mev(v5,vgx = 6)
        v7,e6 = self.mesh.mev(v6,vgx = 5)
        v8,e7 = self.mesh.mev(v7,vgx = 4)
        e8,f2 = self.mesh.mfe(v8,v1)

        oe1 = getoe(v6,f1);oe2 = getoe(v3,f1)
        e9,f3 = self.mesh.mfe(v3,v6,oe1,oe2)
        oe1 = getoe(v7,f1);oe2 = getoe(v2,f1)
        e10,f4 = self.mesh.mfe(v2,v7,oe1,oe2)
        oe1 = getoe(v1,f2);oe2 = getoe(v4,f2)
        e11,f5 = self.mesh.mfe(v4,v1,oe1,oe2)
        oe1 = getoe(v8,f2);oe2 = getoe(v5,f2)
        e12,f6 = self.mesh.mfe(v5,v8,oe1,oe2)

        self.assertEqual(self.mesh.vcnt(),8)
        self.assertEqual(self.mesh.ecnt(),12)
        self.assertEqual(self.mesh.lcnt(),6)
        self.assertEqual(self.mesh.fcnt(),6)

        self.plotmesh()

if __name__ == '__main__':
    unittest.main()









 


