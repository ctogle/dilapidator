import dilap.core.tools as dpr

from dilap.geometry.vec3 import vec3

import dilap.modeling.model as dmo

import dilap.mesh.tools as dtl
import matplotlib.pyplot as plt

import unittest,numpy,math

import pdb

#python3 -m unittest discover -v ./ "*tests.py"

class test_model(unittest.TestCase):

    def plot(self,mesh):
        ax = dtl.plot_axes()
        for f in mesh.faces:
            ps = self.mod.gvps(mesh,f)
            ax = dtl.plot_polygon(ps,ax)
        plt.show()

    def quad(self):
        gmesh = self.mod.agfxmesh()
        self.v1  = gmesh.avert(*self.mod.avert(vec3(-1,-1,-1)))
        self.v2  = gmesh.avert(*self.mod.avert(vec3( 1,-1,-1)))
        self.v3  = gmesh.avert(*self.mod.avert(vec3( 1, 1,-1)))
        self.v4  = gmesh.avert(*self.mod.avert(vec3(-1, 1,-1)))
        self.f1 = gmesh.aface(self.v1,self.v2,self.v3) 
        self.f2 = gmesh.aface(self.v1,self.v3,self.v4) 
        return gmesh

    def dome(self):
        gmesh = self.mod.agfxmesh()
        self.v1  = gmesh.avert(*self.mod.avert(vec3(-1,-1,-1)))
        self.v2  = gmesh.avert(*self.mod.avert(vec3( 1,-1,-1)))
        self.v3  = gmesh.avert(*self.mod.avert(vec3( 1, 1,-1)))
        self.v4  = gmesh.avert(*self.mod.avert(vec3(-1, 1,-1)))
        self.v5  = gmesh.avert(*self.mod.avert(vec3(-1,-1, 1)))
        self.v6  = gmesh.avert(*self.mod.avert(vec3( 1,-1, 1)))
        self.v7  = gmesh.avert(*self.mod.avert(vec3( 1, 1, 1)))
        self.v8  = gmesh.avert(*self.mod.avert(vec3(-1, 1, 1)))
        self.f1  = gmesh.aface(self.v1,self.v3,self.v2) 
        self.f2  = gmesh.aface(self.v1,self.v4,self.v3) 
        self.f3  = gmesh.aface(self.v5,self.v6,self.v7) 
        self.f4  = gmesh.aface(self.v5,self.v7,self.v8) 
        self.f5  = gmesh.aface(self.v1,self.v2,self.v6) 
        self.f6  = gmesh.aface(self.v1,self.v6,self.v5) 
        self.f7  = gmesh.aface(self.v2,self.v3,self.v7) 
        self.f8  = gmesh.aface(self.v2,self.v7,self.v6) 
        self.f9  = gmesh.aface(self.v3,self.v4,self.v8) 
        self.f10 = gmesh.aface(self.v3,self.v8,self.v7) 
        return gmesh

    def cube(self):
        gmesh = self.mod.agfxmesh()
        self.v1  = gmesh.avert(*self.mod.avert(vec3(-1,-1,-1)))
        self.v2  = gmesh.avert(*self.mod.avert(vec3( 1,-1,-1)))
        self.v3  = gmesh.avert(*self.mod.avert(vec3( 1, 1,-1)))
        self.v4  = gmesh.avert(*self.mod.avert(vec3(-1, 1,-1)))
        self.v5  = gmesh.avert(*self.mod.avert(vec3(-1,-1, 1)))
        self.v6  = gmesh.avert(*self.mod.avert(vec3( 1,-1, 1)))
        self.v7  = gmesh.avert(*self.mod.avert(vec3( 1, 1, 1)))
        self.v8  = gmesh.avert(*self.mod.avert(vec3(-1, 1, 1)))
        self.f1  = gmesh.aface(self.v1,self.v3,self.v2) 
        self.f2  = gmesh.aface(self.v1,self.v4,self.v3) 
        self.f3  = gmesh.aface(self.v5,self.v6,self.v7) 
        self.f4  = gmesh.aface(self.v5,self.v7,self.v8) 
        self.f5  = gmesh.aface(self.v1,self.v2,self.v6) 
        self.f6  = gmesh.aface(self.v1,self.v6,self.v5) 
        self.f7  = gmesh.aface(self.v2,self.v3,self.v7) 
        self.f8  = gmesh.aface(self.v2,self.v7,self.v6) 
        self.f9  = gmesh.aface(self.v3,self.v4,self.v8) 
        self.f10 = gmesh.aface(self.v3,self.v8,self.v7) 
        self.f11 = gmesh.aface(self.v4,self.v1,self.v5) 
        self.f12 = gmesh.aface(self.v4,self.v5,self.v8) 
        return gmesh

    def polycube(self):
        pmesh = self.mod.apolymesh()
        gm = self.mod.gfx(pmesh)
        
        self.plot(gm)

    def assert_tmcounts(self,mesh,vcnt,ecnt,fcnt):
        self.assertEqual(mesh.vcnt(),vcnt)
        self.assertEqual(mesh.ecnt(),ecnt)
        self.assertEqual(mesh.fcnt(),fcnt)

    def assert_pmcounts(self,mesh,vcnt,ecnt,lcnt,fcnt):
        self.assert_tmcounts(mesh,vcnt,ecnt,fcnt)
        self.assertEqual(mesh.lcnt(),lcnt)

    def setUp(self):
        self.mod = dmo.model()

    def test_init(self):
        self.assertTrue(self.mod.__str__() == 'model:')
        #self.assert_counts(0,0,0)

    def test_agfxmesh(self):
        gmesh = self.mod.agfxmesh()
        self.assert_tmcounts(gmesh,0,0,0)

    def test_acolmesh(self):
        cmesh = self.mod.acolmesh()
        self.assert_tmcounts(cmesh,0,0,0)

    def test_alodmesh(self):
        lmesh = self.mod.alodmesh()
        self.assert_tmcounts(lmesh,0,0,0)

    def test_apolymesh(self):
        pmesh = self.mod.apolymesh()
        self.assert_pmcounts(pmesh,1,1,1,1)

    #def test_avert(self):
    #def test_gvps(self,vxs):

    def test_subdiv_cube(self):
        gm = self.cube()
        self.assert_tmcounts(gm,8,36,12)
        self.mod.subdiv(gm)
        self.assert_tmcounts(gm,20,108,36)
        self.mod.subdiv(gm)
        self.assert_tmcounts(gm,56,324,108)

    def test_subdiv_dome(self):
        gm = self.dome()
        self.assert_tmcounts(gm,8,30,10)
        self.mod.subdiv(gm)
        self.assert_tmcounts(gm,18,90,30)

    def test_subdiv_quad(self):
        gm = self.quad()
        self.assert_tmcounts(gm,4,6,2)
        self.mod.subdiv(gm)
        self.assert_tmcounts(gm,6,18,6)

    #def test_trn(self):
    #def test_rot(self):
    #def test_scl(self):

    def test_cube(self):
        gm = self.cube()
        self.assert_tmcounts(gm,8,36,12)

if __name__ == '__main__':
    unittest.main()








 

