import dilap.core.tools as dpr
from dilap.geometry.vec3 import vec3

import dilap.topology.polygonmesh as dmsh

import dilap.mesh.tools as dtl
import matplotlib.pyplot as plt

import unittest,numpy,math

import pdb

#python3 -m unittest discover -v ./ "*tests.py"

class test_polygonmesh(unittest.TestCase):

    # could use a plotting method
    # could use a plotting method
    # could use a plotting method
    # could use a plotting method
    # could use a plotting method
    # could use a plotting method
    def plotmesh(self):
        mesh = self.mesh

        ax = dtl.plot_axes()
        for f in mesh.faces:
            ps = mesh.gvps(f)
            ax = dtl.plot_polygon(ps,ax)
        plt.show()

    def cube(self):
        self.v1 = self.mesh.avert(vec3(-1,-1,-1))
        self.v2 = self.mesh.avert(vec3( 1,-1,-1))
        self.v3 = self.mesh.avert(vec3( 1, 1,-1))
        self.v4 = self.mesh.avert(vec3(-1, 1,-1))
        self.v5 = self.mesh.avert(vec3(-1,-1, 1))
        self.v6 = self.mesh.avert(vec3( 1,-1, 1))
        self.v7 = self.mesh.avert(vec3( 1, 1, 1))
        self.v8 = self.mesh.avert(vec3(-1, 1, 1))
        bot = [self.v1,self.v2,self.v3,self.v4]
        top = [self.v5,self.v6,self.v7,self.v8]
        self.l1 = self.mesh.aloop(bot)
        self.l2 = self.mesh.aloop(top)
        self.f1 = self.mesh.aface(self.l1,())
        self.f2 = self.mesh.aface(self.l2,())

    def setUp(self):
        self.mesh = dmsh.polygonmesh()

    def test_init(self):
        self.assertEqual(self.mesh.vcnt(),1)
        self.assertEqual(self.mesh.ecnt(),1)
        self.assertEqual(self.mesh.lcnt(),1)
        self.assertEqual(self.mesh.fcnt(),1)

    def test_avert(self):
        self.mesh.avert(self.mesh.verts[0].p.cp())
        self.assertEqual(self.mesh.vcnt(),2)
        self.assertEqual(self.mesh.ecnt(),1)
        self.assertEqual(self.mesh.lcnt(),1)
        self.assertEqual(self.mesh.fcnt(),1)

    def test_sedge(self):
        v1,e1 = self.mesh.verts[0],self.mesh.edges[0]

        v2 = self.mesh.avert(vec3(2,0,0))
        e1,e2 = self.mesh.sedge(e1,v2)
        v3 = self.mesh.avert(vec3(2,2,0))
        e2,e3 = self.mesh.sedge(e2,v3)
        v4 = self.mesh.avert(vec3(0,2,0))
        e3,e4 = self.mesh.sedge(e3,v4)

        self.assertEqual(self.mesh.vcnt(),4)
        self.assertEqual(self.mesh.ecnt(),4)
        self.assertEqual(self.mesh.lcnt(),1)
        self.assertEqual(self.mesh.fcnt(),1)

        e5,e6,l2,f2 = self.mesh.mefl(v1,v2)
        v5 = self.mesh.avert(vec3(2,0,2))
        e6,e7 = self.mesh.sedge(e6,v5)
        v6 = self.mesh.avert(vec3(0,0,2))
        e7,e8 = self.mesh.sedge(e7,v6)

        self.assertEqual(self.mesh.vcnt(),6)
        self.assertEqual(self.mesh.ecnt(),8)
        self.assertEqual(self.mesh.lcnt(),2)
        self.assertEqual(self.mesh.fcnt(),2)

        self.plotmesh()

    #def test_mefl(self):
    #def test_mask(self):

    #def test_gvps(self):

    def test_cube(self):
        #self.cube()
        pass

if __name__ == '__main__':
    unittest.main()









 


