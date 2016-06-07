from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat
import dilap.geometry.tools as gtl
import dilap.geometry.polymath as pym

import dilap.topology.wiregraph as pgr

import dilap.core.plotting as dtl
import matplotlib.pyplot as plt

import unittest,random,numpy

import pdb

#python3 -m unittest discover -v ./ "*tests.py"

###############################################################################

class test_wiregraph(unittest.TestCase):

    def test_avrvaerese(self):
        wg = pgr.wiregraph()
        i1 = wg.av(**{})
        i2 = wg.av(**{})
        i3 = wg.av(**{})
        i4 = wg.av(**{})
        r1 = wg.ae(i1,i2)
        r2 = wg.ae(i2,i3)
        r3 = wg.ae(i3,i4)
        r4 = wg.ae(i4,i1)
        r2 = wg.re(i2,i3)
        i5 = wg.av(**{})
        r5,r6 = wg.se(i1,i4,i5)
        i1 = wg.rv(i1)

    def test_orings(self):
        wg = pgr.wiregraph()
        i1 = wg.av()
        i2 = wg.av()
        i3 = wg.av()
        i4 = wg.av()
        r1 = wg.ae(i1,i2)
        r2 = wg.ae(i1,i3)
        r3 = wg.ae(i1,i4)
        self.assertEqual(wg.orings[0],[1,3,2])

    def test_mev(self):
        wg = pgr.wiregraph()
        i1 = wg.av(**{})
        i2,r1 = wg.mev(i1,{},{})
        i3,r2 = wg.mev(i2,{},{})

    def atest_loop(self):
        def pl(il):
            ilp = [rg.vs[j][1]['p'] for j in il]
            ilp = pym.contract(ilp,2.0)
            ax = rg.plot()
            ax = dtl.plot_polygon(ilp,ax,col = 'b',lw = 4)
            plt.show()
        rg = pgr.wiregraph()
        i1 = rg.av(p = vec3( 10,-5,0),l = 0)
        i2 = rg.av(p = vec3( 10, 5,0),l = 0)
        i3 = rg.av(p = vec3(-10, 5,0),l = 0)
        i4 = rg.av(p = vec3(-10,-5,0),l = 0)
        r1 = rg.ae(i1,i2)
        r2 = rg.ae(i2,i3)
        r3 = rg.ae(i3,i4)
        r4 = rg.ae(i4,i1)
        i5 = rg.av(p = vec3(2,-10,0),l = 0)
        r5,r6 = rg.se(i1,i4,i5)
        i6 = rg.av(p = vec3(-2,10,0),l = 0)
        r7,r8 = rg.se(i2,i3,i6)
        r9 = rg.ae(i5,i6)

        il = rg.loop(i5,i6,'cw')
        self.assertEqual(il,[i5,i6,i2,i1])
        #pl(il)
        il = rg.loop(i5,i6,'ccw')
        self.assertEqual(il,[i5,i6,i3,i4])
        #pl(il)
        il = rg.loop(i5,i1,'cw')
        self.assertEqual(il,[i5,i1,i2,i6,i3,i4])
        #pl(il)
        il = rg.loop(i5,i1,'ccw')
        self.assertEqual(il,[i5,i1,i2,i6])
        #pl(il)
        il = rg.loop(i1,i5,'cw')
        self.assertEqual(il,[i1,i5,i6,i2])
        #pl(il)

        i7,r10 = rg.mev(i1,{'p':vec3(12,-20,0),'l':0},{})

        il = rg.loop(i5,i6,'cw')
        self.assertEqual(il,[i5,i6,i2,i1])
        #pl(il)
        il = rg.loop(i5,i6,'ccw')
        self.assertEqual(il,[i5,i6,i3,i4])
        #pl(il)
        il = rg.loop(i7,i1,'cw')
        self.assertEqual(il,[i7,i1,i2,i6,i3,i4,i5,i1])
        #pl(il)
        il = rg.loop(i7,i1,'ccw')
        self.assertEqual(il,[i7,i1,i5,i4,i3,i6,i2,i1])
        #pl(il)

        i8,r11 = rg.mev(i3,{'p':vec3(-5,0,0),'l':0},{})

        il = rg.loop(i3,i4,'ccw')
        self.assertEqual(il,[i3,i4,i5,i6,i3,i8])
        #pl(il)
        il = rg.loop(i3,i4,'cw')
        self.assertEqual(il,[i3,i4,i5,i1,i7,i1,i2,i6])
        #pl(il)

    def atest_uloops(self):
        def pl():
            ax = rg.plot()
            for lp in loops:
                lpps = [rg.vs[j][1]['p'] for j in lp]
                lpps = pym.contract(lpps,2)
                ax = dtl.plot_polygon(lpps,ax,lw = 3,col = 'b')
            plt.show()
        rg = pgr.wiregraph()
        i1 = rg.av(p = vec3( 10,-5,0),l = 0)
        i2 = rg.av(p = vec3( 10, 5,0),l = 0)
        i3 = rg.av(p = vec3(-10, 5,0),l = 0)
        i4 = rg.av(p = vec3(-10,-5,0),l = 0)
        r1 = rg.ae(i1,i2)
        r2 = rg.ae(i2,i3)
        r3 = rg.ae(i3,i4)
        r4 = rg.ae(i4,i1)
        i5 = rg.av(p = vec3(2,-10,0),l = 0)
        r5,r6 = rg.se(i1,i4,i5)
        i6 = rg.av(p = vec3(-2,10,0),l = 0)
        r7,r8 = rg.se(i2,i3,i6)
        r9 = rg.ae(i5,i6)
        loops = rg.uloops('ccw')

        #pl()
        self.assertEqual(len(loops),3)

        i7,r10 = rg.mev(i1,{'p':vec3(12,-20,0),'l':0},{})
        loops = rg.uloops('ccw')

        #pl()
        self.assertEqual(len(loops),3)

###############################################################################

if __name__ == '__main__':unittest.main()

###############################################################################



 

