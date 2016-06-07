from dilap.geometry.vec3 import vec3
import dilap.geometry.tools as gtl
import dilap.geometry.polymath as pym

import dilap.topology.partitiongraph as ptg

import dilap.worldly.blockletters as dbl

import dilap.core.plotting as dtl
import matplotlib.pyplot as plt

import unittest,random

import pdb

#python3 -m unittest discover -v ./ "*tests.py"

class test_partitiongraph(unittest.TestCase):

    def test_av(self):
        fp = dbl.block('H',10,25,25)
        rg = ptg.partitiongraph()
        bs = pym.bsegsxy(fp[0],vec3(0,-100,0),vec3(0,100,0))
        for b in bs:
            i = rg.av(b = [b,[]],p = vec3(0,0,0).com(b),l = 0)
        rg.plotxy()
        plt.show()

        pg = rg.bgraph()

        #pdb.set_trace()

        pgpy = pg.polygon(1,'ccw')
        ax = pg.plotxy()
        #ax = dtl.plot_axes_xy(20)
        #ax = dtl.plot_polygon_full_xy(pgpy,ax,lw = 2)
        ax = dtl.plot_polygon_xy(pgpy[0],ax,lw = 2)
        for pgp in pgpy[1]:
            ax = dtl.plot_polygon_xy(pgp,ax,lw = 2)
        plt.show()

    def atest_split(self):
        fp = dbl.block('H',10,25,25)
        rg = ptg.partitiongraph()
        bs = pym.bsegsxy(fp[0],vec3(0,-100,0),vec3(0,100,0))
        b = fp[0]
        i1 = rg.av(b = [b,[]],p = vec3(0,0,0).com(b),l = 0)
        i2 = rg.sv(0,bs[0],bs[1])
        rg.plotxy()
        plt.show()

    def atest_break(self):
        fp = dbl.block('H',10,25,25)
        v1 = {'b':fp,'p':vec3(0,0,0).com(fp[0]),'l':0}
        rg = ptg.partitiongraph()
        ip = vec3(0.5,0.5,0)
        i1 = rg.av(**v1)
        i2 = rg.bv(0,ip,vec3(1,0,0))
        i3 = rg.bv(0,ip,vec3(0,1,0))
        i4 = rg.bv(1,ip,vec3(0,1,0))
        rg.vves()
        rg.plotxy()
        plt.show()

    def atest_bgraph(self):
        fp = dbl.block('H',10,25,25)
        v1 = {'b':fp,'p':vec3(0,0,0).com(fp[0]),'l':0}
        rg = ptg.partitiongraph()
        ip = vec3(0.5,0.5,0)
        i1 = rg.av(**v1)
        #i2 = rg.bv(0,ip,vec3(1,0,0))
        i3 = rg.bv(0,ip,vec3(0,1,0))
        #i4 = rg.bv(1,ip,vec3(0,1,0))
        rg.vves()

        pg = rg.bgraph()

        #pdb.set_trace()

        pgpy = pg.polygon(1,'ccw')
        ax = pg.plotxy()
        #ax = dtl.plot_axes_xy(20)
        #ax = dtl.plot_polygon_full_xy(pgpy,ax,lw = 2)
        ax = dtl.plot_polygon_xy(pgpy[0],ax,lw = 2)
        for pgp in pgpy[1]:
            ax = dtl.plot_polygon_xy(pgp,ax,lw = 2)
        plt.show()

###############################################################################

if __name__ == '__main__':unittest.main()

###############################################################################



