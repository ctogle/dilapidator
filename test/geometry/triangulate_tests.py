from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat

import dilap.geometry.triangulate as dtg

#import dilap.mesh.tools as dtl
import dilap.core.plotting as dtl
import matplotlib.pyplot as plt

import unittest,numpy,math,pdb
import matplotlib.pyplot as plt

#python3 -m unittest discover -v ./ "*tests.py"

class test_triangulate(unittest.TestCase):

    doshow = False

    def show_xy(self,tris,bnds):
        ax = dtl.plot_axes_xy()
        for t in tris:dtl.plot_polygon_xy(t,ax,True,2,'green')
        for b in bnds:dtl.plot_edges_xy(b,ax,True,5,'blue')
        if self.doshow:plt.show()
    def show(self,tris,bnds):
        ax = dtl.plot_axes()
        for t in tris:dtl.plot_polygon(t,ax,True,2,'green')
        for b in bnds:dtl.plot_edges(b,ax,True,5,'blue')
        if self.doshow:plt.show()

    def setUp(self):pass

    ###########################################################################
    ### extra triangulating functions
    ###########################################################################
    def split_nondelauney_edges(self,eb,ibs):
        aps = eb[:]
        for ib in ibs:aps.extend(ib)

        def split_loop(loop):
            x = 0
            while x < len(loop):
                p1,p2 = loop[x-1],loop[x]
                found = False
                for p in aps:
                    if p.isnear(p1) or p.isnear(p2):continue
                    cc = p1.mid(p2)
                    cr = cc.d(p1)
                    if p.inneighborhood(cc,cr):
                        loop.insert(x,cc)
                        found = True
                        break
                if found:continue
                else:x += 1
            return loop

        neb = split_loop(list(eb))
        nibs = [split_loop(list(ib)) for ib in ibs]
        return tuple(neb),tuple(tuple(nib) for nib in nibs)

    def split_nondelauney_edges_chew1(self,eb,ibs):
        els = [eb[x-1].d(eb[x]) for x in range(len(eb))]
        for ib in ibs:els.extend([ib[x-1].d(ib[x]) for x in range(len(ib))])
        hmin = min(els)*math.sqrt(3)

        def split_loop(loop):
            oloop = [loop[0]]
            for x in range(1,len(loop)+1):
                if x == len(loop):x = 0
                p1,p2 = oloop[-1],loop[x]
                el = p1.d(p2)
                m = 1
                while el/m > hmin:m += 1
                divpts = p1.pline(p2,m-1)
                if divpts:
                    for dvp in divpts:
                        oloop.append(dvp)
                if not x == 0:oloop.append(p2)

            return oloop

        neb = split_loop(list(eb))
        nibs = [split_loop(list(ib)) for ib in ibs]
        return hmin,tuple(neb),tuple(tuple(nib) for nib in nibs)

    def tribnd(self,eb,ibs,ref,smo):
        hmin = 1
        eb,ibs = self.split_nondelauney_edges(eb,ibs)
        if ref:hmin,eb,ibs = self.split_nondelauney_edges_chew1(eb,ibs)
        tris,bnds = dtg.triangulate(eb,ibs,hmin,ref,smo)
        return tris,bnds
    ###########################################################################
    ###########################################################################

    def test_squares(self):
        eb,ibs = tuple(vec3(5,1,-1).sq(2,4)),()
        tris,bnds = self.tribnd(eb,ibs,False,False)
        self.show_xy(tris,bnds)
        self.assertTrue(len(tris) == 2)
        self.assertTrue(len(bnds) == 4)
        tris,bnds = self.tribnd(eb,ibs,True,False)
        self.show_xy(tris,bnds)
        self.assertTrue(len(tris) == 4)
        self.assertTrue(len(bnds) == 6)
        eb,ibs = tuple(vec3(1,5,1).pring(2,4)),()
        tris,bnds = self.tribnd(eb,ibs,False,False)
        self.show_xy(tris,bnds)
        self.assertTrue(len(tris) == 2)
        self.assertTrue(len(bnds) == 4)
        tris,bnds = self.tribnd(eb,ibs,True,False)
        self.show_xy(tris,bnds)
        self.assertTrue(len(tris) == 2)
        self.assertTrue(len(bnds) == 4)

    def test_octagons(self):
        eb,ibs = tuple(vec3(1,-2,-1).pring(5,8)),()
        tris,bnds = self.tribnd(eb,ibs,False,False)
        self.show_xy(tris,bnds)
        self.assertTrue(len(tris) == 6)
        self.assertTrue(len(bnds) == 8)
        tris,bnds = self.tribnd(eb,ibs,True,False)
        self.show_xy(tris,bnds)
        self.assertTrue(len(tris) == 6)
        self.assertTrue(len(bnds) == 8)
        eb[0].trn(vec3(-7,2,0))
        tris,bnds = self.tribnd(eb,ibs,True,False)
        self.show_xy(tris,bnds)
        self.assertTrue(len(tris) == 8)
        self.assertTrue(len(bnds) == 10)



    def atest_triangulate(self):
        eb = (vec3(-2,-2,0),vec3(2,-2,0),vec3(2,2,0),vec3(-2,2,0))
        ibs = ()
        hmin,ref,smo = 1,False,False

        tris,bnds = dtg.triangulate(eb,ibs,hmin,ref,smo)
        self.show(tris,bnds)

if __name__ == '__main__':
    unittest.main()





