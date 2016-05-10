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

    doshow = True

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

    def tribnd(self,eb,ibs,ref,smo):
        hmin = 1
        eb,ibs = dtg.split_nondelauney_edges(eb,ibs)
        if ref:hmin,eb,ibs = dtg.split_nondelauney_edges_chew1(eb,ibs)
        tris,bnds = dtg.triangulate(eb,ibs,hmin,ref,smo)
        return tris,bnds

    ###########################################################################
    ###########################################################################

    def test_tinbxy(self):
        a = vec3(6.0, 5.0, 0.0)
        b = vec3(6.75, 4.0, 0.0)
        c = vec3(7.5, 5.0, 0.0)
        bnd = (
            vec3(0.0, 0.0, 0.0),vec3(4.25, 0.0, 0.0),
            vec3(4.25, 4.0, 0.0),vec3(6.75, 4.0, 0.0),
            vec3(6.75, 0.0, 0.0),vec3(12.0, 0.0, 0.0),
            vec3(12.0, 5.0, 0.0),vec3(9.0, 5.0, 0.0),
            vec3(7.5, 5.0, 0.0),vec3(6.0, 5.0, 0.0),
            vec3(4.5, 5.0, 0.0),vec3(3.0, 5.0, 0.0),
            vec3(0.0, 5.0, 0.0))
        #ax = dtl.plot_axes_xy(12)   
        #ax = dtl.plot_polygon_xy(bnd,ax)
        #ax = dtl.plot_polygon_xy((a,b,c),ax,col = 'r',lw = 2)
        #plt.show()
        self.assertTrue(dtg.tinbxy(a,b,c,bnd))

    def test_door(self):
        eb = vec3(6,2.5,0).sq(12,5)
        eb.insert(1,vec3(6.75,0,0))
        eb.insert(1,vec3(6.75,4,0))
        eb.insert(1,vec3(4.25,4,0))
        eb.insert(1,vec3(4.25,0,0))

        eb,ibs = tuple(eb),()
        tris,bnds = self.tribnd(eb,ibs,False,False)
        self.show(tris,bnds)

    def atest_doorwindow(self):
        eb = vec3(0,0,0).sq(8,5)

        ddp = 1.5/(2.0*8.0)
        w1 = eb[0].lerp(eb[1],0.5-ddp)
        w2 = eb[0].lerp(eb[1],0.5+ddp)
        w = (w2,w2.cp().ytrn(3),w1.cp().ytrn(3),w1)
        for wp in w:eb.insert(1,wp)

        ddp = 1.5/(2.0*8.0)
        w1 = eb[0].lerp(eb[5],0.25-ddp).ytrn(1)
        w2 = eb[0].lerp(eb[5],0.25+ddp).ytrn(1)
        w = (w2,w2.cp().ytrn(2),w1.cp().ytrn(2),w1)
        #q = quat(0,0,0,0).uu(vec3(0,0,1),vec3(1,0,0))
        #vec3(0,0,0).fulc(q,eb)

        eb,ibs = tuple(eb),(tuple(w),)
        tris,bnds = self.tribnd(eb,ibs,False,False)
        self.show(tris,bnds)

    def atest_nonxy(self):
        q = quat(0,0,0,0).uu(vec3(0,0,1),vec3(1,0,0))
        eb = vec3(5,1,-1).sq(2,4)
        vec3(0,0,0).fulc(q,eb)
        eb,ibs = tuple(eb),()
        tris,bnds = self.tribnd(eb,ibs,False,False)
        self.show(tris,bnds)

    def atest_squares(self):
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

    def atest_octagons(self):
        eb,ibs = tuple(vec3(1,-2,-1).pring(5,8)),()
        '''#
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
        eb[0].trn(vec3(7,-2,0))
        '''#
        ibs = (tuple(vec3(1,-2,-1).pring(1,8)),)
        #tris,bnds = self.tribnd(eb,ibs,False,False)
        #self.show_xy(tris,bnds)
        tris,bnds = self.tribnd(eb,ibs,False,False)
        self.show_xy(tris,bnds)

    def atest_ushape(self):
        eb = (
            vec3(0,0,0),vec3(1,0,0),vec3(1,1,0),vec3(0.75,1,0),
            vec3(0.75,0.5,0),vec3(0.25,0.5,0),vec3(0.25,1,0),vec3(0,1,0))
        ibs = (tuple(p.cp().uscl(0.35).trn(vec3(0.1,0.1,0)) for p in eb),)

        ax = dtl.plot_polygon_full_xy((eb,ibs))
        plt.show()

        tris,bnds = self.tribnd(eb,(),False,False)
        self.show_xy(tris,bnds)
        tris,bnds = self.tribnd(ibs[0],(),False,False)
        self.show_xy(tris,bnds)

        tris,bnds = self.tribnd(eb,ibs,False,False)
        self.show_xy(tris,bnds)

if __name__ == '__main__':
    unittest.main()





