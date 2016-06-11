from dilap.geometry.quat import quat
from dilap.geometry.vec3 import vec3
import dilap.geometry.tools as gtl
import dilap.geometry.polymath as pym

import dilap.worldly.blockletters as dbl

import dilap.core.plotting as dtl
import matplotlib.pyplot as plt

import unittest,numpy,math

import pdb

#python3 -m unittest discover -v ./ "*tests.py"

class test_polymath(unittest.TestCase):

    #def setUp(self):

    def test_sintsxy(self):
        def perm(s1,s2,s3,s4,c,**kws):
            f = self.assertTrue if c else self.assertFalse
            f(pym.sintsxy(s1,s2,s3,s4,**kws));f(pym.sintsxy(s1,s2,s4,s3,**kws))
            f(pym.sintsxy(s2,s1,s3,s4,**kws));f(pym.sintsxy(s2,s1,s4,s3,**kws))
            f(pym.sintsxy(s3,s4,s1,s2,**kws));f(pym.sintsxy(s3,s4,s2,s1,**kws))
            f(pym.sintsxy(s4,s3,s1,s2,**kws));f(pym.sintsxy(s4,s3,s2,s1,**kws))
        def pl(s11,s12,s21,s22):
            ax = dtl.plot_axes_xy(1.5)
            ax = dtl.plot_point_xy(s11,dtl.plot_point_xy_annotate(s11,ax,'s11'),col = 'b')
            ax = dtl.plot_point_xy(s12,dtl.plot_point_xy_annotate(s12,ax,'s12'),col = 'b')
            ax = dtl.plot_point_xy(s21,dtl.plot_point_xy_annotate(s21,ax,'s21'),col = 'g')
            ax = dtl.plot_point_xy(s22,dtl.plot_point_xy_annotate(s22,ax,'s22'),col = 'g')
            ax = dtl.plot_vector_xy(s11,s11.tov(s21),ax,lw = 5.0,col = 'r')
            ax = dtl.plot_vector_xy(s11,s11.tov(s12),ax,lw = 2.0,col = 'b')
            ax = dtl.plot_vector_xy(s21,s21.tov(s22),ax,lw = 2.0,col = 'g')
            plt.show()

        s1 = vec3(8.0, -3.313708782196045, 0.0)
        s2 = vec3(8.0, 3.313708782196045, 0.0)
        s3 = vec3(5.65685510635376, -3.999999761581421, 0.0)
        s4 = vec3(8.0, -1.6568543910980225, 0.0)
        perm(s1,s2,s3,s4,True)
        perm(s1,s2,s3,s4,False,ie = False)
        #pl(s1,s2,s3,s4)

        s1,s2,s3,s4 = vec3(-1,0,0),vec3(1,1,0),vec3(-1,1,0),vec3(-1,-1,0)
        perm(s1,s2,s3,s4,True)
        perm(s1,s2,s3,s4,False,ie = False)
        #pl(s1,s2,s3,s4)
        s1,s2,s3,s4 = vec3(-1,0.75,0),vec3(1,1,0),vec3(-1,1,0),vec3(-1,-1,0)
        perm(s1,s2,s3,s4,True)
        perm(s1,s2,s3,s4,False,ie = False)
        #pl(s1,s2,s3,s4)
        s1,s2,s3,s4 = vec3(-1,0,0),vec3(1,1,0),vec3(-1,0,0),vec3(-1,-1,0)
        perm(s1,s2,s3,s4,True)
        perm(s1,s2,s3,s4,False,ie = False)
        #pl(s1,s2,s3,s4)
        s1,s2,s3,s4 = vec3(-1,0,0),vec3(1,1,0),vec3(0,1,0),vec3(0,-1,0)
        perm(s1,s2,s3,s4,True)
        #pl(s1,s2,s3,s4)
        s1,s2,s3,s4 = vec3(-1,1,0),vec3(1,1,0),vec3(0,0,0),vec3(0,-1,0)
        perm(s1,s2,s3,s4,False)
        #pl(s1,s2,s3,s4)

        s1,s2,s3,s4 = vec3(-1,1,0),vec3(1,1,0),vec3(0,0,0),vec3(2,0,0)
        perm(s1,s2,s3,s4,False)
        #pl(s1,s2,s3,s4)
        s1,s2,s3,s4 = vec3(-1,1,0),vec3(1,1,0),vec3(3,1,0),vec3(2,1,0)
        perm(s1,s2,s3,s4,False)
        #pl(s1,s2,s3,s4)
        s1,s2,s3,s4 = vec3(-1,1,0),vec3(1,1,0),vec3(0,1,0),vec3(2,1,0)
        perm(s1,s2,s3,s4,True)
        #pl(s1,s2,s3,s4)
        s1,s2,s3,s4 = vec3(-1,1,0),vec3(1,1,0),vec3(1,1,0),vec3(2,1,0)
        perm(s1,s2,s3,s4,True)
        perm(s1,s2,s3,s4,False,ie = False)
        #pl(s1,s2,s3,s4)
        s1,s2,s3,s4 = vec3(-1,1,0),vec3(1,1,0),vec3(-1,1,0),vec3(1,1,0)
        perm(s1,s2,s3,s4,True)
        perm(s1,s2,s3,s4,False,ieb = False)
        #pl(s1,s2,s3,s4)
        s1,s2,s3,s4 = vec3(-1,1,0),vec3(1,1,0),vec3(0,1,0),vec3(0.5,1,0)
        perm(s1,s2,s3,s4,True)
        perm(s1,s2,s3,s4,True,ie = False)
        perm(s1,s2,s3,s4,True,ieb = False)
        perm(s1,s2,s3,s4,True,ie = False,ieb = False)
        #pl(s1,s2,s3,s4)
        s1,s2,s3,s4 = vec3(-1,0,0),vec3(0,0,0),vec3(-1,0,0),vec3(1,0,0)
        perm(s1,s2,s3,s4,True)
        perm(s1,s2,s3,s4,True,ie = False)
        #pl(s1,s2,s3,s4)

    def test_sintsxyp(self):
        def pl(s1,s2,s3,s4):
            ax = dtl.plot_axes_xy(10)
            ax = dtl.plot_edges_xy([s1,s2],ax,lw = 2,col = 'b')
            ax = dtl.plot_edges_xy([s3,s4],ax,lw = 2,col = 'g')
            ip = pym.sintsxyp(s1,s2,s3,s4)
            if type(ip) == type(()):
                ax = dtl.plot_point_xy(ip[0],ax,mk = 'o',col = 'r')
                ax = dtl.plot_point_xy(ip[1],ax,mk = 's',col = 'r')
            elif not ip is None:
                ax = dtl.plot_point_xy(ip,ax,mk = 'o',col = 'r')
            plt.show()
        def ch(ip,rp):
            if ip is None:return rp is None
            elif type(ip) == type(()):
                if not type(rp) == type(()):
                    pdb.set_trace()
                    return False
                b1 = ip[0].isnear(rp[0]) and ip[1].isnear(rp[1])
                b2 = ip[1].isnear(rp[0]) and ip[0].isnear(rp[1])
                return b1 or b2
            else:return ip.isnear(rp)
        def perm(s1,s2,s3,s4,ip,**kws):
            f = self.assertTrue
            f(ch(ip,pym.sintsxyp(s1,s2,s3,s4,**kws)))
            f(ch(ip,pym.sintsxyp(s1,s2,s4,s3,**kws)))
            f(ch(ip,pym.sintsxyp(s2,s1,s3,s4,**kws)))
            f(ch(ip,pym.sintsxyp(s2,s1,s4,s3,**kws)))
            f(ch(ip,pym.sintsxyp(s3,s4,s1,s2,**kws)))
            f(ch(ip,pym.sintsxyp(s3,s4,s2,s1,**kws)))
            f(ch(ip,pym.sintsxyp(s4,s3,s1,s2,**kws)))
            f(ch(ip,pym.sintsxyp(s4,s3,s2,s1,**kws)))
        s1,s2,s3,s4 = vec3(-1,0,0),vec3(1,1,0),vec3(-1,1,0),vec3(-1,-1,0)
        perm(s1,s2,s3,s4,vec3(-1,0,0))
        #pl(s1,s2,s3,s4)
        s1,s2,s3,s4 = vec3(-1,0.75,0),vec3(1,1,0),vec3(-1,1,0),vec3(-1,-1,0)
        perm(s1,s2,s3,s4,vec3(-1,0.75,0))
        perm(s1,s2,s3,s4,None,ie = False)
        #pl(s1,s2,s3,s4)
        s1,s2,s3,s4 = vec3(-1,0,0),vec3(1,1,0),vec3(-1,0,0),vec3(-1,-1,0)
        perm(s1,s2,s3,s4,vec3(-1,0,0))
        perm(s1,s2,s3,s4,None,ie = False)
        #pl(s1,s2,s3,s4)
        s1,s2,s3,s4 = vec3(-1,0,0),vec3(1,1,0),vec3(0,1,0),vec3(0,-1,0)
        perm(s1,s2,s3,s4,vec3(0,0.5,0))
        #pl(s1,s2,s3,s4)
        s1,s2,s3,s4 = vec3(-1,1,0),vec3(1,1,0),vec3(0,0,0),vec3(0,-1,0)
        perm(s1,s2,s3,s4,None)
        #pl(s1,s2,s3,s4)
        s1,s2,s3,s4 = vec3(-1,1,0),vec3(1,1,0),vec3(0,0,0),vec3(2,0,0)
        perm(s1,s2,s3,s4,None)
        #pl(s1,s2,s3,s4)
        s1,s2,s3,s4 = vec3(-1,1,0),vec3(1,1,0),vec3(3,1,0),vec3(2,1,0)
        perm(s1,s2,s3,s4,None)
        #pl(s1,s2,s3,s4)
        s1,s2,s3,s4 = vec3(-1.5,0,0),vec3(1,0,0),vec3(0,0,0),vec3(2.2,0,0)
        perm(s1,s2,s3,s4,(vec3(0,0,0),vec3(1,0,0)))
        #pl(s1,s2,s3,s4)
        s1,s2,s3,s4 = vec3(-1,1,0),vec3(1,1,0),vec3(1,1,0),vec3(2,1,0)
        perm(s1,s2,s3,s4,vec3(1,1,0))
        perm(s1,s2,s3,s4,None,ie = False)
        perm(s1,s2,s3,s4,vec3(1,1,0),ieb = False)
        perm(s1,s2,s3,s4,None,ie = False,ieb = False)
        #pl(s1,s2,s3,s4)
        s1,s2,s3,s4 = vec3(-1,1,0),vec3(1,1,0),vec3(-1,1,0),vec3(1,1,0)
        perm(s1,s2,s3,s4,(vec3(-1,1,0),vec3(1,1,0)))
        perm(s1,s2,s3,s4,None,ieb = False)
        #pl(s1,s2,s3,s4)
        s1,s2,s3,s4 = vec3(-1,1,0),vec3(1,1,0),vec3(0,1,0),vec3(0.5,1,0)
        perm(s1,s2,s3,s4,(vec3(0,1,0),vec3(0.5,1,0)))
        perm(s1,s2,s3,s4,(vec3(0,1,0),vec3(0.5,1,0)),ie = False)
        perm(s1,s2,s3,s4,(vec3(0,1,0),vec3(0.5,1,0)),ieb = False)
        perm(s1,s2,s3,s4,(vec3(0,1,0),vec3(0.5,1,0)),ie = False,ieb = False)
        #pl(s1,s2,s3,s4)
        s1,s2,s3,s4 = vec3(-1,0,0),vec3(0,0,0),vec3(-1,0,0),vec3(1,0,0)
        perm(s1,s2,s3,s4,(vec3(-1,0,0),vec3(0,0,0)))
        #pl(s1,s2,s3,s4)

    def test_sintbxy(self):
        def pl(s1,s2,b):
            ax = dtl.plot_axes_xy()
            ax = dtl.plot_polygon_xy(b,ax,lw = 2.0,col = 'g')
            ax = dtl.plot_edges_xy([s1,s2],ax,lw = 4.0,col = 'b')
            plt.show()
        def perm(s1,s2,b,c,**kws):
            f = self.assertTrue if c else self.assertFalse
            f(pym.sintbxy(s1,s2,b,**kws))
            f(pym.sintbxy(s2,s1,b,**kws))
        b = vec3(0,0,0).sq(8,8)
        s1,s2 = vec3(0,0,0),vec3(8,0,0)
        perm(s1,s2,b,True)
        #pl(s1,s2,b)
        s1,s2 = vec3(0,0,0),vec3(8,8,0)
        perm(s1,s2,b,True)
        #pl(s1,s2,b)
        s1,s2 = vec3(4,4,0),vec3(8,8,0)
        perm(s1,s2,b,True)
        perm(s1,s2,b,False,ie = False)
        #pl(s1,s2,b)
        s1,s2 = vec3(0,0,0),vec3(4,4,0)
        perm(s1,s2,b,True)
        perm(s1,s2,b,False,ie = False)
        #pl(s1,s2,b)
        s1,s2 = vec3(0,0,0),vec3(2,2,0)
        perm(s1,s2,b,False)
        perm(s1,s2,b,False,ie = False)
        #pl(s1,s2,b)
        s1,s2 = vec3(5,5,0),vec3(5,1,0)
        perm(s1,s2,b,False)
        #pl(s1,s2,b)
        s1,s2 = vec3(0,4,0),vec3(8,4,0)
        perm(s1,s2,b,True)
        perm(s1,s2,b,False,ie = False,col = False)
        #pl(s1,s2,b)
        s1,s2 = vec3(-4,4,0),vec3(4,4,0)
        perm(s1,s2,b,True)
        perm(s1,s2,b,False,ieb = False)
        #pl(s1,s2,b)
        s1,s2 = vec3(-5,0,0),vec3(5,0,0)
        perm(s1,s2,b,True)
        #pl(s1,s2,b)

    #def test_sintbxyp(self):

    def test_binbxy(self):
        b1 = vec3(0,0,0).pring(8,8)
        b2 = vec3(0,0,0).pring(4,8)
        b3 = vec3(8,0,0).pring(8,8)
        b4 = vec3(16,0,0).pring(8,8)
        self.assertFalse(pym.binbxy(b1,b2))
        self.assertTrue(pym.binbxy(b2,b1))
        self.assertFalse(pym.binbxy(b1,b3))
        self.assertFalse(pym.binbxy(b1,b4))

    def test_bintbxy(self):
        def pl(b1,b2):
            ax = dtl.plot_axes_xy(10)
            ax = dtl.plot_polygon_xy(b1,ax,lw = 2.0,col = 'g')
            ax = dtl.plot_polygon_xy(b2,ax,lw = 2.0,col = 'b')
            plt.show()
        b1 = vec3(0,0,0).pring(8,8)
        b2 = vec3(4,0,0).pring(4,8)
        self.assertTrue(pym.bintbxy(b1,b2))
        #pl(b1,b2)
        self.assertFalse(pym.bintbxy(b1,b2,col = False))
        b1 = vec3(0,0,0).pring(8,8)
        b2 = vec3(1,0,0).pring(4,8)
        self.assertFalse(pym.bintbxy(b1,b2))
        self.assertFalse(pym.bintbxy(b1,b2,ie = False))
        #pl(b1,b2)
        b1 = vec3(0,0,0).pring(8,8)
        b2 = vec3(5,0,0).pring(4,8)
        self.assertTrue(pym.bintbxy(b1,b2))
        self.assertTrue(pym.bintbxy(b1,b2,ie = False))
        #pl(b1,b2)
        b1 = vec3(-4,0,0).pring(2,8)
        b2 = vec3(4,0,0).pring(4,8)
        self.assertFalse(pym.bintbxy(b1,b2))
        self.assertFalse(pym.bintbxy(b1,b2,ie = False))
        #pl(b1,b2)
        b1 = vec3(-4,0,0).pring(4,8)
        b2 = vec3(4,0,0).pring(4,8)
        self.assertTrue(pym.bintbxy(b1,b2))
        self.assertFalse(pym.bintbxy(b1,b2,col = False))
        #pl(b1,b2)
        b1 = vec3(-4,0,0).pring(4,8)
        b2 = vec3(4,2,0).pring(4,8)
        self.assertTrue(pym.bintbxy(b1,b2))
        self.assertFalse(pym.bintbxy(b1,b2,col = False))
        #pl(b1,b2)
        b1 = vec3(-4,0,0).pring(4,8)
        b2 = vec3(4,4,0).pring(4,8)
        self.assertFalse(pym.bintbxy(b1,b2))
        self.assertFalse(pym.bintbxy(b1,b2,ie = False))
        #pl(b1,b2)
        b1 = vec3(-4,0,0).pring(4,8)
        b2 = vec3(4,2.0*4.0*numpy.tan(numpy.pi/8.0),0).pring(4,8)
        self.assertTrue(pym.bintbxy(b1,b2))
        self.assertFalse(pym.bintbxy(b1,b2,ie = False))
        #pl(b1,b2)

    def test_bsegbxy(self):
        def pl(b1,b2,br):
            ax = dtl.plot_axes_xy()
            ax = dtl.plot_polygon_xy(b1,ax,col = 'b',lw = 5.0)
            ax = dtl.plot_polygon_xy(b2,ax,col = 'g',lw = 5.0)
            for be in br:
                ax = dtl.plot_edges_xy(be,ax,mk = 's',col = 'r',lw = 2.0)
            #ax = dtl.plot_polygon_xy(br,ax,mk = 's',col = 'r',lw = 2.0)
            plt.show()

        boundary = vec3(0,0,0).pring(500,8)
        b1 = pym.contract(boundary,100)
        b2 = vec3(500,0,0).sq(800,300)
        br = pym.bsegbxy(b1,b2)
        #pl(b1,b2,br)
        self.assertTrue(len(br) == 10)

        b1 = vec3(-4,0,0).pring(4,8)
        b2 = vec3(2,2,0).pring(4,8)
        br = pym.bsegbxy(b1,b2)
        #pl(b1,b2,br)
        self.assertTrue(len(br) == 9)
        b1 = vec3(0,0,0).sq(3,3)
        b2 = vec3(2,2,0).sq(3,3)
        br = pym.bsegbxy(b1,b2)
        #pl(b1,b2,br)
        self.assertTrue(len(br) == 6)
        b1 = vec3(0,0,0).sq(3,3)
        b2 = vec3(3,3,0).sq(3,3)
        br = pym.bsegbxy(b1,b2)
        #pl(b1,b2,br)
        self.assertTrue(len(br) == 4)
        b1 = vec3(0,0,0).sq(3,3)
        b2 = vec3(3,0,0).sq(3,3)
        br = pym.bsegbxy(b1,b2)
        #pl(b1,b2,br)
        self.assertTrue(len(br) == 4)
        b1 = vec3(0,0,0).sq(3,3)
        b2 = vec3(3,0,0).sq(3,1)
        br = pym.bsegbxy(b1,b2)
        #pl(b1,b2,br)
        self.assertTrue(len(br) == 6)
        b1 = vec3(0,0,0).sq(3,3)
        b2 = vec3(0,-2,0).sq(3,3)
        quat(0,0,0,0).av(numpy.pi/4.0,vec3(0,0,1)).rotps(b1)
        br = pym.bsegbxy(b1,b2)
        #pl(b1,b2,br)
        self.assertTrue(len(br) == 6)
        b1 = vec3(0,0,0).sq(3,3)
        b2 = vec3(1,1,0).sq(1,1)
        br = pym.bsegbxy(b1,b2)
        #pl(b1,b2,br)
        self.assertTrue(len(br) == 6)
        b1 = vec3(0,0,0).sq(3,3)
        b2 = vec3(2,2,0).sq(3,3)
        quat(0,0,0,0).av(numpy.pi/4.0,vec3(0,0,1)).rotps(b2)
        br = pym.bsegbxy(b1,b2)
        #pl(b1,b2,br)
        self.assertTrue(len(br) == 6)

        b1 = vec3(0,0,0).sq(4,4)
        b2 = vec3(2,0,0).sq(2,6)
        br = pym.bsegbxy(b1,b2)
        #pl(b1,b2,br)
        self.assertTrue(len(br) == 6)

    #def test_sloops(self):

    def test_bsegsxy(self):
        def pl(bs):
            ax = dtl.plot_axes_xy(10)
            ax = dtl.plot_polygon_xy(b,ax,col = 'g',lw = 2)
            ax = dtl.plot_edges_xy((s1,s2),ax,col = 'r',lw = 1)
            for bpy in bs:
                bpy = pym.contract(bpy,0.1)
                ax = dtl.plot_polygon_xy(bpy,ax,col = 'b',lw = 4)
            plt.show()

        b = vec3(0,0,0).pring(7,4)
        s1,s2 = vec3(0,-10,0),vec3(0,10,0)
        bs = pym.bsegsxy(b,s1,s2)
        #pl(bs)
        b = dbl.block('C',1,3,3)[0]
        s1,s2 = vec3(0,-10,0),vec3(0,10,0)
        bs = pym.bsegsxy(b,s1,s2)
        #pl(bs)

        b = dbl.block('H',1,3,3)[0]
        s1,s2 = vec3(0,-100,0),vec3(0,100,0)
        bs = pym.bsegsxy(b,s1,s2)
        pl(bs)

    def test_ebdxy(self):
        def pl(b1,b2,br):
            ax = dtl.plot_axes_xy(10)
            ax = dtl.plot_polygon_xy(b1,ax,col = 'b',lw = 5.0)
            ax = dtl.plot_polygon_xy(b2,ax,col = 'g',lw = 5.0)
            ax = dtl.plot_polygon_xy(br,ax,col = 'r',lw = 2.0)
            plt.show()
        b1 = vec3(-4,0,0).pring(4,8)
        b2 = vec3(2,2,0).pring(4,8)
        br = pym.ebdxy(b1,b2)[0]
        #pl(b1,b2,br)
        self.assertTrue(len(br) == 9)
        b1 = vec3(0,0,0).sq(6,6)
        b2 = vec3(1.5,1.5,0).sq(3,3)
        br = pym.ebdxy(b1,b2)[0]
        #pl(b1,b2,br)
        self.assertTrue(len(br) == 6)
        b1 = vec3(0,0,0).sq(4,4)
        b2 = vec3(2,0,0).sq(2,6)
        br = pym.ebdxy(b1,b2)[0]
        #pl(b1,b2,br)
        #self.assertTrue(len(br) == 4)

        b1 = vec3(0,0,0).pring(500,8)
        b2 = vec3(500,0,0).sq(400,200)
        br = pym.ebdxy(b1,b2)[0]
        #pl(b1,b2,br)
        self.assertTrue(len(br) == 12)

    def test_bnrm(self):
        fp = [vec3(0,0,0),
            vec3(20,0,0),vec3(20,-20,0),vec3(40,-30,0),vec3(60,-30,0),
            vec3(60,40,0),vec3(-60,40,0),vec3(-60,-30,0),vec3(-40,-30,0),
            vec3(-20,-20,0),vec3(-20,0,0)]
        self.assertTrue(pym.bnrm(fp).isnear(vec3(0,0,1)))
        self.assertTrue(pym.bnrm(fp[::-1]).isnear(vec3(0,0,-1)))

        '''#
        # THIS IS AN INVALID BOUNDARY POLYGON... TEST SHOULD FAIL
        b = (
            vec3(-30.200000762939453, 0.20000000298023224, 0.0),
            vec3(-44.25, 0.20000000298023224, 0.0),
            vec3(-44.25, 0.20000000298023224, 3.0),
            vec3(-45.75, 0.20000000298023224, 3.0),
            vec3(-45.75, 0.20000000298023224, 0.0),
            vec3(-44.25, 0.20000000298023224, 0.0),
            vec3(-44.25, 0.20000000298023224, 3.0),
            vec3(-45.75, 0.20000000298023224, 3.0),
            vec3(-45.75, 0.20000000298023224, 0.0),
            vec3(-59.79999923706055, 0.20000000298023224, 0.0),
            vec3(-59.79999923706055, 0.20000000298023224, 5.0),
            vec3(-52.400001525878906, 0.20000000298023224, 5.0),
            vec3(-48.70000076293945, 0.20000000298023224, 5.0),
            vec3(-45.0, 0.20000000298023224, 5.0),
            vec3(-41.29999923706055, 0.20000000298023224, 5.0),
            vec3(-37.599998474121094, 0.20000000298023224, 5.0),
            vec3(-30.200000762939453, 0.20000000298023224, 5.0))

        ax = dtl.plot_axes(50)
        ax = dtl.plot_polygon(b,ax)
        plt.show()

        self.assertTrue(pym.bnrm(b).isnear(vec3(0,1,0)))
        '''#

    def test_contract(self):
        def pl(b):
            ax = dtl.plot_axes_xy(300)
            ax = dtl.plot_polygon_xy(b,ax)
            plt.show()

        b = pym.ebuxy(vec3(0,100,0).sq(500,200),vec3(100,0,0).sq(200,500))[0]
        b = pym.bisectb(b)
        b = pym.bisectb(b)
        b = pym.smoothxy(b,0.2)
        b = pym.contract(b,25)
        #pl(b)

    #def test_smoothxy(self):

    #def test_aggregate(self):

if __name__ == '__main__':
    unittest.main()









 



