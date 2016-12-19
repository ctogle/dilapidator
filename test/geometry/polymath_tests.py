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

    def atest_sintsxy(self):
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

    def atest_sintsxyp(self):
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

    def atest_sintbxy(self):
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

    def atest_binbxy(self):
        b1 = vec3(0,0,0).pring(8,8)
        b2 = vec3(0,0,0).pring(4,8)
        b3 = vec3(8,0,0).pring(8,8)
        b4 = vec3(16,0,0).pring(8,8)
        self.assertFalse(pym.binbxy(b1,b2))
        self.assertTrue(pym.binbxy(b2,b1))
        self.assertFalse(pym.binbxy(b1,b3))
        self.assertFalse(pym.binbxy(b1,b4))

    def atest_bintbxy(self):
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

    def atest_bsegbxy(self):
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

    def atest_bsegsxy(self):
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

        print('H problem!!!')
        #b = dbl.block('H',1,3,3)[0]
        #s1,s2 = vec3(0,-100,0),vec3(0,100,0)
        #bs = pym.bsegsxy(b,s1,s2)
        #pl(bs)

    def atest_ebdxy(self):
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
        self.assertTrue(len(br) == 4)
        b1 = vec3(0,0,0).pring(500,8)
        b2 = vec3(500,0,0).sq(400,200)
        br = pym.ebdxy(b1,b2)[0]
        #pl(b1,b2,br)
        self.assertTrue(len(br) == 12)

    def atest_ebuxy(self):
        def pl(b1,b2,br):
            ax = dtl.plot_axes_xy(20)
            ax = dtl.plot_polygon_xy(b1,ax,col = 'b',lw = 5.0)
            ax = dtl.plot_polygon_xy(b2,ax,col = 'g',lw = 5.0)
            for b in br:ax = dtl.plot_polygon_xy(b,ax,col = 'r',lw = 2.0)
            plt.show()

        vb = vec3(0,0,0).sq(400,400)
        dx,dy,xn,yn = 20,20,4,3
        o = vec3(-dx*xn/2.0,-dy*yn,0).com(vb)
        vgrid = [o.cp().trn(vec3(x*dx,y*dy,0)) for y in range(yn) for x in range(xn)]
        boxes = [p.sq(dx,dy) for p in vgrid]
        boxes = [b for b in boxes if pym.binbxy(b,vb)]
        box = pym.bsuxy(boxes)

        ax = dtl.plot_axes_xy(300)
        ax = dtl.plot_polygon_xy(vb,ax,lw = 2,col = 'b')
        for b in boxes:ax = dtl.plot_polygon_xy(b,ax,lw = 2,col = 'g')
        for b in box:ax = dtl.plot_polygon_xy(b,ax,lw = 2,col = 'r')
        plt.show()

        b1 = vec3(7,0,0).sq(8,8)
        b2 = vec3(0,0,0).sq(8,8)
        br = pym.ebuxy(b1,b2)
        pl(b1,b2,br)
        self.assertTrue(len(br) == 1)
        self.assertTrue(len(br[0]) == 8)
        b1 = vec3(0,0,0).sq(20,8)
        b2 = vec3(0,8,0).sq(20,8)
        br = pym.ebuxy(b1,b2)
        pl(b1,b2,br)
        self.assertTrue(len(br) == 1)
        self.assertTrue(len(br[0]) == 6)
        b1 = vec3(8,2,0).sq(8,8)
        b2 = vec3(0,0,0).sq(8,8)
        br = pym.ebuxy(b1,b2)
        pl(b1,b2,br)
        self.assertTrue(len(br) == 1)
        self.assertTrue(len(br[0]) == 10)

    def test_bnrm(self):
        def pl():
            ax = dtl.plot_axes_xy(700)
            ax = dtl.plot_polygon_xy(fp,ax,lw = 2,col = 'g')
            ax = dtl.plot_points_xy(fp,ax,number = True)
            plt.show()

        fp = [
            vec3(-6.9831085205078125, 121.9827880859375, 0.0),
            vec3(-26.983108520507812, 121.9827880859375, 0.0),
            vec3(-46.98310852050781, 121.9827880859375, 0.0),
            vec3(-66.98310852050781, 121.9827880859375, 0.0),
            vec3(-85.984375, 121.9827880859375, 0.0),
            vec3(-85.98310852050781, 102.9827880859375, 0.0),
            vec3(-85.98310852050781, 82.9827880859375, 0.0),
            vec3(-85.98310852050781, 63.9814453125, 0.0),
            vec3(-66.98310852050781, 63.9827880859375, 0.0),
            vec3(-46.98310852050781, 63.9827880859375, 0.0),
            vec3(-26.983108520507812, 63.9827880859375, 0.0),
            vec3(-6.9831085205078125, 63.9827880859375, 0.0),
            vec3(13.016891479492188, 63.9827880859375, 0.0),
            vec3(33.01689147949219, 63.9827880859375, 0.0),
            vec3(53.01689147949219, 63.9827880859375, 0.0),
            vec3(73.01689147949219, 63.9827880859375, 0.0),
            vec3(93.01689147949219, 63.9827880859375, 0.0),
            vec3(112.017578125, 63.9827880859375, 0.0),
            vec3(112.01689147949219, 82.9827880859375, 0.0),
            vec3(112.01689147949219, 102.9827880859375, 0.0),
            vec3(112.01689147949219, 121.9833984375, 0.0),
            vec3(93.01689147949219, 121.9827880859375, 0.0),
            vec3(73.01689147949219, 121.9827880859375, 0.0),
            vec3(53.01689147949219, 121.9827880859375, 0.0),
            vec3(33.01689147949219, 121.9827880859375, 0.0),
            vec3(13.016891479492188, 121.9827880859375, 0.0)]
        ax = dtl.plot_axes_xy(200)
        ax = dtl.plot_polygon_xy(fp,ax,lw = 3,col = 'b')
        ax = dtl.plot_points_xy(fp,ax,number = True)
        plt.show()
        self.assertTrue(pym.bnrm(fp).isnear(vec3(0,0,1)))

        vb = vec3(-100,100,0).sq(400,400)
        dx,dy,xn,yn = 20,20,10,30
        o = vec3(-dx*xn/2.0,-dy*yn,0).com(vb)
        vgrid = [o.cp().trn(vec3(x*dx,y*dy,0)) for y in range(yn) for x in range(xn)]
        boxes = [p.sq(dx,dy) for p in vgrid]
        boxes = [b for b in boxes if pym.binbxy(b,vb)]
        box = pym.bsuxy(boxes)
        self.assertTrue(pym.bnrm(box[0]).isnear(vec3(0,0,1)))

        fp = [vec3(0,0,0),
            vec3(20,0,0),vec3(20,-20,0),vec3(40,-30,0),vec3(60,-30,0),
            vec3(60,40,0),vec3(-60,40,0),vec3(-60,-30,0),vec3(-40,-30,0),
            vec3(-20,-20,0),vec3(-20,0,0)]
        self.assertTrue(pym.bnrm(fp).isnear(vec3(0,0,1)))
        self.assertTrue(pym.bnrm(fp[::-1]).isnear(vec3(0,0,-1)))

        fp = [
            vec3(-282.1204528808594, 127.64604187011719, 0.0),
            vec3(-334.1018981933594, 179.6274871826172, 0.0),
            vec3(-130.43038940429688, 179.62747192382812, 0.0),
            vec3(-130.43038940429688, 336.434326171875, 0.0),
            vec3(-36.91655731201172, 366.8188171386719, 0.0),
            vec3(-36.916534423828125, 461.72088623046875, 0.0),
            vec3(191.35768127441406, 461.841796875, 0.0),
            vec3(275.3216552734375, 373.85748291015625, 0.0),
            vec3(152.47915649414062, 204.779296875, 0.0),
            vec3(332.7578430175781, -43.35302734375, 0.0),
            vec3(454.96630859375, -3.6450958251953125, 0.0),
            vec3(456.61492919921875, -166.4105224609375, 0.0), 
            vec3(316.5687561035156, -120.90676879882812, 0.0),
            vec3(201.97442626953125, -278.63232421875, 0.0), 
            vec3(277.8673400878906, -383.0899963378906, 0.0),
            vec3(195.58241271972656, -472.19598388671875, 0.0),
            vec3(-10.369781494140625, -468.9277038574219, 0.0),
            vec3(-10.369804382324219, -395.02154541015625, 0.0),
            vec3(-78.82841491699219, -326.56292724609375, 0.0),
            vec3(-175.64352416992188, -326.56292724609375, 0.0),
            vec3(-244.10214233398438, -395.0215759277344, 0.0),
            vec3(-244.10214233398438, -414.89117431640625, 0.0),
            vec3(-465.6239929199219, -192.85736083984375, 0.0),
            vec3(-466.0299072265625, -122.715087890625, 0.0),
            vec3(-385.8233947753906, -122.71507263183594, 0.0),
            vec3(-282.1204528808594, -19.01211929321289, 0.0),
                ]
        bn = pym.bnrm(fp)
        print('bnnn',bn)
        fp.reverse()
        pl()
        bn = pym.bnrm(fp)
        print('bnnnrev',bn)
        self.assertTrue(pym.bnrm(fp).isnear(vec3(0,0,1)))
        self.assertTrue(pym.bnrm(fp[::-1]).isnear(vec3(0,0,-1)))

    def atest_bccw(self):
        fp = [
            vec3(-282.1204528808594, 127.64604187011719, 0.0),
            vec3(-334.1018981933594, 179.6274871826172, 0.0),
            vec3(-130.43038940429688, 179.62747192382812, 0.0),
            vec3(-130.43038940429688, 336.434326171875, 0.0),
            vec3(-36.91655731201172, 366.8188171386719, 0.0),
            vec3(-36.916534423828125, 461.72088623046875, 0.0),
            vec3(191.35768127441406, 461.841796875, 0.0),
            vec3(275.3216552734375, 373.85748291015625, 0.0),
            vec3(152.47915649414062, 204.779296875, 0.0),
            vec3(332.7578430175781, -43.35302734375, 0.0),
            vec3(454.96630859375, -3.6450958251953125, 0.0),
            vec3(456.61492919921875, -166.4105224609375, 0.0), 
            vec3(316.5687561035156, -120.90676879882812, 0.0),
            vec3(201.97442626953125, -278.63232421875, 0.0), 
            vec3(277.8673400878906, -383.0899963378906, 0.0),
            vec3(195.58241271972656, -472.19598388671875, 0.0),
            vec3(-10.369781494140625, -468.9277038574219, 0.0),
            vec3(-10.369804382324219, -395.02154541015625, 0.0),
            vec3(-78.82841491699219, -326.56292724609375, 0.0),
            vec3(-175.64352416992188, -326.56292724609375, 0.0),
            vec3(-244.10214233398438, -395.0215759277344, 0.0),
            vec3(-244.10214233398438, -414.89117431640625, 0.0),
            vec3(-465.6239929199219, -192.85736083984375, 0.0),
            vec3(-466.0299072265625, -122.715087890625, 0.0),
            vec3(-385.8233947753906, -122.71507263183594, 0.0),
            vec3(-282.1204528808594, -19.01211929321289, 0.0),
                ]
        o = pym.bccw(fp)
        self.assertTrue(o == False)
        fp.reverse()
        o = pym.bccw(fp)
        self.assertTrue(o == True)

    def atest_contract(self):
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









 



