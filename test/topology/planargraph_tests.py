from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat
import dilap.geometry.tools as gtl
import dilap.geometry.polymath as pym

import dilap.topology.planargraph as pgr

import dilap.core.plotting as dtl
import matplotlib.pyplot as plt

import unittest,random,numpy

import pdb

#python3 -m unittest discover -v ./ "*tests.py"

###############################################################################

class test_planargraph(unittest.TestCase):

    def test_avrvaerese(self):

        rg = pgr.planargraph()

        i1 = rg.av(p = vec3( 10,-5,0),l = 0)
        i2 = rg.av(p = vec3( 10, 5,0),l = 0)
        i3 = rg.av(p = vec3(-10, 5,0),l = 0)
        i4 = rg.av(p = vec3(-10,-5,0),l = 0)

        r1 = rg.ae(i1,i2)
        r2 = rg.ae(i2,i3)
        r3 = rg.ae(i3,i4)
        r4 = rg.ae(i4,i1)

        r2 = rg.re(i2,i3)

        i5 = rg.av(p = vec3(2,-10,0),l = 0)
        r5,r6 = rg.se(i1,i4,i5)

        #i1 = rg.rv(i1)

        if False:
            ax = rg.plot()
            plt.show()

    def test_orings(self):
        rg = pgr.planargraph()

        i1 = rg.av(p = vec3(  0, 0,0),l = 0)
        i2 = rg.av(p = vec3( 10, 0,0),l = 0)
        i3 = rg.av(p = vec3(-10, 0,0),l = 0)
        i4 = rg.av(p = vec3(  0,10,0),l = 0)

        r1 = rg.ae(i1,i2)
        r2 = rg.ae(i1,i3)
        r3 = rg.ae(i1,i4)

        self.assertEqual(rg.orings[0],[1,3,2])
        if False:
            ax = rg.plot()
            plt.show()

    def test_peninsula(self):
        def pl():
            ax = rg.plot()
            ax = dtl.plot_polygon(fp,ax,col = 'b')
            ax = dtl.plot_point(ex,ax,col = 'r')
            plt.show()

        # given a footprint, a point within the footprint, and a point
        # possibly within the footprint, return a point within the footprint
        # on the line segment between the two points
        def wfp(fp,p1,p2,r = 5):
            ips = pym.sintbxyp(p1,p2,fp)
            if len(ips) > 0:
                for ip in ips:
                    if ip.isnear(p1):continue
                    p2 = p1.lerp(ip,1-r/p1.d(ip))
            return p2

        # given an intersection, make/find another intersection and add a road
        def mr(rg,fp,ix):
            i = rg.vs[ix]
            ip = i[1]['p']
            es = rg.rings[i[0]]
            oes = rg.orings[i[0]]
            oec = len(oes)

            if oec == 0:
                di = ip.tov(vec3(0,0,0).com(fp))
                p2 = ip.cp().trn(di)
                #p2 = wfp(fp,ip,p2)
                res = rg.mev(ix,{'p':p2,'l':0},{})
                
            elif oec == 1:
                op = rg.vs[oes[0]][1]['p']
                a = random.choice((-1,1))*numpy.pi/2.0
                q = quat(1,0,0,0).av(a,vec3(0,0,1))
                p2 = ip.cp().trn(op.tov(ip).rot(q))
                p2 = wfp(fp,ip,p2)
                res = rg.mev(ix,{'p':p2,'l':0},{})

            else:
                for x in range(1,len(oes)):
                    oeas = rg.ea(i[0],oes[x-1],oes[x])

                    print('oeoeoe',oeas)

                    #rg.vs[oes[0]][1]['p']

                #i2,r1 = rg.mev(i1,{'p':p2,'l':0},{})
                pdb.set_trace()

            return res



        # WOULD BE NICE:
        # given a graph, footprint, and position, and radius, 
        #   return an intersection within radius of position if it exists

        fp = vec3(0,0,0).pring(50,8)
        ex = fp[-2].mid(fp[-1])
        rg = pgr.planargraph()

        i1 = rg.av(p = ex,l = 0)
        i2,r1 = mr(rg,fp,i1)
        i3,r2 = mr(rg,fp,i2)
        i4,r3 = mr(rg,fp,i3)
        i5,r4 = mr(rg,fp,i4)
        i6,r5 = mr(rg,fp,i5)

        #i2,r1 = rg.mev(i1,{'p':p2,'l':0},{})

        if False:pl()

    def test_mev(self):
        def pl():
            ax = rg.plot()
            ax = dtl.plot_polygon(fp,ax,col = 'b')
            ax = dtl.plot_point(ex,ax,col = 'r')
            plt.show()

        p1 = vec3(100,0,0)
        p2 = vec3(0,100,0)
        p3 = vec3(100,100,0)

        rg = pgr.planargraph()

        i1 = rg.av(p = p1,l = 0)
        i2,r1 = rg.mev(i1,{'p':p2,'l':0},{})
        i3,r2 = rg.mev(i2,{'p':p3,'l':0},{})

        if False:pl()

    #def test_cw(self):

    def test_ccw(self):
        def pl(rg,u,v,t):
            ax = rg.plotxy(s = 0.02,number = True)
            plt.show()
            ax = rg.plotxy(s = 0.02,number = False)
            #ax = dtl.plot_axes_xy(10)
            up = rg.vs[u][1]['p']
            vp = rg.vs[v][1]['p']
            tp = rg.vs[t][1]['p']
            ax = dtl.plot_point_xy_annotate(up,ax,'u   |')
            ax = dtl.plot_point_xy(up,ax,col = 'r')
            ax = dtl.plot_point_xy_annotate(vp,ax,'v   |')
            ax = dtl.plot_point_xy(vp,ax,col = 'r')
            ax = dtl.plot_point_xy_annotate(tp,ax,'t   |')
            ax = dtl.plot_point_xy(tp,ax,col = 'r')
            uor = rg.orings[u]
            for uox in range(len(uor)):
                if uox == v:continue
                up = rg.vs[uor[uox]][1]['p']
                ax = dtl.plot_point_xy_annotate(up,ax,str(uor[uox]))
                ax = dtl.plot_point_xy(up,ax,col = 'b')
            vor = rg.orings[v]
            for vox in range(len(vor)):
                if vox == u:continue
                vp = rg.vs[vor[vox]][1]['p']
                ax = dtl.plot_point_xy_annotate(vp,ax,str(vor[vox]))
                ax = dtl.plot_point_xy(vp,ax,col = 'b')
            plt.show()
        def perm(rg,u,v,w):
            tip = rg.ccw(u,v)
            #pl(rg,u,v,tip)
            self.assertEqual(tip,w)
        b1 = vec3(-4,0,0).pring(4,8)
        b2 = vec3(2,2,0).pring(4,8)
        b1segs = pym.bsegbxy(b1,b2)
        b2segs = pym.bsegbxy(b2,b1)
        rg = pym.sstopg(b1segs+b2segs)
        #perm(rg,1,3,4)
        perm(rg,0,1,14)
        perm(rg,8,0,1)
        perm(rg,14,3,4)

    def test_loop(self):
        def pl(il):
            ilp = [rg.vs[j][1]['p'] for j in il]
            ilp = pym.contract(ilp,2.0)
            ax = rg.plot()
            ax = dtl.plot_polygon(ilp,ax,col = 'b',lw = 4)
            plt.show()
        rg = pgr.planargraph()
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

        b1 = vec3(-4,0,0).pring(4,8)
        b2 = vec3(2,2,0).pring(4,8)
        b1segs = pym.bsegbxy(b1,b2)
        b2segs = pym.bsegbxy(b2,b1)
        rg = pym.sstopg(b1segs+b2segs)
        il = rg.loop(3,4,'ccw')
        #pl(il)

        il = rg.loop(1,2,'ccw')
        #pl(il)

        #il = rg.loop(1,3,'ccw')
        #pl(il)

        il = rg.loop(2,3,'ccw')
        #pl(il)

        il = rg.loop(4,5,'ccw')
        #pl(il)

        il = rg.loop(5,6,'ccw')
        #pl(il)

        il = rg.loop(6,7,'ccw')
        #pl(il)

        il = rg.loop(7,8,'ccw')
        #pl(il)

        #il = rg.loop(8,9,'ccw')
        #pl(il)

    def test_uloops(self):
        def pl():
            ax = rg.plotxy(l = 500)
            for lp in loops:
                lpps = [rg.vs[j][1]['p'] for j in lp]
                lpps = pym.contract(lpps,0.1)
                ax = dtl.plot_polygon_xy(lpps,ax,lw = 3,col = 'b')
            plt.show()
        rg = pgr.planargraph()
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

        b1 = vec3(-4,0,0).pring(4,8)
        b2 = vec3(2,2,0).pring(4,8)
        b1segs = pym.bsegbxy(b1,b2)
        b2segs = pym.bsegbxy(b2,b1)
        rg = pym.sstopg(b1segs+b2segs)
        loops = rg.uloops('ccw')
        #pl()

        b = [
            vec3(-320.2890319824219, 72.00794982910156, 40.000003814697266),
            vec3(-299.5256042480469, 72.07875061035156, 40.000003814697266),
            vec3(-299.5280456542969, 91.03677368164062, 40.000003814697266),
            vec3(-299.5280456542969, 129.9380645751953, 40.000003814697266),
            vec3(-299.5280456542969, 135.1514129638672, 40.000003814697266),
            vec3(-299.5280456542969, 138.2604522705078, 40.000003814697266),
            vec3(-292.4235534667969, 147.1795654296875, 40.000003814697266),
            vec3(-305.3035583496094, 125.37327575683594, 40.000003814697266),
            vec3(-312.3858947753906, 99.95317077636719, 40.000003814697266),
                ]
        bsegs = [(b[j-1],b[j]) for j in range(len(b))]
        rg = pym.sstopg(bsegs)
        loops = rg.uloops('ccw')
        print('itsher')
        pl()

    def atest_polygon(self):
        def pl():
            ax = dtl.plot_axes_xy(50)
            ax = dtl.plot_polygon_full_xy(py,ax,lw = 2,col = 'b')
            plt.show()

        rg = pgr.planargraph()

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
        i7,r10 = rg.mev(i1,{'p':vec3(12,-20,0),'l':0},{})
        i8,r11 = rg.mev(i3,{'p':vec3(-5,0,0),'l':0},{})

        py = rg.polygon(2,'ccw')
        #pl()

###############################################################################

if __name__ == '__main__':unittest.main()

###############################################################################



 

