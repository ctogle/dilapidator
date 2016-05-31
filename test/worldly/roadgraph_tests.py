from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat
import dilap.geometry.tools as gtl
import dilap.geometry.polymath as pym

import dilap.worldly.roadgraph as rdg

import dilap.core.plotting as dtl
import matplotlib.pyplot as plt

import unittest,random,numpy

import pdb

#python3 -m unittest discover -v ./ "*tests.py"

###############################################################################

class test_roadgraph(unittest.TestCase):

    def test_avrvaerese(self):

        rg = rdg.wgraph()

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

        i1 = rg.rv(i1)

        if False:
            ax = rg.plot()
            plt.show()

    def test_orings(self):
        rg = rdg.wgraph()

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

        def veri(rg,fp,p,r = 25):
            fpnms = pym.bnrmsxy(fp)

            fpds = pym.bdistpxy(fp)
            fpix = fpds.index(min(fpds))
            print('fpds',fpds,fpix)

            ax = dtl.plot_axes_xy(100)
            ax = dtl.plot_polygon_xy(fp,ax)
            ax = dtl.plot_point_xy(p,ax,col = 'r')
            ax = dtl.plot_edges_xy((fp[fpix-1],fp[fpix]),ax,lw = 3,col = 'g')
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
        rg = rdg.wgraph()

        i1 = rg.av(p = ex,l = 0)
        i2,r1 = mr(rg,fp,i1)
        i3,r2 = mr(rg,fp,i2)
        i4,r3 = mr(rg,fp,i3)
        i5,r4 = mr(rg,fp,i4)
        i6,r5 = mr(rg,fp,i5)

        #veri(rg,fp,p2)
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

        rg = rdg.wgraph()

        i1 = rg.av(p = p1,l = 0)
        i2,r1 = rg.mev(i1,{'p':p2,'l':0},{})
        i3,r2 = rg.mev(i2,{'p':p3,'l':0},{})

        if False:pl()

    def test_loop(self):
        def pl(il):
            ilp = [rg.vs[j][1]['p'] for j in il]
            ilp = pym.contract(ilp,2.0)
            ax = rg.plot()
            ax = dtl.plot_polygon(ilp,ax,col = 'b',lw = 4)
            plt.show()
        rg = rdg.wgraph()
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
        pl(il)
        il = rg.loop(i3,i4,'cw')
        self.assertEqual(il,[i3,i4,i5,i1,i7,i1,i2,i6])
        pl(il)

    def atest_uloops(self):
        def pl():
            ax = rg.plot()
            for lp in loops:
                lpps = [rg.vs[j][1]['p'] for j in lp]
                lpps = pym.contract(lpps,2)
                ax = dtl.plot_polygon(lpps,ax,lw = 3,col = 'b')
            plt.show()
        rg = rdg.wgraph()
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

        if True:pl()
        self.assertEqual(len(loops),3)

        i7,r10 = rg.mev(i1,{'p':vec3(12,-20,0),'l':0},{})
        loops = rg.uloops('ccw')

        if True:pl()
        self.assertEqual(len(loops),3)

###############################################################################

if __name__ == '__main__':unittest.main()

###############################################################################







