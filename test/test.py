import dilap.core.vector as dpv
import dilap.core.quaternion as dpq
import dilap.construct as dlc
import dilap.destruct as dld
import dilap.core.tools as dpr
import dilap.mesh.tools as dtl
import dilap.core.profiler as prf

import dilap.core.lsystem as pls
#import dilap.core.tmesh as tms
import dilap.structures.tools as dstl

import dilap.mesh.piecewisecomplex as pwc

import matplotlib.pyplot as plt
import pdb

def cube():
    l = 2.4
    dlc.build(dlc.cube(l).translate_x(3).translate_y(0.5).translate_z(5))
    
def cone():
    l = 2.4
    dlc.build(dlc.cone(r = 2,h = 3,n = 16))

def stonehenge():
    c1 = dlc.cube().scale_z(3).translate_z(1.5)
    c2 = dlc.cylinder(r = 0.5).scale_z(3).translate_z(1.5).translate_x(4)
    c3 = dlc.cube().scale_x(5).translate_z(3.5).translate_x(2)
    c4 = dlc.cone(r = 0.5).translate_z(0.5).scale_z(0.25).translate_z(3).translate_x(4)
    stones = [c1,c2,c3,c4]
    dlc.build(stones,consume = True)

def wall():
    v1 = dpv.zero()
    v2 = dpv.zero().translate_x(5).translate_y(10)
    wa = dlc.wall(v1,v2,3,0.25)
    dlc.build(wa)

def perim():
    vs = dpr.point_ring(10,6)
    pm = dlc.perimeter(vs,3,0.25)
    dlc.build(pm)

def floor():
    gap = (dpv.zero(),2,3)
    f = dlc.floor(gap = gap)
    dlc.build(f)

def pipes():
    p1 = dlc.pipe()
    dlc.build(p1)

def roads():
    r1 = dlc.road()
    dlc.build(r1)

def houselot():
    dds = dld.dilapidors['ivy'](4)
    dlc.realize(dlc.contextualizer['lot'](dilaps = [dds]),7)

def street():
    dds = dld.dilapidors['ivy'](4)
    dlc.realize(dlc.contextualizer['street'](dilaps = [dds]),7)

def cont():
    dds = dld.dilapidors['ivy'](4)
    dlc.realize(dlc.contextualizer['continent'](dilaps = [dds]),7)

def prf_houselot():
    prf.profile_function(houselot)

def prf_street():
    prf.profile_function(street)

def prf_cont():
    prf.profile_function(cont)

def prf_lstest():
    prf.profile_function(pls.test)

def afmtest():
    dlc.build(tms.afmtest())

def tetra():
    pts = dpr.dice_edges(dpr.square(5,5),1)
    pts.extend([d.copy().translate_z(10) for d in pts])

    plc = pwc.piecewise_linear_complex()
    plcxs = plc.add_points(*pts)

    for x in range(8):
        y = x+1 if x < 7 else 0
        plc.add_edge(x,y)

    plc.plot()
    plt.show()
    plc.tetrahedralize()

    pdb.set_trace()

    plc.plot()
    plt.show()

def triang():

    #pts = dpr.dice_edges(dpr.square(50,10),2)
    hpts = [dpv.vector(5,-2,0)]
    hpts.append(hpts[-1].copy().translate(dpv.vector(0,5,0)))
    hpts.append(hpts[-1].copy().translate(dpv.vector(-10,0,0)))
    hpts.append(hpts[-1].copy().translate(dpv.vector(0,-3,0)))
    hpts.append(hpts[-1].copy().translate(dpv.vector(5,0,0)))
    hpts.append(hpts[-1].copy().translate(dpv.vector(0,-2,0)))
    hpts2 = [h.copy().translate_x(-12) for h in hpts]
    pts = dpr.inflate([h.copy().translate_x(-6) for h in hpts],14)

    #pts  = dpv.translate_coords(dpr.square(50,10),dpv.vector(-30,-12,0))
    pts2 = dpv.translate_coords(dpr.square(30,10),dpv.vector(30,20,0))
    pts3 = dpv.translate_coords(dpr.square(20,10),dpv.vector(-30,-20,0))

    pts2.insert(1,dpv.vector(20,15,0))
    pts2.insert(1,dpv.vector(20,20,0))
    pts2.insert(1,dpv.vector(18,20,0))
    pts2.insert(1,dpv.vector(18,15,0))
    #pts2 = dpr.point_ring(100,16)

    points = []
    edges = []
    polygons = [(pts,(hpts,hpts2)),(pts2,()),(pts3,())]
    #polygons = [(pts2,())]
    #polygons = [(pts,(hpts,hpts2))]
    #polygons = [(pts,()),(pts2,())]
    polyhedra = []

    plc = pwc.piecewise_linear_complex()
    plc.add_points(*points)
    plc.add_edges(*edges)
    plc.add_polygons(*polygons)
    #plc.add_polyhedra(*polyhedra)
    #plc.triangulate_xy()

    #plc.translate_polygon(2,dpv.vector(0,0,10))
    plc.extrude_polygon(2,dpv.vector(0,0,10))

    plc.triangulate()

    #ax = plc.plot_xy()
    ax = plc.plot()
    plt.show()

    #pelt = plc.covers['tri'].pelt()

    '''#
    pelt = pwc.model_plc(polygons = polygons)
    dlc.build(pelt)

    '''#

def intriangle_test():

    def show(ins):
        print('inside:',ins)
        ax = dtl.plot_axes_xy()
        ax = dtl.plot_polygon_xy([t1,t2,t3],ax)
        ax = dtl.plot_point_xy(pt,ax)
        plt.show()

    t1 = dpv.vector(0,0,0)
    t2 = dpv.vector(1,0,0)
    t3 = dpv.vector(0,1,0)

    pt = dpv.vector(0,0,0)
    ins = dpr.intriangle(pt,t1,t2,t3)
    if not ins == True:show(ins)

    pt = dpv.vector(1,0,0)
    ins = dpr.intriangle(pt,t1,t2,t3)
    if not ins == True:show(ins)

    pt = dpv.vector(0,1,0)
    ins = dpr.intriangle(pt,t1,t2,t3)
    if not ins == True:show(ins)

    pt = dpv.vector(0.5,0.5,0)
    ins = dpr.intriangle(pt,t1,t2,t3)
    if not ins == True:show(ins)

    pt = dpv.vector(0.5,0,0)
    ins = dpr.intriangle(pt,t1,t2,t3)
    if not ins == True:show(ins)

    pt = dpv.vector(0,0.5,0)
    ins = dpr.intriangle(pt,t1,t2,t3)
    if not ins == True:show(ins)

    pt = dpv.vector(0.75,0.25,0)
    ins = dpr.intriangle(pt,t1,t2,t3)
    if not ins == True:show(ins)

    pt = dpv.vector(1,1,0)
    ins = dpr.intriangle(pt,t1,t2,t3)
    if not ins == False:show(ins)

    pt = dpv.vector(-22.847272872924805, -22.847267150878906, 0.0)
    t1 = dpv.vector(-21.464466094970703, -28.53553009033203, 0.0)
    t2 = dpv.vector(-17.15900993347168, -24.230073928833008, 0.0)
    t3 = dpv.vector(-24.230077743530273, -17.159006118774414, 0.0)

    ins = dpr.intriangle(pt,t1,t2,t3)
    #ins = dpr.inside(pt,[t1,t2,t3])
    #ins = dpv.distance_to_border_xy(pt,[t1,t2,t3])
    if not ins == True:show(ins)

    print('barrry',dpr.barycentric(pt,t1,t2,t3))
    pdb.set_trace()

def sub_v_v1v2_test():
    one = dpv.vector(100,0,0)
    two = dpv.vector(0,100,100)
    
    #for x in range(1000000):
    #    three = one - two
    for x in range(1000000):
        three = dpv.v1_v2(one,two)

def mergetest():
    #left = dpr.square(0.5,8)
    #lefth = dpr.square(0.25,2)
    #right = [x.translate_x(0.5) for x in dpr.square(0.5,5)]
    #righth = [x.translate_x(0.5) for x in dpr.square(0.25,2)]
    #dpr.rotate_polygon((tuple(left),(tuple(lefth),)),dpq.q_from_av(dpr.PI/2.0,dpv.x()))
    #dpr.rotate_polygon((tuple(right),(tuple(righth),)),dpq.q_from_av(dpr.PI/2.0,dpv.x()))
    mpolys = dstl.post(dpv.vector(0,0,0),4,0.5,8.0)
    #mtest = dstl.merge_two_polygons((tuple(left),()),(tuple(right),()))
    mtest = dstl.merge_polygons(mpolys)
    #mtest = dstl.merge_polygons(
    #    ((tuple(left),(tuple(lefth),)),
    #    (tuple(right),(tuple(righth),))))[0]
    print('mtest',mtest)
    ax = dtl.plot_axes()
    for m in mtest:ax = dtl.plot_polygon_full(m,ax)
    #ax = dtl.plot_polygon_xy(left,ax)
    #ax = dtl.plot_polygon_xy(right,ax)
    plt.show()

def isecttest():
    p1,p2 = dpv.vector(-10,0,0),dpv.vector(10,0,0)
    p3,p4 = dpv.vector(-5,0,0),dpv.vector(5,0,0)
    isect1 = dtl.segments_intersect_at(p1,p2,p3,p4)
    isect2 = dtl.segments_intersect_at(p1,p2,p4,p3)
    isect3 = dtl.segments_intersect_at(p2,p1,p4,p3)
    isect4 = dtl.segments_intersect_at(p2,p1,p3,p4)
    print('isect',isect1,isect2,isect3,isect4)
    pdb.set_trace()
    
def csgtest():
    plc1 = dtl.box(5,5,5)
    plc2 = dtl.box(5,3,4).translate(dpv.vector(4,0,-0.5))

    print('union input')
    ax = dtl.plot_axes()
    ax = plc1.plot(ax)
    ax = plc2.plot(ax)
    plt.show()

    #plc3 = pwc.union(plc1,plc2)
    #plc3 = pwc.intersection(plc1,plc2)
    plc3 = pwc.difference(plc1,plc2)

    print('union output')
    ax = dtl.plot_axes()
    ax = plc3.plot(ax)
    plt.show()

def containmenttest():
    p = dpv.vector(2.5,0,0)
    py = tuple(dpr.square(5,5))

    wn = dpr.winding(p,py)
    ins = dpr.inconcave_xy(p,py)
    print('wtfman',wn,ins)

def break_polygontest():
    pts1 = dpv.translate_coords(dpr.square(30,10),dpv.vector( 0,0,0))
    pts2 = dpv.translate_coords(dpr.square(20,10),dpv.vector( 0,0,0))
    dpv.rotate_coords(pts2,dpq.q_from_av(dpr.PI2,dpv.y()))

    broken = pwc.break_polygon((tuple(pts1),()),[(tuple(pts2),())])
    print('brokenarrow',broken)



#prf.profile_function(sub_v_v1v2_test)

#pdb.set_trace()

#cube()
#cone()
#stonehenge()
#wall()
#perim()
#floor()
#pipes()
#roads()
#prf_houselot()
#prf_cont()
#prf_lstest()
#pls.test()
#afmtest()
#tetra()
#triang()
#cont()
csgtest()
#break_polygontest()
#containmenttest()
#mergetest()
#isecttest()
#intriangle_test()
#facade()



