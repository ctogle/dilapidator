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



    print('input plc')
    plc.plot()
    plt.show()
    print('begin tetra')
    plc.tetrahedralize()
    print('end tetra')
    print('output plc')
    plc.plot()
    plt.show()

def triang():

    #pts = dpr.dice_edges(dpr.square(50,10),2)
    hpts = [dpv.vector(5,-2,3)]
    hpts.append(hpts[-1].copy().translate(dpv.vector(0,5,0)))
    hpts.append(hpts[-1].copy().translate(dpv.vector(-10,0,0)))
    hpts.append(hpts[-1].copy().translate(dpv.vector(0,-3,0)))
    hpts.append(hpts[-1].copy().translate(dpv.vector(5,0,0)))
    hpts.append(hpts[-1].copy().translate(dpv.vector(0,-2,0)))
    hpts2 = [h.copy().translate_x(-12).translate_z(-6) for h in hpts]
    pts = dpr.inflate([h.copy().translate_x(-6).translate_z(-3) for h in hpts],14)

    #pts  = dpv.translate_coords(dpr.square(50,10),dpv.vector(-30,-12,0))
    pts2 = dpv.translate_coords(dpr.square(30,10),dpv.vector(30,20,0))
    pts3 = dpv.translate_coords(dpr.square(20,10),dpv.vector(-30,-20,0))

    pts2.insert(1,dpv.vector(25,15,0))
    pts2.insert(1,dpv.vector(25,20,0))
    pts2.insert(1,dpv.vector(20,20,0))
    pts2.insert(1,dpv.vector(20,15,0))
    #pts2 = dpr.point_ring(100,16)

    ax = dtl.plot_axes_xy()
    ax = dtl.plot_points_xy(pts2,ax,number = True)
    plt.show()

    points = []
    edges = []
    polygons = [(pts,(hpts,hpts2)),(pts2,()),(pts3,())]
    #polygons = [(pts2,())]
    #polygons = [(pts,(hpts,hpts2))]
    #polygons = [(pts,()),(pts2,())]
    polyhedra = []

    plc = pwc.piecewise_linear_complex(
        refine = True,smooth = False)
    plc.add_points(*points)
    plc.add_edges(*edges)
    plc.add_polygons(*polygons)
    #plc.add_polyhedra(*polyhedra)
    plc.triangulate()

    #ax = plc.plot_xy()
    ax = plc.plot()
    plt.show()

    '''#
    pelt = pwc.model_plc(polygons = polygons)
    dlc.build(pelt)

    '''#

def csgtest():
    plc1 = dtl.box(5,5,5)
    #plc2 = dtl.icosphere(2,1)
    plc2 = dtl.box(5,3,4).translate(dpv.vector(2,2,-0.5))
    #plc2 = dtl.box(5,2,2).translate(dpv.vector(0,0,1.5))
    #plc2 = dtl.box(5,2,2).translate(dpv.vector(0,0,1.5))
    #plc2.translate(dpv.vector(0,0,5.0))
    #plc2.translate(dpv.vector(0,0,-3.0))
    #plc2.rotate(dpq.q_from_av(dpr.PI/6.0,dpv.x()))
    #plc2.translate(dpv.vector(0,0, 3.0))
    #plc1.triangulate()
    #plc2.triangulate()

    print('union input')
    ax = dtl.plot_axes()
    ax = plc1.plot(ax)
    ax = plc2.plot(ax)
    plt.show()

    plc3 = pwc.union(plc1,plc2)
    #plc3 = pwc.intersection(plc1,plc2)
    #plc3 = pwc.difference(plc1,plc2)

    print('union output')
    ax = dtl.plot_axes()
    ax = plc3.plot(ax)
    plt.show()

cube()
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

#tetra()
#triang()
#cont()
#csgtest()





