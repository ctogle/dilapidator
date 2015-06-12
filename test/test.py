import dilap.construct as dlc
import dilap.destruct as dld
import dilap.core.tools as dpr
import dilap.core.profiler as prf

import dilap.core.lsystem as pls
import dilap.core.tmesh as tms

import dilap.core.ugrid as dgrid

import dilap.core.mesh.piecewisecomplex as pwc

import dp_vector as dpv

import matplotlib.pyplot as plt

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
    pts = dpr.dice_edges(dpr.corners(5,5),1)
    pts.extend([d.copy().translate_z(10) for d in pts])

    plc = pwc.piecewise_linear_complex()
    plcxs = plc.add_points(*pts)

    for x in range(8):
        y = x+1 if x < 7 else 0
        plc.connect(x,y)

    plc.plot()
    plt.show()
    plc.tetrahedralize()

    pdb.set_trace()

    plc.plot()
    plt.show()

def triang():
    pts = dpr.dice_edges(dpr.corners(5,5),1)

    plc = pwc.piecewise_linear_complex()
    plcxs = plc.add_points(*pts)

    for x in range(8):
        y = x+1 if x < 7 else 0
        plc.connect(x,y)

    plc.plot()
    plt.show()
    plc.triangulate_xy()

    pdb.set_trace()

    plc.plot()
    plt.show()

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
triang()



