import dilap.construct as dlc
import dilap.primitive.tools as pdr

import dp_vector as dpv

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
    vs = pdr.point_ring(10,6)
    pm = dlc.perimeter(vs,3,0.25)
    dlc.build(pm)

def floor():
    gap = (dpv.zero(),2,3)
    f = dlc.floor(gap = gap)
    dlc.build(f)

def houselot():
    dcx = dlc.lot()
    dcx.generate()
    dcx.graph()

#cube()
#cone()
#stonehenge()
#wall()
#perim()
#floor()
houselot()


