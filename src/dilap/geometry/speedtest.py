import dilap.core.tools as dpr
import dilap.core.vector as dpv
import dilap.core.quaternion as dpq

import dilap.geometry.vec3 as v3

import numpy



def test():
    v3_1 = v3.vec3(1,1,1)
    v3_2 = v3.vec3(10,8,6)

    dpv_1 = dpv.vector(1,1,1)
    dpv_2 = dpv.vector(10,8,6)

    for x in range(100000000):
        #three = v3_1.crs(v3_2)
        three = dpv.cross(dpv_1,dpv_2)

test()




