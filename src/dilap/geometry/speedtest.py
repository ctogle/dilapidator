import dilap.core.tools as dpr
import dilap.core.vector as dpv
import dilap.core.quaternion as dpq

import dilap.geometry.vec3 as v3

import numpy



def test():
    one = v3.vec3(1,1,0)
    two = dpv.vector(1,1,0)
    a = dpr.PI2
    q = dpq.q_from_av(a,dpv.vector(1,0,0))

    for x in range(100000000):
        #one = (1,1,0)
        #one = v3.vec3(1,1,0)
        #one = numpy.zeros((3),dtype = numpy.float)
        #q = dpq.q_from_av(a,dpv.vector(1,0,0))
        #one.xrot(a)
        #one.qrot(q)
        two.rotate(q)

def test2():
    x,y = 10,13

    for z in range(100000000):
        #w = x-y
        #w2 = w*w

        w = abs(x-y)

test2()




