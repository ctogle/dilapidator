import dilap.core.tools as dpr
import dilap.core.vector as dpv
import dilap.structures.structure as dst

import dilap.mesh.tools as dtl
import matplotlib.pyplot as plt

class house(dst.structure):

    def __init__(self,*args,**kwargs):
        dst.structure.__init__(self,*args,**kwargs)
        self._def('l',20,**kwargs)
        self._def('w',20,**kwargs)
        self._def('wallwidth',0.75,**kwargs)
        self.gen(*args,**kwargs)

    def gen(self,*args,**kwargs):
        l,w,rh = self.l,self.w,10
        base1 = dpr.square(l,w)
        base2 = dpr.square(l,w)
        dpv.translate_coords_x(base2,l)
        room1 = ((base1,rh,(-1,-2,1,-2)),{'style':'empty'})
        room2 = ((base2,rh,(0,-2,-1,-2)),{'style':'empty'})
        rooms = (room1,room2)
        #self.graph.add_room(base1,(0,),**{'height':rh,'style':'empty'})

        r = 20
        self.graph._add_room((
            ((0,0,0),(r,0,0),(r,r,0),(0,r,0)),
            ((0,0,1),(r,0,1),(r,r,1),(0,r,1)),
            #((r,0,0),(r,r,0),(0,r,0)),
            #((r,0,1),(r,r,1),(0,r,1)),
                ))
        self.graph._add_room((((r,0,0),(2*r,0,0),(2*r,r,0),(r,r,0)),))
        self.graph._add_room((((r,r,0),(2*r,r,0),(2*r,2*r,0),(r,2*r,0)),))

        self.graph._connect_rooms(0,1)
        self.graph._connect_rooms(1,2)
        self.graph._connect_rooms(0,-1)

        #self.graph.plot()
        #plt.show()

        #self.plan(rooms)






