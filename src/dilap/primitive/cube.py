import dilap.core.model as dmo
import dilap.primitive.tools as dpr

class cube(dmo.model):

    def __init__(self,*args,**kwargs):
        ps,ns,us,fs = self._geo()
        dmo.model.__init__(self,*args,**kwargs)

    def _geo(self):
        bottom = dpr.point_ring(0.5,4)

        pdb.set_trace()

        print 'cube geometry!'


