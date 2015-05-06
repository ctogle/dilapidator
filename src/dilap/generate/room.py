import dilap.core.base as db
import dilap.core.tools as dpr
import dilap.core.sgraph as dsg
import dilap.core.context as dgc
import dilap.io.io as dio
import dilap.primitive.cube as dcu
import dilap.primitive.cone as dco
import dilap.primitive.cylinder as dcyl
import dilap.primitive.wall as dw
import dilap.primitive.floor as df     

import dp_vector as dpv

# replacing the functionality of floor_sector from make_places
class room(dgc.context):

    def __init__(self,*args,**kwargs):
        dgc.context.__init__(self,*args,**kwargs)
        self._def('x',0,**kwargs)
        self._def('y',0,**kwargs)
        self._def('l',10,**kwargs)
        self._def('w',10,**kwargs)
        self._def('fh',0.25,**kwargs)
        self._def('ch',0.25,**kwargs)
        self._def('fgap',None,**kwargs)
        self._def('cgap',None,**kwargs)
        self._def('shafted',False,**kwargs)

    def generate(self,wallheight = 4.0,worn = 0):
        fgap,cgap = self.fgap,self.cgap
        fl = df.floor(self.l,self.w,h = self.fh,gap = fgap)
        fl.translate_x(self.x).translate_y(self.y).translate_z(self.fh)
        cl = df.floor(self.l,self.w,h = self.ch,gap = cgap)
        cl.translate_x(self.x).translate_y(self.y).translate_z(wallheight)
        #rnode = self._node_consume(self._node_wrap(fl,cl))
        rnode = self._node_wrap(fl,cl)
        self._nodes_to_graph(rnode)


