from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat
from dilap.geometry.tform import tform

import dilap.topology.tree as dtr
import dilap.topology.vert as dvt

class scenevert(dvt.vert):

    def __init__(self,ix,tform,models = None):
        dvt.vert.__init__(self,ix)
        self.tform = tform
        if models is None:models = []
        self.models = models

    def scl(self,s):
        for mod in self.models:mod.scl(s)
        return self

    def rot(self,q):
        for mod in self.models:mod.rot(q)
        return self

    def trn(self,t):
        for mod in self.models:mod.trn(t)
        return self

    def graph(self,io,wtform,**kws):
        s,r,t = wtform.scl,wtform.rot,wtform.pos
        self.scl(s).rot(r).trn(t)
        for mod in self.models:io.build_model2(mod,**kws)
        t = wtform.pos.cpf()
        r = wtform.rot.cpf()
        s = wtform.scl.cpr()
        self.trn(t).rot(r).scl(s)

class scenegraph(dtr.tree):
    
    vertclass = scenevert

    def __init__(self):
        ntf = tform(vec3(0,0,0),quat(1,0,0,0),vec3(1,1,1))
        dtr.tree.__init__(self,ntf)

    def avert(self,p,q,s,*args,parent = None,**kwargs):
        ntf = tform(p,q,s)
        vrt = dtr.tree.avert(self,parent,ntf,*args,**kwargs)
        return vrt

    def graphvert(self,vrt,ptf,io,**kws):
        wtform = vrt.tform.true(ptf)
        print('vrt!\n',vrt.ix,wtform,'\nfrom\n',vrt.tform,'\nfrom\n',ptf)
        vrt.graph(io,wtform,**kws)
        for x in self.below(vrt):self.graphvert(x,wtform,io,**kws)

    def graph(self,io,**kws):
        self.graphvert(self.root,None,io,**kws)
      









