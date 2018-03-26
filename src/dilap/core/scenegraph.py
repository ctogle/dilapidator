from dilap.geometry import *
from dilap.topology import *
from dilap.io import exporters as io
import os


class scenevert(vert):


    def __init__(self,ix,tform,models = None):
        vert.__init__(self,ix)
        self.tform = tform
        if models is None:
            self.models = []
        else:
            self.models = models


    def scl(self,s):
        for mod in self.models:
            mod.scl(s)
        return self


    def rot(self,q):
        for mod in self.models:
            mod.rot(q)
        return self


    def trn(self,t):
        for mod in self.models:
            mod.trn(t)
        return self


    def worldspace(self,wtform):
        s,r,t = wtform.scl,wtform.rot,wtform.pos
        self.scl(s).rot(r).trn(t)
        yield from self.models
        t = wtform.pos.cpf()
        r = wtform.rot.cpf()
        s = wtform.scl.cpr()
        self.trn(t).rot(r).scl(s)


class scenegraph(tree):
    

    vertclass = scenevert


    def __init__(self):
        ntf = tform(vec3(0,0,0),quat(1,0,0,0),vec3(1,1,1))
        tree.__init__(self,ntf)


    def avert(self,p,q,s,*args,parent = None,**kwargs):
        if p is None:p = vec3(0,0,0)
        if q is None:q = quat(1,0,0,0)
        if s is None:s = vec3(1,1,1)
        ntf = tform(p,q,s)
        vrt = tree.avert(self,parent,ntf,*args,**kwargs)
        return vrt


    def worldspace(self,vrt = None,ptf = None):
        '''Walk the graph yielding models in world space
        To walk the entire graph: sg.worldspace(sg.root,None)'''
        if vrt is None:vrt = self.root
        wtform = vrt.tform.true(ptf)
        yield from vrt.worldspace(wtform)
        for x in self.below(vrt):
            yield from self.worldspace(x,wtform)


    def output(self, exporter, output):
        if not os.path.isdir(output):
            os.mkdir(output)
        io[exporter].build_scenegraph(self, output)
