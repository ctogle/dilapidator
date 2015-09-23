from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat
from dilap.geometry.tform import tform

import dilap.topology.tree as dtr
import dilap.topology.vert as dvt










class scenevert(dvt.vert):

    def __init__(self,ix,tform):
        dvt.vert.__init__(self,ix)
        self.tform = tform
        self.models = []

class scenegraph(dtr.tree):
    
    vertclass = scenevert

    def __init__(self):
        ntf = tform(vec3(0,0,0),quat(0,0,0,1),vec3(1,1,1))
        dtr.tree.__init__(self,ntf)

    def avert(self,p,q,s,parent = None):
        ntf = tform(p,q,s)
        vrt = dtr.tree.avert(self,parent,ntf)
        return vrt

    def graphvert(self,vrt,ptf):
        wtform = vrt.tform.true(ptf)
        print('vrt!\n',vrt.ix,wtform,'\nfrom\n',vrt.tform,'\nfrom\n',ptf)
        for x in self.below(vrt):
            self.graphvert(x,wtform)

    def graph(self):
        self.graphvert(self.root,None)
      









