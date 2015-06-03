import dilap.core.base as db

import dp_vector as dpv
import dp_quaternion as dpq

class tform(db.base):

    def __str__(self,b = 0):
        strr  = '\t'*b + 'tform' + '\n'
        strr += '\t'*b + 'pos: ' + str(self.pos) + '\n'
        strr += '\t'*b + 'rot: ' + str(self.rot) + '\n'
        strr += '\t'*b + 'scl: ' + str(self.scl)
        return strr

    def __init__(self,owner,**kwargs):
        self.owner = owner
        self._def('parent',None,**kwargs)
        self._def('children',[],**kwargs)
        self._def('pos',dpv.zero(),**kwargs)
        #self._def('rot',dpv.zero(),**kwargs)
        self._def('rot',dpq.zero(),**kwargs)
        self._def('scl',dpv.one(),**kwargs)

    # should return a tform to world space
    def true(self):
        np = self.pos.copy()
        nr = self.rot.copy()
        ns = self.scl.copy()
        if self.parent:
            tpar = self.parent.true()
            #np.rotate_z(tpar.rot.z)
            np.rotate(tpar.rot)
            np.translate(tpar.pos)
            #nr.translate(tpar.rot)
            nr.rotate(tpar.rot)
            ns.scale(tpar.scl)
        new = tform(self.owner,pos = np,rot = nr,scl = ns)
        return new


