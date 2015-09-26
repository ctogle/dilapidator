import dilap.core.tools as dpr

from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat










__doc__ = '''dilapidator\'s implementation of a topological edge'''
# dilapidators implementation of a topological edge
class edge:

    def __str__(self):return 'edge:'+str(self.ix)

    def connect(self,f):
        self.ring.append(f)

    def disconnect(self,f):
        self.ring.remove(f)

    def snxt(self,nxt):
        self.nxt = nxt

    def gnxt(self,nxt):
        return self.nxt

    # has an index unique per mesh
    def __init__(self,one,two,ix):
        self.one = one
        self.two = two
        self.ix = ix
        self.ring = []
        self.nxt = None
        self.lst = None
        self.loop = None











