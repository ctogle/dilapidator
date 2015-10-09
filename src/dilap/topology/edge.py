import dilap.core.tools as dpr

from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat










__doc__ = '''dilapidator\'s implementation of a topological edge'''
# dilapidators implementation of a topological edge
class edge:

    def __str__(self):return 'edge:'+str(self.ix)

    #def connect(self,l):
    #    self.ring.append(l)

    #def disconnect(self,l):
    #    self.ring.remove(l)

    #def snxt(self,nxt):
    #    self.nxt = nxt

    #def gnxt(self,nxt):
    #    return self.nxt

    # has an index unique per mesh
    def __init__(self,one,two,ix):
        self.ix = ix   # unique index

        self.one = one # start vertex
        self.two = two # end vertex

        self.rl = None # right loop
        self.ll = None # left loop
        self.rcw = None  # right clockwise edge
        self.rccw = None # right counter-clockwise edge 
        self.lcw = None  # left clockwise edge 
        self.lccw = None # left counter-clockwise edge 



        #self.ring = [] # loop-edge links

        #self.nxt = None
        #self.lst = None
        #self.loop = None









