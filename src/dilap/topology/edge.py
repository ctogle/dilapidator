import dilap.core.tools as dpr

from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat










__doc__ = '''dilapidator\'s implementation of a topological edge'''
# dilapidators implementation of a topological edge
class edge:

    def __str__(self):return 'edge:'+str(self.ix)+','+str(self.gx)

    # has an index unique per mesh, and a geometric point index
    def __init__(self,one,two,ix,gx = 0):
        self.one = one
        self.two = two
        self.ix = ix
        self.gx = gx











