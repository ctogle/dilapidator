import dilap.core.tools as dpr

from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat










__doc__ = '''dilapidator\'s implementation of a topological vertex'''
# dilapidators implementation of a topological vertex
class vert:

    def __str__(self):return 'vert:'+str(self.ix)+','+str(self.gx)
    #def __iter__(self):yield self.x;yield self.y;yield self.z
    #def __add__(self,o):return vec3(self.x+o.x,self.y+o.y,self.z+o.z)
    #def __sub__(self,o):return vec3(self.x-o.x,self.y-o.y,self.z-o.z)
    #def __mul__(self,o):return self.cp_c().mul_c(o)
    def __is_equal(self,o):return self.gx == o.gx
    # could use <,> for lexicographic ordering?
    def __richcmp__(x,y,op):
        if op == 2:return x.__is_equal(y)
        else:assert False

    # has an index unique per mesh, and a geometric point index
    def __init__(self,ix,gx = 0):
        self.ix = ix
        self.gx = gx










