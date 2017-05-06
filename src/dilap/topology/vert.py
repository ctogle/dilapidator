from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat



__doc__ = '''dilapidator\'s implementation of a topological vertex'''
# dilapidators implementation of a topological vertex
class vert:

    def __str__(self):return 'vert:'+str(self.ix)
    def __is_equal(self,o):return self.ix == o.ix
    def __richcmp__(x,y,op):
        if op == 2:return x.__is_equal(y)
        else:assert False

    def connect(self,e):
        self.ring.append(e)

    def disconnect(self,e):
        self.ring.remove(e)

    # has an index unique per mesh, and list of 1-ring verts
    def __init__(self,ix,*args,**kwargs):
        self.ix = ix
        self.ring = []
        self.args = args
        self.kwargs = kwargs










