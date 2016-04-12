from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat










__doc__ = '''dilapidator\'s implementation of a topological loop'''
# dilapidators implementation of a topological loop
class loop:

    def __str__(self):return 'loop:'+str(self.ix)

    # has an index unique per mesh
    def __init__(self,edge,ix):
        self.edge = edge
        self.ix = ix
        self.edge.loop = self












