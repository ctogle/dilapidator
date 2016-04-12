from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat










__doc__ = '''dilapidator\'s implementation of a topological face'''
# dilapidators implementation of a topological face
class face:

    def __str__(self):return 'face:'+str(self.ix)

    # has an index unique per mesh
    def __init__(self,ix):
        self.ix = ix













