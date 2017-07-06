from dilap.geometry import *


class edge:


    def __str__(self):
        return 'edge: %d' % self.ix


    # has an index unique per mesh
    def __init__(self,one,two,ix):
        self.ix = ix   # unique index
        self.one = one # start vertex
        self.two = two # end vertex
