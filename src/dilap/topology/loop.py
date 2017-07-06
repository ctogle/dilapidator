class loop:

    def __str__(self):
        return 'loop: %d' % self.ix

    # has an index unique per mesh
    def __init__(self,edge,ix):
        self.edge = edge
        self.ix = ix
        self.edge.loop = self
