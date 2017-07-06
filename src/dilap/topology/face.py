class face:

    def __str__(self):
        return 'face: %d' % self.ix

    # has an index unique per mesh
    def __init__(self,ix):
        self.ix = ix
