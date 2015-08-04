

class pointset:

    def __iter__(self):
        return self.ps.__iter__()

    def __init__(self):
        self.ps = []
        self.pcnt = 0
        self.snapgrid = 0.001

    def lexicographic(self):raise NotImplemented

    #####
    def get_points_copy(self,rng = None):
        if rng is None:rng = range(self.pcnt)
        return [self.ps[x].copy() for x in rng]

    #def get_points(self,rng = None):
    def get_points(self,*rng):
        if not rng:rng = range(self.pcnt)
        #if rng is None:rng = range(self.pcnt)
        return [self.ps[x] for x in rng]

    def add_point(self,np):
        px = self.pcnt
        self.ps.append(np)
        self.pcnt += 1
        return px

    def add_points(self,*nps):
        pxs = self.pcnt
        for np in nps:
            self.ps.append(np)
            self.pcnt += 1
        return range(pxs,self.pcnt)


