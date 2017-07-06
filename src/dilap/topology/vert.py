from dilap.geometry import *


class vert:


    def __str__(self):
        return 'vert: %d' % self.ix


    def __is_equal(self,o):
        return self.ix == o.ix


    def __richcmp__(x,y,op):
        if op == 2:return x.__is_equal(y)
        else:assert False


    def connect(self,e):
        self.ring.append(e)


    def disconnect(self,e):
        self.ring.remove(e)


    def __init__(self,ix,*args,**kwargs):
        self.ix = ix
        self.ring = []
        self.args = args
        for k in kwargs:
            if not hasattr(self,k):
                self.__setattr__(k,kwargs[k])
