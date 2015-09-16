import dilap.core.base as db

import dilap.topology.vertex as dvt

lix = 0
def index():
    global lix
    lix +=1 
    return lix

class loop(db.base):

    def new_vertex(self):
        nv = dvt.vertex()
        self.vertex = nv
        return nv

    def __init__(self):
        self.ix = index()
        self.hedge = None

        self.vertex = None


