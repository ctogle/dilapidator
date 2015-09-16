import dilap.core.base as db

import dilap.topology.loop as dll

fix = 0
def index():
    global fix
    fix +=1 
    return fix

class face(db.base):

    def new_loop(self):
        nl = dll.loop()
        nv = nl.new_vertex()
        self.loops.append(nl)
        self.loopcount += 1
        return nv,nl

    def __init__(self):
        self.ix = index()
        self.loops = []
        self.loopcount = 0
        self.cplex = None



