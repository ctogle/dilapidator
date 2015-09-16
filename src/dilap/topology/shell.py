import dilap.core.base as db

import dilap.topology.face as dtf

six = 0
def index():
    global six
    six +=1 
    return six

class shell(db.base):

    def new_face(self):
        nf = dtf.face()
        nv,nl = nf.new_loop()
        self.faces.append(nf)
        self.facecount += 1
        return nv,nl,nf
        
    def __init__(self):
        self.ix = index()
        self.faces = []
        self.facecount = 0




