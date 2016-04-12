import dilap.core.base as db
import dilap.core.base as db

###############################################################################
###############################################################################

class plan(db.base):

    # a plan contains data used to create 
    # and arrange models in a scenegraph

    def plot(self):
        raise NotImplemented

    def __init__(self,*args,**kwargs):
        self._def('iotype','obj',**kwargs)
        self._def('models',[],**kwargs)
        self._def('manifest',[],**kwargs)

# interpret a plan into the files necessary
# to represent it, such as obj files
def interpret(wplan,fid = None):
    raise NotImplemented

###############################################################################





