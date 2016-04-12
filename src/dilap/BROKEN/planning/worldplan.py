import dilap.core.base as db

###############################################################################
###############################################################################

class wplan(db.base):
    #wplan is a persistent representation of a plan to construct a world
    #it is effectively a graph encoding that should permit deterministic 
    #generation of a world, possibly using mpi and many machines

    def __init__(self,*args,**kwargs):
        self._def('sealevel',0,**kwargs)
        self._def('iotype','obj',**kwargs)

# interpret a wplan into the files necessary
# to represent the world, such as obj files
def interpret(wplan,fid = None):
    raise NotImplemented

class wfactory(db.base):
    #wfactory generates wplans using a possibly nondeterminstic process

    # 1. lay out the ground level infrastructure - represent as a graph
    # 2. create the set of areas separated by this infrastructure
    # 3. in each area, using stylistic sort of choices, 
    #   divide up each area into structures and terrain which falls between
    # 4. paint the terrain with objects such as foliage and trash
    # 5. construct the graph of each building from its footprint
    # 6. paint the interior spaces of these buildings with objects

    def __init__(self,*args,**kwargs):
        pass





