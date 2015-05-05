import dilap.core.base as db

###############################################################################
### dilapidor modifies a context to apply some time-dep effect
###
### this should include things like ivy overgrowing sturctures
### rain wearing away stone or land
### decay and collapse of non-concrete structures
### sinking and deforming of structures and their supports
###
### it should effectively simulate the passage of time on a context
###############################################################################

class dilapidor(db.base):

    def __init__(self,*args,**kwargs):
        pass

    # strength of effect should somehow be proportional to years
    # effects are applied to the nodes in the scenegraph of context
    def wither(self,context,years):
        print('base dilapidor withering',years,'on context',context)


