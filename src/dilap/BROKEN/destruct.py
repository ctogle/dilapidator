import dilap.core.uinfo as di
import dilap.core.sgraph as dsg
import dilap.core.model as dm
import dilap.io.io as dio
import dilap.degenerate.overgrown as dog

import pdb

iotypes = dio.iotypes

###############################################################################
### simple functions which return simple dilapidors
###############################################################################

def overgrown(z_max = 10):
    ivy = dog.ivy(z_max = z_max)
    return ivy

###############################################################################
### convenient collections of functions
###############################################################################

# dilapidors is a dict of funcs which return dilapidor objects
dilapidors = {
    'ivy':overgrown,
}

###############################################################################

                                                                                  
