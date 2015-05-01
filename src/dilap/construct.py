import dilap.io.obj as objio
import dilap.core.sgraph as dsg
import dilap.primitives.cube as dc

import types

iotypes = {
    'obj':objio,
        }

def manufacture(mods,io = 'obj'):
    if not type(mods) is types.ListType:mods = [mods]
    sgr = dsg.sgraph(nodes = mods)
    sgr.graph(iotype[io])

def cube(l = 1.0,io = 'obj'):
    cu = dc.cube()
    cu.scale_u(l)
    manufacture(cu,io = io)
    return cu


