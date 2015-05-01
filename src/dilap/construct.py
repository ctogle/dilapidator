import dilap.core.uinfo as di
import dilap.core.sgraph as dsg
import dilap.io.obj as objio
import dilap.primitive.cube as dc

import types

iotypes = {
    'obj':objio,
        }

def manufacture(mods,io = None):
    if io is None:io = di.fetch_info()['exporter']
    if not type(mods) is types.ListType:mods = [mods]
    sgr = dsg.sgraph(nodes = mods)
    sgr.graph(iotype[io])

def cube(l = 1.0):
    cu = dc.cube()
    cu.scale_u(l)
    manufacture(cu)
    return cu


