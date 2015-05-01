import dilap.core.uinfo as di
import dilap.core.sgraph as dsg
import dilap.core.model as dm
import dilap.io.obj as objio
import dilap.primitive.cube as dc

import types

iotypes = {
    'obj':objio,
        }

# given either one or many models or nodes
# build all models in world space using iotype io
def build(mods,consume = False,io = None):
    if io is None:io = di.fetch_info()['exporter']
    if not type(mods) is types.ListType:mods = [mods]
    for mdx in range(len(mods)):
        m = mods[mdx]
        if issubclass(m.__class__,dm.model): 
            mods[mdx] = dsg.node(models = [m])
    if consume:
        consumenode = dsg.node(children = mods,consumption = True)
        mods = [consumenode]
    sgr = dsg.sgraph(nodes = mods)
    sgr.graph(iotypes[io])

# return a cube model with sidelength l
def cube(l = 1.0):
    cu = dc.cube()
    cu.scale_u(l)
    return cu


