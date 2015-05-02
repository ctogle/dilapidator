import dilap.core.uinfo as di
import dilap.core.sgraph as dsg
import dilap.core.model as dm
import dilap.io.obj as objio
import dilap.primitive.cube as dcu
import dilap.primitive.cone as dco
import dilap.primitive.cylinder as dcyl
import dilap.primitive.wall as dw
import dilap.primitive.floor as df

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
    cu = dcu.cube()
    cu.scale_u(l)
    return cu

# return a cylinder model of radius r, height h, with n sides
def cylinder(r = 1.0,h = 1.0,n = 8):
    cu = dcyl.cylinder(n = n)
    cu.scale_u(r).scale_z(h/float(r))
    return cu

# return a cylinder model of radius r, height h, with n sides
def cone(r = 1.0,h = 1.0,n = 8):
    cu = dco.cone(n = n)
    cu.scale_u(r).scale_z(h/float(r))
    return cu

def wall(v1,v2,h = 1.0,w = 0.5,gs = []):
    wa = dw.wall(v1,v2,h = h,w = w,gaps = gs)
    return wa

def perimeter(vs,h = 1.0,w = 0.5):
    v1,v2 = vs[0],vs[1]
    pwa = wall(v1,v2,h = h,w = w)
    for vdx in range(2,len(vs)):
        v1,v2 = vs[vdx-1],vs[vdx]
        wa = wall(v1,v2,h = h,w = w)
        pwa._consume(wa)
    return pwa

def floor(l = 10.0,w = 10.0,h = 0.5,gap = None,m = 'generic'):
    fl = df.floor(l,w,h = h,gap = gap,m = m)
    return fl

                                                                                  
