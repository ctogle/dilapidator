import dilap.core.uinfo as di
import dilap.core.sgraph as dsg
import dilap.core.model as dm
import dilap.io.io as dio
import dilap.primitive.cube as dcu
import dilap.primitive.cone as dco
import dilap.primitive.cylinder as dcyl
import dilap.primitive.wall as dw
import dilap.primitive.floor as df
import dilap.primitive.pipe as dp
import dilap.primitive.road as dr
import dilap.generate.context as dgc
import dilap.generate.landscape as dls
import dilap.generate.lot as dlot
import dilap.generate.street as dstr
import dilap.generate.continent as dct

import dp_vector as dpv

import pdb

iotypes = dio.iotypes

###############################################################################
### functions to realize models,nodes,contexts
###############################################################################

# given either one or many models or nodes
# build all models in world space using iotype io
def build(mods,consume = False,io = None):
    if io is None:io = di.fetch_info()['exporter']
    if not type(mods) is type([]):mods = [mods]
    for mdx in range(len(mods)):
        m = mods[mdx]
        if issubclass(m.__class__,dm.model): 
            mods[mdx] = dsg.node(models = [m])
    if consume:
        consumenode = dsg.node(children = mods,consumption = True)
        mods = [consumenode]
    sgr = dsg.sgraph(nodes = mods)
    sgr.graph(iotypes[io])

def realize(context,years = 0):
    context.generate(worn = years)
    context.passtime(years)
    context.graph()

###############################################################################

###############################################################################
### simple functions which return simple models
###############################################################################

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

def pipe(curve = None,loop = None,m = 'generic'):
    pi = dp.pipe(curve = curve,loop = loop,m = m)
    return pi

def road(start = None,end = None,tip = None,tail = None):
    if start is None:start = dpv.zero()
    if end is None:end = dpv.vector(100,100,-10)
    if tip is None:tip = dpv.vector(0,1,0)
    if tail is None:tail = dpv.vector(0,1,0)
    rd = dr.road(start,end,tip,tail)
    return rd

###############################################################################

def context(io = 'obj',dilaps = []):
    cx = dgc.context(iotype = io,dilapidors = dilaps)
    return cx

def landscape(io = 'obj',dilaps = []):
    cx = dls.landscape(iotype = io,dilapidors = dilaps)
    return cx

def lot(io = 'obj',dilaps = []):
    cx = dlot.lot(iotype = io,dilapidors = dilaps)
    return cx

def street(io = 'obj',dilaps = []):
    cx = dstr.street(iotype = io,dilapidors = dilaps)
    return cx

def continent(io = 'obj',dilaps = []):
    cx = dct.continent(iotype = io,dilapidors = dilaps)
    return cx

###############################################################################

###############################################################################
### convenient collections of functions
###############################################################################

# generator is a dict of funcs which return model objects
generator = {
    'cube':cube,
    'cylinder':cylinder,
    'cone':cone,
    'wall':wall,
    'perimeter':perimeter,
    'floor':floor,
    'pipe':pipe,
    'road':road,
    #'context':context,
    #'lot':lot,
}

# contextualizer is a dict of funcs which return context objects
contextualizer = {
    'context':context,
    'landscape':landscape,
    'lot':lot,
    'street':street,
    'continent':continent,
}

###############################################################################

                                                                                  
