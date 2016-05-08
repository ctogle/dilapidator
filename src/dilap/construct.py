from dilap.geometry.vec3 import vec3

import dilap.core.uinfo as di
import dilap.core.context as dgc
import dilap.io.io as dio

import dilap.modeling.model as dmo

import dilap.topology.worldly.treeskin as ltr
import dilap.topology.worldly.building as blg

import dilap.core.plotting as dtl
import matplotlib.pyplot as plt

import pdb

iotypes = dio.iotypes

###############################################################################
### functions to realize models,nodes,contexts
###############################################################################

# given a model, output its representation in world space
def build2(mod,io = None):
    if io is None:io = di.fetch_info()['exporter']
    elif type(io) is type(''):io = iotypes[io]
    io.build_model2(mod)

#def realize(context,years = 0):
def realize(context,years = 0,io = None):
    if io is None:io = di.fetch_info()['exporter']
    elif type(io) is type(''):io = iotypes[io]
    context.generate(worn = years)
    #context.passtime(years)
    context.graph(io)

###############################################################################

###############################################################################

def context(io = 'obj',dilaps = []):
    cx = dgc.context(iotype = io,dilapidors = dilaps)
    return cx

###############################################################################

###############################################################################

def tridome(mod):
    gmesh = mod.atridome()
    mod.subdiv(gmesh,True)
    mod.subdiv(gmesh,True)
    mod.subdiv(gmesh,True)
    mod.subdiv(gmesh,True)
    #mod.subdiv(gmesh,True)
    #mod.subdiv(gmesh,False)
    #mod.subdiv(gmesh,True)
    return mod

def house(mod):
    gm = mod.agfxmesh()

    eb = (vec3(-2,-2,0),vec3(2,-2,0),vec3(2,2,0),vec3(-2,2,0))
    ibs = ()
    hmin,ref,smo = 1,False,False

    gm.tripoly(eb,ibs,hmin,ref,smo)

    #v1  = gm.avert(*mod.avert(vec3(-1,-1,-1)))
    #v2  = gm.avert(*mod.avert(vec3( 1,-1,-1)))
    #v3  = gm.avert(*mod.avert(vec3( 1, 1,-1)))
    #v4  = gm.avert(*mod.avert(vec3(-1, 1,-1)))
    #f1  = gm.aface(v1,v2,v3) 
    #f2  = gm.aface(v1,v3,v4) 
    return mod





def teststage(**kwargs):
    #kwargs['years'] = 0
    #p,d = vec3(0,0,0),vec3(0,0,1)
    #ax = dtl.plot_axes()
    #cx = ltr.tree(p,d,ax = ax)
    #realize(cx,**kwargs)
    kwargs['years'] = 0
    ax = dtl.plot_axes()
    cx = blg.building(ax = ax)
    realize(cx,**kwargs)






def teststageold(**kwargs):
    p,d = vec3(0,0,0),vec3(0,0,1)
    ax = dtl.plot_axes()
    mod = ltr.treeskin(p,d,ax = ax).skeleton()

    #mod = dmo.model()

    #mod = tridome(mod)
    #mod = house(mod)

    #p,d = vec3(0,0,0),vec3(0,0,1)
    #ax = dtl.plot_axes()
    #gm = mod.agfxmesh(ltr.treeskinmesh(p,d,ax = ax).skeleton(mod))

    '''#
    ax = dtl.plot_axes()
    for gmesh in mod.gfxmeshes:
        for f in gmesh.faces:
            ps = mod.gvps(gmesh,f)
            ax = dtl.plot_polygon(ps,ax)
    '''#

    plt.show()

    print('build2 cube now')
    build2(mod,**kwargs)

###############################################################################
###############################################################################








                                                                                  

