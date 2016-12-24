from dilap.geometry.vec3 import vec3

import dilap.core.uinfo as di
import dilap.core.context as dgc
import dilap.io.io as dio

import dilap.modeling.model as dmo
import dilap.worldly.world as dwo

import dilap.core.plotting as dtl
import matplotlib.pyplot as plt

import os
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

def realize(context,years = 0,io = None,**kws):
    if io is None:io = di.fetch_info()['exporter']
    elif type(io) is type(''):io = iotypes[io]
    context.generate(worn = years)
    context.graph(io,**kws)

def write_materials(io = None,world_dir = None,**kws):
    if world_dir is None:world_dir = os.getcwd()
    if io is None:io = di.fetch_info()['exporter']
    io = iotypes[io]
    if hasattr(io,'write_materials'):
        io.write_materials(world_dir)

def write_world_script(io = None,world_dir = None,**kws):
    if world_dir is None:world_dir = os.getcwd()
    if io is None:io = di.fetch_info()['exporter']
    io = iotypes[io]
    if hasattr(io,'write_world_script'):
        io.write_world_script(world_dir)

# an example of a minimal world context
def cube_stage(p = None,q = None,s = None):
    m = dmo.model()
    cx = dcx.context()
    sgv = cx.amodel(p,q,s,m,cx.sgraph.root)
    gm = m.atricube()
    m.normals(gm)
    return cx

# stage is a function that returns a context for a world
# wdir is an output directory to associate with this world
def world(stage = cube_stage,wdir = None):
    if wdir is None:wdir = os.getcwd()
    write_materials(world_dir = wdir)
    realize(stage(),world_dir = wdir)
    write_world_script(world_dir = wdir)

###############################################################################

'''#
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
'''#

def teststage(**kws):
    kws['years'] = 0

    #p,d = vec3(0,0,0),vec3(0,0,1)
    #ax = dtl.plot_axes()
    #cx = ltr.tree(p,d,ax = ax)
    #realize(cx,**kws)
    #ax = dtl.plot_axes()
    #cx = blg.building()

    '''#

    bfp = vec3(0,0,0).sq(50,15)
    bsq = test_bseq()
    p,q,s = None,None,None
    bfa = dwo.blg.blgfactory()
    cx = bfa.new(p,q,s,footprint = bfp,sequence = bsq,floorheight = 5)

    '''#

    #cx = dwo.world()
    cx = dwo.worldfactory().new()
    realize(cx,**kws)

def test_bseq():
    seq  = ''
    seq += 'L<0>'
    seq += 'L<1>'
    fseq = 'S<0,0.3,0.5,0,0,1,0>' # make living room on left
    fseq += 'S<0,0.3,0.5,0,0,1,0>' # make living room on left
    fseq += 'S<0,0.3,0.5,0,1,0,0>' # make living room on left
    fseq += 'R<2,rtype,closed>'
    fseq += 'C<0,0.2,0.2,0.2,0.2>'
    fseq += 'C<1,0.2,0.2,0.2,0.2>'
    fseq += 'C<3,0.2,0.2,0.2,0.2>'
    fseq += 'C<2,0.2,2,0.2,2>'
    fseq += 'E<1,0>'
    fseq += 'E<1,2>'
    fseq += 'E<2,3>'
    fseq += 'E<2,0>'
    #fseq += 'X<3>'
    seq += 'I<0,'+fseq+'>'
    seq += 'I<1,'+fseq+'>'
    #seq += 'E<1,0>'
    #seq += 'V<2>'
    #seq += 'V<6>'
    seq += 'X<2>'
    return seq

def atest_bseq():
    seq  = ''
    seq += 'S<0,0.3,0.5,0,0,1,0>' # make living room on left
    seq += 'S<0,0.6,0.5,0,0,1,0>' # make back room on right
    seq += 'S<2,0.5,0.4,0,1,0,0>' # make kitchen on bottom
    seq += 'S<3,0.5,0.3,0,1,0,0>' # make front room on top
    seq += 'S<2,0.6,0.5,0,0,1,0>' # make bathrom on left of kitchen
    seq += 'S<5,0.5,0.7,0,1,0,0>' # make hall closest above kitchen
    seq += 'S<2,0.5,0.4,0,1,0,0>'
    seq += 'S<0,0.5,0.3,0,1,0,0>'
    seq += 'S<4,0.3,0.5,0,0,1,0>'
    seq += 'S<9,0.5,0.5,0,1,0,0>'
    seq += 'E<1,5>E<1,3>E<5,2>E<3,7>E<3,4>E<3,8>E<3,6>E<8,0>E<4,9>E<4,10>X<1>'
    seq += 'R<1,rtype,open>'
    seq += 'R<3,rtype,open>'
    seq += 'R<5,rtype,open>'
    seq += 'R<0,rtype,closed>'
    seq += 'R<10,rtype,closed>'
    return seq

###############################################################################
###############################################################################








                                                                                  

