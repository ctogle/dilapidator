import dilap.core.base as db
import dilap.core.uinfo as di
#import dilap.core.material as dma
import dilap.modeling.material as dma

#import cStringIO as sio
from io import StringIO as sio
import os,shutil,pdb
user_info = di.fetch_info()

#world_dir = user_info['contentdir']

#########################################################################
### materials
#########################################################################

# initialize the materials script
def write_materials(world_dir):
    matfile = os.path.join(world_dir,'materials.mtl')
    matstring = sio()
    matstring.write('\n')
    write_default_materials_mtl(matstring)
    matstring.write('\n')
    with open(matfile,'w') as h:
        h.write(matstring.getvalue())
    for dm in def_mats:
        shutil.copy(db.resource_path(dm.dtexture),world_dir)
    print('new materials file',matfile)

def_mats = []
def write_default_materials_mtl(msio):
    global def_mats
    generic = dma.material('generic',dtexture = 'orangeboxtex.png')
    def_mats.append(generic)
    grass2 = dma.material('grass2',dtexture = 'grass2.jpg')
    def_mats.append(grass2)
    concrete1 = dma.material('concrete1',dtexture = 'concrete1.png')
    def_mats.append(concrete1)
    mcount = len(def_mats)
    msio.write('# Material Count: ')
    msio.write(str(mcount))
    msio.write('\n')
    for dm in def_mats:dm._write('obj',msio)

# write a script based on a scenegraph of models
# this version will write a javascript using scenejs
def write_world_script(world_dir):
    worldfile = os.path.join(world_dir,'world.js')
    js = sio()
    d = 0
    s = lambda l : js.write('\n'+'\t'*d+l)
    s('SceneJS.setConfigs({');d += 1
    s('pluginPath: "scenejs_api/latest/plugins"');d -= 1
    s('});\n')
    s('SceneJS.createScene({');d += 1
    s('nodes: [');d += 1
    s('{');d += 1

    #s('type: "environments/holodeck",')
    #s('type: "environments/modelView",')
    #s('type: "environments/gridRoom",')
    #s('type: "environments/lawn",')

    #s('type: "scenejs_api/latest/plugins/nodes/models/backgrounds/gradient"')
    #depth: -30, (default)
    #colors:[0.05, 0.06, 0.07, 1.0, // top left (R,G,B,A)
    #s('colors: [0.05, 0.06, 0.07, 1.0, 0.05, 0.06, 0.07, 1.0, 0.85, 0.9, 0.98, 1.0, 0.85, 0.9, 0.98, 1.0]')

    #s('nodes: [');d += 1
    #s('{');d += 1

    s('type: "cameras/pickFlyOrbit",')
    s('yaw: -40,')
    s('pitch: -20,')
    s('zoom: 200,')
    s('zoomSensitivity: 10.0,')
    s('nodes: [');d += 1
    s('{');d += 1
    s('type: "rotate",')
    s('x: 1,')
    s('angle: -90,')
    s('nodes: [');d += 1
    textures = {}
    for dfm in def_mats:
        if dfm.name in tobj_filenames:
            textures[dfm.name] = (dfm.dtexture,tobj_filenames[dfm.name])
    modline = lambda t : ',\n'.join([modnode_line % (f,) for f in textures[t][1]])
    matline = ',\n'.join([matnode_line % (textures[t][0],modline(t)) for t in textures])
    s(matline)
    d -= 1;s(']');d -= 1;s('}');
    d -= 1;s(']');d -= 1;s('}');d -= 1;s(']');d -= 1;s('});')
    with open(worldfile,'w') as h:h.write(js.getvalue())
    print('new world file',worldfile)

matnode_line = '''\
{
    type: "texture",
    src: "world0/%s",
    nodes: [\n%s\n]
}'''
modnode_line = '''\
{
    type: "import/obj",
    src: "world0/%s",
}'''

#########################################################################

#########################################################################
### obj file creation
#########################################################################

obj_filenames = []
tobj_filenames = {}
def unique_objfile(ofile):
    if not ofile.endswith('.obj'):ofile += '.obj'
    if ofile in obj_filenames:
        ofile = ofile[:ofile.rfind('mesh.obj')]
        onum = len(obj_filenames) + 1
        ofile += '.'.join([str(onum),'mesh','obj'])
    print('new objfile',ofile)
    obj_filenames.append(ofile)
    return ofile

# represent the model mod as a string and give it a safe filename
def obj_from_model(mod):
    global tobj_filenames
    ofile = unique_objfile(mod.filename)
    faces = mod.face_dict()
    mats = [m for m in faces.keys()]
    mcnt = len(mats)

    #if mcnt > 1:
    #    print('js world doesnt seem to support multi-material obj files')
    #objmat = mod.
    #pdb.set_trace()

    sioio = sio()
    sioio.write('mtllib materials.mtl\n')

    for p,n,u in zip(mod.pset,mod.nset,mod.uset):
        sioio.write( 'v %f %f %f\n'%(p.x,p.y,p.z))
        sioio.write('vn %f %f %f\n'%(n.x,n.y,n.z))
        sioio.write('vt %f %f\n'   %(u.x,u.y))
        
    for mdx in range(mcnt):
        m = mats[mdx]

        if m in tobj_filenames:tobj_filenames[m].append(ofile)
        else:tobj_filenames[m] = [ofile]

        mfaces = faces[m]
        fcnt = len(mfaces)
        sioio.write('usemtl ')
        sioio.write(m)
        sioio.write('\n')
        sioio.write('s off\n')
        for fdx in range(fcnt):
            f = mfaces[fdx]
            f1 = f[0] + 1
            f2 = f[1] + 1
            f3 = f[2] + 1
            sioio.write('f')
            sioio.write(' %i/%i/%i'%(f1,f1,f1))
            sioio.write(' %i/%i/%i'%(f2,f2,f2))
            sioio.write(' %i/%i/%i'%(f3,f3,f3))
            sioio.write('\n')

    '''#
    sioio.write('usemtl Material\n')
    sioio.write('s off\n')
    for face in prim.faces:
        sioio.write('f')
        for vert in face:
            vi = vert + 1
            sioio.write( ' %i/%i/%i' % (vi,vi,vi) )
        sioio.write('\n')
    '''#

    orep = sioio.getvalue()
    return orep,ofile

#########################################################################

#########################################################################
### build functions
#########################################################################

# create models as obj files
def build_models(*mods,**kwargs):
    for m in mods:
        if m is None:return         
        else:build_model(m,**kwargs)
    
# create one model as one obj file
def build_model2(mod,world_dir = None,**kwargs):
    orep,ofile = obj_from_model(mod)
    mod.reps['obj'] = orep
    objpath = os.path.join(world_dir,ofile)
    with open(objpath,'w') as h:h.write(orep)










