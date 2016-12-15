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
    js.write('\nSceneJS.setConfigs({\n\tpluginPath: "scenejs/api/latest/plugins"\n});\n')
    #js.write('\nSceneJS.createScene({\n\tnodes: [\n\t\t{\n\t\t\ttype: "cameras/orbit",')
    js.write('\nSceneJS.createScene({\n\tnodes: [\n\t\t{\n\t\t\ttype: "cameras/pickFlyOrbit",')
    js.write('\n\t\t\tyaw: -40,\n\t\t\tpitch: -20,\n\t\t\tzoom: 200,')
    js.write('\n\t\t\tzoomSensitivity: 10.0,')
    js.write('\n\t\t\tnodes: [\n\t\t\t\t{\n\t\t\t\t\ttype: "rotate",\n\t\t\t\t\tx: 1,')
    js.write('\n\t\t\t\t\tangle: -90,\n\t\t\t\t\tnodes: [')
    js.write('\n\t\t\t\t{\n\t\t\t\ttype: "texture",')
    js.write('\n\t\t\t\tsrc: "world0/orangeboxtex.png",\n\t\t\t\tnodes: [\n')
    nodes = [imp_line % (f,) for f in obj_filenames]
    js.write(',\n'.join(nodes))
    js.write('\n\t\t\t\t]\n\t\t\t\t}\n\t\t\t\t]\n\t\t\t\t}\n\t\t\t]\n\t\t}\n\t]\n});')
    with open(worldfile,'w') as h:h.write(js.getvalue())
    print('new world file',worldfile)

imp_line = '''\
{
    type: "import/obj",
    src: "world0/%s",
}'''

#########################################################################

#########################################################################
### obj file creation
#########################################################################

obj_filenames = []
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
    ofile = unique_objfile(mod.filename)
    faces = mod.face_dict()
    mats = [m for m in faces.keys()]
    mcnt = len(mats)

    sioio = sio()
    sioio.write('mtllib materials.mtl\n')

    for p,n,u in zip(mod.pset,mod.nset,mod.uset):
        sioio.write( 'v %f %f %f\n'%(p.x,p.y,p.z))
        sioio.write('vn %f %f %f\n'%(n.x,n.y,n.z))
        sioio.write('vt %f %f\n'   %(u.x,u.y))
        
    for mdx in range(mcnt):
        m = mats[mdx]
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










