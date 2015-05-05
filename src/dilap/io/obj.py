import dilap.core.base as db
import dilap.core.uinfo as di
import dilap.core.material as dma

#import cStringIO as sio
from io import StringIO as sio
import os,pdb

user_info = di.fetch_info()

world_dir = user_info['contentdir']

#########################################################################
### materials
#########################################################################

# initialize the materials script
def write_materials():
    matfile = os.path.join(world_dir,'materials.mtl')
    matstring = sio()
    #matstring = sio.StringIO()
    matstring.write('\n')
    write_default_materials_mtl(matstring)
    matstring.write('\n')
    with open(matfile,'w') as h:
        h.write(matstring.getvalue())
    print('new materials file',matfile)

def write_default_materials_mtl(msio):
    generic = dma.material('generic',dtexture = 'orangeboxtex.png')
    def_mats = [generic]
    mcount = len(def_mats)
    msio.write('# Material Count: ')
    msio.write(str(mcount))
    msio.write('\n')
    for dm in def_mats:dm._write('obj',msio)

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
    faces = mod._face_dict()
    mats = [m for m in faces.keys()]
    mcnt = len(mats)

    sioio = sio()
    #sioio = sio.StringIO()
    sioio.write('mtllib materials.mtl\n')
    for vdx in range(len(mod.pcoords)):
        pvert = mod.pcoords[vdx]
        nvert = mod.ncoords[vdx]
        uvert = mod.ucoords[vdx]
        sioio.write( 'v %f %f %f\n'%(pvert.x,pvert.y,pvert.z))
        sioio.write('vn %f %f %f\n'%(nvert.x,nvert.y,nvert.z))
        sioio.write('vt %f %f\n'   %(uvert.x,uvert.y))
        
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
def build_model(mod,**kwargs):
    orep,ofile = obj_from_model(mod)
    mod.reps['obj'] = orep
    objpath = os.path.join(world_dir,ofile)
    with open(objpath,'w') as h:h.write(orep)


