import dilap.core.base as db
import dilap.core.uinfo as di
#import dilap.core.material as dma
import dilap.modeling.material as dma

#import cStringIO as sio
from io import StringIO as sio
import os,shutil,pdb
#user_info = di.fetch_info()

#world_dir = user_info['contentdir']

class interface(object):


    def_mats = []
    obj_filenames = []
    tobj_filenames = {}


    @staticmethod
    def build_context(cx,path):
        for ch in cx.children:
            interface.build_context(ch,path)
        write_materials(path)
        sg = cx.sgraph
        for m in sg.graphvert(sg.root,None):
            orep,ofile = obj_from_model(m)
            m.reps['obj'] = orep
            objpath = os.path.join(path,ofile)
            with open(objpath,'w') as h:
                h.write(orep)


# initialize the materials script
def write_materials(world_dir):
    matfile = os.path.join(world_dir,'materials.mtl')
    matstring = sio()
    matstring.write('\n')
    write_default_materials_mtl(matstring)
    matstring.write('\n')
    with open(matfile,'w') as h:
        h.write(matstring.getvalue())
    for dm in interface.def_mats:
        shutil.copy(db.resource_path(dm.dtexture)[0],world_dir)
    print('new materials file',matfile)

#def_mats = []
def write_default_materials_mtl(msio):
    #global def_mats
    generic = dma.material('generic',dtexture = 'orangeboxtex.png')
    interface.def_mats.append(generic)
    grass2 = dma.material('grass2',dtexture = 'grass2.jpg')
    interface.def_mats.append(grass2)
    concrete1 = dma.material('concrete1',dtexture = 'concrete1.png')
    interface.def_mats.append(concrete1)
    mcount = len(interface.def_mats)
    msio.write('# Material Count: ')
    msio.write(str(mcount))
    msio.write('\n')
    for dm in interface.def_mats:dm._write('obj',msio)

obj_filenames = []
tobj_filenames = {}
def unique_objfile(ofile):
    if not ofile.endswith('.obj'):ofile += '.obj'
    if ofile in interface.obj_filenames:
        ofile = ofile[:ofile.rfind('mesh.obj')]
        onum = len(interface.obj_filenames) + 1
        ofile += '.'.join([str(onum),'mesh','obj'])
    print('new objfile',ofile)
    interface.obj_filenames.append(ofile)
    return ofile

# represent the model mod as a string and give it a safe filename
def obj_from_model(mod):
    #global tobj_filenames
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

        if m in interface.tobj_filenames:
            interface.tobj_filenames[m].append(ofile)
        else:
            interface.tobj_filenames[m] = [ofile]

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
