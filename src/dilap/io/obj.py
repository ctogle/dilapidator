from io import StringIO as sio
import appdirs
import shutil
import glob
import os
import pdb


class interface(object):


    res_path = os.path.join(appdirs.user_data_dir(),'dilap_resources')
    def_mats = []
    obj_filenames = []
    tobj_filenames = {}


    @staticmethod
    def resource_path(res = None):
        '''Find a resource file or the directory where it resources are stored'''
        if res is None:
            rpaths = [interface.res_path[:]]
        else:
            rpaths = []
            for x in os.walk(interface.res_path):
                for y in glob.glob(os.path.join(x[0],res)):
                    rpaths.append(y)
                for d in x[1]:
                    z = os.path.join(x[0],d)
                    if res == d:
                        rpaths.append(z)
        return rpaths


    @staticmethod
    def build_scenegraph(sg,path):
        write_materials(path)
        for m in sg.worldspace():
            orep,ofile = obj_from_model(m)
            m.reps['obj'] = orep
            objpath = os.path.join(path,ofile)
            with open(objpath,'w') as h:
                h.write(orep)


def write_materials(world_dir):
    matfile = os.path.join(world_dir,'materials.mtl')
    matstring = sio()
    matstring.write('\n')
    write_default_materials_mtl(matstring)
    matstring.write('\n')
    with open(matfile,'w') as h:
        h.write(matstring.getvalue())
    for name,texture in interface.def_mats:
        shutil.copy(interface.resource_path(texture)[0],world_dir)
    print('new materials file',matfile)


def write_default_materials_mtl(msio):
    generic = ('generic','orangeboxtex.png')
    interface.def_mats.append(generic)
    grass2 = ('grass2','grass2.jpg')
    interface.def_mats.append(grass2)
    concrete1 = ('concrete1','concrete1.png')
    interface.def_mats.append(concrete1)
    mcount = len(interface.def_mats)
    msio.write('# Material Count: ')
    msio.write(str(mcount))
    msio.write('\n')
    for name,texture in interface.def_mats:
        msio.write('\nnewmtl ')
        msio.write(name)
        msio.write('\n')
        msio.write('Ka 1.000 1.000 1.000\n')
        msio.write('Kd 1.000 1.000 1.000\n')
        msio.write('Ks 0.000 0.000 0.000\n')
        msio.write('d 1.0\n')
        msio.write('illum 2\n')
        dtexture = interface.resource_path(texture)[0]
        msio.write('map_Kd '+dtexture+'\n')
        msio.write('\n')


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
