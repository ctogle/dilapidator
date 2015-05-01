import dilap.core.base as db
import dilap.core.uinfo as di

import cStringIO as sio
import os,pdb

user_info = di.fetch_info()

world_dir = user_info['contentdir']

#########################################################################
### materials
#########################################################################

matstring = sio.StringIO()

# initialize the materials script
def reset_materials_script():
    matfile = os.path.join(world_dir, 'materials.mtl')
    matstring.write('\n')
    write_default_materials_mtl(matstring)
    matstring.write('\n')
    with open(matfile, 'w') as handle:
        handle.write(matstring.getvalue())

def write_default_materials_mtl(msio):
    cubemat = material('cubemat',diffuseMap = 'textures/generic/orangeboxtex.png')
    sidewalk1 = material('sidewalk1',diffuseMap = 'textures/concrete/sidewalk1.jpg')
    asphalt = material('asphalt',diffuseMap = 'textures/concrete/asphalt.jpg')
    roadline_y = material('roadline_y',diffuseMap = 'textures/concrete/roadline.png')
    roadline_y_cont = material('roadline_y_cont',diffuseMap = 'textures/concrete/roadline_w_continuous.png')
    roadline_w_cont = material('roadline_w_cont',diffuseMap = 'textures/concrete/roadline_w_continuous.png')
    def_mats = [cubemat,sidewalk1,asphalt,roadline_y,roadline_y_cont,roadline_w_cont]

    mcount = len(def_mats)
    msio.write('# Material Count: ')
    msio.write(str(mcount))
    msio.write('\n')
    for dm in def_mats: dm._write(msio)

class material(db.base):

    def _write_path_property(self,msio,prop,val):
        if val is None:return
        msio.write(prop)
        msio.write(' ')
        msio.write(val)
        msio.write('\n')

    def _write(self,msio):
        msio.write('\nnewmtl ')
        msio.write(self.name)
        msio.write('\n')
        msio.write('Ka 1.000 1.000 1.000\n')
        msio.write('Kd 1.000 1.000 1.000\n')
        msio.write('Ks 0.000 0.000 0.000\n')
        msio.write('d 1.0\n')
        msio.write('illum 2\n')
        self._write_properties(msio)
        msio.write('\n')

    def _write_properties(self,msio):
        self._write_diffuse_properties(msio)

    def _write_diffuse_properties(self,msio):
        #self._write_vector_property(msio,'diffuseColour',self.diffuseColour)
        #self._write_path_property(msio,'diffuseMap',self.diffuseMap)
        #self._write_bool_property(msio,'vertexDiffuse',self.vertexDiffuse)
        self._write_path_property(msio,'map_Kd',self.diffuseMap)

    def __init__(self,name,**kwargs):
        self.name = name
        self._default_('diffuseMap',None,**kwargs)

'''#
newmtl Textured
   Ka 1.000 1.000 1.000
   Kd 1.000 1.000 1.000
   Ks 0.000 0.000 0.000
   d 1.0
   illum 2
   map_Ka lenna.tga           # the ambient texture map
   map_Kd lenna.tga           # the diffuse texture map (most of the time, it will
                              # be the same as the ambient texture map)
   map_Ks lenna.tga           # specular color texture map
   map_Ns lenna_spec.tga      # specular highlight component
   map_d lenna_alpha.tga      # the alpha texture map
   map_bump lenna_bump.tga    # some implementations use 'map_bump' instead of 'bump' below
'''#

#########################################################################

obj_filenames = []
def unique_objfile(ofile):
    if not ofile.endswith('.obj'):ofile += '.obj'
    if ofile in obj_filenames:
        ob = ofile[:ofile.rfind('mesh.obj')]
        onum = len(obj_filenames) + 1
        om += '.'.join([str(onum),'mesh','obj'])
        ofile = ofile.replace(ob,om)
    print 'new objfile',ofile
    obj_filenames.append(ofile)
    return ofile

# represent the model mod as a string and give it a safe filename
def obj_from_model(mod):
    ofile = unique_objfile(mod.filename)
    faces = mod._face_dict()
    mats = faces.keys()
    mcnt = len(mats)

    sioio = sio.StringIO()
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

#########################################################################
#########################################################################










