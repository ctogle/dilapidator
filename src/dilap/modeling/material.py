import dilap.core.base as db

import pdb

__doc__ = '''A material class'''

unused_material_id = 0
class material(db.base):
    '''A material class'''

    def _dpid(self):
        global unused_material_id
        self.dpid = unused_material_id
        unused_material_id += 1

    def __init__(self,name,**kwargs):
        self._dpid()
        self.name = name
        self._def('dtexture',None,**kwargs)

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
    def _write_obj(self,msio):
        msio.write('\nnewmtl ')
        msio.write(self.name)
        msio.write('\n')
        msio.write('Ka 1.000 1.000 1.000\n')
        msio.write('Kd 1.000 1.000 1.000\n')
        msio.write('Ks 0.000 0.000 0.000\n')
        msio.write('d 1.0\n')
        msio.write('illum 2\n')
        dtexture = db.resource_path(self.dtexture)[0]
        msio.write('map_Kd '+dtexture+'\n')
        msio.write('\n')

    def _write(self,exp,msio):
        if exp == 'obj':self._write_obj(msio)
        else:pdb.set_trace()


