from .obj import interface as obj
from io import StringIO as sio
import shutil
import os
import pdb


class interface(obj):


    @staticmethod
    def build_scenegraph(sg,path):
        site = obj.resource_path('site')[0]
        sitepath = os.path.join(path,'site')
        if os.path.isdir(sitepath):
            shutil.rmtree(sitepath)
        shutil.copytree(site,sitepath)
        worldpath = os.path.join(sitepath,'world')
        if os.path.isdir(worldpath):
            shutil.rmtree(worldpath)
        os.makedirs(worldpath)
        obj.build_scenegraph(sg,worldpath)
        write_world_script(sitepath)


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
    for name,texture in obj.def_mats:
        if name in obj.tobj_filenames:
            textures[name] = (texture,obj.tobj_filenames[name])
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
    src: "../world/%s",
    nodes: [\n%s\n]
}'''
modnode_line = '''\
{
    type: "import/obj",
    src: "../world/%s",
}'''
