import dilap.io.obj as dobj
from dilap.worldly import *
from dilap.geometry import vec3
import bpy,sys
from bpy_extras.io_utils import unpack_list
from bpy_extras.io_utils import unpack_face_list


bl_info = {
    'name':'dilapidator',
    'description':'dilapidator procedural mesh generator',
    'category':'Object',
    'author':'Curtis Ogle',
    'version':(1,0),
}


class dilap_panel(bpy.types.Panel):
    bl_label = "dilap Generator Settings"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_context = "scene"

    def draw(self,context):
        box = self.layout
        col = box.column(align=True)
        col.operator('object.dilaprun',icon="SCRIPT")
        col.operator('object.dilappurge',icon="SCRIPT")


def register():
    print('register dilapidator')
    bpy.utils.register_module(__name__)


def unregister():
    print('unregister dilapidator')
    bpy.utils.unregister_module(__name__)


class dilap_run(bpy.types.Operator):
    '''dilap run'''
    # blender will use this as a tooltip for menu items and buttons.

    # unique identifier for buttons and menu items to reference.
    bl_idname = 'object.dilaprun'
    bl_label = 'dilap run'          # display name in the interface.
    bl_options = {'REGISTER','UNDO'} # enable undo for the operator.
    
    # moved assignment from execute() to the body of the class...
    years = bpy.props.IntProperty(name = 'years',default = 10,min = 1,max = 100)
    #dcontext = bpy.props.StringProperty(name = 'context',default = 'lot')
    dcontext = bpy.props.StringProperty(name = 'context',default = 'continent')

    # execute() is called by blender when running the operator.
    def execute(self,context):
        b = vec3(0,0,0).pring(200,8)
        xpj = vec3(1,0,0).prjps(b)
        e = (xpj[1]-xpj[0])/1000.0
        scenegraph = continent(b, e)
        default_materials()
        for m in scenegraph.worldspace():
            build_model(m)
        return {'FINISHED'}


class dilap_purge(bpy.types.Operator):
    '''dilap purge'''
    # blender will use this as a tooltip for menu items and buttons.

    # unique identifier for buttons and menu items to reference.
    bl_idname = 'object.dilappurge'
    bl_label = 'dilap purge'          # display name in the interface.
    bl_options = {'REGISTER','UNDO'} # enable undo for the operator.
    
    # execute() is called by blender when running the operator.
    def execute(self,context):
        global dmodels_bobjs

        #unselect everything
        # <insert code here, this can vary depending on your situation>
        # bpy.ops.object.select_all()
               
        candidate_list = [o.name for o in bpy.data.objects 
                if o.type == "MESH" and o in dmodels_bobjs]
        for object_name in candidate_list:
            bpy.data.objects[object_name].select = True
        # remove all selected.
        bpy.ops.object.delete()
        # remove the meshes, they have no users anymore.
        for item in bpy.data.meshes:
            try:bpy.data.meshes.remove(item) 
            except:print('mesh still in use...')
        return {'FINISHED'}


# create material from texture image
def material_image(name,texture):
    mat = bpy.data.materials.new(name)
    imgpath = dobj.resource_path(texture)[0]
    tex = bpy.data.textures.new(name,type = 'IMAGE')
    tex.image = bpy.data.images.load(imgpath)
    tex.use_alpha = True
    #mat.use_shadeless = True
    mtex = mat.texture_slots.add()
    mtex.texture = tex
    mtex.texture_coords = 'UV'
    mtex.use_map_color_diffuse = True
    return mat


# create material based on colors
def material_solid(name,diffuse,specular,alpha):
    mat = bpy.data.materials.new(name)
    mat.diffuse_color = diffuse
    mat.diffuse_shader = 'LAMBERT'
    mat.diffuse_intensity = 1.0
    mat.specular_color = specular
    mat.specular_shader = 'COOKTORR'
    mat.specular_intensity = 0.5
    mat.alpha = alpha
    mat.ambient = 1
    return mat


materials = {}
# global list of loaded blender materials
def default_materials():
    global materials
    if not materials.keys():
        materials['generic'] = material_image('generic','orangeboxtex.png')
        materials['grass1'] = material_image('grass1','grass1.dds')
        materials['grass2'] = material_image('grass2','grass2.jpg')
        materials['grass3'] = material_image('grass3','grass3.jpg')
        materials['grass4'] = material_image('grass4','grass4.png')
        materials['brick1'] = material_image('brick1','brick1.jpg')
        materials['brick2'] = material_image('brick2','brick2.jpg')
        materials['concrete1'] = material_image('concrete1','concrete1.png')
        materials['concrete2'] = material_image('concrete2','concrete2.png')
        materials['concrete3'] = material_image('concrete3','concrete3.jpg')
    return materials


# put the object into the scene, making it active if make_active
def object_to_scene(obj,make_active = True):
    bpy.context.scene.objects.link(obj)
    if make_active:bpy.context.scene.objects.active = obj
    dmodels_bobjs.append(obj)


# create a blender object from a blender mesh
def object_from_mesh(name,mesh,obj_loc = (0,0,0),mats = None):
    obj = bpy.data.objects.new(name,mesh)
    obj.location = obj_loc
    if not mats is None:
        mats = [materials[mat] for mat in mats]
        [obj.data.materials.append(ma) for ma in mats]
    return obj


# create a blender mesh from model geometry data
def mesh_from_data(name,coords,uvs,faces,face_mats,mats):
    mesh = bpy.data.meshes.new(name)
    if not mats is None:
        [mesh.materials.append(materials[ma]) for ma in mats]
    mesh.vertices.add(len(coords))
    mesh.vertices.foreach_set('co',unpack_list(coords))
    mesh.tessfaces.add(len(faces))
    mesh.tessfaces.foreach_set('vertices_raw',unpack_face_list(faces))
    mesh.tessfaces.foreach_set('material_index',face_mats)
    mesh.tessface_uv_textures.new()
    for fdx in range(len(faces)):
        fa = faces[fdx]
        mesh.tessface_uv_textures[0].data[fdx].uv1 = uvs[fa[0]].to_tuple()
        mesh.tessface_uv_textures[0].data[fdx].uv2 = uvs[fa[1]].to_tuple() 
        mesh.tessface_uv_textures[0].data[fdx].uv3 = uvs[fa[2]].to_tuple() 
    mesh.update()
    return mesh


dmodels_bobjs = []
# build a single model into the blender world
# return the resulting blender space object
#def build_model2(mod,**kwargs):
def build_model(mod):
    oname = mod.filename.replace('.mesh','.000')
    mname = oname+'.'+'mesh'
    ps = mod.pset.ps
    us = mod.uset.ps
    mats = ['generic','concrete1','grass2']
    fs_lookup = {}
    for fmx in range(len(mats)):
        fs_lookup[mats[fmx]] = fmx

    for gfx in mod.gfxmeshes:
        faces = [f for f in gfx.faces if not f is None]
        face_mats = [fs_lookup[gfx.fs_mats[f][1]] for f in faces]
        oloc = (0,0,0)
        mesh = bpy.data.meshes.new(mname)
        if not mats is None:
            [mesh.materials.append(materials[ma]) for ma in mats]
        mesh.vertices.add(len(ps))
        mesh.vertices.foreach_set('co',unpack_list(ps))
        mesh.tessfaces.add(len(faces))
        mesh.tessfaces.foreach_set('vertices_raw',unpack_face_list(faces))
        mesh.tessfaces.foreach_set('material_index',face_mats)
        mesh.tessface_uv_textures.new()
        for fdx in range(len(faces)):
            fa = faces[fdx]
            mesh.tessface_uv_textures[0].data[fdx].uv1 = tuple(us[fa[0]])[:-1]
            mesh.tessface_uv_textures[0].data[fdx].uv2 = tuple(us[fa[1]])[:-1]
            mesh.tessface_uv_textures[0].data[fdx].uv3 = tuple(us[fa[2]])[:-1]
        mesh.update()
    obj = object_from_mesh(oname,mesh,oloc,mats)
    object_to_scene(obj)
    return obj


if __name__ == "__main__":
    register()
