import dilap.topology.tree as dtr
import dilap.topology.vert as dvt
import dilap.geometry.tools as gtl
import dilap.core.base as db
import collections,re,os,pdb
from io import StringIO as sio

class fbxvert(dvt.vert):

    def attribute(self,l):
        if l.startswith(';'):
            self.attributes.append((None,l))
            self.commentcount += 1
        elif ':' in l:
            if l.endswith(':'):
                k,v = l[:-1].strip(),''
            else:
                s = l.find(':')
                k,v = l[:s].strip(),l[s+1:].strip()
            self.attributes.append((k,v))
        else:
            self.attributes.append((None,None))
            self.newlinecount += 1

    def identify(self,l):
        if l.endswith('{'):l = l[:-1].strip()
        if l.endswith(':'):
            self.name = l[:-1].strip()
            self.properties = ''
        else:
            s = l.find(':')
            self.name = l[:s].strip()
            self.properties = l[s+1:].strip()

    def findeach(self,name):
        for j,a in enumerate(self.attributes):
            if a[0] == name:
                yield j,a

    def find(self,name):
        each = [e for e in self.findeach(name)]
        return each

    def __str__(self):
        s = ':'.join([self.name,self.properties])
        return s

    def __init__(self,ix,tag):
        dvt.vert.__init__(self,ix)
        self.identify(tag)
        self.attributes = []
        self.commentcount = 0
        self.newlinecount = 0

class fbxtree(dtr.tree):

    def print(self,v = None,depth = 0):
        if v is None:
            print('\n'+'-'*50)
            v = self.root
        strv = str(v)
        print(strv)
        depth += 1
        for c in self.below(v):
            print(' '*6*depth+'|'+('='*(len(strv)-5))+'> ',end = '')
            self.print(c,depth)

    def ascii(self,v = None,depth = -1):
        if v is None:v = self.root
        depth += 1
        below = self.below(v)
        for a in v.attributes:
            ak,av = a
            if ak is None:
                if av is None:
                    yield ''
                else:
                    yield '%s%s' % ('\t'*depth,av)
            elif ak.startswith('node'):
                c = below[int(ak.replace('node',''))]
                yield '%s%s: %s  {' % ('\t'*depth,c.name,c.properties)
                for l in self.ascii(c,depth):yield l
                yield '%s}' % ('\t'*depth,)
            else:
                yield '%s%s: %s' % ('\t'*depth,ak,av)

    def __repr__(self):
        s = '\n'.join([l for l in self.ascii()])
        return s

    def __init__(self):
        self.vertclass = fbxvert
        dtr.tree.__init__(self,'ROOT:')

    def findeach(self,name,v = None):
        if v is None:v = self.root
        for c in self.below(v):
            for c in self.findeach(name,c):
                if c.name == name:
                    yield c
        if v.name == name:
            yield v

    def find(self,name,v = None):
        each = [e for e in self.findeach(name,v)]
        return each

    attribute = re.compile('^(.+:.*[^{]?|;.*)$')
    continued = re.compile('^[^:{;]+$')
    nodepush = re.compile('^.+:.*{$')
    nodepop = re.compile('^}$')

    def write(self,path):
        with open(path,'w') as f:
            for l in self.ascii():
                f.write(l+os.linesep)

    @classmethod
    def readlines(cls,lineiter):
        tree = cls()
        cursor = tree.root
        for l in lineiter:        
            s = l.strip()
            if cls.nodepush.match(s):
                nodename = 'node%d:' % len(tree.below(cursor))
                cursor.attribute(nodename)
                cursor = tree.avert(cursor,s)
            elif cls.nodepop.match(s):cursor = tree.above(cursor)
            elif cls.attribute.match(s):
                cursor.attribute(s)
            elif cls.continued.match(s):
                a = cursor.attributes[-1]
                cursor.attributes[-1] = (a[0],a[1]+s)
            elif s:
                print('unknown line format: \'%s\'' % s)
                pdb.set_trace()
                raise NotImplementedError('unknown line format: \'%s\'' % s)
            else:cursor.attribute(s)
        return tree

    @classmethod
    def read(cls,path):
        if os.path.exists(path):
            with open(path,'r') as f:
                tree = cls.readlines(f)
        elif os.linesep in path:
            tree = cls.readlines((l.strip() for l in path.split(os.linesep)))
        return tree

    @classmethod
    def new(cls,sg,**kws):
        modstrings = []
        posestring = []
        modelrelationstring = []
        objectconnectionsstring = []

        def graphvert(vrt,ptf,**kws):
            wtform = vrt.tform.true(ptf)
            print('vrt!\n',vrt.ix,wtform,'\nfrom\n',vrt.tform,'\nfrom\n',ptf)
            wtform.scl.uscl(100)
            s,r,t = wtform.scl,wtform.rot,wtform.pos
            for m in vrt.models:
                m.uset.uscl(100)
                for gm in m.gfxmeshes:
                    mname = uniquefilename(m,gm)
                    modstrings.append(mstring % mstringargs(mname,m,gm,s,r,t))
                    posestring.append(posenode % mname)
                    modelrelationstring.append(objectrelationmesh % mname)
                    objectconnectionsstring.append(
                        objectconnectionstring % (mname,mname))
                m.uset.uscl(0.01)
            wtform.scl.uscl(0.01)
            for x in sg.below(vrt):
                graphvert(x,wtform,**kws)
        graphvert(sg.root,None,**kws)

        dt = db.nowdt()
        nmodels = len(modstrings)
        targs = (
            dt.year,dt.month,dt.day,
            dt.hour,dt.minute,dt.second,
            db.timestamp(dt = dt),
            definitions % (nmodels,nmodels,nmodels),
            objectproperties(modstrings,posestring),
            objectrelations(modelrelationstring),
            '\n'.join(objectconnectionsstring),
                )
        tree = cls.read(template % targs)
        return tree

def read(path,*ags,**kws):return fbxtree.read(path,*ags,**kws)
def new(*ags,**kws):return fbxtree.new(*ags,**kws)

worlddir = None
modellookup = {}
def uniquefilename(m,gm):
    name = '%s.%d.%d'
    if m.filename in modellookup:
        modellookup[m.filename] += 1
    else:modellookup[m.filename] = 0
    n = modellookup[m.filename]
    name = name % (m.filename,n,m.gfxmeshes.index(gm))
    return name

def mstringargs(mname,m,gm,s,r,t):
    # i expect issues if vertices are missing (any were deleted)
    # and thus vertex indices wont naturally line up
    #m.normals(gm)
    m.uvs(gm)
    vertices = sio()
    for p in m.gvps_i(gm,range(gm.vertcount)):
        vertices.write('%.6f,%.6f,%.6f,' % (p.x,p.y,p.z))
    vertexcnt = len(vertices.getvalue())
    if vertexcnt > 0:vertices.truncate(vertexcnt-1)
    uvs = sio()
    for u in m.gvus_i(gm,range(gm.vertcount)):
        uvs.write('%.6f,%.6f,' % (u.x,u.y))
    uvcnt = len(uvs.getvalue())
    if uvcnt > 0:uvs.truncate(uvcnt-1)
    #uvs.truncate(len(uvs.getvalue())-1)
    polygonvertexindex = sio()
    normals = sio()
    uvindex = sio()
    for f in gm.faces:
        if f is None:continue
        if m.isneedle(gm,f):continue
        v1,v2,v3 = f
        polygonvertexindex.write('%d,%d,%d,' % (v1,v2,-(v3+1)))
        u1,u2,u3 = gm.verts[v1][2],gm.verts[v2][2],gm.verts[v3][2]
        uvindex.write('%d,%d,%d,' % (u1,u2,u3))
        for n in m.facenormals(gm,f):
            normals.write('%.6f,%.6f,%.6f,' % (n.x,n.y,n.z))
    polygonvertexindexcnt = len(polygonvertexindex.getvalue())
    if polygonvertexindexcnt > 0:polygonvertexindex.truncate(polygonvertexindexcnt-1)
    #polygonvertexindex.truncate(len(polygonvertexindex.getvalue())-1)
    normalcnt = len(normals.getvalue())
    if normalcnt > 0:normals.truncate(normalcnt-1)
    #normals.truncate(len(normals.getvalue())-1)
    uvindexcnt = len(uvindex.getvalue())
    if uvindexcnt > 0:uvindex.truncate(uvindexcnt-1)
    #uvindex.truncate(len(uvindex.getvalue())-1)
    a = (mname,
        '%.12f,%.12f,%.12f' % (t.x,t.y,t.z),
        '%.12f,%.12f,%.12f' % (-90,0,0),
        '%.12f,%.12f,%.12f' % (s.x,s.y,s.z),
        vertices.getvalue(),polygonvertexindex.getvalue(),
        normals.getvalue(),uvs.getvalue(),uvindex.getvalue())
    return a

template = '''\
; FBX 6.1.0 project file
; Created by dilapidator FBX exporter
; for support mail: cogle@vt.edu
; --------------------------------------------------

FBXHeaderExtension :  {
  FBXHeaderVersion : 1003
  FBXVersion : 6100
  CreationTimeStamp :  {
    Version : 1000
    Year : %d
    Month : %02d
    Day : %02d
    Hour : %d
    Minute : %d
    Second : %d
    Millisecond : 0
  }
  Creator : "dilapidator"
  Otherflags :  {
    FlagPLE : 0
  }
}
CreationTime : "%s"
Creator : "dilapidator"

; Object definitions
;------------------------------------------------------------------

Definitions :  {
  %s
}

; Object properties
;------------------------------------------------------------------

Objects :  {
  %s
}

; Object relations
;------------------------------------------------------------------

Relations :  {
  %s
}

; Object connections
;------------------------------------------------------------------

Connections :  {
  %s
}

; Takes and animation section
;----------------------------------------------------

Takes :  {
  Current : ""
}

; Version 5 settings
;------------------------------------------------------------------

Version5 :  {
  AmbientRenderSettings :  {
    Version : 101
    AmbientLightColor : 0.0,0.0,0.0,0
  }
  FogOptions :  {
    FogEnable : 0
    FogMode : 0
    FogDensity : 0.000
    FogStart : 5.000
    FogEnd : 25.000
    FogColor : 0.1,0.1,0.1,1
  }
  Settings :  {
    FrameRate : "24"
    TimeFormat : 1
    SnapOnFrames : 0
    ReferenceTimeIndex : -1
    TimeLineStartTime : 0
    TimeLineStopTime : 479181389250
  }
  RendererSetting :  {
    DefaultCamera : "Producer Perspective"
    DefaultViewingMode : 0
  }
}'''

definitions = '''\
  Version : 100
  Count : 3
  ObjectType : "Model" {
    Count : %d
  }
  ObjectType : "Geometry" {
    Count : %d
  }
  ObjectType : "Material" {
    Count : %d
  }
  ObjectType : "Pose" {
    Count : 1
  }
  ObjectType : "GlobalSettings" {
    Count : 1
  }'''
  
def objectproperties(modstrings,posestring): 
    matstrings = [material]
    posestring = [pose % (len(modstrings),'\n'.join(posestring))]
    settstring = [globalsettings]
    return '\n'.join(modstrings+matstrings+posestring+settstring)

objectrelationmesh = '''\
  Model : "Model::%s", "Mesh" {
  }'''
def objectrelations(modelrelationstring): 
    modelrelationstring = '\n'.join(modelrelationstring)
    relations = '''\
  %s
  Model : "Model::Producer Perspective", "Camera" {
  }
  Model : "Model::Producer Top", "Camera" {
  }
  Model : "Model::Producer Bottom", "Camera" {
  }
  Model : "Model::Producer Front", "Camera" {
  }
  Model : "Model::Producer Back", "Camera" {
  }
  Model : "Model::Producer Right", "Camera" {
  }
  Model : "Model::Producer Left", "Camera" {
  }
  Model : "Model::Camera Switcher", "CameraSwitcher" {
  }
  Material : "Material::Material", "" {
  }''' % modelrelationstring
    return relations

objectconnectionstring = '''\
    Connect : "OO", "Model::%s", "Model::Scene"
    Connect : "OO", "Material::Material", "Model::%s"'''

material = '''\
  Material: "Material::Material", "" {
    Version: 102
    ShadingModel: "lambert"
    MultiLayer: 0
    Properties60:  {
      Property: "ShadingModel", "KString", "", "Lambert"
      Property: "MultiLayer", "bool", "",0
      Property: "EmissiveColor", "ColorRGB", "",0.8000,0.8000,0.8000
      Property: "EmissiveFactor", "double", "",0.0000
      Property: "AmbientColor", "ColorRGB", "",1.0000,1.0000,1.0000
      Property: "AmbientFactor", "double", "",1.0000
      Property: "DiffuseColor", "ColorRGB", "",0.8000,0.8000,0.8000
      Property: "DiffuseFactor", "double", "",0.8000
      Property: "Bump", "Vector3D", "",0,0,0
      Property: "TransparentColor", "ColorRGB", "",1,1,1
      Property: "TransparencyFactor", "double", "",0.0000
      Property: "SpecularColor", "ColorRGB", "",1.0000,1.0000,1.0000
      Property: "SpecularFactor", "double", "",0.5000
      Property: "ShininessExponent", "double", "",12.3
      Property: "ReflectionColor", "ColorRGB", "",0,0,0
      Property: "ReflectionFactor", "double", "",1
      Property: "Emissive", "ColorRGB", "",0,0,0
      Property: "Ambient", "ColorRGB", "",1.0,1.0,1.0
      Property: "Diffuse", "ColorRGB", "",0.8,0.8,0.8
      Property: "Specular", "ColorRGB", "",1.0,1.0,1.0
      Property: "Shininess", "double", "",12.3
      Property: "Opacity", "double", "",1.0
      Property: "Reflectivity", "double", "",0
    }
  }'''

pose = '''\
  Pose: "Pose::BIND_POSES", "BindPose" {
    Type: "BindPose"
    Version: 100
    Properties60:  {
    }
    NbPoseNodes: %d
    %s
  }'''

posenode = '''\
    PoseNode:  {
      Node: "Model::%s"
      Matrix: 0.000000075497901,0.000000000000000,-1.000000000000000,0.000000000000000,-1.000000000000000,0.000000000000000,-0.000000075497901,0.000000000000000,0.000000000000000,1.000000000000000,0.000000000000000,0.000000000000000,0.000000000000000,0.000000000000000,0.000000000000000,1.000000000000000
    }'''

globalsettings = '''\
  GlobalSettings:  {
    Version: 1000
    Properties60:  {
      Property: "UpAxis", "int", "",1
      Property: "UpAxisSign", "int", "",1
      Property: "FrontAxis", "int", "",2
      Property: "FrontAxisSign", "int", "",1
      Property: "CoordAxis", "int", "",0
      Property: "CoordAxisSign", "int", "",1
      Property: "UnitScaleFactor", "double", "",1
    }
  }'''

mstring = '''\
  Model : "Model::%s", "Mesh" {
    Version : 232
    Properties60 :  {
      Property : "QuaternionInterpolate", "bool", "",0
      Property : "Visibility", "Visibility", "A+",1
      Property : "Lcl Translation", "Lcl Translation", "A+",%s
      Property : "Lcl Rotation", "Lcl Rotation", "A+",%s
      Property : "Lcl Scaling", "Lcl Scaling", "A+",%s
      Property : "RotationOffset", "Vector3D", "",0,0,0
      Property : "RotationPivot", "Vector3D", "",0,0,0
      Property : "ScalingOffset", "Vector3D", "",0,0,0
      Property : "ScalingPivot", "Vector3D", "",0,0,0
      Property : "TranslationActive", "bool", "",0
      Property : "TranslationMin", "Vector3D", "",0,0,0
      Property : "TranslationMax", "Vector3D", "",0,0,0
      Property : "TranslationMinX", "bool", "",0
      Property : "TranslationMinY", "bool", "",0
      Property : "TranslationMinZ", "bool", "",0
      Property : "TranslationMaxX", "bool", "",0
      Property : "TranslationMaxY", "bool", "",0
      Property : "TranslationMaxZ", "bool", "",0
      Property : "RotationOrder", "enum", "",0
      Property : "RotationSpaceForLimitOnly", "bool", "",0
      Property : "AxisLen", "double", "",10
      Property : "PreRotation", "Vector3D", "",0,0,0
      Property : "PostRotation", "Vector3D", "",0,0,0
      Property : "RotationActive", "bool", "",0
      Property : "RotationMin", "Vector3D", "",0,0,0
      Property : "RotationMax", "Vector3D", "",0,0,0
      Property : "RotationMinX", "bool", "",0
      Property : "RotationMinY", "bool", "",0
      Property : "RotationMinZ", "bool", "",0
      Property : "RotationMaxX", "bool", "",0
      Property : "RotationMaxY", "bool", "",0
      Property : "RotationMaxZ", "bool", "",0
      Property : "RotationStiffnessX", "double", "",0
      Property : "RotationStiffnessY", "double", "",0
      Property : "RotationStiffnessZ", "double", "",0
      Property : "MinDampRangeX", "double", "",0
      Property : "MinDampRangeY", "double", "",0
      Property : "MinDampRangeZ", "double", "",0
      Property : "MaxDampRangeX", "double", "",0
      Property : "MaxDampRangeY", "double", "",0
      Property : "MaxDampRangeZ", "double", "",0
      Property : "MinDampStrengthX", "double", "",0
      Property : "MinDampStrengthY", "double", "",0
      Property : "MinDampStrengthZ", "double", "",0
      Property : "MaxDampStrengthX", "double", "",0
      Property : "MaxDampStrengthY", "double", "",0
      Property : "MaxDampStrengthZ", "double", "",0
      Property : "PreferedAngleX", "double", "",0
      Property : "PreferedAngleY", "double", "",0
      Property : "PreferedAngleZ", "double", "",0
      Property : "InheritType", "enum", "",0
      Property : "ScalingActive", "bool", "",0
      Property : "ScalingMin", "Vector3D", "",1,1,1
      Property : "ScalingMax", "Vector3D", "",1,1,1
      Property : "ScalingMinX", "bool", "",0
      Property : "ScalingMinY", "bool", "",0
      Property : "ScalingMinZ", "bool", "",0
      Property : "ScalingMaxX", "bool", "",0
      Property : "ScalingMaxY", "bool", "",0
      Property : "ScalingMaxZ", "bool", "",0
      Property : "GeometricTranslation", "Vector3D", "",0,0,0
      Property : "GeometricRotation", "Vector3D", "",0,0,0
      Property : "GeometricScaling", "Vector3D", "",1,1,1
      Property : "LookAtProperty", "object", ""
      Property : "UpVectorProperty", "object", ""
      Property : "Show", "bool", "",1
      Property : "NegativePercentShapeSupport", "bool", "",1
      Property : "DefaultAttributeIndex", "int", "",0
      Property : "Color", "Color", "A",0.8,0.8,0.8
      Property : "Size", "double", "",100
      Property : "Look", "enum", "",1
    }
    MultiLayer : 0
    MultiTake : 1
    Shading : Y
    Culling : "CullingOff"
    Vertices : %s
    PolygonVertexIndex : %s
    GeometryVersion : 124
    LayerElementNormal : 0 {
      Version : 101
      Name : ""
      MappingInformationType : "ByPolygonVertex"
      ReferenceInformationType : "Direct"
      Normals : %s
    }
    LayerElementSmoothing : 0 {
      Version : 102
      Name : ""
      MappingInformationType : "ByPolygon"
      ReferenceInformationType : "Direct"
      Smoothing : 0,0,0,0,0,0
    }
		LayerElementUV: 0 {
			Version: 101
			Name: "UVMap"
			MappingInformationType: "ByPolygonVertex"
			ReferenceInformationType: "IndexToDirect"
			UV: %s
			UVIndex: %s
		}
		LayerElementTexture: 0 {
			Version: 101
			Name: ""
			MappingInformationType: "NoMappingInformation"
			ReferenceInformationType: "IndexToDirect"
			BlendMode: "Translucent"
			TextureAlpha: 1
			TextureId: 
		}
    LayerElementMaterial : 0 {
      Version : 101
      Name : ""
      MappingInformationType : "AllSame"
      ReferenceInformationType : "IndexToDirect"
      Materials : 0
    }
    Layer : 0 {
      Version : 100
      LayerElement :  {
        Type : "LayerElementNormal"
        TypedIndex : 0
      }
      LayerElement :  {
        Type : "LayerElementSmoothing"
        TypedIndex : 0
      }
			LayerElement:  {
				Type: "LayerElementUV"
				TypedIndex: 0
			}
      LayerElement :  {
        Type : "LayerElementMaterial"
        TypedIndex : 0
      }
    }
  }'''

