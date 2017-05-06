import dilap.core.base as db
from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat
import dilap.modeling.model as dmo
import dilap.worldly.world as dwo
import dilap.construct as dlc

import dilap.io.fbx as fbx

import dilap.core.plotting as dtl
import matplotlib.pyplot as plt

import unittest,numpy,os,pdb

#python3 -m unittest discover -v ./ "*tests.py"          

class test_fbx(unittest.TestCase):

    example = os.path.join(os.getcwd(),'Rock_Medium_SPR.fbx')
    cubefile  = os.path.join(os.getcwd(),'cube.fbx')
    cubefile2 = os.path.join(os.getcwd(),'cube.2.fbx')
    newcubefile  = os.path.join(os.getcwd(),'newcube.fbx')
    thecubefile  = os.path.join(os.getcwd(),'thecube.fbx')

    def aatest_countstuff(self):
        tree = fbx.read(self.thecubefile)
        model = tree.find('Model')[0]
        vertices = model.find('Vertices')[0][1][1]
        polygonvertexindex = model.find('PolygonVertexIndex')
        layer = tree.find('LayerElementNormal')[0]
        normals = layer.find('Normals')[0][1][1]
        layer = tree.find('LayerElementUV')[0]
        uvs = layer.find('UV')[0][1][1]
        uvindex = layer.find('UVIndex')[0][1][1]

    def test_readwrite(self):
        tree1 = fbx.read(self.cubefile)
        tree1.print()
        objects = tree1.find('Objects')
        models = tree1.find('Model',v = objects[0])
        cube = models[0]
        tree1.write(self.cubefile2)
        tree2 = fbx.read(self.cubefile2)
        #self.assertEqual(tree1.__repr__(),tree2.__repr__())

    def atest_writenewblg(self):
        def test_bseq():
            seq  = ''
            seq += 'L<0>'
            seq += 'L<1>'
            fseq = 'S<0,0.3,0.5,0,0,1,0>' # make living room on left
            fseq += 'S<0,0.3,0.5,0,0,1,0>' # make living room on left
            fseq += 'S<0,0.3,0.5,0,1,0,0>' # make living room on left
            fseq += 'R<2,rtype,closed>'
            fseq += 'C<0,0.2,0.2,0.2,0.2>'
            fseq += 'C<1,0.2,0.2,0.2,0.2>'
            fseq += 'C<3,0.2,0.2,0.2,0.2>'
            fseq += 'C<2,0.2,2,0.2,2>'
            fseq += 'E<1,0>'
            fseq += 'E<1,2>'
            fseq += 'E<2,3>'
            fseq += 'E<2,0>'
            #fseq += 'X<3>'
            seq += 'I<0,'+fseq+'>'
            seq += 'I<1,'+fseq+'>'
            #seq += 'E<1,0>'
            #seq += 'V<2>'
            #seq += 'V<6>'
            seq += 'X<2>'
            return seq
        bfp = vec3(0,0,0).sq(50,15)
        bsq = test_bseq()
        p,q,s = None,None,None
        bfa = dwo.blg.blgfactory()
        cx = bfa.new(p,q,s,footprint = bfp,sequence = bsq,floorheight = 5)
        dlc.realize(cx,io = 'fbx')

    def test_writenewworld(self):
        #tree = fbx.new(cx)
        #tree.write(self.newcubefile)
        teststage(io = 'fbx')

    def atest_writetemplate(self):
        tree = fbx.read(self.cubefile)
        del tree.root.attributes['comment0']
        del tree.root.attributes['comment1']
        del tree.root.attributes['comment2']
        tree.root['CreationTime'] = '"'+db.timestamp()+'"'
        tree.root['Creator'] = '"Dilapidator"'
        
        pdb.set_trace()

def teststage(**kws):
    kws['years'] = 0

    s = 736
    s = 682
    s = 189
    s = 916
    s = 286
    #s = random.randint(0,1000)

    fkws = {
        'seed' : s,
            }
    wkws = {
        'boundary' : vec3(0,0,0).pring(500,8),
        'landmasses' : [vec3(0,0,0).pring(250,8)],
            }
    cx = dwo.worldfactory(**fkws).new(**wkws)
    dlc.realize(cx,**kws)

if __name__ == '__main__':
    unittest.main()
