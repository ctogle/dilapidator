#import dilap.core.tools as dpr
import dilap.geometry.tools as dpr

import dilap.topology.tree as dtr

#import dilap.mesh.tools as dtl
import dilap.core.plotting as dtl
import matplotlib.pyplot as plt

import unittest,numpy,math

import pdb

#python3 -m unittest discover -v ./ "*tests.py"

class test_tree(unittest.TestCase):

    def makeex(self):
        self.one = self.tr.avert(self.tr.root)
        self.two = self.tr.avert(self.tr.root)
        self.three = self.tr.avert(self.two)
        self.four = self.tr.avert(self.two)
        self.five = self.tr.avert(self.one)
        self.six = self.tr.avert(self.three)

    def setUp(self):
        self.tr = dtr.tree()

    def test_mask0(self):
        self.makeex()
        self.assertFalse(self.tr.root in self.tr.mask(0,self.tr.root,None))
        self.assertTrue(self.one in self.tr.mask(0,self.tr.root,None))
        self.assertTrue(self.two in self.tr.mask(0,self.tr.root,None))
        self.assertFalse(self.three in self.tr.mask(0,self.tr.root,None))
        self.assertFalse(self.four in self.tr.mask(0,self.tr.root,None))
        self.assertFalse(self.five in self.tr.mask(0,self.tr.root,None))
        self.assertFalse(self.six in self.tr.mask(0,self.tr.root,None))

        self.assertTrue(self.tr.root in self.tr.mask(0,self.one,None))
        self.assertFalse(self.one in self.tr.mask(0,self.one,None))
        self.assertFalse(self.two in self.tr.mask(0,self.one,None))
        self.assertFalse(self.three in self.tr.mask(0,self.one,None))
        self.assertFalse(self.four in self.tr.mask(0,self.one,None))
        self.assertTrue(self.five in self.tr.mask(0,self.one,None))
        self.assertFalse(self.six in self.tr.mask(0,self.one,None))

        self.assertTrue(self.tr.root in self.tr.mask(0,self.two,None))
        self.assertFalse(self.one in self.tr.mask(0,self.two,None))
        self.assertFalse(self.two in self.tr.mask(0,self.two,None))
        self.assertTrue(self.three in self.tr.mask(0,self.two,None))
        self.assertTrue(self.four in self.tr.mask(0,self.two,None))
        self.assertFalse(self.five in self.tr.mask(0,self.two,None))
        self.assertFalse(self.six in self.tr.mask(0,self.two,None))

        self.assertFalse(self.tr.root in self.tr.mask(0,self.three,None))
        self.assertFalse(self.one in self.tr.mask(0,self.three,None))
        self.assertTrue(self.two in self.tr.mask(0,self.three,None))
        self.assertFalse(self.three in self.tr.mask(0,self.three,None))
        self.assertFalse(self.four in self.tr.mask(0,self.three,None))
        self.assertFalse(self.five in self.tr.mask(0,self.three,None))
        self.assertTrue(self.six in self.tr.mask(0,self.three,None))

        self.assertFalse(self.tr.root in self.tr.mask(0,self.four,None))
        self.assertFalse(self.one in self.tr.mask(0,self.four,None))
        self.assertTrue(self.two in self.tr.mask(0,self.four,None))
        self.assertFalse(self.three in self.tr.mask(0,self.four,None))
        self.assertFalse(self.four in self.tr.mask(0,self.four,None))
        self.assertFalse(self.five in self.tr.mask(0,self.four,None))
        self.assertFalse(self.six in self.tr.mask(0,self.four,None))

        self.assertFalse(self.tr.root in self.tr.mask(0,self.five,None))
        self.assertTrue(self.one in self.tr.mask(0,self.five,None))
        self.assertFalse(self.two in self.tr.mask(0,self.five,None))
        self.assertFalse(self.three in self.tr.mask(0,self.five,None))
        self.assertFalse(self.four in self.tr.mask(0,self.five,None))
        self.assertFalse(self.five in self.tr.mask(0,self.five,None))
        self.assertFalse(self.six in self.tr.mask(0,self.five,None))

        self.assertFalse(self.tr.root in self.tr.mask(0,self.six,None))
        self.assertFalse(self.one in self.tr.mask(0,self.six,None))
        self.assertFalse(self.two in self.tr.mask(0,self.six,None))
        self.assertTrue(self.three in self.tr.mask(0,self.six,None))
        self.assertFalse(self.four in self.tr.mask(0,self.six,None))
        self.assertFalse(self.five in self.tr.mask(0,self.six,None))
        self.assertFalse(self.six in self.tr.mask(0,self.six,None))

    def test_mask1(self):
        self.makeex()
        rootes = self.tr.mask(1,self.tr.root,None)
        self.assertTrue(self.tr.edges[0] in rootes)
        self.assertTrue(self.tr.edges[1] in rootes)
        for x in range(2,len(self.tr.edges)):
            self.assertFalse(self.tr.edges[x] in rootes)
        rootes = self.tr.mask(1,self.six,None)
        self.assertTrue(self.tr.edges[-1] in rootes)
        rootes = self.tr.mask(1,self.three,None)
        self.assertTrue(self.tr.edges[-1] in rootes)

    def test_ebove(self):
        tr = dtr.tree()
        one = tr.avert(tr.root)
        two = tr.avert(tr.root)
        three = tr.avert(two)
        re = tr.ebove(three)
        self.assertTrue(re.one is two)
        self.assertTrue(re.two is three)
        re = tr.ebove(one)
        self.assertTrue(re.one is tr.root)
        self.assertTrue(re.two is one)

    def test_above(self):
        self.makeex()
        self.assertTrue(self.tr.above(self.one) is self.tr.root)
        self.assertTrue(self.tr.above(self.two) is self.tr.root)
        self.assertTrue(self.tr.above(self.three) is self.two)
        self.assertTrue(self.tr.above(self.four) is self.two)
        self.assertTrue(self.tr.above(self.five) is self.one)
        self.assertTrue(self.tr.above(self.six) is self.three)

    def test_below(self):
        self.makeex()
        self.assertTrue(self.one in self.tr.below(self.tr.root))
        self.assertTrue(self.two in self.tr.below(self.tr.root))
        self.assertTrue(self.three in self.tr.below(self.two))
        self.assertTrue(self.four in self.tr.below(self.two))
        self.assertTrue(self.five in self.tr.below(self.one))
        self.assertTrue(self.six in self.tr.below(self.three))

    def test_aroot(self):
        self.assertEqual(self.tr.above(self.tr.root),None)
        self.assertEqual(self.tr.below(self.tr.root),[])
        self.assertEqual(self.tr.vertcount,1)
        self.assertEqual(self.tr.edgecount,0)

    def test_avert(self):
        tr = dtr.tree()
        one = tr.avert(tr.root)
        two = tr.avert(tr.root)
        three = tr.avert(two)
        self.assertEqual(tr.below(tr.root),[one,two])
        self.assertEqual(tr.above(one),tr.root)
        self.assertEqual(tr.above(two),tr.root)
        self.assertEqual(tr.below(one),[])
        self.assertEqual(tr.below(two),[three])
        self.assertEqual(tr.above(three),two)

    def test_rvert(self):
        tr = dtr.tree()
        one = tr.avert(tr.root)
        two = tr.avert(tr.root)
        three = tr.avert(two)
        tr.rvert(two)
        self.assertFalse(two in tr.mask(0,tr.root,None))
        self.assertEqual(tr.vcnt(),2)
        self.assertEqual(tr.ecnt(),1)
        self.assertTrue(2 in tr.vistack)
        self.assertTrue(3 in tr.vistack)

    def test_aedge(self):
        tr = dtr.tree()
        one = tr.avert(tr.root)
        two = tr.avert(tr.root)
        three = tr.avert(two)
        for e in tr.edges:
            if e.one is tr.root:
                self.assertTrue(e.two is one or e.two is two)
            elif e.one is two:
                self.assertTrue(e.two is three)
            else:self.assertTrue(False)

    def test_redge(self):
        tr = dtr.tree()
        one = tr.avert(tr.root)
        two = tr.avert(tr.root)
        three = tr.avert(two)
        re = tr.ebove(three)
        tr.redge(re)
        self.assertTrue(not three in tr.verts)

    def test_sedge(self):
        tr = dtr.tree()
        one = tr.avert(tr.root)
        two = tr.avert(tr.root)
        three = tr.avert(two)
        self.assertTrue(tr.above(three) is two)
        four = tr.sedge(tr.ebove(three))
        self.assertTrue(tr.above(three) is four)
        self.assertTrue(tr.above(four) is two)

    def test_medge(self):
        tr = dtr.tree()
        one = tr.avert(tr.root)
        two = tr.avert(tr.root)
        three = tr.avert(two)
        self.assertTrue(tr.above(three) is two)
        tr.medge(two)
        self.assertTrue(tr.above(three) is tr.root)
        self.assertTrue(three in tr.below(tr.root))

if __name__ == '__main__':
    unittest.main()









 

