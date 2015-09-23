import dilap.core.tools as dpr

import dilap.topology.tree as dtr

import dilap.mesh.tools as dtl
import matplotlib.pyplot as plt

import unittest,numpy,math

import pdb

#python3 -m unittest discover -v ./ "*tests.py"

class test_tree(unittest.TestCase):

    def test_init(self):
        tr = dtr.tree()
        self.assertEqual(tr.above(tr.root),None)
        self.assertEqual(tr.below(tr.root),[])
        self.assertEqual(tr.vertcount,1)
        self.assertEqual(tr.edgecount,0)

    #def test_mask(self):
    #def test_aroot(self):

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

    #def test_aedge(self):
    #def test_above(self):
    #def test_below(self):

if __name__ == '__main__':
    unittest.main()









 

