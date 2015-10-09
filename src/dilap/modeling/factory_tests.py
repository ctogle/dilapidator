import dilap.core.tools as dpr

import dilap.modeling.factory as dfc

import dilap.mesh.tools as dtl
import matplotlib.pyplot as plt

import unittest,numpy,math

import pdb

#python3 -m unittest discover -v ./ "*tests.py"

class test_factory(unittest.TestCase):

    def setUp(self):
        self.f = dfc.factory()

    def test_new(self):
        mcube = self.f.new()

if __name__ == '__main__':
    unittest.main()









 
