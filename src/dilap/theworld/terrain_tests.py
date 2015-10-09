import dilap.core.tools as dpr

from dilap.geometry.vec3 import vec3

import dilap.modeling.terrain as dmt

import dilap.mesh.tools as dtl
import matplotlib.pyplot as plt

import unittest,numpy,math

import pdb

#python3 -m unittest discover -v ./ "*tests.py"

class test_terrain(unittest.TestCase):

    def setUp(self):pass

    def test_init(self):
        tm = dmt.terrain()

        pdb.set_trace()

if __name__ == '__main__':
    unittest.main()








 


