import dilap.core.tools as dpr

from dilap.geometry.vec3 import vec3

import dilap.planning.worldplan as wp

import matplotlib.pyplot as plt

import unittest,numpy,math

import pdb

#python3 -m unittest discover -v ./ "*tests.py"

class test_worldplan(unittest.TestCase):

    def test_factory(self):
        fy = wp.wfactory()
        print('amen')

if __name__ == '__main__':
    unittest.main()








 
 
