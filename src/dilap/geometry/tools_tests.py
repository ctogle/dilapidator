from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat

import dilap.geometry.tools as gtl

import dilap.mesh.tools as dtl
import matplotlib.pyplot as plt

import unittest,numpy

#python3 -m unittest discover -v ./ "*tests.py"

class test_tools(unittest.TestCase):

    def setUp(self):
        self.pslist = [
            vec3(-1,-1,0),vec3( 1,-1,0),
            vec3( 1, 1,0),vec3(-1, 1,0)]
        self.pstupl = (
            vec3(-1,-1,0),vec3( 1,-1,0),
            vec3( 1, 1,0),vec3(-1, 1,0))

    def test_com(self):
        self.assertEqual(gtl.com(self.pslist),vec3(0,0,0))
        self.assertEqual(gtl.com(self.pstupl),vec3(0,0,0))

if __name__ == '__main__':
    unittest.main()










