from dilap.geometry.quat import quat
from dilap.geometry.vec3 import vec3
import dilap.geometry.polymath as pym

import dilap.geometry.tools as gtl

import dilap.worldly.blockletters as dbl

import dilap.core.plotting as dtl
import matplotlib.pyplot as plt

import unittest,numpy,math

import pdb

#python3 -m unittest discover -v ./ "*tests.py"

class test_blockletters(unittest.TestCase):

    #def setUp(self):

    def plot(self,l):
        ax = dtl.plot_axes_xy(10)
        ax = dtl.plot_polygon_full_xy(l,ax,lw = 2)
        plt.show()

    def test_bnrm(self):
        c = dbl.block('C',1,3,3)
        self.plot(c)

if __name__ == '__main__':
    unittest.main()



