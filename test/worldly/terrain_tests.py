
from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat
import dilap.geometry.tools as gtl
import dilap.geometry.polymath as pym

import dilap.worldly.roadgraph as rdg

import dilap.core.plotting as dtl
import matplotlib.pyplot as plt

import unittest,random,numpy

import pdb

#python3 -m unittest discover -v ./ "*tests.py"

###############################################################################

class test_roadgraph(unittest.TestCase):

    def test_checkseq(self):

        fp = vec3(0,0,0).pring(100,5)
        seq = 'S<>G<>L<>'
        
        rg = rdg.wgraph().checkseq(fp,seq,True)

###############################################################################

if __name__ == '__main__':unittest.main()

###############################################################################







