from dilap.geometry.vec3 import vec3
import dilap.geometry.tools as gtl
import dilap.geometry.polymath as pym

import dilap.worldly.partitiongraph as ptg

import dilap.worldly.world as dwo

import dilap.core.plotting as dtl
import matplotlib.pyplot as plt

import unittest,random

import pdb

#python3 -m unittest discover -v ./ "*tests.py"

###############################################################################

class test_partitiongraph(unittest.TestCase):

    def test_pgraph_split(self):
        # NOTE: THIS FAILS BECAUSE SINTSXYP HAS TOO MUCH ERROR!!
        #fp = vec3(0,0,0).pring(200,8)

        fp = vec3(0,0,0).pring(100,8)
        seq = 'S<0,0.5,0.5,0,1,root>S<0,0.75,0.5,0,1,root>'
        ptg.checkseq(fp,100,seq)
        seq = 'S<0,0.5,0.5,0,1,root>S<0,0.5,0.5,0,1,vertex>'
        ptg.checkseq(fp,100,seq)

    def test_pgraph_insert(self):
        fp = vec3(0,0,0).pring(100,8)
        seq = 'R<0,inside>S<0,0.5,0.5,0,1,root>E<0,1>X<0,1.0,0.5,0.0,vertex>'
        seq = 'S<0,0.5,0.5,0,1,root>I<1,'+seq+'>'
        ptg.checkseq(fp,100,seq)

    def test_pgraph_addedge(self):
        fp = vec3(0,0,0).pring(100,8)
        seq = 'S<0,0.5,0.5,0,1,root>E<0,1>'
        ptg.checkseq(fp,100,seq)

    def test_pgraph_settype(self):
        fp = vec3(0,0,0).pring(100,8)
        seq = 'S<0,0.5,0.5,0,1,root>E<0,1>R<1,inside>'
        ptg.checkseq(fp,100,seq)

    def test_pgraph_addexit(self):
        fp = vec3(0,0,0).pring(100,8)
        seq = 'S<0,0.5,0.5,0,1,root>E<0,1>X<1,0.5,0.5,0.0,vertex>'
        ptg.checkseq(fp,100,seq)

    def test_pgraph_continent(self):
        fp = vec3(0,0,0).pring(100,8)
        #fp = vec3(0,0,0).sq(200,200)
        seq = 'C<0>'
        eg = {
          'C':dwo.continent,
          'M':dwo.metropolis,
          'T':dwo.stitch,
              }
        ptg.checkseq(fp,100,seq,True,grammer = eg)

###############################################################################

if __name__ == '__main__':unittest.main()

###############################################################################






