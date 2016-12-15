from dilap.geometry.vec3 import vec3
import dilap.geometry.tools as gtl

import dilap.core.lsystem as lsy
import dilap.worldly.world as dwo

import dilap.core.plotting as dtl

import matplotlib.pyplot as plt
import itertools as it
import unittest,numpy,math

import pdb

#python3 -m unittest discover -v ./ "*tests.py"

###############################################################################

class test_blg(unittest.TestCase):

    def checkseq(self,fp,seq):
        kws = {
            'footprint':fp,'sequence':seq,
            'floorheight':5,
                }
        blgg = dwo.blg.blggraph(**kws)
        blgg.graph()

        blgg.plotxy(fp)
        plt.show()

    def test_ltapt(self):
        fp = vec3(0,0,0).sq(100,50)
        seq = dwo.blg.bseq.myltapt()
        self.checkseq(fp,seq)

    def test_rtapt(self):
        fp = vec3(0,0,0).sq(100,50)
        seq = dwo.blg.bseq.myrtapt()
        self.checkseq(fp,seq)

    def test_lbapt(self):
        fp = vec3(0,0,0).sq(100,50)
        seq = dwo.blg.bseq.mylbapt()
        self.checkseq(fp,seq)

    def atest_graph(self):
        fp = vec3(0,0,0).sq(100,50)

        seq = dwo.blg.bseq.simplebuilding()

        print('seqq',seq)

        kws = {
            'footprint':fp,'sequence':seq,
            'floorheight':10,
                }

        blgg = dwo.blg.blggraph(**kws)
        blgg.graph()

        blgg.plotxy(fp)
        plt.show()
        blgg.plot(fp)
        plt.show()

###############################################################################

if __name__ == '__main__':unittest.main()

###############################################################################





