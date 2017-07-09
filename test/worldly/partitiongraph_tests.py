from dilap.geometry import vec3
from dilap.worldly import partition
import matplotlib.pyplot as plt
import unittest
import pdb

#python3 -m unittest discover -v ./ "*tests.py"

###############################################################################

class test_partitiongraph(unittest.TestCase):


    def test_split(self):
        sq = vec3(0,0,0).sq(256,256)
        peg = vec3(128,0,0).pring(64,8)
        p = partition(loop = sq,holes = [],style = 'A')
        l,r = p.split(p.root,peg)
        p.plot()
        plt.show()

        sq = vec3(0,0,0).sq(256,256)
        peg = vec3(0,0,0).pring(64,8)
        p = partition(loop = sq,holes = [],style = 'A')
        l,r = p.split(p.root,peg)
        p.plot()
        plt.show()

        sq = vec3(0,0,0).sq(256,256)
        peg = vec3(0,0,0).pring(300,8)
        p = partition(loop = sq,holes = [],style = 'A')
        l,r = p.split(p.root,peg)
        p.plot()
        plt.show()

        sq = vec3(0,0,0).sq(256,256)
        peg = vec3(0,0,0).pring(168,8)
        p = partition(loop = sq,holes = [],style = 'A')
        l,r = p.split(p.root,peg)
        p.plot()
        plt.show()

        sq = vec3(0,0,0).sq(256,256)
        peg = vec3(128,0,0).sq(256,256)
        p = partition(loop = sq,holes = [],style = 'A')
        l,r = p.split(p.root,peg)
        p.plot()
        plt.show()

        sq = vec3(0,0,0).sq(256,256)
        peg = vec3(256,0,0).sq(256,256)
        p = partition(loop = sq,holes = [],style = 'A')
        l,r = p.split(p.root,peg)
        p.plot()
        plt.show()

###############################################################################

if __name__ == '__main__':unittest.main()

###############################################################################






