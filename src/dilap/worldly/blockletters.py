import dilap.core.base as db

from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat
import dilap.geometry.polymath as pym

import dilap.topology.planargraph as pgr

import dilap.core.plotting as dtl
import matplotlib.pyplot as plt

import pdb

# produce a graph for a particular letter within a square 
#   the letter is fitted to a rectangle of length/width sx/sy
def letter(l,w,sx,sy):
    if w > sx or w > sy:raise ValueError
    sq = vec3(0,0,0).sq(sx-w,sy-w)
    g = pgr.graph()

    if l == 'C':
        i1 = g.av(p = sq[1],l = 0)
        i2,e1 = g.mev(i1,{'p':sq[0],'l':0},{})
        i3,e2 = g.mev(i2,{'p':sq[3],'l':0},{})
        i4,e3 = g.mev(i3,{'p':sq[2],'l':0},{})

    elif l == 'H':
        i1 = g.av(p = sq[1],l = 0)
        i2,e1 = g.mev(i1,{'p':sq[0],'l':0},{})
        i3,e2 = g.mev(i2,{'p':sq[3],'l':0},{})
        i4,e3 = g.mev(i3,{'p':sq[2],'l':0},{})

    return g

# produce a polygon representing letter l, with edges of width w, 
#   the letter is fitted to a rectangle of length/width sx/sy
def block(l,w,sx,sy):
    lg = letter(l,w,sx,sy)
    py = lg.polygon(w/2.0,'ccw')
    return py
    
    



