from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat
import dilap.geometry.polymath as pym

import dilap.topology.planargraph as pgr

import dilap.core.plotting as dtl
import matplotlib.pyplot as plt

import pdb

def C(sq,g):
    i1 = g.av(p = sq[1],l = 0)
    i2,e1 = g.mev(i1,{'p':sq[0],'l':0},{})
    i3,e2 = g.mev(i2,{'p':sq[3],'l':0},{})
    i4,e3 = g.mev(i3,{'p':sq[2],'l':0},{})

def H(sq,g):
    sq = pym.bisectb(sq)
    i1 = g.av(p = sq[0],l = 0)
    i2 = g.av(p = sq[7],l = 0)
    i3 = g.av(p = sq[6],l = 0)
    i4 = g.av(p = sq[2],l = 0)
    i5 = g.av(p = sq[3],l = 0)
    i6 = g.av(p = sq[4],l = 0)
    e1 = g.ae(i1,i2,**{})
    e2 = g.ae(i2,i3,**{})
    e3 = g.ae(i4,i5,**{})
    e4 = g.ae(i5,i6,**{})
    e5 = g.ae(i2,i5,**{})

def I(sq,g):
    sq = pym.bisectb(sq)
    i1 = g.av(p = sq[0],l = 0)
    i2 = g.av(p = sq[1],l = 0)
    i3 = g.av(p = sq[2],l = 0)
    i4 = g.av(p = sq[6],l = 0)
    i5 = g.av(p = sq[5],l = 0)
    i6 = g.av(p = sq[4],l = 0)
    e1 = g.ae(i1,i2,**{})
    e2 = g.ae(i2,i3,**{})
    e3 = g.ae(i4,i5,**{})
    e4 = g.ae(i5,i6,**{})
    e5 = g.ae(i2,i5,**{})

# produce a graph for a particular letter within a square 
#   the letter is fitted to a rectangle of length/width sx/sy
def letter(l,w,sx,sy):
    if w > sx or w > sy:raise ValueError
    sq = vec3(0,0,0).sq(sx-w,sy-w)
    g = pgr.planargraph()
    if   l == 'C':C(sq,g)
    elif l == 'H':H(sq,g)
    elif l == 'I':I(sq,g)
    return g

# produce a polygon representing letter l, with edges of width w, 
#   the letter is fitted to a rectangle of length/width sx/sy
def block(l,w,sx,sy):
    lg = letter(l,w,sx,sy)
    py = pym.pgtopy(lg,w/2.0)
    return py
    
    



