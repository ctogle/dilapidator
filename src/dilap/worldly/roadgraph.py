import dilap.core.base as db
import dilap.core.context as cx
import dilap.modeling.model as dmo
import dilap.modeling.factory as dfa

from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat
import dilap.geometry.tools as gtl
import dilap.geometry.polymath as pym

import dilap.worldly.treeskin as ltr
import dilap.worldly.building as blg
import dilap.worldly.blgsequencing as bseq
import dilap.worldly.partitiongraph as ptg

import dilap.core.plotting as dtl
import matplotlib.pyplot as plt

import math,numpy,random,pdb

###############################################################################

def checkseq(pg,fp,seq,show = False):
    print('check-pseq:',seq)
    easement = 2

    def trimtofp(p1,p2,r = 5):
        ips = pym.sintbxyp(p1,p2,fp)
        if len(ips) > 0:
            for ip in ips:
                if ip.isnear(p1):continue
                p2 = p1.lerp(ip,1-r/p1.d(ip))
        return p2

    # place a pair of vertices on the edge of the footprint
    # pointing inward, constituting an exit from the region
    def seed(rg,ss):
        print('SEED',ss)
        exp = fp[-1].lerp(fp[0],0.5)
        exn = vec3(0,0,1).crs(fp[-1].tov(fp[0])).nrm()
        ex1 = exp.cp()
        ex2 = ex1.cp().trn(exn.cp().uscl(easement))
        i1 = rg.av(p = ex1,l = 0)
        i2,r1 = rg.mev(i1,{'p':ex2,'l':0},{})
        return [i2]

    # make an edge and a vertex attached to an existing vertex
    def grow(rg,ss):
        print('GROW',ss)
        tip = nvs[-1]
        tv = rg.vs[tip]
        avs = [rg.vs[k] for k in rg.orings[tip]]
        if len(avs) == 1:
            op = tv[1]['p']
            nptn = avs[0][1]['p'].tov(op).nrm().uscl(100)
            np = op.cp().trn(nptn)
            np = trimtofp(op,np,easement+10)
            nv,nr = rg.mev(tip,{'p':np,'l':0},{})
        return [nv]

    def loop(rg,ss):
        print('LOOP',ss)

    grammer = {
        'S':seed,
        'G':grow,
        'L':loop,
            }

    scnt = len(seq)
    sx = 0
    while sx < scnt:
        c = seq[sx]
        if c in grammer:
            ex = db.seqread(seq,sx)
            nvs = grammer[c](pg,seq[sx+2:ex])
        else:ex = sx+1
        sx = ex
    if show:
        ax = pg.plot()
        ax = dtl.plot_polygon(fp,ax,col = 'r')
        plt.show()
    return pg

###############################################################################





