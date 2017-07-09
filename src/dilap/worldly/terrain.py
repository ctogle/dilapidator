from dilap.core import *
from dilap.geometry import *
import dilap.geometry.polymath as pym
import pdb


class terrain(topography):


    '''Nesting the loops with respect to containment and offsetting them 
    in the z direction can represent a topographical map for terrain'''


    def raise_earth(self,tips,e,mingrad = 1.0,mindelz = -5.0,
                    depth = 0,maxdepth = None):
        print('... raising earth ... (depth: %d)' % depth)
        if maxdepth and depth == maxdepth:
            print('... reached max depth (%d)! ...' % maxdepth)
        newtips = []
        for tip in tips:
            dr = abs(mindelz/mingrad)
            newloop = [p.cp().ztrn(mindelz) for p in tip.loop]
            for j in range(int(dr)):
                newloop = pym.contract(newloop,1)
                newloop = pym.aggregate(newloop,5)
                if newloop is None:
                    break
                elif len(newloop) <= 3:
                    if len(newloop) < 3:
                        newloop = None
                    break
            if newloop is None:
                continue
            uloops = pym.pinchb(newloop,10)
            uloops = [u for u in uloops if u]
            for u in uloops:
                if not pym.bccw(u):
                    u.reverse()
            for u in uloops:
                ua = pym.bareaxy(u,True)
                if ua > 100:
                    newtip = self.avert(tip,u)
                    newtips.append(newtip)
        if newtips:
            newtips.extend(self.raise_earth(newtips,e,mingrad,mindelz,depth+1))
        return newtips


    def __init__(self, b, e, steps, seald, **kws):
        '''Create a topographical map by recursively operating on polygons
        b     : is the boundary of terrain point location
        e     : is the scale of error to use during generation
        steps : is a sequence of instructions for building topography
        seald : the generation depth which should be marked as at sealevel'''
        topography.__init__(self, loop=b, style=None)
        self.boundary = b
        self.e = e
        for j,s in enumerate(steps):
            svs = [self.avert(self.root, loop=l) for l in s[0]]
            steps[j] = (svs,)+steps[j][1:]
        sealevel = 0
        for tips, dz, raises in steps:
            totaldepth = 0
            for grad, depth in raises:
                totaldepth += 0 if depth is None else depth
                tips = self.raise_earth(tips, e, 
                    mingrad=grad, mindelz=dz, maxdepth=depth)
                if totaldepth >= seald:
                    for l in tips:
                        l.style = 'sealevel'


    def interpolate(self,x,y):
        '''Given the loops of [v]+chs, return an interpolated z value based 
        on proximity of tp to the loops and their respective z positions'''
        tp = vec3(x,y,0)
        v,chs = self.loc(tp)
        vex,vtpd = pym.bnearpxy(v.loop,tp)
        z1 = v.loop[vex].z
        if chs:
            chexs,chpds = zip(*[pym.bnearpxy(c.loop,tp) for c in chs])
            chx = chpds.index(min(chpds))
            z2 = chs[chx].loop[chexs[chx]].z
            if z1 is None:
                z = z2
            else:
                dmax = max([vtpd]+list(chpds))
                dtot = (dmax-vtpd)+sum([dmax-chpd for chpd in chpds])
                if dtot == 0:dtot = 1
                z = ((dmax-vtpd)*z1+sum([(dmax-chpd)*z2 for chpd in chpds]))/dtot
        else:
            z2 = z1
            z = 0 if z1 is None else z1 
        return z
