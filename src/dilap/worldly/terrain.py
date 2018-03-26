from dilap.core import *
from dilap.geometry import *
import dilap.geometry.polymath as pym
from .lsysgen import waterway
import pdb


class terrain(topography):


    '''Nesting the loops with respect to containment and offsetting them 
    in the z direction can represent a topographical map for terrain'''


    def raise_earth(self,tips,e,increment = 1.0,
            mingrad = 1.0,mindelz = -5.0,depth = 0,maxdepth = None):

        print('... raising earth ... (depth: %d)' % depth)
        newtips = []

        if maxdepth and depth == maxdepth:
            print('... reached max depth (%d)! ...' % maxdepth)
            return newtips

        print('... tipcount (%d) ...' % len(tips))
        for tip in tips:
            newloop = [p.cp().ztrn(mindelz) for p in tip.loop]

            # NEED A ROUTINE FOR WHEN TIP HAS A HOLE REPRESENTED
            # WHEN CREATING LOOPS, CREAT FROM HOLES SO AS TO MEET HALFWAY
            useless = lambda l : None if (l is None or len(l) < 3) else l

            dr = abs(mindelz/mingrad)/increment
            #print('... dr (%f) ...' % dr)
            for j in range(int(dr)):

                #newloop = pym.smart_contract(newloop,increment)
                newloop = pym.contract(newloop,increment)
                newloop = pym.aggregate(newloop,5)
                newloop = useless(newloop)

                if newloop is None or len(newloop) == 3:
                    break

            #print('... pinchhh () ...')
            if newloop:
                #uloops = pym.pinchb2(newloop,10)
                uloops = pym.pinchb(newloop,10)
                uloops = [u for u in uloops if u]
                for u in uloops:
                    if not pym.bccw(u):
                        u.reverse()
                    ua = pym.bareaxy(u,True)
                    if ua > 100:
                        newtip = self.avert(tip,u)
                        newtips.append(newtip)

        if newtips:
            newtips.extend(self.raise_earth(
                newtips,e,increment,mingrad,mindelz,depth+1))

        return newtips


    @classmethod
    def from_boundary(cls, b, e):
        wide = pym.contract(b, -e*100)

        high = pym.splitb(wide, e*20)
        high = pym.smoothxyi(high, w = 0.5, epsilon = 0.1, i = 10, constraint = 0)
        water = waterway(wide)

        water.plotxy(l=400)
        plt.show()

        scale = 10.0
        steps = [
            ([[p.cp() for p in pym.splitb(wide, e*20)]], scale, ((0.2 , None), ))
                ]
            #(lows , -scale, ((0.08, 5), (0.02, 3), (0.1, 5), (0.2, None))), ]
            #([l[1] for l in lows], -scale, ((0.2, None), )), ]
        sealevel = 9
        return cls(wide, e, steps, sealevel)


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
        for tips, dz, raises in steps:
            print(len(tips))
            print(tips)
            print(len(tips))
            increment = abs(dz)
            totaldepth = 0
            for grad, depth in raises:
                totaldepth += 0 if depth is None else depth
                tips = self.raise_earth(tips, e, increment=increment,
                            mingrad=grad, mindelz=dz, maxdepth=depth)
                if totaldepth >= seald:
                    for l in tips:
                        l.style = 'sealevel'
        ax = self.plot()
        plt.show()


    def ___interpolate(self,x,y):
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
                '''
                compute an interpolated normal vector
                compute distance to each polygon along this vector
                weight average of two z values based on this distance
                '''
                dmax = max([vtpd]+list(chpds))
                dtot = (dmax-vtpd)+sum([dmax-chpd for chpd in chpds])
                if dtot == 0:
                    dtot = 1
                    print('dtot==0')
                z = ((dmax-vtpd)*z1+sum([(dmax-chpd)*z2 for chpd in chpds]))/dtot
                z = z2
        else:
            z2 = z1
            z = 0 if z1 is None else z1 
        return z


    def interpolate(self, x, y):
        '''Given the loops of [v]+chs, return an interpolated z value based 
        on proximity of tp to the loops and their respective z positions'''
        tp = vec3(x, y, 0)
        v, chs = self.loc(tp)
        #vex, vtpd = pym.bnearpxy_smooth(v.loop, tp)
        vex, vtpd = pym.bnearpxy(v.loop, tp)
        z1 = v.loop[vex].z
        if chs:
            #chexs, chpds = zip(*[pym.bnearpxy_smooth(c.loop,tp) for c in chs])
            chexs, chpds = zip(*[pym.bnearpxy(c.loop,tp) for c in chs])
            chx = chpds.index(min(chpds))
            z2 = chs[chx].loop[chexs[chx]].z
            if z1 is None:
                z = z2
            else:
                '''
                compute an interpolated normal vector
                compute distance to each polygon along this vector
                weight average of two z values based on this distance
                '''
                dmax = max([vtpd]+list(chpds))
                dtot = (dmax-vtpd)+sum([dmax-chpd for chpd in chpds])
                if dtot == 0:
                    dtot = 1
                    print('dtot==0')
                z = ((dmax-vtpd)*z1+sum([(dmax-chpd)*z2 for chpd in chpds]))/dtot
                z = z2
        else:
            z2 = z1
            z = 0 if z1 is None else z1 
        return z
