import dilap.core.context as dgc
import dilap.core.tools as dpr
import dilap.core.tmesh as dtm
import dilap.core.lsystem as dls

import dilap.primitive.cube as dcu
import dilap.primitive.road as dr

import dp_vector as dpv
import dp_quaternion as dpq

import pdb,numpy

class rmesh(dls.lsystem):

    def _realize(self,p,d):
        self.edges = []
        self.nodes = []
        return dls.lsystem._realize(self,p,d)

    loadouts = []
    loadouts.append(('FX',[
        ('X','X/YF/'),
        ('Y','\FX\Y')]))
    loadouts.append(('A',[
        ('F','F~[~FA]F'),
        ('A','[--////FQA]-<F[>>\\\\FQA]>+F[-\>>QA\[X]][</++QA/[Y]]Q'),
        ('X','F+Q'),('Y','F-Q')]))

    #def __init__(self,*args,**kwargs):
    #    self._def('nodes',[],**kwargs)
    #    self._def('edges',[],**kwargs)

    def __init__(self,ldx,*args,**kwargs):
        axiom,rules = self.loadouts[ldx]
        self._def('axiom',axiom,**kwargs)
        self._def('rules',rules,**kwargs)

        self._def('seed',1,**kwargs) # seed used for random numbers
        self._def('iterations',5,**kwargs) # number of iterations

        self._def('angle',numpy.pi/2.0,**kwargs) # angle used for rotations
        self._def('minangle',numpy.pi/12.0,**kwargs) # angle used for rotations
        self._def('maxangle',numpy.pi*2.0,**kwargs) # angle used for rotations

        self._def('branchdraw',self.draw_branch,**kwargs)
        self._def('leafdraw',self.draw_leaf,**kwargs)
        self._def('finaldraw',self.draw,**kwargs)
        dls.lsystem.__init__(self,*args,**kwargs)

    def draw_branch(self,p,n):
        self.edges.append((p.copy(),n.copy()))
        dls.draw_branch(p,n)

    def draw_leaf(self,p):
        self.nodes.append((p.copy()))
        dls.draw_leaf(p)

    def draw(self):
        es,ls = self.edges,self.nodes

        nps = []
        nes = []
        for ex in range(len(es)):
            e1,e2 = es[ex]
            e1x,e2x = None,None
            for pdx in range(len(nps)):
                if e1.near(nps[pdx]):e1x = pdx
                if e2.near(nps[pdx]):e2x = pdx
                if not e1x is None and not e2x is None:break
            if e1x is None:
                e1x = len(nps)
                nps.append(e1)
            if e2x is None:
                e2x = len(nps)
                nps.append(e2)
            nes.append((e1x,e2x))

        m = dtm.meshme(nps,None,None,None,nes,[])
        smod = m.skeleton()
        self.model = smod
        dls.draw()

class infrastructure(dgc.context):

    def _terrain_points(self):
        self.tpts = []
        [self.tpts.extend(r.tpts) for r in self.roads]
        return self.tpts

    def _hole_points(self):
        self.hpts = []
        [self.hpts.extend(r.hpts) for r in self.roads]
        return self.hpts

    def generate(self,seed,boundary,worn = 0):
        self.intersections = []
        self.roads = []
        rcs,ris = dr.circle(dr.highway)
        self.intersections.extend([i.translate(seed) for i in ris])
        self.roads.extend([r.translate(seed) for r in rcs])

        p = dpv.zero()
        d = dpv.xhat.copy()
        rmeshmodel = rmesh(-1)._realize(p,d).model
        rmeshmodel.translate(seed)
        self._nodes_to_graph(self._node_wrap(rmeshmodel))

        [self._nodes_to_graph(self._node_wrap(i)) for i in self.intersections]
        [self._nodes_to_graph(self._node_wrap(r)) for r in self.roads]
        return self



