import dilap.core.base as db
import dilap.core.tools as dpr
import dilap.core.model as dmo

import dilap.core.tmesh as dtm

import dp_vector as dpv
import dp_quaternion as dpq

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import numpy,random,pdb

class lsystem(db.base):

    # apply all rules to char 
    # if any rule is applied, return the new character
    def _rules(self,char):
        if char in self.variables:
            for r in self.rules:
                if r[0] == char:return r[1]
        else:return char

    # apply the rules to state
    def _iterate(self,instate):
        outstate = ''.join([self._rules(s) for s in instate])
        print('iteration',instate,'to',outstate)
        return outstate

    # iterate over axiom applying rules and return
    def _produce(self):
        state = self.axiom[:]
        self.variables = [r[0] for r in self.rules]
        for i in range(self.iterations):
            state = self._iterate(state)
            if len(state) > self.truncate:
                state = state[:self.truncate]
                break
        print('final product:',state)
        return state

    # given a start position and start direction, realize this lsystem
    def _realize(self,p,d):
        random.seed(self.seed)
        p,d = p.copy(),d.copy()
        lstate = self._produce()
        stack = []
        for ls in lstate:
            if   ls == '[':stack.append((p.copy(),d.copy()))
            elif ls == ']':p,d = stack.pop(-1)
            elif ls in self.grammers:self.grammer[ls](p,d)
            else:pass
        self.finaldraw()
        return self

    def __init__(self,*args,**kwargs):
        # general parameters
        self._def('iterations',5,**kwargs) # number of iterations
        self._def('seed',0,**kwargs) # seed used for random numbers

        # grammer interpretation parameters
        self._def('angle',numpy.pi/4.0,**kwargs) # angle used for rotations
        self._def('minangle',numpy.pi/32.0,**kwargs) # angle used for rotations
        self._def('maxangle',numpy.pi,**kwargs) # angle used for rotations
        self._def('xangle',numpy.pi/4.0,**kwargs) # angle used for rotations
        self._def('minxangle',numpy.pi/32.0,**kwargs) # angle used for rotations
        self._def('maxxangle',2.0*numpy.pi,**kwargs) # angle used for rotations
        self._def('yangle',numpy.pi/4.0,**kwargs) # angle used for rotations
        self._def('minyangle',numpy.pi/32.0,**kwargs) # angle used for rotations
        self._def('maxyangle',2.0*numpy.pi,**kwargs) # angle used for rotations

        self._def('width',1.0,**kwargs) # width of new edges
        self._def('rho',1.0,**kwargs) # the length of a new edge
        self._def('minrho',0.1,**kwargs)
        self._def('maxrho',10.0,**kwargs)
        self._def('polar',numpy.pi/12.0,**kwargs)
        self._def('minpolar',0.0,**kwargs)
        self._def('maxpolar',numpy.pi,**kwargs)
        self._def('azimuthal',numpy.pi/24.0,**kwargs)
        self._def('minazimuthal',0.0,**kwargs)
        self._def('maxazimuthal',2.0*numpy.pi,**kwargs)

        # grammer dependent rules/axiom
        self._def('axiom','',**kwargs)
        self._def('rules',[],**kwargs)

        # hard limit on total number of final steps
        self._def('truncate',5000,**kwargs)

        # grammer for the lsystem
        self.grammer = {
            '(':self.polar_up,')':self.polar_down,
            '{':self.azimuthal_up,'}':self.azimuthal_down,

            '+':self.pitch_up,'-':self.pitch_down,
            '/':self.yaw_up,'\\':self.yaw_down,
            '<':self.roll_up,'>':self.roll_down,

            '[':None,']':None,'&':self.randrot,'~':self.randdirrot,
            '^':self.wobblerot,'$':self.azimuthal_flip,

            '!':self.rho_up,'@':self.rho_down,
            '#':self.fatter,'%':self.thinner,

            'F':self.edge,'Q':self.term,'O':self.orient,
                }
        self.grammers = list(self.grammer.keys())

        # functions to call for F,Q, and at the end
        self._def('branchdraw',draw_branch,**kwargs)
        self._def('leafdraw',draw_leaf,**kwargs)
        self._def('finaldraw',draw,**kwargs)
                                                                    
    def rho_up(self,p,d):self.rho = dpr.clamp(2.0*self.rho,self.minrho,self.maxrho)
    def rho_down(self,p,d):self.rho = dpr.clamp(0.5*self.rho,self.minrho,self.maxrho) 
    def polar_up(self,p,d):
        if d.near(dpv.zhat) or d.near(dpv.nzhat):qv = dpv.yhat
        else:qv = d.cross(dpv.zhat)
        d.rotate(dpq.q_from_av(-self.polar,qv))
    def polar_down(self,p,d):
        if d.near(dpv.zhat) or d.near(dpv.nzhat):qv = dpv.yhat
        else:qv = d.cross(dpv.zhat)
        d.rotate(dpq.q_from_av(-self.polar,qv))
    def azimuthal_up(self,p,d):d.rotate(dpq.q_from_av(self.azimuthal,dpv.zhat))
    def azimuthal_down(self,p,d):d.rotate(dpq.q_from_av(-self.azimuthal,dpv.zhat))

    def azimuthal_flip(self,p,d):d.rotate(dpq.q_from_av(numpy.pi,dpv.zhat))

    def pitch_up(self,p,d):d.rotate(dpq.q_from_av(self.angle,dpv.xhat))
    def pitch_down(self,p,d):d.rotate(dpq.q_from_av(self.angle,dpv.nxhat))
    def yaw_up(self,p,d):d.rotate(dpq.q_from_av(self.angle,dpv.zhat))
    def yaw_down(self,p,d):d.rotate(dpq.q_from_av(self.angle,dpv.nzhat))
    def roll_up(self,p,d):d.rotate(dpq.q_from_av(self.angle,dpv.yhat))
    def roll_down(self,p,d):d.rotate(dpq.q_from_av(self.angle,dpv.nyhat))

    def randdirrot(self,p,d):
        which = random.choice(['-','<','/','+','>','\\'])
        self.grammer[which](p,d)
    def randrot(self,p,d):
        if random.random() < 0.5:newangle = 0.5*self.angle
        else:newangle = 2.0*self.angle
        self.angle = dpr.clamp(newangle,self.minangle,self.maxangle)
    # wobble d such that it stays with tolerance of a radial xy projection
    def wobblerot(self,p,d):
        polar = numpy.arcsin(d.cross(dpv.zhat).magnitude())
        azimuthal = dpv.angle_from_xaxis(d)
        if random.random() < 0.1:self.azimuthal_flip(p,d)
        elif random.random() < 0.5:self.azimuthal_up(p,d)
        else:self.azimuthal_down(p,d)
        if abs(polar) < numpy.pi/12.0:self.polar_down(p,d)
        elif polar < -numpy.pi/2.0:self.polar_down(p,d)
        elif polar > numpy.pi/2.0:self.polar_up(p,d)

    def fatter(self,p,d):self.width *= 2.0
    def thinner(self,p,d):self.width *= 0.5
    def edge(self,p,d):self.branchdraw(p.copy(),p.translate(d.copy().scale_u(self.rho)))
    def term(self,p,d):self.leafdraw(p)
    def orient(self,p,d):
        # oient direction so that it points cylindrically away from p
        d.rotate(dpq.q_from_uu(d,p)).normalize()
        self.wobblerot(p,d)

class pythagoras_tree(lsystem):

    def __init__(self,*args,**kwargs):
        self._def('axiom','Q',**kwargs)
        self._def('rules',[('F','FF'),('Q','F[<Q]>Q')],**kwargs)
        lsystem.__init__(self,*args,**kwargs)

class dragon_curve(lsystem):

    def __init__(self,*args,**kwargs):
        self._def('axiom','FX',**kwargs)
        self._def('rules',[('X','X+YF+'),('Y','-FX-Y')],**kwargs)
        self._def('iterations',9,**kwargs) # number of iterations
        self._def('angle',numpy.pi/2.0,**kwargs) # angle used for rotations
        self._def('seed',4,**kwargs) # seed used for random numbers
        lsystem.__init__(self,*args,**kwargs)

class axial_tree(lsystem):

    def __init__(self,*args,**kwargs):
        #self._def('axiom','F',**kwargs)
        #self._def('rules',[('F','F[+F]F[-F]F')],**kwargs)
        self._def('axiom','X',**kwargs)
        self._def('rules',[('X','F[+X][-X]FX'),('F','FF')],**kwargs)
        self._def('iterations',5,**kwargs) # number of iterations
        self._def('angle',dpr.rad(25.7),**kwargs) # angle used for rotations
        self._def('seed',0,**kwargs) # seed used for random numbers
        lsystem.__init__(self,*args,**kwargs)

class plant(lsystem):

    def __init__(self,*args,**kwargs):
        self._def('axiom','X',**kwargs)
        self._def('rules',[('X','F-[[X]+X]+F[+FX]-X'),('F','FF')],**kwargs)
        self._def('iterations',3,**kwargs) # number of iterations
        self._def('angle',dpr.rad(25),**kwargs) # angle used for rotations
        self._def('seed',0,**kwargs) # seed used for random numbers
        lsystem.__init__(self,*args,**kwargs)

# an ltree is an lsystem that produces structured information about a tree
# such that it can produce a mesh from that information representing the tree
class ltree(lsystem):
                                                    
    def _realize(self,p,d):
        self.branches = []
        self.leaves = []
        return lsystem._realize(self,p,d)

    loadouts = []
    loadouts.append(('A',[
        ('F','F~[~FA]F'),
        ('A','[--////FQA]-<F[>>\\\\FQA]>+F[-\>>QA\[X]][</++QA/[Y]]Q'),
        ('X','F+Q'),('Y','F-Q')]))

    def __init__(self,ldx,*args,**kwargs):
        axiom,rules = self.loadouts[ldx]
        self._def('axiom',axiom,**kwargs)
        self._def('rules',rules,**kwargs)

        self._def('seed',1,**kwargs) # seed used for random numbers
        self._def('iterations',5,**kwargs) # number of iterations
        self._def('angle',numpy.pi/12.0,**kwargs) # angle used for rotations
        self._def('minangle',numpy.pi/12.0,**kwargs) # angle used for rotations
        self._def('maxangle',numpy.pi/6,**kwargs) # angle used for rotations

        self._def('branchdraw',self.draw_branch,**kwargs)
        self._def('leafdraw',self.draw_leaf,**kwargs)
        self._def('finaldraw',self.draw,**kwargs)
        lsystem.__init__(self,*args,**kwargs)

    def draw_branch(self,p,n):
        self.branches.append((p.copy(),n.copy()))
        draw_branch(p,n)

    def draw_leaf(self,p):
        self.leaves.append((p.copy()))
        draw_leaf(p)

    def draw(self):
        es,ls = self.branches,self.leaves

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
        draw()

draws = []
def draw_branch(p,n):draws.append((([p.x,n.x],[p.y,n.y]),{'zs':[p.z,n.z]}))
def draw_leaf(p):draws.append((([p.x],[p.y]),{'zs':[p.z],'marker':'o'}))
def draw():
    global draws
    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    for d in draws:ax.plot(*d[0],**d[1])
    plt.show()
    draws = []

def test():
    import dilap.construct as dlc
    p = dpv.zero()
    d = dpv.zhat.copy()

    #pythagoras_tree()._realize(p,d)
    #dragon_curve()._realize(p,d)

    total = dmo.model()

    #for l in range(5):tree(l)._realize(p,d)
    ltree(-1)._realize(p,d)
    tps = [(x,y) for x in range(3) for y in range(3)]
    for i,xy in enumerate(tps):
        x,y = xy
        kws = {'seed':i}
        lmod = ltree(-1,**kws)._realize(p,d).model
        lmod.translate(dpv.vector(10*x,10*y,0))
        total._consume(lmod)
    
    dlc.build(total)

    #plant()._realize(p,d)
    #axial_tree()._realize(p,d)

    #ls = rtree()
    #for x in range(5):
    #    ls._randomize()
    #    ls._realize(p,d)



