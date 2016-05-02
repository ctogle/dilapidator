import dilap.core.base as db
#import dilap.core.vector as dpv
#import dilap.core.quaternion as dpq
#import dilap.core.tools as dpr
from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat
import dilap.geometry.tools as dpr

#import dilap.mesh.tools as dtl
import dilap.core.plotting as dtl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import numpy,random

import pdb

###############################################################################
### simple plotting for branches/leaves
###############################################################################

draws = []
def draw_branch(p,n):draws.append((([p.x,n.x],[p.y,n.y]),{'zs':[p.z,n.z]}))
def draw_leaf(p):draws.append((([p.x],[p.y]),{'zs':[p.z],'marker':'o'}))
def draw(ax = None):
    global draws
    if ax is None:ax = dtl.plot_axes()
    for d in draws:ax.plot(*d[0],**d[1])
    draws = []

###############################################################################
###############################################################################
###############################################################################

###############################################################################
### base lsystem class
###############################################################################

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
        #print('iteration',instate,'to',outstate)
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
        #print('final product:',state)
        return state

    # given a start position and start direction, realize this lsystem
    def _realize(self,p,d,ax = None):
        random.seed(self.seed)
        p,d = p.cp(),d.cp()
        lstate = self._produce()
        stack = []
        for ls in lstate:
            if   ls == '[':stack.append((p.cp(),d.cp()))
            elif ls == ']':p,d = stack.pop(-1)
            elif ls in self.grammers:self.grammer[ls](p,d)
            else:pass
        self.finaldraw(ax)
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
                                                                    
    def rho_up(self,p,d):
        self.rho = dpr.clamp(2.0*self.rho,self.minrho,self.maxrho)
    def rho_down(self,p,d):
        self.rho = dpr.clamp(0.5*self.rho,self.minrho,self.maxrho) 
    def polar_up(self,p,d):
        if d.isnear(vec3(0,0,1)) or d.isnear(vec3(0,0,-1)):qv = vec3(0,1,0)
        else:qv = d.crs(vec3(0,0,1))
        #d.rot(quat(0,0,0,0).av(-self.polar,qv))
        d.rot(quat(0,0,0,0).av(self.polar,qv))
    def polar_down(self,p,d):
        if d.isnear(vec3(0,0,1)) or d.isnear(vec3(0,0,-1)):qv = vec3(0,1,0)
        else:qv = d.crs(vec3(0,0,1))
        d.rot(quat(0,0,0,0).av(-self.polar,qv))
    def azimuthal_up(self,p,d):
        d.rot(quat(0,0,0,0).av(self.azimuthal,vec3(0,0,1)))
    def azimuthal_down(self,p,d):
        d.rot(quat(0,0,0,0).av(-self.azimuthal,vec3(0,0,1)))
    def azimuthal_flip(self,p,d):
        d.rot(quat(0,0,0,0).av(numpy.pi,vec3(0,0,1)))

    def pitch_up(self,p,d):
        d.rot(quat(0,0,0,0).av(self.angle,vec3(1,0,0)))
    def pitch_down(self,p,d):
        d.rot(quat(0,0,0,0).av(-self.angle,vec3(1,0,0)))
    def yaw_up(self,p,d):
        d.rot(quat(0,0,0,0).av(self.angle,vec3(0,0,1)))
    def yaw_down(self,p,d):
        d.rot(quat(0,0,0,0).av(-self.angle,vec3(0,0,1)))
    def roll_up(self,p,d):
        d.rot(quat(0,0,0,0).av(self.angle,vec3(0,1,0)))
    def roll_down(self,p,d):
        d.rot(quat(0,0,0,0).av(-self.angle,vec3(0,1,0)))

    def randdirrot(self,p,d):
        which = random.choice(['-','<','/','+','>','\\'])
        self.grammer[which](p,d)
    def randrot(self,p,d):
        if random.random() < 0.5:newangle = 0.5*self.angle
        else:newangle = 2.0*self.angle
        self.angle = dpr.clamp(newangle,self.minangle,self.maxangle)
    # wobble d such that it stays with tolerance of a radial xy projection
    def wobblerot(self,p,d):
        polar = numpy.arcsin(d.crs(vec3(0,0,1)).mag())
        print('need wobblerot work')
        raise NotImplemented
        azimuthal = dpv.angle_from_xaxis(d)
        if random.random() < 0.1:self.azimuthal_flip(p,d)
        elif random.random() < 0.5:self.azimuthal_up(p,d)
        else:self.azimuthal_down(p,d)
        if abs(polar) < numpy.pi/12.0:self.polar_down(p,d)
        elif polar < -numpy.pi/2.0:self.polar_down(p,d)
        elif polar > numpy.pi/2.0:self.polar_up(p,d)

    def fatter(self,p,d):self.width *= 2.0
    def thinner(self,p,d):self.width *= 0.5
    def edge(self,p,d):self.branchdraw(p.cp(),p.trn(d.cp().uscl(self.rho)))
    def term(self,p,d):self.leafdraw(p)
    def orient(self,p,d):
        # oient direction so that it points cylindrically away from p
        d.rot(quat(0,0,0,0).uu(d,p)).nrm()
        self.wobblerot(p,d)
        
###############################################################################
###############################################################################
###############################################################################





###############################################################################
### some specific lsystems from the internet
###############################################################################

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

###############################################################################
###############################################################################
###############################################################################





