import dilap.core.base as db
from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat
import dilap.geometry.tools as dpr
import dilap.topology.tree as dtr

import dilap.core.plotting as dtl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from io import StringIO as sio

import numpy,random

import pdb

def lgen(p,d,axiom,rules,i,**kws):
    ls = lstate(i,p,d,axiom,rules,**kws)
    for piece in ls:yield piece

class lstate(dtr.tree):
    rho = 1.0
    dpolar = numpy.pi/2
    dazimuthal = numpy.pi/2

    def avert(self):
        if not hasattr(self,'tip'):self.tip = self.root
        oldtip = self.tip.p.cp(),self.tip.d.cp()
        self.tip = dtr.tree.avert(self,self.tip)
        self.tip.p,self.tip.d = oldtip

    def __init__(self,i,p,d,axiom,rules,**kws):
        dtr.tree.__init__(self)
        for k in kws:self.__setattr__(k,kws[k])
        self.tip = self.root
        self.tip.p,self.tip.d = p,d
        self.i = i
        self.axiom = axiom
        self.rules = rules

    def __call__(self,**kws):
        for s in self.__iter__(**kws):pass
        return self

    def __iter__(self,**kws):
        for k in kws:self.__setattr__(k,kws[k])
        for s in lstring(self.axiom,self.rules,self.i):
            if s in lgrammer.dic:
                piece = lgrammer.dic[s](self)
                if piece:yield piece

class lstring:
    def __init__(self,axiom,rules,i = 3):
        self.axiom = axiom
        self.rules = rules
        self.i = i
    def produce(self):
        state = self.axiom[:]
        for i in range(self.i):
            out = sio()
            for c in state:
                if c in self.rules:
                    out.write(self.rules[c])
                else:out.write(c)
            state = out.getvalue()
        return state
    def __iter__(self):
        nonf = False
        for c in self.produce():
            if nonf:
                yield c
            elif c != 'F':nonf = True

x ,y ,z  = vec3( 1,0,0),vec3(0, 1,0),vec3(0,0, 1)
nx,ny,nz = vec3(-1,0,0),vec3(0,-1,0),vec3(0,0,-1)
class lgrammer:
    push = lambda ls : ls.avert()
    def pop(ls):ls.tip = ls.above(ls.tip)
    edge = lambda ls : (ls.tip.p.cp(),ls.tip.p.trn(ls.tip.d.cp().uscl(ls.rho)))
    term = lambda ls : ls.tip.p.cp()

    getqv = lambda d : y if (d.isnear(z) or d.isnear(nz)) else z
    def polar(ls,f = 1.0):
        ls.tip.d.rot(quat(0,0,0,0).av(f*ls.dpolar,lgrammer.getqv(ls.tip.d)))
    polar_u = lambda ls : lgrammer.polar(ls, 1)
    polar_d = lambda ls : lgrammer.polar(ls,-1)
    def azimuthal(ls,f = 1.0):
        ls.tip.d.rot(quat(0,0,0,0).av(f*ls.dazimuthal,lgrammer.getqv(ls.tip.d)))
    azimuthal_u = lambda ls : lgrammer.azimuthal(ls, 1)
    azimuthal_d = lambda ls : lgrammer.azimuthal(ls,-1)
    azimuthal_f = lambda ls : lgrammer.azimuthal(ls,numpy.pi/ls.dazimuthal)

    dic = {
        '{':push,'}':pop,'F':edge,'Q':term,
        '(':polar_u,'[':azimuthal_u,
        ')':polar_d,']':azimuthal_d,
                    '$':azimuthal_f,
            }
        #'+':self.pitch_up,'-':self.pitch_down,
        #'/':self.yaw_up,'\\':self.yaw_down,
        #'<':self.roll_up,'>':self.roll_down,
        #'&':self.randrot,'~':self.randdirrot,
        #'^':self.wobblerot,
        #'!':self.rho_up,'@':self.rho_down,
        #'#':self.fatter,'%':self.thinner,
        #'O':self.orient,









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
        self._def('maxpolar',numpy.pi/2.0,**kwargs)
        self._def('azimuthal',numpy.pi/24.0,**kwargs)
        self._def('minazimuthal',0.0,**kwargs)
        self._def('maxazimuthal',2.0*numpy.pi,**kwargs)

        # grammer dependent rules/axiom
        self._def('axiom','',**kwargs)
        self._def('rules',[],**kwargs)
        self.rules = [r for r in self.rules if not r is None]

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





