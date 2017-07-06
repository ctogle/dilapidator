import dilap.geometry.tools as dpr
from dilap.geometry import *
from dilap.topology import *
import dilap.core.plotting as dtl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from io import StringIO as sio
import numpy
import random
import pdb


def lgen(p,d,axiom,rules,i,**kws):
    yield from lstate(i,p,d,axiom,rules,**kws)


class lstate(tree):


    drho = 1
    dpolar     = numpy.pi/2
    dazimuthal = numpy.pi/2
    dpitch     = numpy.pi/2
    dyaw       = numpy.pi/2
    droll      = numpy.pi/2


    def avert(self):
        if not hasattr(self,'tip'):self.tip = self.root
        oldtip = self.tip.p.cp(),self.tip.d.cp()
        self.tip = tree.avert(self,self.tip)
        self.tip.p,self.tip.d = oldtip


    def __init__(self,i,p,d,axiom,rules,**kws):
        tree.__init__(self)
        for k in kws:self.__setattr__(k,kws[k])
        self.tip = self.root
        self.tip.p,self.tip.d = p,d
        self.i = i
        self.axiom = axiom
        self.rules = rules


    def __call__(self,**kws):
        for s in self.__iter__(**kws):
            pass
        return self


    def __iter__(self,**kws):
        for k in kws:
            self.__setattr__(k,kws[k])
        for s in lstring(self.axiom,self.rules,self.i).produce():
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


    #def __iter__(self):
    #    yield from self.produce()


x ,y ,z  = vec3( 1,0,0),vec3(0, 1,0),vec3(0,0, 1)
nx,ny,nz = vec3(-1,0,0),vec3(0,-1,0),vec3(0,0,-1)
class lgrammer:

    push = lambda ls : ls.avert()
    def pop(ls):ls.tip = ls.above(ls.tip)
    edge = lambda ls : (ls.tip.p.cp(),ls.tip.p.trn(ls.tip.d.cp().uscl(ls.drho)))
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
    
    pitch_u = lambda ls : d.rot(quat(0,0,0,0).av( ls.dpitch,vec3(1,0,0)))
    pitch_d = lambda ls : d.rot(quat(0,0,0,0).av(-ls.dpitch,vec3(1,0,0)))
    yaw_u   = lambda ls : d.rot(quat(0,0,0,0).av( ls.dyaw,  vec3(0,0,1)))
    yaw_d   = lambda ls : d.rot(quat(0,0,0,0).av(-ls.dyaw,  vec3(0,0,1)))
    roll_u  = lambda ls : d.rot(quat(0,0,0,0).av( ls.droll, vec3(0,1,0)))
    roll_d  = lambda ls : d.rot(quat(0,0,0,0).av(-ls.droll, vec3(0,1,0)))

    dic = {
        '{':push,'}':pop,'F':edge,'Q':term,
        '(':polar_u,'[':azimuthal_u,'<':pitch_u,'+':yaw_u,'=':roll_u,
        ')':polar_d,']':azimuthal_d,'>':pitch_d,'-':yaw_d,'|':roll_d,
                    '$':azimuthal_f,
            }
        #'&':self.randrot,'~':self.randdirrot,
        #'^':self.wobblerot,
        #'O':self.orient,

    '''#
    self._def('truncate',5000,**kwargs)
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
def orient(self,p,d):
    # oient direction so that it points cylindrically away from p
    d.rot(quat(0,0,0,0).uu(d,p)).nrm()
    self.wobblerot(p,d)
    '''#
        

