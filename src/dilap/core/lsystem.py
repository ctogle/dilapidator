from dilap.geometry import *
from dilap.topology import *
from io import StringIO as sio
from numpy import pi


class lstate(tree):


    drho = 1
    dpolar     = pi/2
    dazimuthal = pi/2
    dpitch     = pi/2
    dyaw       = pi/2
    droll      = pi/2


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
                if piece:
                    yield piece


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
    azimuthal_f = lambda ls : lgrammer.azimuthal(ls,pi/ls.dazimuthal)
    
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
