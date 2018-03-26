from dilap.geometry import *
from dilap.topology import *
from .plotting import *
from io import StringIO as sio
from numpy import pi, cos, sin, arctan


class lstate(tree):

    drho = 1
    dpolar     = pi/2
    dazimuthal = pi/2
    dtwist     = pi/2
    dpitch     = pi/2
    dyaw       = pi/2
    droll      = pi/2

    def avert(self):
        if not hasattr(self,'tip'):self.tip = self.root
        oldtip = self.tip.p.cp(),self.tip.d.cp(),self.tip.t,self.tip.ld
        self.tip = tree.avert(self,self.tip)
        self.tip.p,self.tip.d,self.tip.t,self.tip.ld = oldtip
        self.tip.term = False

    def __init__(self,i,p,d,axiom,rules,grammer=None,**kws):
        tree.__init__(self)
        for k in kws:self.__setattr__(k,kws[k])
        self.rootp = p.cp()
        self.tip = self.root
        self.tip.p,self.tip.d,self.tip.t,self.tip.ld = p,d,0,d.cp()
        #print('lsystem ROOT',self.tip.p)
        self.tip.term = False
        self.i = i
        self.axiom = axiom
        self.rules = rules
        if grammer is None:
            grammer = lgrammer
        self.grammer = grammer

    def __call__(self,**kws):
        for s in self.__iter__(**kws):
            pass
        return self

    def __iter__(self,**kws):
        for k in kws:
            self.__setattr__(k,kws[k])
        for s in lstring(self.axiom,self.rules,self.i).produce():
            if s in self.grammer.dic:
                piece = self.grammer.dic[s](self)
                if piece:
                    yield piece

    def plot(self, v=None, ax=None):
        if ax is None:
            ax = plot_axes(200)
        if v is None:
            v = self.root
            plot_edges((self.rootp,v.p), ax)
        a = self.above(v)
        if a:
            plot_edges((v.p,a.p), ax)
        for b in self.below(v):
            self.plot(b, ax=ax)


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
    def edge(ls):
        st = ls.tip.p.cp()
        tw = quat(0,0,0,0).av(ls.tip.t, ls.tip.ld)
        d = ls.tip.d.cp().rot(tw)
        ls.tip.ld = d.cp()
        ed = ls.tip.p.trn(d.uscl(ls.drho))
        return st,ed
    def term(ls):
        ls.tip.term = True
        return ls.tip.p.cp()

    def twist(ls, f=1.0):
        ls.tip.t += f*ls.dtwist
    twist_u = lambda ls : lgrammer.twist(ls, 1.0)
    twist_d = lambda ls : lgrammer.twist(ls,-1.0)

    theta = lambda p : arctan(p.y/p.x) if not p.x == 0 else (pi/2 if p.y > 0 else -pi/2)
    def polar(ls, f=1.0):
        th = lgrammer.theta(ls.tip.p)
        if ls.tip.d.isnear(z) or ls.tip.d.isnear(nz):
            qv = y
        else:
            qv = vec3(-sin(th), cos(th), 0)
        ls.tip.d.rot(quat(0,0,0,0).av(f*ls.dpolar,qv))
    polar_u = lambda ls : lgrammer.polar(ls, 1)
    polar_d = lambda ls : lgrammer.polar(ls,-1)

    getqv = lambda d : y if (d.isnear(z) or d.isnear(nz)) else z

    def azimuthal(ls, f=1.0):
        ls.tip.d.rot(quat(0,0,0,0).av(f*ls.dazimuthal,lgrammer.getqv(ls.tip.d)))
    azimuthal_u = lambda ls : lgrammer.azimuthal(ls, 1)
    azimuthal_d = lambda ls : lgrammer.azimuthal(ls,-1)
    azimuthal_f = lambda ls : lgrammer.azimuthal(ls,pi/ls.dazimuthal)

    pitch_u = lambda ls : ls.tip.d.rot(quat(0,0,0,0).av( ls.dpitch,vec3(1,0,0)))
    pitch_d = lambda ls : ls.tip.d.rot(quat(0,0,0,0).av(-ls.dpitch,vec3(1,0,0)))
    yaw_u   = lambda ls : ls.tip.d.rot(quat(0,0,0,0).av( ls.dyaw,  vec3(0,0,1)))
    yaw_d   = lambda ls : ls.tip.d.rot(quat(0,0,0,0).av(-ls.dyaw,  vec3(0,0,1)))
    roll_u  = lambda ls : ls.tip.d.rot(quat(0,0,0,0).av( ls.droll, vec3(0,1,0)))
    roll_d  = lambda ls : ls.tip.d.rot(quat(0,0,0,0).av(-ls.droll, vec3(0,1,0)))

    dic = {
        '{':push,'}':pop,'F':edge,'Q':term,
        '^':twist_u,'(':polar_u,'[':azimuthal_u,'<':pitch_u,'+':yaw_u,'=':roll_u,
        '_':twist_d,')':polar_d,']':azimuthal_d,'>':pitch_d,'-':yaw_d,'|':roll_d,
                               '$':azimuthal_f,
            }
        #'&':self.randrot,'~':self.randdirrot,
        #'^':self.wobblerot,
        #'O':self.orient,
