import dilap.core.base as db
import dilap.core.tools as dpr

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
        draw()

    def __init__(self,*args,**kwargs):
        # general parameters
        self._def('iterations',5,**kwargs) # number of iterations
        self._def('angle',numpy.pi/4.0,**kwargs) # angle used for rotations
        self._def('minangle',numpy.pi/32.0,**kwargs) # angle used for rotations
        self._def('maxangle',numpy.pi,**kwargs) # angle used for rotations

        self._def('seed',0,**kwargs) # seed used for random numbers
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

        self._def('axiom','',**kwargs)
        self._def('rules',[],**kwargs)

        self._def('truncate',5000,**kwargs)

        self.grammer = {
            '(':self.polar_up,')':self.polar_down,
            '{':self.azimuthal_up,'}':self.azimuthal_down,

            '+':self.pitch_up,'-':self.pitch_down,
            '/':self.yaw_up,'\\':self.yaw_down,
            '<':self.roll_up,'>':self.roll_down,

            '[':None,']':None,'&':self.randrot,
            '^':self.wobblerot,'$':self.azimuthal_flip,

            '!':self.rho_up,'@':self.rho_down,
            '#':self.fatter,'%':self.thinner,

            'F':self.edge,'Q':self.term,'O':self.orient,
                }
        self.grammers = list(self.grammer.keys())
                                                                    
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
    def edge(self,p,d):draw_branch(p.copy(),p.translate(d.copy().scale_u(self.rho)))
    def term(self,p,d):draw_leaf(p)
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

class tree(lsystem):

    def __init__(self,ldx,*args,**kwargs):
        loadouts = []
        loadouts.append(('+++A',[('A','F[<+A][<-A]')]))
        loadouts.append(('F',[('F','FF-[-F+F+F]+[+F-F-F]')]))
        loadouts.append(('VZFFF',[
            ('I','VZFFF'),('V','[+++W][---W]YV'),('W','+X[-W]Z'),
            ('X','-W[+X]Z'),('Y','YZ'),('Z','[-FFF][+FFF]F'),('F','FF')]))
        loadouts.append(('F[+Q]',[
            ('F','F[//&//@F!+<//@F!&+//+Q][-->\\@F!\\&\\@F!]'),
            ('Q','F[@F!<//&//Q][\\&\\@F!&->Q][&\\&+\\&\\&++&<@F!Q]')]))
        loadouts.append(('Q',[
            ('F','F<[++/F]+[-\\F]->'),('Q','F[X][Y]'),
            ('X','+>\\@F!Q'),('Y','<-@F!Q//@F!')]))
        axiom,rules = loadouts[ldx]

        self._def('axiom',axiom,**kwargs)
        self._def('rules',rules,**kwargs)

        self._def('iterations',4,**kwargs) # number of iterations
        self._def('angle',numpy.pi/12.0,**kwargs) # angle used for rotations
        self._def('minangle',numpy.pi/12.0,**kwargs) # angle used for rotations
        self._def('maxangle',numpy.pi/6,**kwargs) # angle used for rotations

        self._def('seed',0,**kwargs) # seed used for random numbers
        lsystem.__init__(self,*args,**kwargs)
                                                    
class axial_tree(lsystem):

    def __init__(self,*args,**kwargs):
        #self._def('axiom','F',**kwargs)
        #self._def('rules',[('F','F[+F]F[-F]F')],**kwargs)
        self._def('axiom','X',**kwargs)
        #self._def('rules',[('X','F[+X][-X]FX'),('F','FF')],**kwargs)
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

class rtree(lsystem):

    def __init__(self,*args,**kwargs):
        lsystem.__init__(self,*args,**kwargs)
        self._randomize()



        self.grammer_variety = {
            'rotations':['+','-','/','\\','<','>'],
            'edgehow':['!','@','#','%'],
            'grow':['F','Q'],
                }
        self.grammer_varieties = list(self.grammer_variety.keys())

    def _random_grammer(self):
        #variety = random.choice(['rotations','edgehow','grow'])
        #variety = random.choice(['rotations','grow'])
        variety = random.choice(['grow'])
        rg = random.choice(self.grammer_variety[variety])
        return rg

    def _random_sequence(self,start = 1,stop = 5):
        slen = random.randint(start,stop)
        seqs = []
        for x in range(slen):
            rg = self._random_grammer()
            #rg = random.choice(self.grammers)
            seqs.append(rg)
        return seqs

    def _pairop(self,seqs,op = ('[',']')):
        x = random.randint(0,len(seqs)-2)
        y = random.randint(x+2,len(seqs))
        seqs.insert(x,op[0])
        seqs.insert(y,op[1])

    def _random_rule(self):
        rin  = self._random_sequence(1,1)
        rout = self._random_sequence(2,10)
        for pdx in range(3):
            if random.random() < 0.25:continue
            op = random.choice([('[',']'),('<','>'),('+','-'),('/','\\')])
            self._pairop(rout,op)
        return (''.join(rin),''.join(rout))

    def _random_rules(self,rulemin = 1,rulemax = 5):
        rcnt = random.randint(rulemin,rulemax)
        rules = [self._random_rule() for x in range(rcnt)]
        for r in rules:print('random rule:',r)  
        return rules


    # a stemrule replaces a segment of stem with more stem
    def _stemrule(self):
        return ('F','FF')

    # a branchrule replaces a segment of stem with more stems
    # and uses push/pop to make branches
    def _branchrule(self):
        return ('Q','F[+[<FQ]>FQ]-')




        branch = ['F']

        openers = ['[','<','+','/']
        closers = {'[':']','<':'>','+':'-','/':'\\'}
        later = []
        for i in range(5):
            if branch[-1] == 'F':
                branch.append('[')
                later.append(']')
                b = random.choice(openers)
                later.append(closers[b])
            elif branch[-1] in openers:
                b = 'F'*random.randint(0,2)+'Q'
            branch.append(b)
            if later and random.random() < 0.5:
                branch.append(later.pop(0))
            if branch[-1] == ']':
                branch.extend(later)
                later = []
        for l in later:branch.append(l)
        branch = ''.join(branch)

        print('branchrule',branch)

        return ('Q',branch)
        #return ('Q','F[<Q]>Q')

    def _axiom(self):
        return 'Q'

    def _randomize(self):
        self.rules = []
        self.rules.append(self._stemrule())
        self.rules.append(self._branchrule())
        self.axiom = self._axiom()
        self.iterations = 5
        self.angle = dpr.rad(25)

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
    p = dpv.zero()
    d = dpv.zhat.copy()

    #ls = pythagoras_tree()
    #ls._realize(p,d)

    #ls = dragon_curve()

    ls = tree(3)
    ls._realize(p,d)

    #ls = plant()
    ls = axial_tree()
    ls._realize(p,d)

    #ls = rtree()
    #for x in range(5):
    #    ls._randomize()
    #    ls._realize(p,d)



