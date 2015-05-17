import dilap.core.base as db

import dp_vector as dpv
import dp_quaternion as dpq

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import numpy,random,pdb

class lsystem(db.base):

    def __init__(self,*args,**kwargs):
        self._def('iterations',5,**kwargs)
        # the true definition of the lsystem
        self._def('variables',[],**kwargs)
        self._def('constants',[],**kwargs)
        self._def('axiom','',**kwargs)
        self._def('rules',[],**kwargs)

        link = self._link # function return a link length
        push = self._push # function called on p,d after pushing
        pop = self._pop # function called on p,d after poping
        extd = self._extend # function called which "extends"
        term = self._terminate # function called which "terminates"
        self.grammer_operations = [
            self._link,self._push,self._pop,self._extend,self._terminate]

    # apply rule r to character char
    # if a change was applied, return the new character
    # otherwise return None
    def _rule(self,r,char):
        if char == r[0]:return r[1]

    # apply all rules to char 
    # if any rule is applied, return the new character
    def _rules(self,char):
        if char in self.constants:return char
        else:
            for r in self.rules:
                rc = self._rule(r,char)
                if rc:return rc

    # apply the rules to state
    def _iterate(self,instate):
        outstate = ''.join([self._rules(s) for s in instate])
        print('iteration',instate,'to',outstate)
        return outstate

    # iterate over axiom applying rules and return
    def _produce(self):
        state = self.axiom[:]
        for i in range(self.iterations):
            state = self._iterate(state)
        print('final product:',state)
        return state

    def _link(self):return 1.0
    def _push(self,p,d):pass
    def _pop(self,p,d):pass
    def _extend(self,p,d,l):
        n = p.copy().translate(d.copy().scale_u(l))
        draw_segment(p,n)
        return n
    def _terminate(self,p):
        draw_leaf(p)

    '''
    lsystem defines a grammer implicitly with _realize:
      0 : extend and terminate
      1 : extend
      [ : push - push to stack and perform some operation
      ] : pop - pop from stack and perform some operation
    '''
    # given a start position and start direction, realize this lsystem
    # return the locus of points in the order generated...
    def _realize(self,p,d):
        lstate = self._produce()
        output = [p.copy()]
        link,push,pop,extd,term = self.grammer_operations

        stack = []
        for ls in lstate:
            if ls == '0':term(extd(p,d,link()))
            elif ls == '1':p = extd(p,d,link())
            elif ls == '[':
                stack.append((p.copy(),d.copy()))
                push(p,d)
            elif ls == ']':
                p,d = stack.pop(-1)
                pop(p,d)
        draw()

class pythagoras_tree(lsystem):

    def __init__(self,*args,**kwargs):
        self._def('variables',['1','0'],**kwargs)
        self._def('constants',['[',']'],**kwargs)
        self._def('axiom','0',**kwargs)
        self._def('rules',[('1','11'),('0','1[0]0')],**kwargs)
        lsystem.__init__(self,*args,**kwargs)
        self._def('angle',numpy.pi/4.0,**kwargs)

        #####
        self._def('tropism',45,**kwargs)#direction of drift for F opertor (gravity?)
        self._def('tropism_size',5,**kwargs)#size of drift for F opertor (gravity?)

    def _push(self,p,d):
        #self.lastangle = (2*random.random()-1.0)*self.angle
        self.lastangle = self.angle
        q1 = dpq.q_from_av(self.lastangle,dpv.zhat)
        d.rotate(q1)
    def _pop(self,p,d):
        q2 = dpq.q_from_av(self.lastangle,dpv.nzhat)
        d.rotate(q2)

class tree(pythagoras_tree):

    def __init__(self,*args,**kwargs):
        pythagoras_tree.__init__(self,*args,**kwargs)
        self._def('angle2',numpy.pi,**kwargs)

    def _push(self,p,d):
        #self.lastangle = (2*random.random()-1.0)*self.angle
        q1 = dpq.q_from_av(self.angle, dpv.xhat)
        q2 = dpq.q_from_av(self.angle2,dpv.zhat)
        d.rotate(q1)
        d.rotate(q2)
    def _pop(self,p,d):
        q1 = dpq.q_from_av(self.angle, dpv.nxhat)
        q2 = dpq.q_from_av(self.angle2,dpv.nzhat)
        d.rotate(q2)
        d.rotate(q1)

draws = []
def draw_segment(p,n):
    draws.append((([p.x,n.x],[p.y,n.y]),{'zs':[p.z,n.z]}))

def draw_leaf(p):
    draws.append((([p.x],[p.y]),{'zs':[p.z],'marker':'o'}))

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

    ls = tree()
    ls._realize(p,d)



