import dilap.core.lsystem as lsy
import dilap.core.plotting as dtl
import dilap.core.context as cx
import dilap.geometry.tools as gtl
from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat
from dilap.geometry.pointset import pointset
import dilap.modeling.model as dmo
import dilap.topology.trimesh as dtm
import dilap.topology.tree as dtr
import dilap.topology.vert as dvt
import dilap.topology.edge as deg
import dilap.topology.loop as dlp
import dilap.topology.face as dfc

#import dilap.topology.tools.triangulate as dtg
import dilap.geometry.triangulate as dtg

import numpy
import pdb



###############################################################################
### ltree creates a set of points/edges for leaves/branches and has loadouts
###############################################################################

# an ltree is an lsystem that produces structured information about a tree
# such that it can produce a mesh from that information representing the tree
class ltree(lsy.lsystem):
                                                    
    def _realize(self,p,d,ax = None):
        self.branches = []
        self.leaves = []
        return lsy.lsystem._realize(self,p,d,ax)

    loadouts = []
    loadouts.append(('A',[
        ('F','F~[~FA]F'),
        ('A','[--////FQA]-<F[>>\\\\FQA]>+F[-\>>QA\[X]][</++QA/[Y]]Q'),
        ('X','F+Q'),('Y','F-Q')]))
    loadouts.append(('Q',[
        ('Q','<FF[)}+FQ][){Q]')]))

    def __init__(self,ldx,*args,**kwargs):
        axiom,rules = self.loadouts[ldx]
        self._def('axiom',axiom,**kwargs)
        self._def('rules',rules,**kwargs)

        self._def('seed',1,**kwargs) # seed used for random numbers
        self._def('iterations',4,**kwargs) # number of iterations
        self._def('angle',numpy.pi/12.0,**kwargs) # angle used for rotations
        self._def('minangle',numpy.pi/12.0,**kwargs) # angle used for rotations
        self._def('maxangle',numpy.pi/6,**kwargs) # angle used for rotations

        self._def('branchdraw',self.draw_branch,**kwargs)
        self._def('leafdraw',self.draw_leaf,**kwargs)
        self._def('finaldraw',self.draw,**kwargs)
        lsy.lsystem.__init__(self,*args,**kwargs)

    def draw_branch(self,p,n):
        self.branches.append((p.cp(),n.cp()))
        lsy.draw_branch(p,n)

    def draw_leaf(self,p):
        self.leaves.append((p.cp()))
        lsy.draw_leaf(p)

    def draw(self,ax = None):
        es,ls = self.branches,self.leaves
        nps,nes = [],[]
        for ex in range(len(es)):
            e1,e2 = es[ex]
            e1x,e2x = None,None
            for pdx in range(len(nps)):
                if e1.isnear(nps[pdx]):e1x = pdx
                if e2.isnear(nps[pdx]):e2x = pdx
                if not e1x is None and not e2x is None:break
            if e1x is None:
                e1x = len(nps)
                nps.append(e1)
                tr = dtr.tree()
                tr.root.MARK = (e1x,e1)
            if e2x is None:
                e2x = len(nps)
                nps.append(e2)
                nvrt = tr.avert(tr.verts[e1x])
                nvrt.MARK = (e2x,e2)
            nes.append((e1x,e2x))
        lsy.draw(ax)
        self.topology = (tr,nps,nes)

###############################################################################
###############################################################################
###############################################################################

class knuckle:
                                                    
    def __init__(self,base,corner,l,ml,mw,**kwargs):
        self.base = base
        self.corner = corner
        self.baseprops = (l,ml,mw)
        self.chprops = []

    def digit(self,c,l,ml,mw):
        self.chprops.append((c,(l,ml,mw)))

    def tribridge(self,m,gm,bt,tp,fm):
        for x in range(len(bt)):
            p1,p2,p3,p4 = bt[x-1],bt[x],tp[x],tp[x-1]
            v1 = gm.avert(*m.avert(p1.cp()))
            v2 = gm.avert(*m.avert(p2.cp()))
            v3 = gm.avert(*m.avert(p3.cp()))
            v4 = gm.avert(*m.avert(p4.cp()))
            f1  = gm.aface(v1,v2,v3,fm) 
            f2  = gm.aface(v1,v3,v4,fm) 

    def calc_w(self,x = None):
        if x is None:chp = self.baseprops
        else:chp = self.chprops[x][1]
        return 0.1*max((chp[2]-(chp[0]/chp[1])**(0.25)),0.2)

    def gen(self,m,gm,r = False):
        if r:bt = vec3(0,0,-1).pring(self.calc_w(),8)
        else:bt = self.base.MARK[1].pring(self.calc_w(),8)
        mid = self.corner.MARK[1].pring(self.calc_w(),8)
        tps = [chp[0].MARK[1].pring(self.calc_w(x),8) 
                  for x,chp in enumerate(self.chprops)]
        u2s = [self.corner.MARK[1].tov(chp[0].MARK[1]) for chp in self.chprops]
        u1 = self.base.MARK[1].tov(self.corner.MARK[1])
        for u2,tp in zip(u2s,tps):
            if u1.mag() == 0:q1 = quat(0,0,0,1)
            else:q1 = quat(0,0,0,1).uu(vec3(0,0,1),u1)
            q2 = quat(0,0,0,1).uu(vec3(0,0,1),u2)
            vec3(0,0,0).com(mid).fulc(q1,mid)
            vec3(0,0,0).com(tp).fulc(q2,tp)
            self.tribridge(m,gm,mid,tp,'concrete1')
            if r:self.tribridge(m,gm,bt,mid,'concrete1')
            vec3(0,0,0).com(mid).fulc(q1.flp(),mid)

###############################################################################
###############################################################################
###############################################################################

class tree(cx.context):

    def __init__(self,p,d,*args,ax = None,**kwargs):
        self._def('name','treecontext',**kwargs)
        self._def('params',(('iterations',(1,10,1)),),**kwargs)
        cx.context.__init__(self,*args,**kwargs)
        self.lsys = ltree(1)._realize(p,d,ax)

    # do something which fills the scenegraph
    def generate(self,worn = 0):
        print('generate with worn',worn)
        self.skeleton()
        return self

    def leaf(self,s,e):
        m = dmo.model()
        gm = m.atricube('generic')
        p = s.MARK[1].cp()
        tn = s.MARK[1].tov(e.MARK[1])
        q = quat(0,0,0,0).uu(vec3(0,0,1),tn)
        #w = 0.1*max((mw-(l/ml)**(0.25)),0.2)
        #s = vec3(w,w,0.5*tn.mag())
        #q = quat(0,0,0,0).av(0.0,vec3(0,0,1))
        s = vec3(0.2,0.2,0.5*tn.mag())
        m.trn(vec3(0,0,1))
        sgv = self.amodel(p,q,s,m,None)

    def branch(self,s,e,l,ml,mw):
        m = dmo.model()
        gm = m.atricube('generic')
        m.trn(vec3(0,0,1))
        p = s.MARK[1].cp()
        tn = s.MARK[1].tov(e.MARK[1])
        q = quat(0,0,0,0).uu(vec3(0,0,1),tn)
        w = 0.1*max((mw-(l/ml)**(0.25)),0.2)
        s = vec3(w,w,0.5*tn.mag())
        sgv = self.amodel(p,q,s,m,None)

    def skeleton(self):
        m = dmo.model()
        gm = m.agfxmesh()
        tr,nps,nes = self.lsys.topology
        mw = 2.0
        knuckles = []
        def walk(v,l,ml):
            ch = tr.below(v)
            if not ch:self.leaf(tr.above(v),v)
            #if not ch:self.branch(tr.above(v),v,l,ml,mw)
            else:
                kn = knuckle(tr.above(v),v,l,ml,mw)
                knuckles.append(kn)
                l += 1.0
                ml = len(tr.allbelow(v))
                for c in ch:
                    #self.branch(v,c,l,ml,mw)
                    kn.digit(c,l,ml,mw)
                    walk(c,l,ml)
                if len(knuckles) == 1:kn.gen(m,gm,True)
                else:kn.gen(m,gm,False)
        walk(tr.root,0,len(tr.allbelow(tr.root)))
        p,q,s = vec3(0,0,0),quat(0,0,0,0),vec3(1,1,1)
        sgv = self.amodel(p,q,s,m,None)
        return self

        '''#
        m = dmo.model()
        gm = m.atricube()
        p = vec3(0,0,1)
        q = quat(0,0,0,0).av(numpy.pi/3.0,vec3(0,0,1))
        #q = quat(0,0,0,0).av(0.0,vec3(0,0,1))
        s = vec3(1,2,1)
        self.amodel(p,q,s,m,None)
        '''#

        for v in nps:
            m = dmo.model()
            gm = m.atricube()
            p = v.cp()
            q = quat(0,0,0,0).av(0.0,vec3(0,0,1))
            s = vec3(0.1,0.1,0.2)
            sgv = self.amodel(p,q,s,m,None)
        for e in nes:
            m = dmo.model()
            gm = m.atricube()
            m.trn(vec3(0,0,1))
            p = nps[e[0]].cp()
            tn = nps[e[0]].tov(nps[e[1]])
            q = quat(0,0,0,0).uu(vec3(0,0,1),tn)
            s = vec3(0.05,0.05,0.5*tn.mag())
            sgv = self.amodel(p,q,s,m,None)
        return self





