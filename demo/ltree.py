#!/usr/bin/python3
import dilap.core.qtgui as qg
import dilap.core.lsystem as lsy

import dilap.geometry.tools as gtl
from dilap.geometry.vec3 import vec3

import dilap.core.plotting as dtl
import matplotlib.pyplot as plt

import pdb

class lsyswindow(qg.mwindow):

    def __init__(self,**kws):
        self.lsykws = {
            'axiom' : 'Q',
            'rules' : [('Q','<FF[)}+FQ][){Q]')],
            'seed' : 0,
            'iterations' : 4,
            #'polar' : gtl.rad(10),
            #'azimuthal' : gtl.rad(55),
            'polar' : 10,
            'azimuthal' : 55,
                }
        #self.lsys = lsy.lsystem(**self.lsykws)
        kws['title'] = 'LSystem Window'
        qg.mwindow.__init__(self,**kws)

    def pcb(self,ax):
        print('plot callback!',len([r for r in self.lsykws['rules'] if not r is None]))
        #self.lsykws['iterations'] -= 1
        self.lsykws['polar'] = gtl.rad(self.lsykws['polar'])
        self.lsykws['azimuthal'] = gtl.rad(self.lsykws['azimuthal'])
        self.lsys = lsy.lsystem(**self.lsykws)
        self.lsys._realize(vec3(0,0,0),vec3(0,0,1),ax)
        self.lsykws['polar'] = gtl.deg(self.lsykws['polar'])
        self.lsykws['azimuthal'] = gtl.deg(self.lsykws['azimuthal'])
        return ax

    def lsycb(self,k,f):
        def cb(v):
            self.lsykws[k] = f(v)
            self.pltwidget.update()
        return cb

    def lsyslider(self,k,f,minv,maxv,step,t,o):
        iv = self.lsykws[k]
        s = qg.slider(iv,self.lsycb(k,f),minv,maxv,step,t,o)
        return s
  
    def lsyspin(self,k,f,minv,maxv,step,t,o):
        iv = self.lsykws[k]
        s = qg.spin(minv,maxv,iv,step,self.lsycb(k,f),t,o)
        return s
  
    def lsytextbox(self,k,t):
        iv = self.lsykws[k]
        tb = qg.textbox(iv,self.lsycb(k,str),t)
        return tb

    def arule(self,i = '',f = ''):
        rcnt = len(self.lsykws['rules'])
        self.lsykws['rules'].append((i,f))
        def cbl(v):
            r = self.lsykws['rules'][rcnt]
            self.lsykws['rules'][rcnt] = (str(v),r[1])
            print(self.lsykws['rules'][rcnt])
            self.pltwidget.update()
        def cbr(v):
            r = self.lsykws['rules'][rcnt]
            self.lsykws['rules'][rcnt] = (r[0],str(v))
            print(self.lsykws['rules'][rcnt])
            self.pltwidget.update()
        def cbb():
            self.lsykws['rules'][rcnt] = None
            print([r for r in self.lsykws['rules'] if not r is None])
            print(self.lsykws['rules'],rcnt)
            self.rwidgs[rcnt].hide()
            self.pltwidget.update()
        rb = qg.buttons((cbb,),('clicked',),('Remove',),None,None,'')
        tbl = qg.textbox(i,cbl,'')
        tbr = qg.textbox(f,cbr,'')
        rw = qg.mwidget(qg.layout((tbl,tbr,rb),'h'),'')
        self.rulewidget.awidg(rw)
        self.rwidgs.append(rw)
        return self

    def adefrule(self):
        self.arule('F','FQ')
        self.pltwidget.update()
        return self

    def content(self,**kws):
        qg.init_figure()
        kws['plot_callback'] = self.pcb
        self.pltwidget = qg.pltwidget(self,**kws)

        si = self.lsyspin('iterations',int,1,10,1,'Iterations','v')
        ss = self.lsyspin('seed',int,0,10000,1,'Seed','v')
        sp = self.lsyspin('polar',float,0,180,1,'Polar Angle','v')
        sa = self.lsyspin('azimuthal',float,0,360,1,'Azimuthal Angle','v')
        ax = self.lsytextbox('axiom','Axiom')
        
        rules = self.lsykws['rules']
        self.rcnt = len(rules)
        self.lsykws['rules'] = []
        self.rwidgs = []
        self.rulewidget = qg.mwidget(qg.layout((),'v'))
        arb = qg.buttons((self.adefrule,),('clicked',),('Add Rule',),None,None,'')
        ub = qg.buttons((self.pltwidget.update,),('clicked',),('Update',),None,None,'')
        self.rulewidget.awidg(arb)
        for r in rules:self.arule(*r)
        self.cntrlwidget = qg.mwidget(
            qg.layout((ub,ax,self.rulewidget,sp,sa,si,ss,),'v'))
        content = qg.layout((self.cntrlwidget,self.pltwidget,),'h')
        self.pltwidget.update()
        return content

###############################################################################

if __name__ == '__main__':
    kws = {}
    qg.runapp(lsyswindow,kws)

###############################################################################





