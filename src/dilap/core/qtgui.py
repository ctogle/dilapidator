#import modular4.base as mb
import dilap.geometry.tools as gtl
import dilap.core.base as mb
import dilap.core.plotting as dtl

import os,sys,numpy,matplotlib,six,multiprocessing
matplotlib.rcParams['backend.qt4'] = 'PySide'
matplotlib.rcParams['pdf.fonttype'] = 42
#matplotlib.use('Qt4Agg')
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as figure_canvas
from matplotlib.backend_bases import NavigationToolbar2
import matplotlib.pyplot as plt

from PySide import QtGui,QtCore

import pdb





###############################################################################
### utility functions
###############################################################################

def runapp(windowclass,kws):
    print(5,'opening window',windowclass)
    mapp(windowclass,**kws).exec_()
    print(5,'closed window',windowclass)

def convert_pixel_space(w,h):
    good_w,good_h = 1920.0,1080.0
    screensize = QtGui.QApplication.desktop().availableGeometry()
    runtime_w,runtime_h = screensize.width(),screensize.height()
    w_conversion = runtime_w / good_w
    h_conversion = runtime_h / good_h
    w_size,h_size = w * w_conversion,h * h_conversion
    return w_size,h_size

def bound(label,lay):
    group = QtGui.QGroupBox(title = label)
    group.setLayout(lay)
    return layout((group,))

def sortto(l1,l2):
    new = []
    for l in l1:
        if not l in l2:continue
        else:new.append(l)
    for x in range(len(l2)):l2.pop(0)
    for n in new:l2.append(n)

###############################################################################
### functions to create layouts/widgets
###############################################################################

def layout(widgets = [],orientation = 'v'):
    if orientation == 'g':
        l = QtGui.QGridLayout()
        l.setSpacing(0.1)
        for w,p,s in widgets:
            wgs = (w,)+p+s
            l.addWidget(*wgs)
    else:
        if orientation == 'v':l = QtGui.QVBoxLayout()
        elif orientation == 'h':l = QtGui.QHBoxLayout()
        for w in widgets:l.addWidget(w)
    return l

def splitter(widgets = [],orientation = 'v',boxlabel = 'mwidgetsplitter'):
    if orientation == 'v':o = QtCore.Qt.Vertical
    elif orientation == 'h':o = QtCore.Qt.Horizontal
    split = QtGui.QSplitter(o)
    for w in widgets:split.addWidget(w)
    #if not callback is None:split.splitterMoved.connect(f)
    return bound(boxlabel,layout((split,)))

def buttons(funcs,events,labels,fw = None,fh = None,
            boxlabel = 'mwidgetbuttons',ori = 'v'):
    bs = []
    for bx in range(len(funcs)):
        f,e,l = funcs[bx],events[bx],labels[bx]
        b = QtGui.QPushButton(l)
        if not fw is None:b.setFixedWidth(fw)
        if not fh is None:b.setFixedHeight(fh)
        b.__getattribute__(e).connect(f)
        bs.append(b)
    return mwidget(layout(bs,ori),boxlabel)

def check(label,initial,callback,boxlabel = 'mwidgetcheck',ori = 'v'):
    '''
    create a widget containing a single check box which calls a function when toggled
    '''
    c = QtGui.QCheckBox(label)
    if initial:c.setCheckState(QtCore.Qt.Checked)
    else:c.setCheckState(QtCore.Qt.Unchecked)
    togg = lambda i : callback(False if i == 0 else True)
    c.stateChanged.connect(togg)
    return mwidget(layout((c,),ori),boxlabel)

def checks(tlist,labels,master = True,callback = None,boxlabel = None,ori = 'v'):
    '''
    create a widget containing a set of check boxes which add/remove items from a list
    '''
    qck,quck = QtCore.Qt.CheckState.Checked,QtCore.Qt.CheckState.Unchecked
    def togg(c,t):
        if not t in tlist:
            tlist.append(t)
            c.setCheckState(qck)
        elif t in tlist:
            tlist.remove(t)
            c.setCheckState(quck)
        sortto(tlisto,tlist)
        if tlisto == tlist:m.setCheckState(qck)
        else:m.setCheckState(quck)
    def flipall():
        def f():
            s = m.checkState()
            for lx in range(len(labels)):
                c,t = cs[lx+1],tlisto[lx]
                if not c.checkState() is s:
                    togg(c,t)
        return f
    def toggle(c,t):
        def f():togg(c,t)
        return f
    tlisto = tlist[:]
    if labels is tlist:labels = tlist[:]
    cs = [QtGui.QCheckBox(l) for l in labels]
    for c,l in zip(cs,labels):
        c.setCheckState(qck if l in tlist else quck)
        c.clicked.connect(toggle(c,l))
        if callback:c.clicked.connect(callback)
    if master:
        m = QtGui.QCheckBox('All')
        m.setCheckState(qck)
        for l in labels:
            if not l in tlist:
                m.setCheckState(quck)
                break
        m.clicked.connect(flipall())
        if callback:m.clicked.connect(callback)
        cs.insert(0,m)
    if type(boxlabel) == type(''):
        return mwidget(layout(cs,ori),boxlabel)
    else:return cs

def selector(labels,initial,callback,boxlabel = 'mwidgetselector'):
    def pick():
        c = sel.currentIndex()
        callback(lcopy[c])
    lcopy = labels[:]
    sel = QtGui.QComboBox()
    for l in labels:sel.addItem(l)
    sel.currentIndexChanged.connect(pick)
    sel.setCurrentIndex(labels.index(initial))
    return mwidget(layout((sel,)),boxlabel)

def radios(labels,initial,callback,boxlabel = 'mwidgetradios',ori = 'v'):
    def pick(x):
        f = lambda : callback(lcopy[x])
        return f
    lcopy = labels[:]
    rads = [QtGui.QRadioButton(l) for l in labels]
    rads[labels.index(initial)].setChecked(True)
    for rx in range(len(rads)):rads[rx].clicked.connect(pick(rx))
    return mwidget(layout(rads,ori),boxlabel)

def spin(minv,maxv,init = None,step = None,callback = None,boxlabel = None,ori = 'v'): 
    bind = lambda : callback(spin.value())
    spin = QtGui.QDoubleSpinBox()
    spin.setDecimals(10)
    spin.setMinimum(minv)
    spin.setMaximum(maxv)
    if not step is None:spin.setSingleStep(step)
    if not init is None:spin.setValue(init)
    if not callback is None:spin.valueChanged.connect(bind)
    if type(boxlabel) == type(''):
        return mwidget(layout((spin,)),boxlabel)
    else:return spin

def textbox(initial,callback = None,boxlabel = 'mwidgettextbox',ori = 'v'):
    def f():
        btx = str(box.text())
        callback(btx)
    box = QtGui.QLineEdit()
    box.setText(str(initial))
    box.editingFinished.connect(f)
    return mwidget(layout((box,),ori),boxlabel)

def slider(initial,callback = None,minv = 1,maxv = 10,intv = 1,
                        boxlabel = 'mwidgettextbox',ori = 'v'):
    o = QtCore.Qt.Vertical if ori == 'v' else QtCore.Qt.Horizontal
    def sf():callback(slide.value())
    slide = QtGui.QSlider(o)
    slide.setMinimum(minv)
    slide.setMaximum(maxv)
    slide.setSingleStep(intv)
    slide.setTickPosition(slide.TicksBothSides)
    slide.setTickInterval(intv)
    slide.setValue(initial)
    slide.sliderReleased.connect(sf)
    slide.setGeometry(30,40,100,30)
    return mwidget(layout((slide,),ori),boxlabel)

###############################################################################
### classes useful for making applications
###############################################################################

class mapp(QtGui.QApplication):

    def __init__(self,main_window_class = None,**kws):
        QtGui.QApplication.__init__(self,sys.argv)
        QtGui.QApplication.setStyle(QtGui.QStyleFactory.create('Plastique'))
        if main_window_class is None:main_window_class = mwindow
        self.main_window = main_window_class(**kws)

class mwindow(QtGui.QMainWindow):

    def _standards(self,**st):
        if 'title' in st:wt = st['title']
        else:wt = 'mwindow'
        if 'geometry' in st:geo = st['geometry']
        else:
            x,y = convert_pixel_space(300,300)
            x_size,y_size = convert_pixel_space(512,512)
            geo = (x,y,x_size,y_size)
        gearicon = QtGui.QIcon(mb.resource_path('gear.png'))
        self.setWindowIcon(gearicon)
        self.setWindowTitle(wt)
        self.setGeometry(*geo)

    def __init__(self,**kws):
        QtGui.QMainWindow.__init__(self)
        self._standards(**kws)
        w = QtGui.QWidget()
        w.setLayout(self.content(**kws))
        self.setCentralWidget(w)
        self.show()

    def content(self,**kws):
        content = QtGui.QVBoxLayout()
        return content

class mwidget(QtGui.QWidget,mb.base):

    def awidg(self,w):
        self.lay.addWidget(w)
        self.update()

    def __init__(self,lay = None,lab = '',scroll = False,**kws):
        QtGui.QWidget.__init__(self)
        if not lay is None:
            if lab:lay = bound(lab,lay)
            if scroll:
                scroll = QtGui.QScrollArea()
                scroll.setBackgroundRole(QtGui.QPalette.Window)
                scroll.setWidget(mwidget(lay))
                self.lay = layout((scroll,))
            else:self.lay = lay
            self.setLayout(self.lay)

###############################################################################
### custom widget classes
###############################################################################

figure = None
def init_figure():
    global figure
    figure = plt.figure()

class pltwidget(mwidget):

    def __init__(self,parent,**kws):
        mwidget.__init__(self)
        self.parent = parent
        self.fig = figure
        self.canvas = figure_canvas(self.fig)
        self.toolbar = plttoolbar(self.canvas)
        self.ax = None
        self._def('plot_callback',None,**kws)
        self.setLayout(layout([self.canvas,self.toolbar],'v'))
        self.setBackgroundRole(QtGui.QPalette.Window)

    def show(self):
        mwidget.show(self.update())
        return self

    def update(self):
        ax = self.clear_ax()
        if self.plot_callback:
            self.plot_callback(ax)
        self.canvas.draw()
        mwidget.update(self)
        return self

    def clear_ax(self,proj = '3d'):
        #self.fig.clf()
        #ax = dtl.plot_axes(25,self.fig)
        ax = self.fig.gca(projection = proj)
        x = 25
        ax.set_xlim([-x,x])
        ax.set_ylim([-x,x])
        ax.set_zlim([-(9.0/16.0)*x,(9.0/16.0)*x])
        ax.cla()
        ax.grid(False)
        return ax
    
class plttoolbar(NavigationToolbar2,QtGui.QToolBar):

    message = QtCore.Signal(str)
    if hasattr(NavigationToolbar2,'toolitems'):
        titems = NavigationToolbar2.toolitems 
        toolitems = [t for t in titems if t[0] in ('Pan','Zoom','Save')]
    else:toolitems = []

    def dynamic_update(self):
        self.canvas.draw()

    def pan(self,*ags):
        super(plttoolbar,self).pan(*ags)
        self._update_buttons_checked()

    def zoom(self,*ags):
        super(plttoolbar,self).zoom(*ags)
        self._update_buttons_checked()

    def _update_buttons_checked(self):
        self._actions['pan'].setChecked(self._active == 'PAN')
        self._actions['zoom'].setChecked(self._active == 'ZOOM')

    def _init_toolbar(self):
        for text,tooltip_text,image_file,callback in self.toolitems:
            if text is None:self.addSeparator()
            else:
                i = QtGui.QIcon()
                a = self.addAction(i,text,getattr(self,callback))
                self._actions[callback] = a
            if callback in ('zoom','pan'):a.setCheckable(True)
            if tooltip_text is not None:a.setToolTip(tooltip_text)
        self.locLabel = QtGui.QLabel("", self)
        self.locLabel.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignTop)
        self.locLabel.setSizePolicy(QtGui.QSizePolicy(
            QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Ignored))
        labelAction = self.addWidget(self.locLabel)
        labelAction.setVisible(True)

    def __init__(self,canvas):
        self.canvas = canvas
        self.img_extensions = 'Image (*.png,*.pdf)'
        self._actions = {}
        QtGui.QToolBar.__init__(self)
        NavigationToolbar2.__init__(self,canvas)

    def draw_rubberband(self,event,x0,y0,x1,y1):
        height = self.canvas.figure.bbox.height
        y1 = height - y1
        y0 = height - y0
        w = abs(x1 - x0)
        h = abs(y1 - y0)
        rect = [int(val)for val in (min(x0,x1),min(y0,y1),w,h)]
        self.canvas.drawRectangle(rect)

    def save_figure(self,*ags):
        fname = QtGui.QFileDialog.getSaveFileName(self,
            'Choose Filename','aplot.pdf',self.img_extensions) 
        if fname:
            try:
                self.canvas.print_figure(six.text_type(fname[0]))
                print(5,'saved figure at',fname)
            except Exception as e:
                QtGui.QMessageBox.critical(
                    self,'error saving file',str(e),
                    QtGui.QMessageBox.Ok,QtGui.QMessageBox.NoButton)






usemp = False
def displaycontext(cx,**kws):
    if usemp:
        p = multiprocessing.Process(target = runapp,args = (cwindow,kws))
        p.start()
        return p
    else:
        kws['context'] = cx
        runapp(cwindow,kws)
        return None



'''#
class plttree_book(mwidget):

    def set_page(self,pgx):
        self.tree_pages[self.page].hide()
        self.tree_pages[pgx].show()
        self.page = pgx

    def _change_page(self,cpg,ppg):
        for ix in range(len(self.tree_items)):
            if self.tree_items[ix] is cpg:
                self.set_page(ix)
                return

    def _header(self,header):
        if not type(header) is type(''):header = ''
        self.tree.setHeaderLabel(header)

    def _pages(self,pages):
        self.tree.setColumnCount(1)
        self.tree.clear()
        self.aux = {}
        self._targets = []
        self._targetlabels = []
        titems,tpages,tops,bottoms = [],[],[],[]
        for x in range(len(pages)):
            if not pages[x]:continue
            pgd,pgt,pge = pages[x]
            if 'header' in pge:h = str(pge['header'])
            else:h = ''
            toplabel = 'pspace page: %i : %s' % (x,h)
            top = QtGui.QTreeWidgetItem(None,[toplabel])
            tops.append(top)
            main_page = mwidget()
            titems.append(top)
            tpages.append(main_page)
            if len(pgd.shape) == 2:
                subs = (((pgd,pgt),'single'),)
            elif len(pgd.shape) == 3:
                subs = tuple(((pgd[x],pgt),'trajectory: %i' % x) for x in range(pgd.shape[0]))
            else:
                print(5,'unknown tree widget scenario')
                raise ValueError
            for subpg,subh in subs:
                for t in subpg[1]:
                    if not t in self._targets:
                        self._targets.append(t)
                        self._targetlabels.append(t)
                bottom = QtGui.QTreeWidgetItem(top,[subh])
                bottoms.append(bottom)
                titems.append(bottom)
                t1 = subpg[1][0]

                if 'extra_trajectory' in pge:
                    extras = pge['extra_trajectory']
                else:extras = None
                if 'pspaceaxes' in pge:
                    self.aux['pspaceaxes'] = pge['pspaceaxes']
                self.aux['entry'] = subpg

                sub_page = mpltwidget(self,subpg,extras = extras)
                tpages.append(sub_page)
            self.tree.addTopLevelItem(top)
        if self.xdomain is None:self.xdomain = t1
        if self.ydomain is None:self.ydomain = t1
        if self.zdomain is None:self.zdomain = t1
        if self.linestyles is None:self.linestyles = ['-' for t in self._targets]
        if self.linewidths is None:self.linewidths = [1 for t in self._targets]
        if self.linemarkers is None:self.linemarkers = ['' for t in self._targets]
        if self.linecolors is None:
            linsp = numpy.linspace(0,0.9,len(self._targets))
            self.linecolors = [self.colormap(i) for i in linsp]
        for page in tpages:
            self.hsplit.addWidget(page)
            page.hide()
        self.tree_items = titems
        self.tree_pages = tpages
        self.tree_tops = tops
        self.tree_bottoms = bottoms
        self.set_page(self.page)

    def _axisslice_widgets(self):
        def slice_callback(a):
            def f(c):
                v = float(c)
                nearest = [abs(x-v) for x in self.axisvalues[a]]
                nearest = nearest.index(min(nearest))
                self.axisdefaults[a] = self.axisvalues[a][nearest]
                self.update()
            return f
        axslices = []
        for axx in range(len(self.axisnames)):
            ls = tuple(str(v) for v in mb.uniq(self.axisvalues[axx]))
            si = str(self.axisdefaults[axx])
            cb = slice_callback(axx)
            axslices.append(selector(ls,si,cb,boxlabel = self.axisnames[axx]))
        return mwidget(layout(axslices),'Parameter Space Axes')

    def _domain_widgets(self,dom):
        i = bool(self.__getattribute__(dom+'log'))
        lab = textbox(self.xlabel,self._defbind(dom+'label'),dom+'-label')
        sel = selector(self._targets,self._targets[0],
            self._defbind(dom+'domain'),dom+'-domain')
        lg = check('Use log('+str(dom)+')',i,self._defbind(dom+'log'),'')
        i = bool(self.__getattribute__(dom+'bnd'))
        usebnd = check('Use Bounds',i,self._defbind(dom+'bnd'),'')
        lowbnd = textbox(
            str(self.xmin),self._defbind(dom+'min'),
            boxlabel = dom+'-minimum',ori = 'v')
        highbnd = textbox(
            str(self.xmax),self._defbind(dom+'max'),
            boxlabel = dom+'-maximum',ori = 'v')
        sel = mwidget(layout((
            mwidget(layout((sel,lab),'h')), 
            mwidget(layout((lg,usebnd,lowbnd,highbnd),'h'))),'v'),dom+'-axis')
        sel.setFixedWidth(self._panelwidth)
        return sel

    def _target_widgets(self):
        def tnbind(x):
            def f(tname):
                self._targetlabels.pop(x)
                self._targetlabels.insert(x,tname)
                self.update()
            return f
        def lwbind(x):
            def f(lw):
                self.linewidths[x] = int(lw)
                self.update()
            return f
        def lsbind(x):
            def f(ls):
                self.linestyles[x] = str(ls)
                self.update()
            return f
        def lmbind(x):
            def f(lm):
                self.linemarkers[x] = str(lm)
                self.update()
            return f
        def clbind(x):
            def f():
                col = QtGui.QColorDialog.getColor()
                if col.isValid():self.linecolors[x] = col.getRgbF()
                self.update()
            return f
        tcs = checks(self._targets,self._targets,True,self.update,None)
        lwidgs = []
        for tx in range(len(self._targets)):
            tnamebox = textbox(self._targets[tx],tnbind(tx),None)
            lws = [str(x) for x in range(10)]
            lwsel = selector(lws,str(self.linewidths[tx]),lwbind(tx),None)
            lss = ['','-','--','-.',':']
            lssel = selector(lss,str(self.linestyles[tx]),lsbind(tx),None)
            lms = ['','o','v','^','<','>','s','p','*','h','H','D','d','x','+']
            lmsel = selector(lms,str(self.linemarkers[tx]),lmbind(tx),None)
            lcbtn = buttons((clbind(tx),),('clicked',),('Col',),30,30,None)
            lwidgs.append((tcs[tx+1],(tx+1,0),(1,1)))
            lwidgs.append((tnamebox,(tx+1,1),(1,1)))
            lwidgs.append((lwsel,(tx+1,2),(1,1)))
            lwidgs.append((lssel,(tx+1,3),(1,1)))
            lwidgs.append((lmsel,(tx+1,4),(1,1)))
            lwidgs.append((lcbtn,(tx+1,5),(1,1)))
        lwidgs.insert(0,(tcs[0],(0,0),(1,1)))
        sls = mwidget(layout(lwidgs,'g'),'Plot Targets')
        return sls

    def _domain_target_ptype_widgets(self):
        plab = textbox(self.plottitle,self._defbind('plottitle'),'Plot Title')
        pleg = check('Show Legend',self.legend,self._defbind('legend'),'')
        popt = mwidget(layout((plab,pleg),'v'),'')
        plab.setFixedWidth(self._panelwidth)
        xaxis = self._domain_widgets('x')
        yaxis = self._domain_widgets('y')
        zaxis = self._domain_widgets('z')
        rds = radios(self.plottypes,self.plottype,self._defbind('plottype'),'Plot Type')
        tcs = self._target_widgets()
        if self.axisnames:axs = self._axisslice_widgets()
        else:axs = mwidget()
        bot = mwidget(layout((rds,axs),'h'),'')
        return mwidget(splitter((popt,xaxis,yaxis,zaxis,tcs,bot),'v',''),'Plot Filter',True)

    def _widgets(self):
        self.vsplit = QtGui.QSplitter(QtCore.Qt.Vertical)
        self.hsplit = QtGui.QSplitter(QtCore.Qt.Horizontal)
        self.tree = QtGui.QTreeWidget()
        self.vsplit.addWidget(self.tree)
        self.hsplit.addWidget(self.vsplit)
        self.tree.currentItemChanged.connect(self._change_page)
        self._header(self.header)
        self._pages(self.pages)
        self._set_axis_info()
        updatebutton = buttons((self.update,),('clicked',),('Update Plot',),600,100,'')
        self.plt_controls = self._domain_target_ptype_widgets()
        self.vsplit.addWidget(updatebutton)
        self.vsplit.addWidget(self.plt_controls)
        return (self.hsplit,)

    def _set_axis_info(self):
        if 'pspaceaxes' in self.aux and self.axisdefaults is None:
            self.reduce_lines = True
            self.axisnames = self.aux['pspaceaxes']
            d,t = self.aux['entry']
            if not t == self._targets:
                print('probably a serious problem!')
                pdb.set_trace()
            axxs = tuple(t.index(a) for a in self.axisnames)
            self.axisvalues = [d[a] for a in axxs]
            self.axisdefaults = [vs[0] for vs in self.axisvalues] 

    def _defbind(self,k):
        def f(c):
            self.__setattr__(k,c)
            print(5,'set pltwidget attribute: %s : %s' % (k,str(c)))
            self.update()
        return f

    def update(self):
        self.tree_pages[self.page].update()
        mwidget.update(self)

    def calc_lines_callback(self,pwidg,ax,d,t,x,ys):
        print('calc_lines_callback!')
        #ax.plot([500,500],[-1,100],linewidth = 5.0,marker = 'o',color = 'b')
        return ax

    def calc_color_callback(self,pgwidg,ax,d,t,x,y,z):
        print('calc_color_callback!')
        return ax

    _panelwidth = 500

    def __init__(self,**kws):
        mwidget.__init__(self,**kws)
        self.kws = kws
        self._def('line_callbacks',[],**kws)
        self._def('pages',[],**kws)
        self._def('page',0,**kws)
        self._def('header','Data Selection',**kws)
        self._def('_targets',[],**kws)
        self._def('xdomain',None,**kws)
        self._def('ydomain',None,**kws)
        self._def('zdomain',None,**kws)
        self._def('linestyles',None,**kws)
        self._def('linewidths',None,**kws)
        self._def('linemarkers',None,**kws)
        self._def('linecolors',None,**kws)
        self._def('xlabel','',**kws)
        self._def('xlabelsize',20,**kws)
        self._def('ylabel','',**kws)
        self._def('ylabelsize',20,**kws)
        self._def('zlabel','',**kws)
        self._def('zlabelsize',20,**kws)
        self._def('plottitle','',**kws)
        self._def('plottitlesize',18,**kws)
        self._def('legend',True,**kws)
        self._def('xlog',False,**kws)
        self._def('ylog',False,**kws)
        self._def('zlog',False,**kws)
        self._def('xbnd',False,**kws)
        self._def('xmin','',**kws)
        self._def('xmax','',**kws)
        self._def('xticksize',20,**kws)
        self._def('ybnd',False,**kws)
        self._def('ymin','',**kws)
        self._def('ymax','',**kws)
        self._def('yticksize',20,**kws)
        self._def('zbnd',False,**kws)
        self._def('zmin','',**kws)
        self._def('zmax','',**kws)
        self._def('zlabsize',20,**kws)
        self._def('axisnames',[],**kws)
        self._def('axisvalues',[],**kws)
        self._def('axisdefaults',None,**kws)
        self._def('reduce_lines',False,**kws)
        self._def('maxlinecount',20,**kws)
        self._def('colorplot_interpolation','nearest',**kws)
        self._def('colormap',plt.get_cmap('jet'),**kws)
        self._def('plottypes',('lines','color'),**kws)
        self._def('plottype','lines',**kws)
        wgs = self._widgets()
        self._layout = layout(wgs,'h')
        self.setLayout(self._layout)
        for top in self.tree_tops:self.tree.expandItem(top)
'''#




