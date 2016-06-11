import dilap.core.base as db
import dilap.core.context as cx
import dilap.modeling.model as dmo
import dilap.modeling.factory as dfa

from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat
import dilap.geometry.polymath as pym

import dilap.topology.planargraph as pgr
import dilap.topology.partitiongraph as ptg

import dilap.worldly.terrain as ter
import dilap.worldly.treeskin as ltr
import dilap.worldly.building as blg
import dilap.worldly.blgsequencing as bseq
#import dilap.worldly.partitiongraph as ptg
import dilap.worldly.roadgraph as rdg
import dilap.worldly.blockletters as dbl

import dilap.core.plotting as dtl
import matplotlib.pyplot as plt

import math,numpy,random,pdb

###############################################################################
### world factory generates full world contexts from a footprint
###############################################################################

class worldfactory(dfa.factory):

    # process of making a world:
    #   generate an xy boundary polygon for the boundary of the world
    #   generate topology in the form of a deterministic function with feedback
    #   partition the interior of the boundary into subsections with type information
    #
    #   place waterways?
    #
    #   create road networks based on region topology
    #       each region has a finite nonzero number of exits
    #           create a network of roads which respects these exits
    #   partition natural and developed regions with 
    #       corridor regions based on road graphs
    #   in developed regions, use the road system to place structures, 
    #       further partitioning the region into developed, buildings, and easements
    #   in natural regions, road graphs are simple and no buildings are present
    #       instead, vegetation and perhaps more extreme topology are added



    # generate the topographical data for the world
    def genregions(self,t):
        pg = ptg.partitiongraph()
        ### create the ocean vertex
        b = t.root.loop
        pv = pg.av(b = [b,[]],p = vec3(0,0,0).com(b),t = ['ocean'])
        ### create landmass vertices
        ls = t.looptree.below(t.root)
        lmvs = []
        for l in ls:
            pv = pg.sv(pv,b,l.loop)
            pg.vs[pv][1]['t'] = ['natural']
            lmvs.append(pv)
        ### create the details of each landmass vertex
        for lm in lmvs:self.genregion(lm,pg)
        ### attach the terrain mesh to the partition vertices...
        for vx in range(pg.vcnt):
            pv = pg.vs[vx]
            if pv is None:continue
            pv[1]['tmesh'] = t
        
        ###
        ax = pg.plotxy()
        plt.show()
        ###

        ### return partition graph
        return pg

    def genregion(self,v,pg):
            
        #last1 = v
        #new = pg.bv(last1,vec3(0.25,0.25,0),vec3(1,0,0))
        #pg.vs[last1][1]['t'] = ['developed']

        #last2 = new[0]
        #new = pg.bv(last2,vec3(0.75,0.75,0),vec3(0,1,0))
        #last3 = new[0]
        #new = pg.bv(last3,vec3(0.66,0.66,0),vec3(1,0,0))
        #last4 = new[0]
        #new = pg.bv(last3,vec3(0.33,0.33,0),vec3(0,1,0))
        #pg.vs[last3][1]['t'] = ['developed']

        ### add infrastructural connectivity for regions
        # create a road network graph...
        #pv = pg.av(b = fp,p = vec3(0,0,0).com(fp[0]),l = 0)
        return pg

    # generate the models of the world and add them to the context object
    # NOTE: this should eventually be parallelizable
    def gen(self,w,pg):
        random.seed(0)
        for vx in range(pg.vcnt):
            v = pg.vs[vx]
            if v is None:continue
            self.vgen(w,v,pg)

    ###########################################################################
    ### methods for generating models from descriptive world data
    ###########################################################################

    # vgen takes a context, vertex, and partitiongraph and from that
    #   it generates all models within the space of the vertex
    # this consists of two fundamental steps:
    #   based on the nature of this and adjacent vertices, 
    #       the boundary of this vertex is constrained, 
    #       and meshes which cross it must be forced to appear continuous
    #   based on the nature of this vertex, the interior of the space 
    #       associated with the vertex must be populated with models
    #        - these models must geometrically cover the entire ground surface
    # types:
    #   "natural" - cover with terrain vegetation only edges should only require a terrain seam
    #   "developed" - cover with structures, filling gaps with nature models
    #       edges should only require terrain seams and road seams
    #       where roads may intersect non-colinearly with the boundary only
    #   "corridor" - cover with infrastructure, filling gaps with concrete?
    #       edges should only require terrain seams and road seams
    #       road seems may include colinear intersections with the boundary
    # this supports the idea that all space is either infrastructure, structure, or nature
    def vgen(self,w,v,pg):
        l = 10
        if 't' in v[1]:
            vtypes = v[1]['t']
            if 'ocean' in vtypes:self.vgen_ocean(w,v,pg,l)
            elif 'natural' in vtypes:self.vgen_natural(w,v,pg,l)
            elif 'developed' in vtypes:self.vgen_developed(w,v,pg,l)
        else:self.vgen_ocean(w,v,pg,l)
        return
    
    def vgen_ocean(self,w,v,pg,l,depth = 2):
        print('ocean vertex',v[0])
        ### create an unrefined flat surface for the bottom of the ocean
        ###  and an unrefined flat surface with contracted holes for the 
        ###  surface of the ocean
        m = dmo.model()
        sgv = w.amodel(None,None,None,m,w.sgraph.root)
        tm = m.agfxmesh()
        gb = v[1]['b']
        wb = [[p.cp().ztrn(depth) for p in gb[0]],
            [pym.contract([p.cp().ztrn(depth) for p in ib],20) for ib in gb[1]]]
        ngvs = m.asurf(gb,tm,fm = 'concrete1',ref = False,hmin = 100,zfunc = None)
        ngvs = m.asurf(wb,tm,fm = 'concrete1',ref = False,hmin = 100,zfunc = None)

    def vgen_natural(self,w,v,pg,l):
        print('natural vertex',v[0])
        ### create terrain without holes
        vb = self.vstitch(v,pg,l)
        self.vgen_terrain(w,v,pg,[],l)

    def vgen_developed(self,w,v,pg,l):
        print('developed vertex',v[0])
        ### create a child context for a set of buildings within the vertex
        blgs,blgfps = self.vgen_buildings(v,pg,l)
        for blg in blgs:w.achild(blg.generate(0))
        ### create terrain with holes where buildings meet the ground
        self.vgen_terrain(w,v,pg,[[p.cpxy() for p in f] for f in blgfps],l)

    ###########################################################################

    def __str__(self):return 'world factory:'
    def __init__(self,*ags,**kws):
        self._def('bclass',cx.context,**kws)
        dfa.factory.__init__(self,*ags,**kws)

        s = random.randint(0,1000)
        print('landmass seed:',s)
        random.seed(s)

    # create a context representing an entire world
    def new(self,*ags,**kws):
        ### create the boundary of the world
        boundary = vec3(0,0,0).pring(500,8)
        ### generate the topographical structure of the world
        t = ter.continent(boundary)
        ### generate the region partitions of each landmass
        pg = self.genregions(t)
        ### create a world context from the partition graph
        w = self.bclass(*ags,**kws)
        self.gen(w,pg)
        return w


        '''#
        worldh = 50
        #worldseq = 'J<0,2>C<1>'

        toposeq = 'FF'

        contseq =  'R<0,natural>'
        contseq += 'S<0,0.5,0.5,0,1,vertex>'
        contseq += 'S<1,0.5,0.5,1,0,vertex>'
        contseq += 'S<0,0.5,0.5,1,0,vertex>'
        #contseq += 'E<0,1>E<1,2>E<2,3>E<3,0>'
        #contseq += 'M<0,S<>G<>L<>,>'
        contseq += 'T<2>'
        #contseq += 'M<0>T<2>'


        # USE SPLOTCH TO ESTABLISH A COASTLINE!!

        # WHAT IS THE EXACT INTERPLAY OF THE 
        # TERRAINMESH OBJECT AND THE CONTINENT SPLIT SEQUENCE....?


        worldseq = 'J<0,2,AV>V<0,3>R<1,natural>L<1,50,-1,'+toposeq+'>I<1,'+contseq+'>'

        eg = {
          'L':cartographer,
          #'C':continent, # partitions an island (coasts, forests, mountains, etc)
          'M':metropolis, # partition an area into a city...
          'T':stitch, # ensure the terrain seam of each region is functional
              }
        pg = ptg.checkseq(self.boundary,worldh,worldseq,show = True,grammer = eg)
        self.gen(w,pg)
        return w
        '''#





    # segment the boundary of a vertex based on edges of other vertices
    def vstitch(self,v,pg,l = 10):
        vb = [b.cp() for b in v[1]['b'][0]]
        fnd = True
        while fnd:
            fnd = False
            for ox in range(pg.vcnt):
                o = pg.vs[ox]
                if o is None or ox == v[0]:continue
                ob = o[1]['b'][0]
                adjs = pym.badjbxy(vb,ob,0.1)
                for adj in adjs:
                    vbx,obx = adj
                    vb1,vb2,ob1,ob2 = vb[vbx-1],vb[vbx],ob[obx-1],ob[obx]
                    ips = pym.sintsxyp(vb1,vb2,ob1,ob2,ieb = 0,skew = 0,ie = 0)
                    if ips is None:continue
                    ip1,ip2 = ips
                    if not ip1 in vb or not ip2 in vb:
                        if ip1.onsxy(vb1,vb2,0):vb.insert(vbx,ip1);fnd = True
                        if ip2.onsxy(vb1,vb2,0):vb.insert(vbx,ip2);fnd = True
                if fnd:break
        vbx = 0
        while vbx < len(vb):
            vb1,vb2 = vb[vbx-1],vb[vbx]
            if vb1.d(vb2) < l:vbx += 1
            else:vb.insert(vbx,vb1.mid(vb2))
        return vb

    def vgen_terrain(self,w,v,pg,vhs,l):
        m = dmo.model()
        sgv = w.amodel(None,None,None,m,w.sgraph.root)
        tm = m.agfxmesh()
        vb = self.vstitch(v,pg,l)
        ngvs = m.asurf((tuple(vb),tuple(tuple(h) for h in vhs)),
            tm,fm = 'generic',ref = True,hmin = 100,zfunc = v[1]['tmesh'])
        return m

    def vblgsplotch(self,v,pg,l):
        blgfps = []
        bx,by,sx,sy = -30,30,30,30
        blgfps.append(vec3(bx,by,v[1]['tmesh'](bx,by)).sq(sx,sy))
        bx,by,sx,sy = 30,30,50,30
        blgfps.append(vec3(bx,by,v[1]['tmesh'](bx,by)).sq(sx,sy))
        bx,by,sx,sy = 30,-30,40,60
        blgfps.append(vec3(bx,by,v[1]['tmesh'](bx,by)).sq(sx,sy))
        return blgfps

    def vgen_buildings(self,v,pg,l):
        bfa = blg.blgfactory()
        blgfps = self.vblgsplotch(v,pg,l)
        blgs = []
        for blgfp in blgfps:
            #v[1]['tmesh'].al(blgfp)??
            #blgloop = t.al(blgfp,t.root)
            #seq = bseq.simplebuilding()
            seq = ''
            nblg = bfa.new(None,None,None,footprint = blgfp,sequence = seq)
            blgs.append(nblg)
        return blgs,blgfps



###############################################################################
###############################################################################



###############################################################################
### sequences
###############################################################################

# create a topography for use within a vertex 
def cartographer(g,subseq):
    ss = subseq.split(',')
    tv = g.vs[int(ss[0])]
    fp = [p.cp().ztrn(float(ss[2])) for p in tv[1]['fp']]
    h = float(ss[1])
    tseq = ss[3]
    tm = ter.checkseq(fp,h,tseq,False)
    tv[1]['info']['terrainmesh'] = tm
    return g

# given a "natural" space, 
# create a sequence which generates a "developed" space
# 
# this must include laying of roads through an iterative process
# in the spaces adjacent to the roads, produce hubs of buildings with parking
# if necessary, lay roads which exit the space
#   this means laying a road seam on the boundary
#   this should include the vertex where this road should find
#       another metropolitan area
def metropolis(g,subseq):
    irx,rseq,bseq = subseq.split(',')
    irx = int(irx)
    iv = g.vs[irx]
    print('METROPOLIS!',subseq)

    # the goal is to contrive a sequence representing an urban area
    #   surrounded by rural areas so that it can be connected to other areas
    #   found on a continent including other areas generated by this method
    # this means splitting a region based on geometric information...

    easement = 2

    iv[1]['type'] = ['developed']
    fp = iv[1]['fp']

    rg = pgr.graph()
    rg = rdg.checkseq(rg,fp,rseq,True)

    rgpy = pym.pgtopy(rg,1)

    rgpyi = list(rgpy[0])
    #rgpyi = pym.ebixy(fp,rgpy[0])
    if pym.bnrm(rgpyi).z < 0:
        rgpyi.reverse()
        rgpy = (tuple(rgpyi),rgpy[1])

    #easement_py = pym.ebdxy(fp,rgpyi)
    easement_py = fp

    #rgpyi = pym.ebixy(fp,rgpy[0])

    print('eas',easement_py)

    ax = rg.plot()
    ax = dtl.plot_polygon(easement_py,ax,col = 'b',lw = 2)
    #ax = dtl.plot_polygon(rgpyi,ax,col = 'g',lw = 2)
    #ax = dtl.plot_polygon_full(rgpy,ax,col = 'r',lw = 2)
    plt.show()
    
    new = g.sv(irx,easement_py,rgpyi)
    g.vs[new][1]['type'] = ['corridor']



    '''#
    mseq = subseq+','
    mseq += 'R<0,natural>'
    mseq += 'S<0,0.5,0.5,0,1,root>'
    mseq += 'S<1,0.5,0.5,1,0,root>'
    mseq += 'S<0,0.5,0.5,1,0,root>'
    mseq += 'E<0,1>E<1,2>E<2,3>E<3,0>'
    mseq += 'R<2,developed>'

    ptg.checkseq(g.vs[irx][1]['fp'],20,mseq,True,grammer = g.grammer)

    seq = subseq+','+mseq
    ptg.insert(g,seq)
    '''#






###############################################################################





