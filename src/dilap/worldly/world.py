import dilap.core.base as db
import dilap.core.context as cx
import dilap.modeling.model as dmo
import dilap.modeling.factory as dfa

from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat
import dilap.geometry.polymath as pym

import dilap.topology.planargraph as pgr

import dilap.worldly.terrain as ter
import dilap.worldly.treeskin as ltr
import dilap.worldly.building as blg
import dilap.worldly.blgsequencing as bseq
import dilap.worldly.partitiongraph as ptg
import dilap.worldly.roadgraph as rdg

import dilap.core.plotting as dtl
import matplotlib.pyplot as plt

import math,numpy,random,pdb

###############################################################################
### world factory generates full world contexts from a footprint
###############################################################################

class worldfactory(dfa.factory):

    # process of making a world:
    #   generate an xy boundary polygon for the boundary of the world
    #   partition the interior of this polygon into subsections with type information
    #   generate topology in the form of a deterministic function with feedback
    #       maybe use interpolation over a grid?
    #   place waterways?
    #   create road networks based on region topology
    #       each region has a finite nonzero number of exits
    #           create a network of roads which respects these exits
    #   partition natural and developed regions with 
    #       corridor regions based on road graphs
    #   in developed regions, use the road system to place structures, 
    #       further partitioning the region into developed, buildings, and easements
    #   in natural regions, road graphs are simple and no buildings are present
    #       instead, vegetation and perhaps more extreme topology are added

    def __str__(self):return 'world factory:'
    def __init__(self,*ags,**kws):
        self._def('bclass',cx.context,**kws)
        self._def('boundary',vec3(0,0,0).pring(50,8),**kws)
        dfa.factory.__init__(self,*ags,**kws)
    # generate the models of the world and add them to the context object
    def gen(self,w,pg):
        for vx in range(pg.vcnt):
            v = pg.vs[vx]
            if v is None:continue
            self.vgen(w,v,pg)

    # create a context representing an entire world
    def new(self,*ags,**kws):
        w = self.bclass(*ags,**kws)
        worldh = 50
        #worldseq = 'J<0,2>C<1>'

        toposeq = 'FF'

        contseq =  'R<0,natural>'
        contseq += 'S<0,0.5,0.5,0,1,vertex>'
        contseq += 'S<1,0.5,0.5,1,0,vertex>'
        contseq += 'S<0,0.5,0.5,1,0,vertex>'
        contseq += 'E<0,1>E<1,2>E<2,3>E<3,0>'
        contseq += 'M<0,S<>G<>L<>,>'
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
        vtypes = v[1]['type']
        if 'ocean' in vtypes:
            print('ocean vertex',v[0])
            self.vgen_ocean(w,v,pg)
        elif 'natural' in vtypes:
            print('natural vertex',v[0])
            self.vgen_natural(w,v,pg)
        elif 'corridor' in vtypes:
            print('corridor vertex',v[0])
            self.vgen_corridor(w,v,pg)
        elif 'developed' in vtypes:
            print('developed vertex',v[0])
            self.vgen_developed(w,v,pg)
        return


    def vgen_ocean(self,w,v,pg):
        m = dmo.model()
        sgv = w.amodel(None,None,None,m,w.sgraph.root)
        tm = m.agfxmesh()
        fp = v[1]['fp']
        voids = tuple(tuple(j) for j in v[1]['voids'])
        ngvs = m.asurf((fp,voids),tm,fm = 'concrete1',ref = True,hmin = 100)

    def vgen_natural(self,w,v,pg):
        m = dmo.model()
        sgv = w.amodel(None,None,None,m,w.sgraph.root)
        tm = m.agfxmesh()

        fp = v[1]['fp']
        if 'terrainseam' in v[1]['info']:
            tseam = v[1]['info']['terrainseam']
        else:tseam = [j.cp() for j in fp]

        xpj = vec3(1,0,0).prjps(fp)
        fpl = xpj[1]-xpj[0]
        #ifp = vec3(0,0,0).com(fp).sq(fpl-pg.stitchsize,fpl-pg.stitchsize)
        print('pg.stitchsize',pg.stitchsize)
        ifp = pym.contract(fp,pg.stitchsize)

        if False:ngvs = m.asurf((tseam,(ifp,)),tm,ref = True,hmin = pg.stitchsize)
        else:ngvs = m.asurf((tseam,()),tm,ref = True,hmin = pg.stitchsize)
        #ngvs = m.asurf((tseam,()),tm,ref = True,hmin = pg.stitchsize)

        if 'terrainmesh' in v[1]['info']:
            for ngv in ngvs:
                p = m.pset.ps[tm.verts[ngv][0]]
                p.ztrn(v[1]['info']['terrainmesh'](p.x,p.y))

        #for ngv in ngvs:
        #    p = m.pset.ps[tm.verts[ngv][0]]
        #    p.ztrn(terrain_zfunc(p.x,p.y))

    def vgen_developed(self,w,v,pg):
        #bfa = blg.blgfactory()

        m = dmo.model()
        sgv = w.amodel(None,None,None,m,w.sgraph.root)
        tm = m.agfxmesh()

        fp = v[1]['fp']
        tseam = v[1]['info']['terrainseam']
        voids = tuple(tuple(j) for j in v[1]['voids'])

        xpj = vec3(1,0,0).prjps(fp)
        fpl = xpj[1]-xpj[0]
        #ifp = vec3(0,0,0).com(fp).sq(fpl-pg.stitchsize,fpl-pg.stitchsize)
        print('pg.stitchsize',pg.stitchsize)
        ifp = pym.contract(fp,pg.stitchsize)

        ngvs = m.asurf((tseam,voids),tm,ref = True,hmin = pg.stitchsize)

        if 'terrainmesh' in v[1]['info']:
            for ngv in ngvs:
                p = m.pset.ps[tm.verts[ngv][0]]
                p.ztrn(v[1]['info']['terrainmesh'](p.x,p.y))

        '''#
        m = dmo.model()
        sgv = w.amodel(None,None,None,m,w.sgraph.root)
        tm = m.agfxmesh()
        for x in range(len(ifp)):
            p1,p2 = ifp[x-1],ifp[x]
            p4,p3 = p1.cp().ztrn(10),p2.cp().ztrn(10)
            m.asurf(((p1,p2,p3,p4),()),tm,ref = False)
        '''#
        #m.asurf(([p.cp().ztrn(10) for p in ifp],()),tm,ref = False)

    def vgen_corridor(self,w,v,pg):
        m = dmo.model()
        sgv = w.amodel(None,None,None,m,w.sgraph.root)
        tm = m.agfxmesh()

        fp = v[1]['fp']
        tseam = v[1]['info']['terrainseam']

        xpj = vec3(1,0,0).prjps(fp)
        fpl = xpj[1]-xpj[0]

        print('pg.stitchsize',pg.stitchsize)

        ngvs = m.asurf((tseam,()),tm,fm = 'concrete1',ref = True,hmin = pg.stitchsize)

        if 'terrainmesh' in v[1]['info']:
            for ngv in ngvs:
                p = m.pset.ps[tm.verts[ngv][0]]
                p.ztrn(v[1]['info']['terrainmesh'](p.x,p.y))

        return

        '''#
        m = dmo.model()
        sgv = w.amodel(None,None,None,m,w.sgraph.root)
        tm = m.agfxmesh()
        for x in range(len(ifp)):
            p1,p2 = ifp[x-1],ifp[x]
            p4,p3 = p1.cp().ztrn(10),p2.cp().ztrn(10)
            m.asurf(((p1,p2,p3,p4),()),tm,ref = False)
        m.asurf(([p.cp().ztrn(10) for p in ifp],()),tm,ref = False)
        '''#



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

    rgpy = rg.polygon(1,'ccw')

    rgpyi = list(rgpy[0])
    #rgpyi = pym.ebixy(fp,rgpy[0])
    if pym.bnrm(rgpyi).z < 0:
        rgpyi.reverse()
        rgpy = (tuple(rgpyi),rgpy[1])

    easement_py = pym.ebdxy(fp,rgpyi)
    #easement_py = fp

    #rgpyi = pym.ebixy(fp,rgpy[0])

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






# stitch ensures that the seams between vertices can smoothly meet at their
#   boundary (assuming deterministic z position determination)
def stitch(g,subseq):
    l = float(subseq)
    g.stitchsize = l

    # first create the terrain seam for each vertex starting with the footprint
    #   wherever adjacent vertex boundaries intersect the footprint at points 
    #   which are not currently in the terrain seam, add those points
    for vx in range(g.vcnt):
        v = g.vs[vx]
        if v is None:continue

        # the terrain seam is always the footprint plus additional 
        # points which lie along the edges of the footprint
        vb = [p.cp() for p in v[1]['fp']]
        v[1]['info']['terrainseam'] = vb

        fnd = True
        while fnd:
            fnd = False
            for ox in range(g.vcnt):
                o = g.vs[ox]
                if o is None or ox == vx:continue
                ob = o[1]['fp']
                adjs = pym.badjbxy(vb,ob,0.1)
                acnt = len(adjs)
                if acnt == 0:continue
                else:

                    if acnt > 1:print('multiadjacency!!')

                    for adj in adjs:
                        vbx,obx = adj
                        vb1,vb2 = vb[vbx-1],vb[vbx]
                        ob1,ob2 = ob[obx-1],ob[obx]
                        ips = pym.sintsxyp(vb1,vb2,ob1,ob2,ieb = 0,skew = 0,ie = 0)
                        if ips is None:continue
                        ip1,ip2 = ips
                        if not ip1 in vb or not ip2 in vb:
                            if ip1.onsxy(vb1,vb2,0):
                                vb.insert(vbx,ip1);fnd = True
                            if ip2.onsxy(vb1,vb2,0):
                                vb.insert(vbx,ip2);fnd = True
                    if fnd:break

    for vx in range(g.vcnt):
        v = g.vs[vx]
        if v is None:continue
        vb = v[1]['info']['terrainseam']
        vbx = 0
        while vbx < len(vb):
            vb1,vb2 = vb[vbx-1],vb[vbx]
            if vb1.d(vb2) < l:vbx += 1
            else:vb.insert(vbx,vb1.mid(vb2))

###############################################################################





