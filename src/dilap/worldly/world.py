import dilap.core.base as db
import dilap.core.context as cx
import dilap.modeling.model as dmo
import dilap.modeling.factory as dfa

from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat
import dilap.geometry.polymath as pym

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

    def __str__(self):return 'world factory:'
    def __init__(self,*ags,**kws):
        #self._def('bclass',world,**kws)
        self._def('bclass',cx.context,**kws)
        self._def('boundary',vec3(0,0,0).pring(50,8),**kws)
        #self._def('boundary',vec3(0,0,0).sq(100,100),**kws)
        dfa.factory.__init__(self,*ags,**kws)
    def new(self,*ags,**kws):
        w = self.bclass(*ags,**kws)
        worldh = 50
        worldseq = 'C<0>'
        eg = {
          'C':continent,
          'M':metropolis,
          'T':stitch,
              }
        pg = ptg.checkseq(self.boundary,worldh,worldseq,show = True,grammer = eg)
        self.gen(w,pg)
        return w

    # generate the models of the world and add them to the context object
    def gen(self,w,pg):
        for vx in range(pg.vcnt):
            v = pg.vs[vx]
            if v is None:continue
            self.vgen(w,v,pg)




    # vgen takes a context, vertex, and partitiongraph and from that
    #   it generates all models within the space of the vertex
    # this consists of two fundamental steps:
    #   based on the nature of this and adjacent vertices, 
    #       the boundary of this vertex is constrained, 
    #       and meshes which cross it must be forced to appear continuous
    #   based on the nature of this vertex, the interior of the space 
    #       associated with the vertex must be populated with models
    #        - these models must geometrically cover the entire ground surface
    def vgen(self,w,v,pg):

        # types:
        #   "natural" - cover with terrain vegetation only
        #       edges should only require a terrain seam
        #   "developed" - cover with structures, filling gaps with nature models
        #       edges should only require terrain seams and road seams
        #       where roads may intersect non-colinearly with the boundary only
        #   "corridor" - cover with infrastructure, filling gaps with concrete?
        #       edges should only require terrain seams and road seams
        #       road seems may include colinear intersections with the boundary
        #
        # this supports the idea that all space is 
        #   either infrastructure, structure, or nature

        vtypes = v[1]['type']
        if 'natural' in vtypes:
            print('natural vertex',v[0])
            self.vgen_natural(w,v,pg)
        elif 'corridor' in vtypes:
            print('corridor vertex',v[0])
            self.vgen_corridor(w,v,pg)
        elif 'developed' in vtypes:
            print('developed vertex',v[0])
            self.vgen_developed(w,v,pg)

        return


    def vgen_natural(self,w,v,pg):
        m = dmo.model()
        sgv = w.amodel(None,None,None,m,w.sgraph.root)
        tm = m.agfxmesh()

        fp = v[1]['fp']
        tseam = v[1]['info']['terrainseam']

        xpj = vec3(1,0,0).prjps(fp)
        fpl = xpj[1]-xpj[0]
        #ifp = vec3(0,0,0).com(fp).sq(fpl-pg.stitchsize,fpl-pg.stitchsize)
        print('pg.stitchsize',pg.stitchsize)
        ifp = pym.contract(fp,pg.stitchsize)

        if False:ngvs = m.asurf((tseam,(ifp,)),tm,ref = True,hmin = pg.stitchsize)
        else:ngvs = m.asurf((tseam,()),tm,ref = True,hmin = pg.stitchsize)
        #ngvs = m.asurf((tseam,()),tm,ref = True,hmin = pg.stitchsize)

        for ngv in ngvs:
            p = m.pset.ps[tm.verts[ngv][0]]
            p.ztrn(terrain_zfunc(p.x,p.y))

    def vgen_developed(self,w,v,pg):
        #bfa = blg.blgfactory()

        m = dmo.model()
        sgv = w.amodel(None,None,None,m,w.sgraph.root)
        tm = m.agfxmesh()

        fp = v[1]['fp']
        tseam = v[1]['info']['terrainseam']

        xpj = vec3(1,0,0).prjps(fp)
        fpl = xpj[1]-xpj[0]
        #ifp = vec3(0,0,0).com(fp).sq(fpl-pg.stitchsize,fpl-pg.stitchsize)
        print('pg.stitchsize',pg.stitchsize)
        ifp = pym.contract(fp,pg.stitchsize)

        if False:ngvs = m.asurf((tseam,(ifp,)),tm,ref = True,hmin = pg.stitchsize)
        else:ngvs = m.asurf((tseam,()),tm,ref = True,hmin = pg.stitchsize)
        #ngvs = m.asurf((tseam,()),tm,ref = True,hmin = pg.stitchsize)

        for ngv in ngvs:
            p = m.pset.ps[tm.verts[ngv][0]]
            p.ztrn(terrain_zfunc(p.x,p.y))

        m = dmo.model()
        sgv = w.amodel(None,None,None,m,w.sgraph.root)
        tm = m.agfxmesh()
        for x in range(len(ifp)):
            p1,p2 = ifp[x-1],ifp[x]
            p4,p3 = p1.cp().ztrn(10),p2.cp().ztrn(10)
            m.asurf(((p1,p2,p3,p4),()),tm,ref = False)
        m.asurf(([p.cp().ztrn(10) for p in ifp],()),tm,ref = False)

    def vgen_corridor(self,w,v,pg):
        m = dmo.model()
        sgv = w.amodel(None,None,None,m,w.sgraph.root)
        tm = m.agfxmesh()

        fp = v[1]['fp']
        tseam = v[1]['info']['terrainseam']

        xpj = vec3(1,0,0).prjps(fp)
        fpl = xpj[1]-xpj[0]
        #ifp = vec3(0,0,0).com(fp).sq(fpl-pg.stitchsize,fpl-pg.stitchsize)
        print('pg.stitchsize',pg.stitchsize)
        ifp = pym.contract(fp,pg.stitchsize)

        #if False:ngvs = m.asurf((tseam,(ifp,)),tm,fm = 'concrete1',ref = True,hmin = pg.stitchsize)
        #else:ngvs = m.asurf((tseam,()),tm,ref = True,hmin = pg.stitchsize)
        ngvs = m.asurf((tseam,()),tm,fm = 'concrete1',ref = True,hmin = pg.stitchsize)

        for ngv in ngvs:
            p = m.pset.ps[tm.verts[ngv][0]]
            p.ztrn(terrain_zfunc(p.x,p.y))

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

# begin with self.boundary - the shoreline of a continent
# establish a port region - all roads can trace back to this region
# each segment of road specified potentially splits the boundary
#   into additional subregions with distinct characteristics
# regions grow iteratively
#   region growth means additional roads/buildings/decorative meshes
# when roads/buildings/etc are generated, 
#   cover remaining void with terrain meshes


# continent must partition an entire continent of "natural" space
#   into "natural", "developed", and "corridor" regions
def continent(g,subseq):
    irx = int(subseq)
    print('CONTINENT!',subseq)
    seq = subseq+','
    seq += 'R<0,natural>'
    seq += 'S<0,0.5,0.5,0,1,vertex>'
    seq += 'S<1,0.5,0.5,1,0,vertex>'
    seq += 'S<0,0.5,0.5,1,0,vertex>'
    seq += 'E<0,1>E<1,2>E<2,3>E<3,0>'
    seq += 'M<0>T<2>'

    #ptg.checkseq(g.vs[irx][1]['fp'],20,seq,True,grammer = g.grammer)

    ptg.insert(g,seq)

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
    irx = int(subseq)
    iv = g.vs[irx]
    print('METROPOLIS!',subseq)

    # the goal is to contrive a sequence representing an urban area
    #   surrounded by rural areas so that it can be connected to other areas
    #   found on a continent including other areas generated by this method
    # this means splitting a region based on geometric information...

    iv[1]['type'] = ['developed']
    fp = iv[1]['fp']
    mc = vec3(0,0,0).com(fp)

    exits = [fp[-1].lerp(fp[0],0.5)]


    rg = rdg.wgraph()
    i1 = rg.av(p = exits[0],l = 0)
    i2,r1 = rg.mev(i1,{'p':mc,'l':0},{})

    ax = rg.plot()
    ax = dtl.plot_polygon(fp,ax,col = 'r',lw = 2)
    plt.show()

    


    raise NotImplemented


    ixs,rds = rdg.road_graph(fp,exits)
    print('made roads...',ixs,rds)

    ibbs,rbbs = rdg.bboxes(fp,ixs,rds)

    ebd = fp
    rks = list(rbbs.keys())

    stack = [ixs[0]]
    tip = 0
    while True:
        if stack:i = stack.pop(tip)
        else:break

        ebd = pym.ebdxy(ebd,ibbs[i[0]])
        ebi = pym.ebixy(fp,ibbs[i[0]])

        ax = dtl.plot_axes_xy(100)
        ax = dtl.plot_polygon_xy(pym.contract(ebd,1),ax,col = 'b',lw = 3.0)
        ax = dtl.plot_polygon_xy(pym.contract(ebi,1),ax,col = 'g',lw = 3.0)
        plt.show()

        new = g.sv(irx,ebd,ebi)
        g.vs[new][1]['type'] = ['corridor']

        #ring = [ixs[j] for j in i[2]]
        for r in i[2]:
            k = (i[0],r)
            if k in rks:
                rks.remove(k)

                ebd = pym.ebdxy(ebd,rbbs[(i[0],r)])
                ebi = pym.ebixy(fp,rbbs[(i[0],r)])
                #rdb = pym.ebdxy(ebi,ibbs[i[0]])
                #rdb = pym.ebdxy(rdb,ibbs[r])

                ax = dtl.plot_axes_xy(100)
                ax = dtl.plot_polygon_xy(pym.contract(ebd,1),ax,col = 'b',lw = 3.0)
                #ax = dtl.plot_polygon_xy(pym.contract(rdb,1),ax,col = 'g',lw = 3.0)
                ax = dtl.plot_polygon_xy(pym.contract(ebi,1),ax,col = 'g',lw = 3.0)
                plt.show()

                new = g.sv(irx,ebd,ebi)
                g.vs[new][1]['type'] = ['corridor']

                stack.append(ixs[r])

    return

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

                # PROBABLY SHOULD CALL BADJBXY IN WHILE LOOP UNTIL IT FINDS NOTHING
                # SINCE IT CHANGES THE LOOP TO FIX AN ADJACENCY, IT HAS TO BE 
                # RECALLED FOR EACH POTENTIAL ADJACENCY ANYWAY...

                adjs = pym.badjbxy(vb,ob,0.1)
                acnt = len(adjs)
                if acnt == 0:continue
                else:
                    if acnt > 1:print('multiadjacency!!')
                    for adj in adjs:
                        vbx,obx = adj
                        vb1,vb2 = vb[vbx-1],vb[vbx]
                        ob1,ob2 = ob[obx-1],ob[obx]

                        #print('urrgh',vbx,obx,vb1,vb2,ob1,ob2)

                        ips = pym.sintsxyp(vb1,vb2,ob1,ob2,ieb = 0,skew = 0,ie = 0)
                        #if ips is None:
                        if False:
                            ax = dtl.plot_axes(100)
                            ax = dtl.plot_polygon(vb,ax,lw = 2,col = 'b')
                            ax = dtl.plot_polygon(ob,ax,lw = 2,col = 'g')
                            ax = dtl.plot_edges((vb1,vb2),ax,lw = 3,col = 'b')
                            ax = dtl.plot_edges((ob1,ob2),ax,lw = 3,col = 'g')
                            plt.show()
                        if ips is None:
                            continue

                        ip1,ip2 = ips
                        #ip1,ip2 = pym.sintsxyp(vb1,vb2,ob1,ob2)
                        if not ip1 in vb:
                            if ip1.onsxy(vb1,vb2,0):
                                vb.insert(vbx,ip1)
                                fnd = True
                        if not ip2 in vb:
                            if ip2.onsxy(vb1,vb2,0):
                                vb.insert(vbx,ip2)
                                fnd = True
                        #if ip1.onsxy(vb1,vb2,0):vb.insert(vbx,ip1)
                        #if ip2.onsxy(vb1,vb2,0):vb.insert(vbx,ip2)
                        #if ip1.onsxy(ob1,ob2,0):ob.insert(obx,ip1)
                        #if ip2.onsxy(ob1,ob2,0):ob.insert(obx,ip2)

                    if fnd:break

    for vx in range(g.vcnt):
        v = g.vs[vx]
        if v is None:continue
        #vb = v[1]['fp']
        vb = v[1]['info']['terrainseam']

        vbx = 0
        while vbx < len(vb):
            vb1,vb2 = vb[vbx-1],vb[vbx]
            if vb1.d(vb2) < l:vbx += 1
            else:vb.insert(vbx,vb1.mid(vb2))

###############################################################################
###############################################################################

# use this function to provide deterministic z location of terrain point 
# as a function of x and y
#   each partition graph vertex will use this function 
#   to locate its bounding polygons
def terrain_zfunc(x,y):
    if x == 0 and y == 0:t = 0
    elif x == 0:t = numpy.pi/2.0 if y > 0 else -numpy.pi/2.0
    else:
        t = math.atan(y/x)
    r = x**2 + y**2
    z = 10*r**(1.0/8.0)*math.sin(t)

    #z = (10+r)**0.125 + math.cos(x/10.0) + math.sin(y/10.0)
    z = (10+r)**0.25 + math.cos(x/10.0) + math.sin(y/10.0)
    return z

###############################################################################





