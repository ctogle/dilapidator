import dilap.core.base as db
import dilap.core.context as cx
import dilap.geometry.tools as gtl
from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat
from dilap.geometry.pointset import pointset
import dilap.geometry.polymath as pym
import dilap.modeling.model as dmo
import dilap.topology.trimesh as dtm
import dilap.topology.tree as dtr
import dilap.topology.vert as dvt
import dilap.topology.edge as deg
import dilap.topology.loop as dlp
import dilap.topology.face as dfc

import dilap.modeling.factory as dfa
import dilap.geometry.triangulate as dtg

import dilap.core.plotting as dtl
import matplotlib.pyplot as plt

import numpy,math,random
import pdb



###############################################################################
### functions which create geometry in the form of polygons
###############################################################################

# return geometry describing a wall from p1 to p2 of height wh with holes hs
#   return polygons constituting a trimesh for the wall
#   return boundary polygons where portals will connect to the wall
def awall(p1,p2,wh,hs,hp1 = None,hp2 = None):
    if hp1 is None:hp1 = p1
    if hp2 is None:hp2 = p2
    wn = p1.tov(p2).nrm().crs(vec3(0,0,1))
    hp1.prj(p1,wn);hp2.prj(p2,wn)
    polys,portals = [],[]

    eb = [p1.cp(),p2.cp(),p2.cp().ztrn(wh),p1.cp().ztrn(wh)]
    ibs = []

    for h in hs:
        # which side of eb,where along it, how wide, how tall, dist from side
        dx,dp,dw,dh,dz = h
        ddp = dw/(2.0*hp1.d(hp2))
        d1,d2 = hp1.lerp(hp2,dp-ddp),hp1.lerp(hp2,dp+ddp)
        h = [d2.cp().ztrn(dz),d2.cp().ztrn(dz+dh),
            d1.cp().ztrn(dz+dh),d1.cp().ztrn(dz)]
        portals.append(h)
        if dz == 0:
            for hp in h:eb.insert(1,hp)
        else:ibs.append(h)

    wpy = (tuple(eb),tuple(tuple(x) for x in ibs))
    polys.append(wpy)

    return polys,portals

# return polygons constituting a trimesh of a portal 
def aportal(oloop,p1,p2,ww):
    polys = []
    pn = vec3(0,0,1).crs(p1.tov(p2).nrm()).uscl(ww)
    iloop = [p.cp().trn(pn) for p in oloop]
    for lx in range(len(oloop)):
        polys.append(((iloop[lx-1],oloop[lx-1],oloop[lx],iloop[lx]),()))
    return polys

###############################################################################
### functions which modify room instructions operationally 
###############################################################################

# given two rooms, determine which walls if any are shared by both rooms
def adjacentwalls(r1,r2,minovlp = 1):
    adjs = []
    for r1x in range(len(r1)):
        r1p1,r1p2 = r1[r1x-1],r1[r1x]
        for r2x in range(len(r2)):
            r2p1,r2p2 = r2[r2x-1],r2[r2x]
            r1p1nr = r1p1.isnear(r2p1) or r1p1.isnear(r2p2)
            r1p2nr = r1p2.isnear(r2p1) or r1p2.isnear(r2p2)
            ips = pym.sintsxyp(r1p1,r1p2,r2p1,r2p2,ie = False,skew = False)
            if type(ips) == type(()):
                if ips[0].d(ips[1]) > minovlp:
                    adjs.append((r1x,r2x))
    return adjs

###############################################################################
### functions which add trimesh geometry to models
###############################################################################

def genroof(m,rmkws):
    tm = m.agfxmesh()

    return

    for lx in range(self.floors):
        lvl = self.rims[lx]
        if lx == self.floors-1:
            tlvl = lvl
            print('make doors if the level below is not fully eclipsed')
        else:tlvl = pym.ebdxy(lvl,self.rims[lx+1])

        #ax = dtl.plot_axes_xy(10)
        #ax = dtl.plot_polygon_xy(tlvl,ax,lw = 2)
        #plt.show()

        rbnd = (tuple(tlvl),())
        m.asurf(rbnd,tm)

        for x in range(len(tlvl)):
            p1,p2 = tlvl[x-1],tlvl[x]

            wpys,portals = awall(p1,p2,2,[])
            for wpy in wpys:m.asurf(wpy,tm)

        itlvl = pym.contract(tlvl,0.125)
        for x in range(len(itlvl)):
            p1,p2 = itlvl[x-1],itlvl[x]

            wpys,portals = awall(p2,p1,2,[])
            for wpy in wpys:m.asurf(wpy,tm)


        # must zoffset tlvl and itlvl?
        #m = self.asurf((tlvl,(itlvl,)),m = m,gm = gm)

###############################################################################
###############################################################################
###############################################################################





###############################################################################
### context representing any building
###############################################################################

class building(cx.context):

    def __init__(self,p,q,s,*ags,**kws):
        self._def('name','buildingcontext',**kws)
        self._def('bgraph',None,**kws)
        self.p = p
        self.q = q
        self.s = s
        cx.context.__init__(self,*ags,**kws)

    # insert a hole specification for a particular wall
    def wallhole(self,r,w,h):
        hp = h[1]
        whs = self.bgraph.vs[r][2]['wholes'][w]
        for wx in range(len(whs)):
            if hp > whs[wx][1]:
                whs.insert(wx,h)
                return
            elif hp == whs[wx][1]:
                print('cannot add another hole here...',hp)
                return
        whs.append(h)

    # add information to make doors for this room based on complete topology
    def holes(self,rx):
        sl,sw = 8,8
        dw,dh,dz = 1.5,3.0,0.0
        rv = self.bgraph.vs[rx]
        rb = rv[2]['bound']
        rexits = rv[2]['exits']
        for adj in self.bgraph.vs[rx][1]:
            ov = self.bgraph.vs[adj]
            if ov is None:continue
            elif ov[2]['level'] > rv[2]['level']:rv[2]['shaft'] = True
            elif ov[2]['level'] < rv[2]['level']:rv[2]['shaft'] = True
            else:
                ob = ov[2]['bound']
                aws = adjacentwalls(rb,ob,minovlp = dw+2)
                for aw in aws:
                    rwx,owx = aw
                    ips = pym.sintsxyp(rb[rwx-1],rb[rwx],ob[owx-1],ob[owx])
                    imp = ips[0].mid(ips[1])
                    ridp = 1-rb[rwx-1].d(imp)/rb[rwx-1].d(rb[rwx])
                    self.wallhole(rx,rwx-1,(adj,ridp,dw,dh,dz))
                    rv[2]['wtypes'][rwx-1] = 'i'
        gadjs = self.bgraph.geoadj(rx)
        for adj in gadjs:
            ov = self.bgraph.vs[adj[0]]
            if ov is None:continue
            ob = ov[2]['bound']
            aws = adjacentwalls(rb,ob,minovlp = 0)
            for aw in aws:
                rwx,owx = aw
                if rv[2]['wtypes'][rwx-1] == 'e':
                    rv[2]['wtypes'][rwx-1] = 'i'
        for wx in range(len(rb)):
            wp1,wp2 = rb[wx-1],rb[wx]
            wt = rv[2]['wtypes'][wx-1]
            if wt == 'e':
                isexit = False
                if rv[2]['level'] == 0:
                    for rex in rexits:
                        if gtl.onseg_xy(rex,wp1,wp2) and not wp2.isnear(rex):
                            print('external exit!',rex)
                            self.wallhole(rx,wx-1,('exit',0.5,dw,dh,dz))
                            isexit = True
                            break
                if not isexit:
                    self.wallhole(rx,wx-1,('window',0.5,dw,dh-1,dz+1))

    def genshafts(self):

        def genwrap(shb,shx,skirt = True,crown = True):
            vkw = self.bgraph.vs[shvx][2]
            for x in range(len(shb)):
                b2,b1,w2,w1 = shb[x-1],shb[x],fbnd[x-1],fbnd[x]
                hs,wt = vkw['wholes'][x-1],vkw['wtypes'][x-1]
                wh,ww = vkw['wheights'][x-1],vkw['wwidths'][x-1]
                sh,ch = vkw['skirt'],vkw['crown']
                wpys,portals = awall(w1,w2,wh,hs,b1.cp(),b2.cp())
                for wpy in wpys:m.asurf(wpy,tm)
                for hps in portals:
                    wpys = aportal(hps,w1,w2,ww)
                    for wpy in wpys:m.asurf(wpy,tm)
                if skirt:
                    sp1,sp2 = w1.cp().ztrn(-sh),w2.cp().ztrn(-sh)
                    wpys,portals = awall(sp1,sp2,sh,[])
                    for wpy in wpys:m.asurf(wpy,tm)
                if crown:
                    cp1,cp2 = w1.cp().ztrn(wh),w2.cp().ztrn(wh)
                    wpys,portals = awall(cp1,cp2,ch,[])
                    for wpy in wpys:m.asurf(wpy,tm)

        def genplatform(shb,shx):
            platform = [p.cp() for p in shb]

            buff,rw = 4,4

            r1,r2 = platform[2].cp(),platform[3].cp()

            rtn = r1.tov(r2).nrm().uscl(buff)
            rnm = vec3(0,0,1).crs(rtn).nrm().uscl(rw)
            r2,r1 = r2.trn(rtn.flp()),r1.trn(rtn.flp())
            r3,r4 = r2.cp().trn(rnm),r1.cp().trn(rnm)

            belv = self.bgraph.vs[sh[shx-1]][2]
            rh = max(belv['wheights'])+belv['skirt']+belv['crown']
            ramp = [r4,r1,r2,r3]
            platform = pym.ebdxy(platform,ramp)
            ramp[0].ztrn(-rh)
            ramp[1].ztrn(-rh)

            m.asurf((ramp,()),tm)
            m.asurf((platform,()),tm)

        def genbottom(shb,abvx):
            shaftfloor = shb
            m.asurf((shaftfloor,()),tm)
            genwrap(shb,shx-1,skirt = False,crown = shcnt > 1)

        def gentop(shb,belx):
            shaftceiling = [p.cp().ztrn(wht) for p in shb]
            m.asurf((shaftceiling,()),tm,rv = True)
            genplatform(shb,shx)
            genwrap(shb,shx+1,crown = False)

        def genlink(shb,shx):
            genplatform(shb,shx)
            genwrap(shb,shx)

        for sh in self.bgraph.shafts:
            m = dmo.model()
            sgv = self.amodel(self.p,self.q,self.s,m,self.sgraph.root)
            tm = m.agfxmesh()

            shcnt = len(sh)
            for shx in range(shcnt):
                shvx = sh[shx]
                shv = self.bgraph.vs[shvx]
                shb = shv[2]['bound']
                wht = max(shv[2]['wheights'])
                wwi = min(shv[2]['wwidths'])
                fbnd = pym.contract(shb,wwi)

                if shx == 0:genbottom(fbnd,shx+1)
                elif shx == shcnt-1:gentop(fbnd,shx-1)
                else:genlink(fbnd,shx)

    # add the trimeshes of a room interior to a model
    def genrooms(self):
        for v in self.bgraph.vs:
            if v is None:continue
            vx,ves,vkw = v

            # insert doors based on topology
            self.holes(vx)

            if vkw['shaft']:continue

            # add a new trimesh to the model
            m = dmo.model()
            sgv = self.amodel(self.p,self.q,self.s,m,self.sgraph.root)
            tm = m.agfxmesh()

            # create the floor and ceiling surfaces
            fbnd = pym.contract(vkw['bound'],min(vkw['wwidths']))

            floor,ceiling = vkw['floor'],vkw['ceiling']
            if floor is None:
                floor = (tuple([p.cp() for p in fbnd]),
                    tuple(tuple(fh) for fh in vkw['fholes']))
            if floor:m.asurf(floor,tm)
            if ceiling is None:
                wh = min(vkw['wheights'])
                ceiling = (tuple([p.cp().ztrn(wh) for p in fbnd]),
                    tuple(tuple(ch) for ch in vkw['choles']))
            if ceiling:m.asurf(ceiling,tm,rv = True)

            # create the wall surfaces and portal surfaces
            for x in range(len(vkw['bound'])):
                b2,b1 = vkw['bound'][x-1],vkw['bound'][x]
                w2,w1 = fbnd[x-1],fbnd[x]
                hs,wt = vkw['wholes'][x-1],vkw['wtypes'][x-1]
                wh,ww = vkw['wheights'][x-1],vkw['wwidths'][x-1]

                # create wall surfaces
                wpys,portals = awall(w1,w2,wh,hs,b1.cp(),b2.cp())
                for wpy in wpys:m.asurf(wpy,tm)

                # create portal surfaces
                for hps in portals:
                    wpys = aportal(hps,w1,w2,ww)
                    for wpy in wpys:m.asurf(wpy,tm)

    def genshell(self):
        m = dmo.model()
        sgv = self.amodel(self.p,self.q,self.s,m,self.sgraph.root)
        tm = m.agfxmesh()
        z = 0
        for lvx in range(self.bgraph.levels):
            fp = [p.cp().ztrn(z) for p in self.bgraph.footprint]
            z += 7
            rim = []
            for fpx in range(len(fp)):
                fp1,fp2 = fp[fpx-1],fp[fpx]
                for v in self.bgraph.vs:
                    if v is None:continue
                    vx,ves,vkw = v
                    vb = vkw['bound']
                    for vbx in range(len(vb)):
                        vb1,vb2 = vb[vbx-1],vb[vbx]
                        ips = pym.sintsxyp(fp1,fp2,vb1,vb2,ie = False,skew = False)
                        if not type(ips) is type(()):continue
                        hs = vkw['wholes'][vbx-1]
                        if not hs:hs.append((0,0.5,1.5,2,1))
                        ww,wh = vkw['wwidths'][vbx-1],vkw['wheights'][vbx-1]
                        wtn = vb1.tov(vb2).crs(vec3(0,0,1)).nrm().uscl(ww)
                        wp1,wp2 = vb1.cp().trn(wtn),vb2.cp().trn(wtn)
                        wpys,portals = awall(wp1,wp2,wh,hs)
                        for wpy in wpys:m.asurf(wpy,tm)
                        for hps in portals:
                            wpys = aportal(hps,wp1,wp2,ww)
                            for wpy in wpys:m.asurf(wpy,tm)
                        sh,ch = vkw['skirt'],vkw['crown']
                        sp1,sp2 = wp1.cp().ztrn(-sh),wp2.cp().ztrn(-sh)
                        cp1,cp2 = wp1.cp().ztrn(wh),wp2.cp().ztrn(wh)
                        wpys,portals = awall(sp1,sp2,sh,[])
                        for wpy in wpys:m.asurf(wpy,tm)
                        wpys,portals = awall(cp1,cp2,ch,[])
                        for wpy in wpys:m.asurf(wpy,tm)

                        # need last w2 and current w1
                        if not rim:rim.append(wp1)
                        elif not rim[-1].isnear(wp1):rim.append(wp1)
                        if not rim:rim.append(wp2)
                        elif not rim[-1].isnear(wp2):rim.append(wp2)

            #if lvx == 0:
            #    ax = dtl.plot_axes_xy(50)
            #    ax = dtl.plot_polygon_xy(rim,ax)
            #    plt.show()

            #wpys,portals = awall(cp1,cp2,sh+wh+ch,[])
            #for wpy in wpys:m.asurf(wpy,tm)

    # do something which fills the scenegraph
    def generate(self,worn = 0):
        self.floors = len(list(self.bgraph.fl_look.keys()))

        self.genrooms()
        self.genshafts()
        self.genshell()

        #self.genroof()

        # create the models which cover the roof
        #m = dmo.model()
        #sgv = self.amodel(p,q,s,m,None)
        #genroof(m,self.rmkws)

        return self

###############################################################################
### graph structure to create building contexts from
###############################################################################

def splitb(b,pline = None,mindx = 5,mindy = 5,pos = None,ori = None):
    if pline is None:

        x,y = vec3(1,0,0),vec3(0,1,0)
        bx,by = x.prjps(b),y.prjps(b)
        bxd = bx[1]-bx[0]
        byd = by[1]-by[0]

        if pos is None:
            dx = random.uniform(0,bxd/2.0-mindx)
            dy = random.uniform(0,byd/2.0-mindy)
            spos = vec3(dx,dy,0).com(b)
        elif type(pos) is type(0.0):
            #minbx,maxbx = bx[0]+min(mindx,bxy/2.0),bx[1]-
            minbx,maxbx = bx[0]+mindx,bx[1]-mindx
            minby,maxby = by[0]+mindy,by[1]-mindy
            dx = minbx + (maxbx - minbx)*pos
            dy = minby + (maxby - minby)*pos
            spos = vec3(dx,dy,0)

        if ori is None:bdir = y.cp() if bxd > byd else x.cp()
        elif ori == 'v':bdir = y.cp()
        elif ori == 'h':bdir = x.cp()
        else:raise ValueError

        bdir.uscl(100)

        s1,s2 = spos.cp().trn(bdir.flp()),spos.cp().trn(bdir.flp())

        '''#
        ax = dtl.plot_axes_xy(50)
        ax = dtl.plot_polygon_xy(b,ax)
        ax = dtl.plot_edges_xy((s1,s2),ax,col = 'r')
        plt.show()
        '''#

        l,r = pym.bsegsxy(b,s1,s2)
    else:
        s1,s2 = pline
        l,r = pym.bsegsxy(b,s1,s2)
    return l,r

class blggraph(db.base):

    def av(self,os,kws):
        vx = self.vcnt
        self.vcnt += 1
        self.vs.append([vx,os,kws])
        lvl = kws['level']
        if lvl in self.fl_look:self.fl_look[lvl].append(vx)
        else:self.fl_look[lvl] = [vx]
        return vx

    # connect can be True, False, or a list of points supposedly along the footprint
    def rv(self,rx,nvs = [],connect = False):
        if self.vs[rx] is None:
            print('replacing none vertex!',self.vs[rx])
            return
        ev = self.vs[rx]
        for nvx in nvs:
            nv = self.vs[nvx]
            nv[2]['level'] = ev[2]['level']
        for res in ev[1]:
            if self.vs[res] is None:continue
            for nvx in nvs:
                nv = self.vs[nvx]
                if connect is True:
                    if not res in nv[1]:nv[1].append(res)
                    if not nv[0] in self.vs[res][1]:
                        self.vs[res][1].append(nv[0])
                elif type(connect) is type(()):
                    offset = min(nvs)
                    if res == connect[0] and nvx-offset == connect[1]:
                        if not res in nv[1]:nv[1].append(res)
                        if not nv[0] in self.vs[res][1]:
                            self.vs[res][1].append(nv[0])
                    print('insert connection!',connect)
        self.vs[rx] = None
        
    # merge two vertices into a single vertex of their combined volume
    def mv(self,r1x,r2x):
        r1,r2 = self.vs[r1x],self.vs[r2x]
        if r1[2]['level'] == r2[2]['level']:
            print('merge vertices on the same floor',r1x,r2x)
        else:
            print('merge vertices on different floors',r1x,r2x)

        vx = self.vcnt
        self.vcnt += 1

        pdb.set_trace()

        self.vs.append([vx,os,kws])

        lvl = kws['level']
        if lvl in self.fl_look:self.fl_look[lvl].append(vx)
        else:self.fl_look[lvl] = [vx]

        return vx

    # add an edge between two vertices if it does not exists
    def ae(self,v1,v2):
        vv1,vv2 = self.vs[v1][1],self.vs[v2][1]
        if not v2 in vv1:vv1.append(v2)
        if not v1 in vv2:vv2.append(v1)

    # remove the edge between two vertices if it exists
    def re(self,v1,v2):
        self.vs[v1][1].remove(v2)
        self.vs[v2][1].remove(v1)

    # given the index of a vertex, return the indices of vertices 
    #   who are geometrically adjacent within the building
    def geoadj(self,rx,minovlp = 0):
        rv = self.vs[rx]
        bd = rv[2]['bound']
        fnd = []
        for ra in range(self.vcnt):
            av = self.vs[ra]
            if av is None or av is rv:continue
            ab = av[2]['bound']
            adjwalls = adjacentwalls(bd,ab)
            for adjw in adjwalls:fnd.append((ra,adjw))
        return fnd

    # given the index of a vertex, return the indices of vertices 
    #   who are topologically and geometrically adjacent within the building
    def anyadj(self,rx,minovlp = 0):
        rv = self.vs[rx]
        bd = rv[2]['bound']
        es = rv[1]
        fnd = []
        if not es:return fnd
        for ra in rv[1]:
            av = self.vs[ra]
            adjwalls = adjacentwalls(rv[2]['bound'],av[2]['bound'])
            for adjw in adjwalls:
                fnd.append((ra,adjw))
        return fnd

    # consider all rooms that room rx should be connected to based on edges
    def verifyedges(self,rx):
        dw,dh,dz = 1.5,3.0,0.0
        rv = self.vs[rx]
        if rv is None:return
        rb = rv[2]['bound']
        for adj in self.vs[rx][1]:
            ov = self.vs[adj]
            if ov is None:continue
            ob = ov[2]['bound']
            aws = adjacentwalls(rb,ob)
            if not aws:
                if adj in rv[1]:rv[1].remove(adj)
                if rx in ov[1]:ov[1].remove(rx)
        # remove any exits not along this wall
    
    # attach another graph to this one based on geometry
    # srx/orx are the indices of the adjacent rooms in each graph
    def attach(self,srx,orx,o):
        vc = self.vcnt
        nvs = []
        for ox in range(o.vcnt):
            if ox == orx:
                print('bridge the graphs')
                pdb.set_trace()
            else:
                ov = o.vs[ox]
                okws = ov[2]
                new = self.av([ovx+vc for ovx in ov[1]],okws)
                nvs.append(new)
        return nvs

    ###########################################################################

    def plotxy(self,fp = None,ax = None):
        if ax is None:ax = dtl.plot_axes_xy(80)
        for rmv in self.vs:
            if rmv is None:continue
            ax = dtl.plot_polygon_xy(rmv[2]['bound'],ax,lw = 4)
            rc = vec3(0,0,0).com(rmv[2]['bound'])
            rs = str(rmv[0])+','+str(rmv[1])
            ax = dtl.plot_point_xy(rc,dtl.plot_point_xy_annotate(rc,ax,rs))
            if rmv[1]:
                for re in rmv[1]:
                    if self.vs[re] is None:continue
                    rt = (rc,vec3(0,0,0).com(self.vs[re][2]['bound']))
                    ax = dtl.plot_edges_xy(rt,ax,col = 'g')
            for exit in rmv[2]['exits']:
                ax = dtl.plot_point_xy(exit,dtl.plot_point_xy_annotate(exit,ax,'exit'))
        if not fp is None:
            ax = dtl.plot_polygon_xy(fp,ax,lw = 2,col = 'r')
        return ax

    def plot(self,fp = None,ax = None):
        if ax is None:ax = dtl.plot_axes(80)
        for rmv in self.vs:
            if rmv is None:continue
            ax = dtl.plot_polygon(rmv[2]['bound'],ax,lw = 4)
            rc = vec3(0,0,0).com(rmv[2]['bound'])
            rs = str(rmv[0])+','+str(rmv[1])
            ax = dtl.plot_point(rc,dtl.plot_point_xy_annotate(rc,ax,rs))
            if rmv[1]:
                for re in rmv[1]:
                    if self.vs[re] is None:continue
                    rt = (rc,vec3(0,0,0).com(self.vs[re][2]['bound']))
                    ax = dtl.plot_edges(rt,ax,col = 'g')
            for exit in rmv[2]['exits']:
                ax = dtl.plot_point(exit,dtl.plot_point_xy_annotate(exit,ax,'exit'))
        if not fp is None:
            ax = dtl.plot_polygon(fp,ax,lw = 2,col = 'r')
        return ax

    ###########################################################################

    def defroom(self,b,es,**kws):
        def defkw(k,d):
            if not k in kws:kws[k] = d
        wcnt = len(b)
        defkw('bound',tuple(b));defkw('exits',[])
        defkw('floor',None);defkw('ceiling',None);defkw('shaft',False)
        defkw('level',0);defkw('skirt',0.5);defkw('crown',1.0)
        defkw('fholes',[]);defkw('choles',[])
        defkw('wholes',[[] for x in range(wcnt)])
        defkw('wtypes',['e' for x in range(wcnt)])
        defkw('wheights',[5 for x in range(wcnt)])
        defkw('wwidths',[0.2 for x in range(wcnt)])
        return self.av(es,kws)

    ###########################################################################

    def __init__(self,*ags,**kws):
        self.vs = []
        self.vcnt = 0
        self.fl_look = {}

        self._def('stack',[],**kws)
        self._def('levels',1,**kws)
        self._def('footprint',None,**kws)
        self._def('floorheight',5,**kws)
        self._def('shafts',[],**kws)

    ###########################################################################

    # sup is a room - split it into 2 rooms and return the new one
    def splitr(self,sup,connect = True,**kws):
        sv = self.vs[sup]
        bd = sv[2]['bound']
        l,r = splitb(bd,**kws)
        sv[2]['bound'] = tuple(l)
        newes = sv[1][:]+([sup] if connect else [])
        new = self.defroom(r,newes,
            exits = sv[2]['exits'][:],level = sv[2]['level'],
            skirt = sv[2]['skirt'],crown = sv[2]['crown'])
        # add edge for existing room
        for vx in sv[1]:self.vs[vx][1].append(new)
        if connect:
            sv[1].append(new)
        return new
    
    ###########################################################################

    # nest another graph within self by replacing vertex rix
    def insert(self,rix,o,connect = True):
        nvs = self.attach(None,None,o)
        self.rv(rix,nvs,connect = connect)
        for vx in range(self.vcnt):self.verifyedges(vx)
        return nvs

    ###########################################################################

    def rgraph(self,fprint,ax = None,**kws):
        self.footprint = fprint
        rooms = [self.defroom(fprint,[],**kws)]
        while self.stack:
            rx,line = self.stack.pop(0)
            new = self.splitr(rx,pline = line)
            rooms.append(new)
        for vx in range(self.vcnt):self.verifyedges(vx)
        for vx in range(self.vcnt):self.verifyedges(vx)
        return rooms       

###############################################################################
###############################################################################



###############################################################################
### building factory to generate full building contexts from a footprint
###############################################################################

class blgfactory(dfa.factory):

    def __str__(self):return 'building factory:'
    def __init__(self,*ags,**kws):
        self._def('bclass',building,**kws)
        dfa.factory.__init__(self,*ags,**kws)
    def new(self,exits,lvls,footprints,
            fplanstacks,topocuts,topoadds,shafts,*ags,**kws):

        skirt,crown,wheight,wwidth = 1,1,5,0.25

        floorheight = skirt+crown+wheight
        floorheights = [floorheight for x in range(lvls)]

        blgg = blggraph(stack = [],levels = lvls,
            floorheight = floorheight,footprint = footprints[0])

        shafttopo = []
        floorslots = []
        for lvx in range(lvls):
            lvbelow = [lvx-1] if lvx > 0 else []

            fp,fplan,fshafts = footprints[lvx],fplanstacks[lvx],shafts[lvx]
            ftopocuts,ftopoadds = topocuts[lvx],topoadds[lvx]
            wcnt = len(fp)
            whs = [wheight for x in range(wcnt)]
            wws = [wwidth for x in range(wcnt)]

            floorslots.append(blgg.defroom(fp,lvbelow,
                exits = exits,level = lvx,wheights = whs,wwidths = wws,
                skirt = skirt,crown = crown))

            floor = blggraph(stack = fplan)
            apts = floor.rgraph(fp,exits = exits,level = lvx,skirt = 1,crown = 1)
            for ftc in ftopocuts:floor.re(*ftc)
            for fta in ftopoadds:floor.ae(*fta)

            vc = blgg.vcnt
            blgg.insert(floorslots[-1],floor)

            # HERE I ONLY KEEP THE FIRST SHAFT... FIX THIS LATER
            shafttopo.append(fshafts[0]+vc)
            if lvx == 0:
                floor.plotxy(floor.vs[fshafts[0]][2]['bound'])
                plt.show()

        for lvx in range(1,lvls):
            blgg.ae(shafttopo[lvx-1],shafttopo[lvx])
            #blgg.mv(shafttopo[lvx-1],shafttopo[lvx])
        blgg.shafts.append(shafttopo)

        blgg.plot()
        plt.show()

        kws['bgraph'] = blgg
        n = self.bclass(*ags,**kws)
        return n

###############################################################################
###############################################################################
###############################################################################





