import dilap.core.base as db
import dilap.core.lsystem as lsy
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

# create a hole in the floor of room1 and the ceiling of room2
def roofhole(r1,r2):
    f1 = r1['floor']
    fh = vec3(0,0,0).com(f1[0]).sq(4,4)
    r1['floor'] = (f1[0],(fh,))
    c2 = r2['ceiling']
    ch = vec3(0,0,0).com(c2[0]).sq(4,4)
    r2['ceiling'] = (c2[0],(ch,))

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



#def boundlevel(self,m,lvl):
def genshell(m,rmkws,floors):

    rims = []
    for lvl in range(floors):

        tm = m.agfxmesh()

        bnd,mbnd = [],[]

        for rmkws in rmkws:
            #if not rmkws['level'] == lvl:continue
            for wx in range(len(rmkws['wtypes'])):
                wt,wh,ww = rmkws['wtypes'][wx-1],rmkws['wheight'],rmkws['wwidth']
                wp1,wp2 = rmkws['bound'][wx-1],rmkws['bound'][wx]
                sh,ch = rmkws['skirt'],rmkws['crown']
                hs = rmkws['wholes'][wx-1]
                if wt == 'e':
                    if not hs:hs.append((0,0.5,1.5,2,1))

                    bnd.append((wp1,wp2))

                    wtn = wp1.tov(wp2).crs(vec3(0,0,1)).nrm().uscl(ww)
                    wp1 = wp1.cp().trn(wtn)
                    wp2 = wp2.cp().trn(wtn)

                    mbnd.append((wp1,wp2))

                    wpys,portals = awall(wp1,wp2,wh,hs)
                    for wpy in wpys:m.asurf(wpy,tm)

                    rmkws['wholepoints'][wx-1].extend(portals)
                    for hps in portals:
                        wpys = aportal(hps,wp1,wp2,2*ww)
                        wpys = aportal(hps,wp1.cp(),wp2.cp(),2*ww)
                        for wpy in wpys:m.asurf(wpy,tm)
                    sp1,sp2 = wp1.cp().ztrn(-sh),wp2.cp().ztrn(-sh)
                    cp1,cp2 = wp1.cp().ztrn(wh),wp2.cp().ztrn(wh)

                    wpys,portals = awall(sp1,sp2,sh,[])
                    for wpy in wpys:m.asurf(wpy,tm)
                    wpys,portals = awall(cp1,cp2,ch,[])
                    for wpy in wpys:m.asurf(wpy,tm)

        tbnd = [(bnd[0],mbnd[0])]

        #pdb.set_trace()
        
        lll = 0
        while not len(tbnd) == len(bnd):
            lst = tbnd[-1][0][1]
            for mbe,be in zip(mbnd,bnd):
                lll += 1
                if lll > 1000:
                    print('oh dear god')
                    break
                if (be,mbe) in tbnd:continue
                if lst.isnear(be[0]):
                    tbnd.append((be,mbe))
                    lst = tbnd[-1][0][1]
            if lll > 1000:
                print('oh dear god')
                break

        rims.append([])
        for x in range(len(tbnd)):
            p1,p2 = tbnd[x-1][1]
            p3,p4 = tbnd[x][1]
            if p2.isnear(p3):continue

            cpy = (
                (p2.cp().ztrn(-sh),p3.cp().ztrn(-sh),
                p3.cp().ztrn(wh+ch),p2.cp().ztrn(wh+ch)),())
            m.asurf(cpy,tm)

            rims[-1].append(cpy[0][3].cp())
            rims[-1].append(cpy[0][2].cp())

    return rims

###############################################################################
###############################################################################
###############################################################################





###############################################################################
### context representing any building
###############################################################################

class building(cx.context):

    def __init__(self,*ags,**kws):
        self._def('name','buildingcontext',**kws)
        self._def('bgraph',None,**kws)
        cx.context.__init__(self,*ags,**kws)

    # insert a hole specification for a particular wall
    def wallhole(self,r,w,h):
        hp = h[1]
        whs = self.bgraph.vs[r][2]['wholes'][w]
        for wx in range(len(whs)):
            if hp > whs[wx][1]:
                whs.insert(wx,h)
                return
        whs.append(h)

    # add information to make doors for this room based on complete topology
    def wallholes(self,rx):
        dw,dh,dz = 1.5,3.0,0.0
        rv = self.bgraph.vs[rx]
        rb = rv[2]['bound']

        rexits = rv[2]['exits']

        # find exits which connect interior rooms and fix topology accordingly
        # exits appearing on exterior walls just become doors

        for adj in self.bgraph.vs[rx][1]:
            ov = self.bgraph.vs[adj]
            if ov is None:continue
            ob = ov[2]['bound']
            aws = adjacentwalls(rb,ob,minovlp = dw+2)
            for aw in aws:
                rwx,owx = aw
                ips = pym.sintsxyp(rb[rwx-1],rb[rwx],ob[owx-1],ob[owx])

                '''#
                for rex in rexits:
                    if gtl.onseg_xy(rex,ips[0],ips[1]):
                        print('internal exit!',rex)


                        bcom = vec3(0,0,0).com(rb+ob)
                        brad = 25
                        ax = dtl.plot_axes_xy(brad,(bcom.x,bcom.y))
                        ax = self.bgraph.plot(ax = ax)
                        ax = dtl.plot_polygon_xy(rb,ax,lw = 3,col = 'r')
                        ax = dtl.plot_polygon_xy(ob,ax,lw = 3,col = 'b')
                        ax = dtl.plot_edges_xy(ips,ax,lw = 5,col = 'g')
                        ax = dtl.plot_point_xy(rex,ax,col = 'r')
                        plt.show()
                '''#

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

        for wx in range(len(rv[2]['wtypes'])):
            wp1,wp2 = rb[wx-1],rb[wx]
            wt = rv[2]['wtypes'][wx-1]
            if wt == 'e':
                isexit = False
                for rex in rexits:
                    if gtl.onseg_xy(rex,wp1,wp2) and not wp2.isnear(rex):
                        print('external exit!',rex)
                        self.wallhole(rx,wx-1,('exit',0.5,dw,dh,dz))
                        isexit = True
                        break
                if not isexit:
                    self.wallhole(rx,wx-1,('window',0.5,dw,dh-1,dz+1))

    # add the trimeshes of a room interior to a model
    def genrooms(self):
        for v in self.bgraph.vs:
            if v is None:continue
            vx,ves,vkw = v

            # insert doors based on topology
            self.wallholes(vx)

            # add a new trimesh to the model
            m = dmo.model()
            sgv = self.amodel(None,None,None,m,None)
            tm = m.agfxmesh()

            # create the floor and ceiling surfaces
            floor,ceiling = vkw['floor'],vkw['ceiling']
            if floor is None or ceiling is None:
                fbnd = pym.contract(vkw['bound'],min(vkw['wwidths']))
                floor   = (tuple([p.cp() for p in fbnd]),())
                wh = min(vkw['wheights'])
                ceiling = (tuple([p.cp().ztrn(wh) for p in fbnd]),())
            m.asurf(floor,tm)
            m.asurf(ceiling,tm,rv = True)

            # create the wall surfaces and portal surfaces
            for x in range(len(floor[0])):
                b2,b1 = vkw['bound'][x-1],vkw['bound'][x]
                w2,w1 = floor[0][x-1],floor[0][x]
                hs,wt = vkw['wholes'][x-1],vkw['wtypes'][x-1]
                wh,ww = vkw['wheights'][x-1],vkw['wwidths'][x-1]

                # create wall surfaces
                wpys,portals = awall(w1,w2,wh,hs,b1.cp(),b2.cp())
                for wpy in wpys:m.asurf(wpy,tm)

                # create portal surfaces
                for hps in portals:
                    wpys = aportal(hps,w1,w2,ww)
                    for wpy in wpys:m.asurf(wpy,tm)

    # do something which fills the scenegraph
    def generate(self,worn = 0):
        self.floors = len(list(self.bgraph.fl_look.keys()))

        self.genrooms()

        #self.bgraph.plot(self.bgraph.vs[-1][2]['bound'])

        #self.genshell()
        #self.genroof()

        # create the models which cover the exterior walls
        #m = dmo.model()
        #sgv = self.amodel(p,q,s,m,None)
        #self.rims = genshell(m,self.rmkws,self.floors)
                  
        # create the models which cover the roof
        #m = dmo.model()
        #sgv = self.amodel(p,q,s,m,None)
        #genroof(m,self.rmkws)

        return self

###############################################################################
### graph structure to create building contexts from
###############################################################################

'''#
# split a boundary polygon into two boundary polygons using a single line segment
def splitline(b,mindx = 8,mindy = 8,pos = None,ori = None):
    x,y = vec3(1,0,0),vec3(0,1,0)
    bx,by = x.prjps(b),y.prjps(b)
    bxd = bx[1]-bx[0]
    byd = by[1]-by[0]

    if pos is None:
        dx = random.uniform(0,bxd/2.0-mindx)
        dy = random.uniform(0,byd/2.0-mindy)
        dx = 0
        dy = 0
        spos = vec3(dx,dy,0).com(b)
    elif type(pos) is type(0.0):
        minbx,maxbx = bx[0]+mindx,bx[1]-mindx
        dx = minbx + (maxbx - minbx)*pos
        dy = minby + (maxby - minby)*pos
        spos = vec3(dx,dy,0)

    if ori is None:bdir = y.cp() if bxd > byd else x.cp()
    elif ori == 'v':bdir = y.cp()
    elif ori == 'h':bdir = x.cp()
    else:raise ValueError

    bdir.uscl(100)

    s1,s2 = spos.cp().trn(bdir.flp()),spos.cp().trn(bdir.flp())

    ax = dtl.plot_axes_xy(50)
    ax = dtl.plot_polygon_xy(b,ax)
    ax = dtl.plot_edges_xy((s1,s2),ax,col = 'r')
    plt.show()

    l,r = pym.bsegsxy(b,s1,s2)
    return l,r
'''#

def splitb(b,pline = None,mindx = 8,mindy = 8,pos = None,ori = None):
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
        #return splitline(b)
    else:
        print('split polyline?')
        raise NotImplemented
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
        for res in ev[1]:
            if self.vs[res] is None:continue
            for nvx in nvs:
                nv = self.vs[nvx]
                if connect is True:
                    if not res in nv[1]:nv[1].append(res)
                    if not nv[0] in self.vs[res][1]:self.vs[res][1].append(nv[0])
                elif type(connect) is type(()) and connect[0] == nvx:
                    print('insert connection!',connect)
                #elif False:nv[2]['exits'] = connect[:]
                else:
                    nv[2]['exits'] = self.vs[rx][2]['exits'][:]
        self.vs[rx] = None
        
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

    def __init__(self,*ags,**kws):
        self.vs = []
        self.vcnt = 0

        # fl_look is a lookup of vertices per floor of the building
        self.fl_look = {}
        self._def('rstack',[],**kws)

    def plot(self,fp = None,ax = None):
        if ax is None:ax = dtl.plot_axes_xy(50)
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

    # consider all rooms that room rx should be connected to based on edges
    #   if it is no longer tangential to a connected room, disconnect them
    #   if it is still tangential to a connected room, confirm their portals 
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

    # sup is a room - split it into 2 rooms and return the new one
    def splitr(self,sup,**kws):
        sv = self.vs[sup]
        bd = sv[2]['bound']

        # compute new boundaries
        #l,r = splitb(bd,pline,8,8,None,None)
        l,r = splitb(bd,**kws)

        # set existing boundary
        sv[2]['bound'] = tuple(l)

        # create new room
        new = self.defroom(r,sv[1][:]+[sup],exits = sv[2]['exits'][:])

        # add edge for existing room
        for vx in sv[1]:self.vs[vx][1].append(new)
        sv[1].append(new)

        # return the index of the new resulting room
        return new
    
    # attach another graph to this one based on geometry
    # srx/orx are the indices of the adjacent rooms in each graph
    def attach(self,srx,orx,o):
        # all edges indices in o have to be updated by +vc
        # each room in o must be added to self
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

        #for vx in range(self.vcnt):
        #    self.verifyedges(vx)

        return nvs

    # nest another graph within self by replacing vertex rix
    def insert(self,rix,o,connect = True):
        #rifp = self.vs[rix][2]['bound']
        #exits = self.vs[rix][2]['exits']
        #o.rgraph(rifp,exits)
        #if type(connect) is type([]):
        #    exits += connect
        #    connect = False
        nvs = self.attach(None,None,o)
        self.rv(rix,nvs,connect = connect)
        for vx in range(self.vcnt):self.verifyedges(vx)
        return nvs

    def defroom(self,b,es,**kws):
        def defkw(k,d):
            if not k in kws:kws[k] = d
        wcnt = len(b)
        defkw('bound',tuple(b));defkw('exits',[])
        defkw('floor',None);defkw('ceiling',None)
        defkw('level',0);defkw('skirt',0.5);defkw('crown',1.0)
        defkw('wholes',[[] for x in range(wcnt)])
        defkw('wtypes',['e' for x in range(wcnt)])
        defkw('wheights',[5 for x in range(wcnt)])
        defkw('wwidths',[0.2 for x in range(wcnt)])
        return self.av(es,kws)

    grammer = {
        'S':self.seqsplit,
            }

    def seqsplit(self,tip,kws):
        line = {'pos':kws['pos'],'ori':kws['ori']}
        news = []
        for tr in tip:
            news.append(self.splitr(tr,**line))
        rooms.extend(news)
        if nxtsch == 'ext':tip.extend(news)
        elif nxtsch == 'new':tip = news
        return tip

    # given a footprint boundary polygon in the xy plane
    # construct a graph of rooms which fit within the footprint
    #   exits are locations on the boundary of footprint
    #   where a door will result in whichever room ends up there
    def rgraph(self,footprint = None,exits = None,rseq = '[0]SS'):
        if footprint is None:footprint = vec3(0,0,0).sq(25,25)
        if not type(exits) == type([]):exits = []
        exits = [e for e in exits if e.onbxy(footprint)]
        rooms = [self.defroom(footprint,[],exits = exits)]
        #lst = rooms[-1]

        def readto(c):
            closex = rseq[cx:].find(c)+cx
            read = rseq[cx+1:closex]
            print('read',read)
            return read,closex

        rlen = len(rseq)
        nxtsch = 'ext'
        connectinserted = False
        tip = []                    # tip holds the set of rooms being operated on
        cx = 0                      # cx is the position of the cursor within rseq
        p = 0.5                     # p is the position used for the split operation
        o = 'v'
        while cx < rlen:
            c = rseq[cx]

            if c == '[':            # set the set of rooms being operated on
                tip,cx = readto(']')
                tip = [int(v) for v in tip.split(',')]

            elif c == '|':          # toggle the nxt scheme for how to set the operating set after an operation
                nxtsch = 'ext' if nxtsch == 'new' else 'new'
                tip = []

            elif c == '-':          # toggle the cnt scheme for how to connect inserted graphs to existing topology
                connectinserted = True if connectinserted is False else False

            elif c == 'I':          # insert a subgraph
                subrseq,cx = readto('>')
                news = []
                for tr in tip:
                    subg = blggraph()
                    print('trrr',tr)
                    rifp = self.vs[tr][2]['bound']
                    exits = self.vs[tr][2]['exits']
                    subg.rgraph(rifp,exits,rseq = subrseq)
                    news.extend(self.insert(tr,subg,connectinserted))
                if nxtsch == 'ext':tip.extend(news)
                elif nxtsch == 'new':tip = news
                
            elif c == 'S':          # split each room in the operating set
                line = {'pos':p,'ori':o}
                news = []
                for tr in tip:
                    news.append(self.splitr(tr,**line))
                rooms.extend(news)
                if nxtsch == 'ext':tip.extend(news)
                elif nxtsch == 'new':tip = news

            elif c == 'R':          # update the cut position randomly
                p += random.uniform(-0.5,0.5)
                print('randp',p)

            elif c == '{':          # update the cut position/orientation
                subrseq,cx = readto('}')
                subrseq = subrseq.split(',')
                p = float(subrseq[0])
                o = subrseq[1]
                if o == 'a':o = None
                print('setpo',p,o)

            cx += 1
            
            #lst = new

        for vx in range(self.vcnt):self.verifyedges(vx)
        for vx in range(self.vcnt):self.verifyedges(vx)
        return rooms





        '''#
        rstack = self.rstack[:]
        while rstack:

            #self.plot(self.vs[lst][2]['bound'])
            #plt.show()

            tos,line = rstack.pop(0)
            new = self.splitr(tos,**line)
            rooms.append(new)
            lst = new

        for vx in range(self.vcnt):
            self.verifyedges(vx)

        #self.plot(self.vs[new][2]['bound'])
        #plt.show()

        return rooms
        '''#

###############################################################################
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

    # 
    def test(self,footprint,exits,*ags,**kws):

        blgg = blggraph(rstack = [(0,{'pos':0.5})]).rgraph(footprint,exits)

        blg1 = blggraph(rstack = [
            (0,{'pos':0.5,'ori':'h'}),
            (1,{'pos':0.5,'ori':'v'}),
            (2,{'pos':0.5,'ori':'h'}),
            (0,{'pos':0.5,'ori':'v'})])
        blgg.insert(0,blg1,[])

        blg2 = blggraph(rstack = [
            (0,{'pos':0.5,'ori':'v'}),
            (1,{'pos':0.5,'ori':'h'}),
            (2,{'pos':0.5,'ori':'v'}),
            (0,{'pos':0.5,'ori':'h'})])
        blgg.insert(1,blg2,[])

        kws['bgraph'] = blgg
        n = self.bclass(*ags,**kws)
        return n

    # 
    def hotel(self,footprint,exits,*ags,**kws):

        aptl = '{0.5,h}[0]S{0.5,v}S[0]S'
        aptr = '{0.5,h}[0]S{0.5,v}S[0]S'

        aptl = '[0]S'
        aptr = '[0]S'
        rseq = '[0]SS-[0]I<'+aptl+'>'+'[3]I<'+aptr+'>'
        rseq = '{0.5,h}[0]S[0]I<'+rseq+'>[1]I<'+rseq+'>'

        print('rseq',rseq)

        hall = blggraph()
        hallvs = hall.rgraph(footprint,exits,rseq)

        kws['bgraph'] = hall

        n = self.bclass(*ags,**kws)
        return n

    # 
    def new(self,footprint,exits,*ags,**kws):
        n = self.hotel(footprint,exits,*ags,**kws)
        n.bgraph.plot()
        plt.show()
        return n

###############################################################################
###############################################################################
###############################################################################





