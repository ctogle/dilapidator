from dilap.core import *
from dilap.geometry import *
import dilap.geometry.tools as gtl
import dilap.geometry.triangulate as dtg
import dilap.geometry.polymath as pym
import dilap.topology.trimesh as dtm
import dilap.worldly.blgsequencing as bseq
import dilap.worldly.polygen as pyg
import dilap.core.plotting as dtl
import matplotlib.pyplot as plt
import numpy
import math
import random
import pdb


dkw = lambda kws,k,d : kws[k] if k in kws else d

class building(scenegraph):

    def __init__(self,p,q,s,**kws):
        self.name = kws.get('name','building')
        self.bgraph = kws.get('bgraph',None)
        self.p = p
        self.q = q
        self.s = s
        scenegraph.__init__(self)

    def windowholes(self,r,w,wp1,wp2,**wopts):
        wt = dkw(wopts,'type','e')
        ww = dkw(wopts,'ww',2.5)
        wh = dkw(wopts,'wh',2.0)
        wz = dkw(wopts,'wz',1.0)
        wm = dkw(wopts,'wm',1.0)
        ws = dkw(wopts,'ws',4)
        wpd = wp1.d(wp2)
        wcnt = int((wpd-1)/(ww+ws*wm))
        if wcnt < 1:return
        wpoff = 1.0/wcnt
        for wx in range(wcnt):
            wp = (wx+1.0)/wcnt-wpoff/2
            self.wallhole(r,w,('window',wp,ww,wh,wz))

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
        dw,dh,dz = 1.5,3,0.0
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
                aws = pym.badjbxy(rb,ob,minovlp = dw+2)
                for aw in aws:
                    rwx,owx = aw
                    ips = pym.sintsxyp(rb[rwx-1],rb[rwx],ob[owx-1],ob[owx])
                    imp = ips[0].mid(ips[1])
                    ridp = 1-rb[rwx-1].d(imp)/rb[rwx-1].d(rb[rwx])
                    if rv[2]['rtype'] == 'open' and ov[2]['rtype'] == 'open':
                        tdh = dh*100
                        tww = max(rv[2]['wwidths'][rwx-1],ov[2]['wwidths'][rwx-1])
                        tdw = ips[0].d(ips[1])-2*tww
                        # hack door that consumes the wall
                        self.wallhole(rx,rwx-1,(adj,ridp,tdw,tdh,dz)) 
                    else:
                        # typical door that does not consume the wall
                        self.wallhole(rx,rwx-1,(adj,ridp,dw,dh,dz))
                    #rv[2]['wtypes'][rwx-1] = 'i'
                    rv[2]['wmetas'][rwx-1]['type'] = 'i'
        gadjs = self.bgraph.geoadj(rx)
        for adj in gadjs:
            ov = self.bgraph.vs[adj[0]]
            if ov is None:continue
            ob = ov[2]['bound']
            aws = pym.badjbxy(rb,ob,minovlp = 0)
            for aw in aws:
                rwx,owx = aw
                #if rv[2]['wtypes'][rwx-1] == 'e':
                #    rv[2]['wtypes'][rwx-1] = 'i'
                if rv[2]['wmetas'][rwx-1]['type'] == 'e':
                    rv[2]['wmetas'][rwx-1]['type'] = 'i'
        hasexit = False
        for wx in range(len(rb)):
            wp1,wp2 = rb[wx-1],rb[wx]
            #wt = rv[2]['wtypes'][wx-1]
            wt = rv[2]['wmetas'][wx-1]['type']
            if wt == 'e':
                isexit = False
                if rv[2]['level'] == 0:
                    for rex in rexits:
                        if rex is True:
                            if not hasexit:
                                self.wallhole(rx,wx-1,('exit',0.5,dw,dh,dz))
                                hasexit = True
                                isexit = True
                                break
                        elif rex.onsxy(wp1,wp2) and not wp2.isnear(rex):
                        #elif gtl.onseg_xy(rex,wp1,wp2) and not wp2.isnear(rex):
                            print('external exit!',rex)
                            self.wallhole(rx,wx-1,('exit',0.5,dw,dh,dz))
                            isexit = True
                            break
                wpd = wp1.d(wp2)
                #if not isexit and wpd > dw+2 and not rv[2]['rtype'] == 'closed':
                if wpd > dw+2 and not rv[2]['rtype'] == 'closed':
                    self.windowholes(rx,wx-1,wp1,wp2,**rv[2]['wmetas'][wx])

    def genshafts(self):

        def genwrap(shb,shx,skirt = True,crown = True):
            vkw = self.bgraph.vs[shvx][2]
            fbnd = pym.contract(shb,wwi)
            for x in range(len(shb)):
                b2,b1,w2,w1 = shb[x-1],shb[x],fbnd[x-1],fbnd[x]
                #hs,wt = vkw['wholes'][x-1],vkw['wtypes'][x-1]
                hs,wt = vkw['wholes'][x-1],vkw['wmetas'][x-1]['type']
                wh,ww = vkw['wheights'][x-1],vkw['wwidths'][x-1]
                sh,ch = vkw['skirt'],vkw['crown']
                wpys,portals = pyg.awall(w1,w2,wh,hs,b1.cp(),b2.cp())
                for wpy in wpys:m.asurf(wpy,tm)
                for hps in portals:
                    wpys = pyg.aportal(hps,w1,w2,ww)
                    for wpy in wpys:m.asurf(wpy,tm)
                if skirt:
                    sp1,sp2 = w1.cp().ztrn(-sh),w2.cp().ztrn(-sh)
                    wpys,portals = pyg.awall(sp1,sp2,sh,[])
                    for wpy in wpys:m.asurf(wpy,tm)
                if crown:
                    cp1,cp2 = w1.cp().ztrn(wh),w2.cp().ztrn(wh)
                    wpys,portals = pyg.awall(cp1,cp2,ch,[])
                    for wpy in wpys:m.asurf(wpy,tm)

        def genplatform(shb,shx):
            buff,rw = 6,6
            fbnd = pym.contract(shb,wwi)
            platform = [p.cp() for p in fbnd]
            r1,r2 = platform[2].cp(),platform[3].cp()
            rtn = r1.tov(r2).nrm().uscl(buff)
            #rnm = vec3(0,0,1).crs(rtn).nrm().uscl(rw)
            rnm = pym.bnrm(fbnd).crs(rtn).nrm().uscl(rw)
            r2,r1 = r2.trn(rtn.flp()),r1.trn(rtn.flp())
            r3,r4 = r2.cp().trn(rnm),r1.cp().trn(rnm)

            belv = self.bgraph.vs[sh[shx-1]][2]
            rh = max(belv['wheights'])+belv['skirt']+belv['crown']
            ramp = [r4,r1,r2,r3]

            platform = pym.ebdxy(platform,ramp)
            #print(len(platform))
            platform = platform[len(platform)//2] # HACK TO GET CORRECT POLYGON...

            ramp[0].ztrn(-rh)
            ramp[1].ztrn(-rh)

            m.asurf((ramp,()),tm)
            m.asurf((platform,()),tm)

        def genbottom(shb):
            shaftfloor = pym.contract(shb,wwi)
            m.asurf((shaftfloor,()),tm)
            genwrap(shb,shx-1,skirt = False,crown = shcnt > 1)

        def gentop(shb):
            shaftceiling = [p.cp().ztrn(wht) for p in pym.contract(shb,wwi)]
            m.asurf((shaftceiling,()),tm,rv = True)
            genplatform(shb,shx)
            genwrap(shb,shx+1,crown = False)

        def genlink(shb,shx):
            genplatform(shb,shx)
            genwrap(shb,shx)

        for sh in self.bgraph.shafts:
            sg,m = self.sgraph,model()
            sgv = sg.avert(None,None,None,models = [m],parent = sg.root)
            tm = m.agfxmesh()

            shcnt = len(sh)
            for shx in range(shcnt):
                shvx = sh[shx]
                shv = self.bgraph.vs[shvx]
                shb = shv[2]['bound']
                wht = max(shv[2]['wheights'])
                wwi = min(shv[2]['wwidths'])

                if shx == 0:genbottom(shb)
                elif shx == shcnt-1:gentop(shb)
                else:genlink(shb,shx)

    # add the trimeshes of a room interior to a model
    def genrooms(self):
        for v in self.bgraph.vs:
            if v is None:continue
            vx,ves,vkw = v

            # insert doors based on topology
            self.holes(vx)

            # dont build the room here if it contains a shaft
            if vkw['shaft']:continue

            # add a new trimesh to the model
            sg,m = self.sgraph,model()
            sgv = sg.avert(None,None,None,models = [m],parent = sg.root)
            tm = m.agfxmesh()

            # create the floor and ceiling surfaces
            #fbnd = pym.contract(vkw['bound'],max(vkw['wwidths']))
            fbnd = pym.contract(vkw['bound'],vkw['wwidths'])
            # need to place fbnd based on differing wwidths!!

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
                hs,wt = vkw['wholes'][x-1],vkw['wmetas'][x-1]['type']
                wh,ww = vkw['wheights'][x-1],vkw['wwidths'][x-1]

                # create wall surfaces
                wpys,portals = pyg.awall(w1,w2,wh,hs,b1.cp(),b2.cp())
                for wpy in wpys:m.asurf(wpy,tm)

                # create portal surfaces
                for hps in portals:
                    wpys = pyg.aportal(hps,w1,w2,ww)
                    for wpy in wpys:m.asurf(wpy,tm)

    def genshell(self):
        sg,m = self.sgraph,model()
        sgv = sg.avert(None,None,None,models = [m],parent = sg.root)
        tm = m.agfxmesh()
        z,eww,ech,fh = 0,0.2,0.5,self.bgraph.floorheight

        self.rims = []
        for lvx in range(self.floors):
            fp = [p.cp().ztrn(z) for p in self.bgraph.footprint]
            z += self.bgraph.floorheight
            self.rims.append([])
            for fpx in range(len(fp)):
                fp0,fp1,fp2 = fp[fpx-2],fp[fpx-1],fp[fpx]
                wtn1 = fp0.tov(fp1).nrm()
                wtn2 = fp1.tov(fp2).nrm()
                tcrstz = gtl.near(wtn1.crs(wtn2).z,0)
                if tcrstz > 0:
                    wnm1 = vec3(0,0,1).crs(wtn1).nrm().uscl(-eww)
                    wnm2 = vec3(0,0,1).crs(wtn2).nrm().uscl(-eww)
                    crnp1 = fp1.cp().trn(wnm1).ztrn(-ech)
                    crnp2 = fp1.cp().trn(wnm2).ztrn(-ech)
                    wpys,portals = pyg.awall(crnp1,crnp2,fh,[])
                    for wpy in wpys:m.asurf(wpy,tm)
                    self.rims[-1].extend([crnp1.cp().ztrn(fh),crnp2.cp().ztrn(fh)])
                elif tcrstz < 0:
                    wnm1 = vec3(0,0,1).crs(wtn1).nrm().uscl(-eww)
                    wnm2 = vec3(0,0,1).crs(wtn2).nrm().uscl(-eww)
                    crnp0 = fp0.cp().trn(wnm1).ztrn(-ech)
                    crnp1 = fp1.cp().trn(wnm1).ztrn(-ech)
                    crnp2 = fp1.cp().trn(wnm2).ztrn(-ech)
                    crnp3 = fp2.cp().trn(wnm2).ztrn(-ech)
                    ip = pym.sintsxyp(crnp0,crnp1,crnp2,crnp3)
                    self.rims[-1].append(ip.cp().ztrn(fh))
                for v in self.bgraph.vs:
                    if v is None:continue
                    vx,ves,vkw = v
                    vb = vkw['bound']
                    for vbx in range(len(vb)):
                        vb1,vb2 = vb[vbx-1],vb[vbx]
                        ips = pym.sintsxyp(fp1,fp2,vb1,vb2,ie = False,skew = False)
                        if not type(ips) is type(()):continue
                        hs = vkw['wholes'][vbx-1]
                        sh,ch = vkw['skirt'],vkw['crown']
                        wh = self.bgraph.floorheight-sh-ch
                        wtn = vb1.tov(vb2).crs(vec3(0,0,1)).nrm().uscl(eww)
                        wp1,wp2 = vb1.cp().trn(wtn),vb2.cp().trn(wtn)
                        wpys,portals = pyg.awall(wp1,wp2,wh,hs)
                        for wpy in wpys:m.asurf(wpy,tm)
                        for hps in portals:
                            wpys = pyg.aportal(hps,wp1,wp2,eww)
                            for wpy in wpys:m.asurf(wpy,tm)
                        sp1,sp2 = wp1.cp().ztrn(-sh),wp2.cp().ztrn(-sh)
                        cp1,cp2 = wp1.cp().ztrn( wh),wp2.cp().ztrn( wh)
                        wpys,portals = pyg.awall(sp1,sp2,sh,[])
                        for wpy in wpys:m.asurf(wpy,tm)
                        wpys,portals = pyg.awall(cp1,cp2,ch,[])
                        for wpy in wpys:m.asurf(wpy,tm)

    # generate a cover of the top of the top of the building
    def genroof(self):
        sg,m = self.sgraph,model()
        sgv = sg.avert(None,None,None,models = [m],parent = sg.root)
        tm = m.agfxmesh()
        m.asurf((tuple(self.rims[-1]),()),tm)

    # do something which fills the scenegraph
    def generate(self,worn = 0):
        self.floors = len(list(self.bgraph.fl_look.keys()))
        skirt = min([v[2]['skirt'] for v in self.bgraph.vs if v])
        if self.p is None:self.p = vec3(0,0,skirt)
        else:self.p.ztrn(skirt)
        self.genrooms()
        self.genshafts()
        self.genshell()
        self.genroof()
        return self

###############################################################################
###############################################################################



###############################################################################
### graph structure to create building contexts from
###############################################################################

class blggraph:

    def av(self,os,kws,null = False):
        vx = self.vcnt
        self.vcnt += 1
        if null:self.vs.append(None)
        else:
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
            ringv = self.vs[res]
            if ringv is None:continue
            for nvx in nvs:
                nv = self.vs[nvx]
                if True in nv[2]['exits']:
                    self.ae(nvx,res)
                    nv[2]['exits'].remove(True)
                    break
        self.vs[rx] = None
        
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
            adjwalls = pym.badjbxy(bd,ab)
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
            adjwalls = pym.badjbxy(rv[2]['bound'],av[2]['bound'])
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
            if ov[2]['level'] != rv[2]['level']:
                if not self.shafts:
                    self.shafts.append([adj if ov[2]['level'] < rv[2]['level'] else rx])
                for sh in self.shafts:
                    if adj in sh and not rx in sh:
                        sh.append(rx)
                        break
                continue
            ob = ov[2]['bound']
            aws = pym.badjbxy(rb,ob)
            if not aws:
                if adj in rv[1]:rv[1].remove(adj)
                if rx in ov[1]:ov[1].remove(rx)
        if rv[2]['shaft']:
            #print('found shaft',rx)
            if not self.shafts:
                self.shafts.append([rx])
            else:
                #print('need to id correct shaft...')
                self.shafts[-1].append(rx)
    
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
                if ov is None:
                    self.av(None,None,null = True)
                    continue
                okws = ov[2]
                newring = [ovx+vc for ovx in ov[1]]
                new = self.av(newring,okws)
                nvs.append(new)
        return nvs

    ###########################################################################

    def plotxy(self,fp = None,ax = None):
        if ax is None:ax = dtl.plot_axes_xy(50)
        for rmv in self.vs:
            if rmv is None:continue
            if rmv[2]['level'] > 0:continue
            bcol = 'b' if rmv[2]['shaft'] else None
            rconbs = 0.5 if pym.bccw(rmv[2]['bound']) else -0.5
            rconb = pym.contract(rmv[2]['bound'],rconbs)
            ax = dtl.plot_polygon_xy(rconb,ax,lw = 4,col = bcol)
            rc = vec3(0,0,0).com(rmv[2]['bound'])
            rs = str(rmv[0])+','+str(rmv[1])
            ax = dtl.plot_point_xy(rc,dtl.plot_point_xy_annotate(rc,ax,rs))
            if rmv[1]:
                for re in rmv[1]:
                    if self.vs[re] is None:continue
                    rt = (rc,vec3(0,0,0).com(self.vs[re][2]['bound']))
                    ax = dtl.plot_edges_xy(rt,ax,col = 'g')
            for exit in rmv[2]['exits']:
                if exit is True:
                    exitp = rc
                    ax = dtl.plot_point_xy(exitp,
                        dtl.plot_point_xy_annotate(exitp,ax,'exit'))
                else:
                    ax = dtl.plot_point_xy(exit,
                        dtl.plot_point_xy_annotate(exit,ax,'exit'))
        if not fp is None:
            ax = dtl.plot_polygon_xy(fp,ax,lw = 2,col = 'r')
        return ax

    def plot(self,fp = None,ax = None):
        if ax is None:ax = dtl.plot_axes(50)
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
                if exit is True:
                    exitp = rc
                    ax = dtl.plot_point(exitp,
                        dtl.plot_point_xy_annotate(exitp,ax,'exit'))
                else:
                    ax = dtl.plot_point(exit,
                        dtl.plot_point_xy_annotate(exit,ax,'exit'))
        if not fp is None:
            ax = dtl.plot_polygon(fp,ax,lw = 2,col = 'r')
        return ax

    ###########################################################################

    def defroom(self,b,es,**kws):
        def defkw(k,d):
            if not k in kws:kws[k] = d
        wcnt = len(b)
        dwh = self.floorheight-1.0
        if pym.bnrm(b).z < 0:b.reverse()
        defkw('bound',tuple(b));defkw('exits',[]);defkw('rtype','room')
        defkw('floor',None);defkw('ceiling',None);defkw('shaft',False)
        defkw('level',0);defkw('skirt',0.5);defkw('crown',0.5)
        defkw('fholes',[]);defkw('choles',[])
        defkw('wholes',[[] for x in range(wcnt)])
        defkw('wmetas',[{'type':'e'} for x in range(wcnt)])
        #defkw('wtypes',['e' for x in range(wcnt)])
        defkw('wheights',[dwh for x in range(wcnt)])
        defkw('wwidths',[1.0 for x in range(wcnt)])
        return self.av(es,kws)

    ###########################################################################

    def __init__(self,*ags,**kws):
        self.vs = []
        self.vcnt = 0
        self.fl_look = {}

        self.sequence = kws.get('sequence', None)
        self.footprint = kws.get('footprint', None)
        self.floorheight = kws.get('floorheight', 8)
        self.shafts = kws.get('shafts', [])

    ###########################################################################

    # sup is a room - split it into 2 rooms and return the new one
    def splitr(self,sup,line,connect = True):
        sv = self.vs[sup]
        bd = sv[2]['bound']
        bprjx = vec3(1,0,0).prjps(bd)
        bprjy = vec3(0,1,0).prjps(bd)
        bprjz = vec3(0,0,1).prjps(bd)
        lp = vec3(
            bprjx[0]+line[0].x*(bprjx[1]-bprjx[0]),
            bprjy[0]+line[0].y*(bprjy[1]-bprjy[0]),
            bprjz[0]+line[0].z*(bprjz[1]-bprjz[0]))
        sp1 = lp.cp().trn(line[1].cp().uscl( 1000))
        sp2 = lp.cp().trn(line[1].cp().uscl(-1000))

        l,r = pym.bsegsxy(bd,sp1,sp2)
        sv[2]['bound'] = tuple(l if pym.bccw(l) else l[::-1])
        newes = sv[1][:]+([sup] if connect else [])
        new = self.defroom(r,newes,
            exits = sv[2]['exits'][:],level = sv[2]['level'],
            skirt = sv[2]['skirt'],crown = sv[2]['crown'])
        # add edge for existing room
        for vx in sv[1]:self.vs[vx][1].append(new)
        if connect:sv[1].append(new)
        return new

    # nest another graph within self by replacing vertex rix
    def insert(self,rix,o):
        nvs = self.attach(None,None,o)
        self.rv(rix,nvs)
        self.shafts = []
        for vx in range(self.vcnt):self.verifyedges(vx)
        return nvs

    ###########################################################################

    def graph(self,fp = None,**kws):
        if not fp is None:self.footprint = fp
        rooms = [self.defroom(self.footprint,[],**kws)]

        def sqrd(seq,sx):
            score = 1
            sx += 1
            while score > 0:
                sx += 1
                if seq[sx] == '<':
                    score += 1
                elif seq[sx] == '>':
                    if score > 0:
                        score -= 1
            return sx

        def level(subseq):
            lrx = int(subseq)
            lvkw = self.vs[lrx][2]
            fh = self.floorheight
            b = [p.cp().ztrn(fh) for p in lvkw['bound']]
            nr = self.defroom(b,[],level = lvkw['level']+1)
            self.ae(nr,lrx)

        def split(subseq):
            srx,spx,spy,spz,sdx,sdy,sdz = [float(v) for v in subseq.split(',')]
            new = self.splitr(int(srx),(vec3(spx,spy,spz),vec3(sdx,sdy,sdz)),False)
            rooms.append(new)

        def rexit(subseq):
            rx = int(subseq)
            self.vs[rx][2]['exits'].append(True)

        def edge(subseq):
            ss = subseq.split(',')
            self.ae(int(ss[0]),int(ss[1]))

        def rkset(subseq):
            ss = subseq.split(',')
            self.vs[int(ss[0])][2][ss[1]] = ss[2]

        def rins(subseq):
            irx = int(subseq[:subseq.find(',')])
            subgseq = subseq[subseq.find(',')+1:]
            kws = {
                'sequence':subgseq,
                'floorheight':self.floorheight,
                    }
            subg = blggraph(**kws)
            subg.graph(fp = self.vs[irx][2]['bound'])
            for sgv in subg.vs:
                if sgv is None:continue
                sgv[2]['level'] = self.vs[irx][2]['level']
            self.insert(irx,subg)

        def shaft(subseq):
            srx = int(subseq)
            self.vs[srx][2]['shaft'] = True

        def contract(subseq):
            ss = subseq.split(',')
            srx = int(ss[0])
            cd = [float(v) for v in ss[1:]]
            for wx in range(len(cd)):
                self.vs[srx][2]['wwidths'][wx] *= cd[wx]

        grammer = {
            'L':level,'S':split,'I':rins,'X':rexit,'E':edge,'V':shaft,'R':rkset,
            'C':contract,
                }
        seq = self.sequence[:]
        scnt = len(seq)
        sx = 0
        while sx < scnt:
            c = seq[sx]
            if c in grammer:
                ex = sqrd(seq,sx)
                grammer[c](seq[sx+2:ex])
            else:ex = sx+1
            sx = ex

        self.shafts = []
        for vx in range(self.vcnt):self.verifyedges(vx)
        return rooms       


def new(*ags,**kws):
    blgg = blggraph(**kws)
    blgg.graph()
    #blgg.plotxy()
    #plt.show()
    kws['bgraph'] = blgg
    n = building(*ags,**kws)
    return n
