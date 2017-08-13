from dilap.core import *
from dilap.geometry import *
from .polygen import awall, aportal
import dilap.geometry.tools as gtl
import dilap.geometry.triangulate as dtg
import dilap.geometry.polymath as pym
import dilap.topology.trimesh as dtm
from dilap.topology import wiregraph
from .partitiongraph import partition
import numpy
import math
import random
import pdb


class blggraph:


    def plot(self, l=100):
        ax = plot_axes(100)
        for v in self.topology.vs:
            vl = v[1]['pv'].loop
            vc = vec3(0,0,0).com(vl)
            ax = plot_polygon(vl, ax)
        esseen = []
        for u, v in self.topology.elook:
            if (u, v) in esseen or (v, u) in esseen:
                continue
            vl = self.topology.vs[u][1]['pv'].loop
            vc = vec3(0,0,0).com(vl)
            ol = self.topology.vs[v][1]['pv'].loop
            oc = vec3(0,0,0).com(ol)
            ax = plot_polygon(vl, ax, col='g')
            ax = plot_polygon(ol, ax, col='g')
            if vc.isnearxy(oc) and not gtl.isnear(vc.z,oc.z,0.1):
                ax = plot_polygon(pym.contract(vl, 1), ax, col='r')
                ax = plot_polygon(pym.contract(ol, 1), ax, col='r')
            ax = plot_edges([vc, oc], ax)
            esseen.append((u, v))
        return ax


    def ashaft(self,start,end,origin):
        '''Add a set of edges from floor start, at room origin, to floor end'''
        print('update wiregraph with shaft')
        bottom = self.floors[start]['partition'].verts[origin]
        below = bottom
        for l in range(start+1,end+1):
            above = self.floors[l]
            candidates,which = [],None
            for room,d in above['partition'].enum():
                if vec3(0,0,0).com(room.loop).isnearxy(vec3(0,0,0).com(below.loop)):
                    candidates.append(room)
            if candidates:
                if len(candidates) == 1:
                    which = candidates[0]
                else:
                    area = pym.bareaxy(below.loop)
                    areas = [pym.bareaxy(r.loop) for r in candidates]
                    adifs = [abs(a-area) for a in areas]
                    which = candidates[adifs.index(min(adifs))]
            else:
                raise ValueError('no candidates for ashaft found')
            if which:
                ovlp = pym.ebuxy(which.loop,below.loop)
                rxs,rms = zip(*self.topology.vs)
                rms = [r['pv'] for r in rms]
                bmx = rms.index(below)
                amx = rms.index(which)
                if not (bmx, amx) in self.topology.elook:
                    self.topology.ae(bmx, amx, overlap=ovlp)
                    below = which
                    print('shaft connected', amx, bmx)


    def afloor(self,pg,floorheight):
        '''Add a floor to the building using a partition'''
        self.floors.append({'partition':pg})
        print('update wiregraph with floor')
        rooms = list(pg.enum(filter = lambda v,d : not v.style == 'skip'))
        room,depth = rooms[0]
        rx = self.topology.av(
            pv=room, style=room.style, skirt=1, crown=1, 
            floor=len(self.floors)-1, 
            wallheight=floorheight-1-1, wallwidth=0.5, 
            oholes=[[] for x in range(len(room.loop))], 
            iholes=[[[] for x in range(len(h))] for h in room.holes], 
            otypes=['e' for x in range(len(room.loop))], 
            itypes=[['e' for x in range(len(h))] for h in room.holes], )
        for room,depth in rooms[1:]:
            rx = self.topology.mev(rx, 
                {'pv':room, 'style':room.style, 'skirt':1, 'crown':1, 
                 'floor':(len(self.floors)-1), 
                 'wallheight':floorheight-1-1, 'wallwidth':0.5, 
                 'oholes':[[] for x in range(len(room.loop))], 
                 'iholes':[[[] for x in range(len(h))] for h in room.holes], 
                 'otypes':['e' for x in range(len(room.loop))], 
                 'itypes':[['e' for x in range(len(h))] for h in room.holes], }, 
                {})
            print('floor added')


    def __init__(self,b):
        self.footprint = b
        self.floors = []
        self.topology = wiregraph()


    @staticmethod
    def awallhole(whs,h):
        for wx in range(len(whs)):
            if h[1] > whs[wx][1]:
                whs.insert(wx,h)
                break
            elif h[1] == whs[wx][1]:
                print('cannot add another hole here...',h[1])
                pdb.set_trace()
                break
        else:
            whs.append(h)


    def wallholes(self,vx):
        vx,v = self.topology.vs[vx]

        rb,oholes,iholes = v['pv'].loop,v['oholes'],v['iholes']
        ww,wh = v['wallwidth'],v['wallheight']

        dw,dh,dz = 1.5,3,0.0
        winw,winh,winz,winm,wins = 2.5,2,1,1,4

        for adj in self.topology.rings[vx]:
            ox,ov = self.topology.vs[adj]
            if ov is None:
                continue
            else:
                # consider adjacency of v.loop and ov.loop
                ob = ov['pv'].loop
                aws = pym.badjbxy(rb,ob,minovlp = dw+2)
                for aw in aws:
                    rwx,owx = aw
                    ips = pym.sintsxyp(rb[rwx-1],rb[rwx],ob[owx-1],ob[owx])
                    imp = ips[0].mid(ips[1])
                    ridp = 1-rb[rwx-1].d(imp)/rb[rwx-1].d(rb[rwx])
                    v['otypes'][rwx-1] = 'i'
                    whs = oholes[rwx-1]
                    self.awallhole(whs,(adj,ridp,dw,dh,dz))
                    #if v[2]['rtype'] == 'open' and ov[2]['rtype'] == 'open':
                    #    tdh = dh*100
                    #    tww = max(v[2]['wwidths'][rwx-1],ov[2]['wwidths'][rwx-1])
                    #    tdw = ips[0].d(ips[1])-2*tww
                    #    # hack door that consumes the wall
                    #    self.wallhole(rx,rwx-1,(adj,ridp,tdw,tdh,dz)) 
                    #else:
                    #    # typical door that does not consume the wall
                    #    self.wallhole(rx,rwx-1,(adj,ridp,dw,dh,dz))

                # consider adjacency of v.loop and each of ov.holes
                for oh in ov['pv'].holes:
                    aws = pym.badjbxy(rb,oh,minovlp = dw+2)
                    for aw in aws:
                        rwx,owx = aw
                        ips = pym.sintsxyp(rb[rwx-1],rb[rwx],oh[owx-1],oh[owx])
                        imp = ips[0].mid(ips[1])
                        ridp = 1-rb[rwx-1].d(imp)/rb[rwx-1].d(rb[rwx])
                        v['otypes'][rwx-1] = 'i'
                        whs = oholes[rwx-1]
                        self.awallhole(whs,(adj,ridp,dw,dh,dz))

                # consider adjacency of each of v.holes and ov.loop
                for hx,rh in enumerate(v['pv'].holes):
                    aws = pym.badjbxy(rh,ob,minovlp = dw+2)
                    for aw in aws:
                        rwx,owx = aw
                        ips = pym.sintsxyp(rh[rwx-1],rh[rwx],ob[owx-1],ob[owx])
                        imp = ips[0].mid(ips[1])
                        ridp = 1-rh[rwx-1].d(imp)/rh[rwx-1].d(rh[rwx])
                        v['itypes'][hx][rwx-1] = 'i'
                        whs = iholes[hx][rwx-1]
                        self.awallhole(whs,(adj,ridp,dw,dh,dz))
        
        # consider windows and external walls
        for wx in range(len(rb)):
            wp1,wp2 = rb[wx-1],rb[wx]
            wt = v['otypes'][wx-1]
            whs = oholes[wx-1]
            if wt == 'e':
                '''
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
                '''
                wpd = wp1.d(wp2)
                #if not isexit and wpd > dw+2 and not rv[2]['rtype'] == 'closed':
                if wpd > dw+2:
                    wcnt = int((wpd-1)/(winw+wins*winm))
                    if wcnt > 0:
                        wpoff = 1.0/wcnt
                        for wx in range(wcnt):
                            winp = (wx+1.0)/wcnt-wpoff/2
                            self.awallhole(whs,('window',winp,winw,winh,winz))


    @staticmethod
    def awall(rmmodel,tm,w1,w2,sh,wh,ww,hs,b1 = None,b2 = None):
        wpys,portals = awall(w1,w2,wh,hs,b1,b2)
        [rmmodel.asurf(wpy,tm) for wpy in wpys]
        for hps in portals:
            wpys = aportal(hps,w1,w2,ww)
            rv = True if ww < 0 else False
            [rmmodel.asurf(wpy,tm,rv = rv) for wpy in wpys]


    @staticmethod
    def room(pv,ohs,ihs,sh,wh,ch,ww,style):
        rmmodel = model()
        tm = rmmodel.agfxmesh()

        fb  = pym.contract([p.cp().ztrn(sh) for p in pv.loop],ww)
        fbh = [pym.contract([p.cp().ztrn(sh) for p in h],-ww) for h in pv.holes]
        cb  = pym.contract([p.cp().ztrn(sh+wh) for p in pv.loop],ww)
        cbh = [pym.contract([p.cp().ztrn(sh+wh) for p in h],-ww) for h in pv.holes]
        rmmodel.asurf((fb, fbh), tm)
        rmmodel.asurf((cb, cbh), tm, rv=True)

        for x in range(len(fb)):
            w2,w1 = fb[x-1],fb[x]
            b2,b1 = pv.loop[x-1],pv.loop[x]
            hs = ohs[x-1]
            blggraph.awall(rmmodel,tm,
                w1,w2,sh,wh,ww,hs,b1.cp().ztrn(sh),b2.cp().ztrn(sh))

        for y in range(len(fbh)):
            for x in range(len(fbh[y])):
                w2,w1 = fbh[y][x-1],fbh[y][x]
                b2,b1 = pv.holes[y][x-1],pv.holes[y][x]
                hs = ihs[y][x-1]
                blggraph.awall(rmmodel,tm,
                    w1,w2,sh,wh,-ww,hs,b1.cp().ztrn(sh),b2.cp().ztrn(sh))

        rmmodel.normals(tm)
        return rmmodel


    lipheight = 1
    @staticmethod
    def roof(pv,ohs,ihs,sh,wh,ch,ww,style):
        rfmodel = model()
        tm = rfmodel.agfxmesh()

        fb  = pym.contract([p.cp().ztrn(sh) for p in pv.loop],ww)
        fbh = [pym.contract([p.cp().ztrn(sh) for p in h],-ww) for h in pv.holes]
        rfmodel.asurf((fb, fbh), tm)

        for x in range(len(fb)):
            w1,w2 = fb[x-1],fb[x]
            b2,b1 = pv.loop[x-1],pv.loop[x]
            blggraph.awall(rfmodel,tm,
                w1,w2,sh,blggraph.lipheight,ww,[],b1.cp().ztrn(sh),b2.cp().ztrn(sh))

        for y in range(len(fbh)):
            for x in range(len(fb)):
                w2,w1 = fbh[y][x-1],fbh[y][x]
                b2,b1 = pv.holes[y][x-1],pv.holes[y][x]
                hs = ihs[y][x-1]
                blggraph.awall(rfmodel,tm,
                    w1,w2,sh,wh,-ww,hs,b1.cp().ztrn(sh),b2.cp().ztrn(sh))

        rfmodel.normals(tm)
        return rfmodel


    @staticmethod
    def shell(bfp,vs,floorheights,eww = 0.5):

        def handleturn(turn,p1,p2,p3,p4,r = 2):
            if turn == 'parallel':
                cp = p2.cp()
            else:
                ltan1,ltan2 = p1.tov(p2),p3.tov(p4)
                l1far = p1.cp().trn(ltan1.cp().uscl(-r))
                l2far = p2.cp().trn(ltan1.cp().uscl( r))
                l3far = p3.cp().trn(ltan2.cp().uscl(-r))
                l4far = p4.cp().trn(ltan2.cp().uscl( r))
                cp = pym.sintsxyp(l1far,l2far,l3far,l4far)
            return cp

        shmodel = model()
        tm = shmodel.agfxmesh()
        for vx,v in vs:
            sty,vb,sh,ch = v['style'],v['pv'].loop,v['skirt'],v['crown'],
            wh,ohs,ihs = v['wallheight'],v['oholes'],v['iholes']

            # need to loop around holes too!!
            lastturn = None
            for vbx in range(-1,len(vb)):

                wt = v['otypes'][vbx-2]
                if wt == 'i':
                    continue

                vb0,vb1,vb2,vb3 = vb[vbx-3],vb[vbx-2],vb[vbx-1],vb[vbx]
                wnm0 = vb0.tov(vb1).crs(vec3(0,0,1)).nrm().uscl(eww)
                wnm1 = vb1.tov(vb2).crs(vec3(0,0,1)).nrm().uscl(eww)
                wnm2 = vb2.tov(vb3).crs(vec3(0,0,1)).nrm().uscl(eww)
                wp0 = vb0.cp().trn(wnm0).ztrn(sh)
                wp1 = vb1.cp().trn(wnm0).ztrn(sh)
                wp2 = vb1.cp().trn(wnm1).ztrn(sh)
                wp3 = vb2.cp().trn(wnm1).ztrn(sh)
                wp4 = vb2.cp().trn(wnm2).ztrn(sh)
                wp5 = vb3.cp().trn(wnm2).ztrn(sh)

                wtn1,wtn2 = vb1.tov(vb2).nrm(),vb2.tov(vb3).nrm()
                tcrstz = gtl.near(wtn1.crs(wtn2).z,0,gtl.epsilon)
                if tcrstz > 0:nextturn = 'convex'
                elif tcrstz < 0:nextturn = 'concave'
                else:nextturn = 'parallel'
                #print('this turn %s  |  last turn %s  | building shell' % (nextturn,lastturn))

                if lastturn:
                    wpp1 = handleturn(lastturn,wp0,wp1,wp2,wp3)
                    wpp2 = handleturn(nextturn,wp2,wp3,wp4,wp5)
                    lastturn = nextturn
                else:
                    lastturn = nextturn
                    continue

                sp1,sp2 = wpp1.cp().ztrn(-sh),wpp2.cp().ztrn(-sh)
                cp1,cp2 = wpp1.cp().ztrn( wh),wpp2.cp().ztrn( wh)
                [shmodel.asurf(wpy,tm) for wpy in awall(sp1,sp2,sh,[])[0]]
                if 'roof' in sty:
                    blggraph.awall(shmodel,tm,
                        wpp1,wpp2,sh,blggraph.lipheight,eww,[])
                    walltop = [wpp1.cp(), wpp2.cp(), 
                        vb2.cp().ztrn(sh).lerp(wpp2, -1), 
                        vb1.cp().ztrn(sh).lerp(wpp1, -1)]
                    vec3(0,0,blggraph.lipheight).trnps(walltop)
                    shmodel.asurf((walltop,()),tm)
                else:
                    blggraph.awall(shmodel,tm,
                        wpp1,wpp2,sh,wh,eww,ohs[vbx-2],wp2,wp3)
                    [shmodel.asurf(wpy,tm) for wpy in awall(cp1,cp2,ch,[])[0]]
        shmodel.normals(tm)
        return shmodel


def amodel(sg,m):
    p, q, s = vec3(0, 0, 0), quat(1, 0, 0, 0), vec3(1, 1, 1)
    sgv = sg.avert(p, q, s, models=[m], parent=sg.root)


def building(fp,e):
    b,hs = fp

    rb = [p.cp() for p in b]
    #inner = pym.contract(rb, 8)
    pg1 = partition(loop=rb, holes=[], style='outer')
    #nsvs,nv = pg1.split(pg1.root, inner, style='inner')
                          
    rb = [p.cp().ztrn(8) for p in rb]
    #inner = pym.contract(rb, 8)
    pg2 = partition(loop=rb, holes=[], style='outer')
    #nsvs,nv = pg2.split(pg2.root, inner, style='inner')

    #rb = [p.cp().ztrn(5) for p in nv.loop]
    #rb = [p.cp().ztrn(5) for p in rb]
    #pg3 = partition(loop=rb, holes=[], style='roof')

    #floors = ((pg1,8),(pg2,6),(pg3,5))
    #shafts = ((0,1,pg1.root.ix),)
    floors = ((pg1,8),(pg2,6))
    shafts = ()



    bg = blggraph(b)
    [bg.afloor(*f) for f in floors]
    [bg.ashaft(*s) for s in shafts]
    #ax = bg.plot()
    #plt.show()

    sg = scenegraph()
    for vx,v in bg.topology.vs:
        bg.wallholes(vx)

        pv,sh,ch = v['pv'],v['skirt'],v['crown']
        wh,ww = v['wallheight'],v['wallwidth']
        ohs,ihs,sty = v['oholes'],v['iholes'],v['style']

        f = bg.roof if 'roof' in sty else bg.room
        amodel(sg,f(pv,ohs,ihs,sh,wh,ch,ww,sty))

    fhs = list(zip(*floors))[1]
    amodel(sg,bg.shell(bg.footprint,bg.topology.vs,fhs))

    # need exits, roofs, atriums, shafts, basements

    return sg















dkw = lambda kws,k,d : kws[k] if k in kws else d

class Abuilding(scenegraph):

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

                if isinstance(b1,vec3) and b1.isnear(b2):
                    pdb.set_trace()

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

                if isinstance(b1,vec3) and b1.isnear(b2):
                    pdb.set_trace()

                # create wall surfaces
                wpys,portals = awall(w1,w2,wh,hs,b1.cp(),b2.cp())
                for wpy in wpys:m.asurf(wpy,tm)

                # create portal surfaces
                for hps in portals:
                    wpys = aportal(hps,w1,w2,ww)
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
                    wpys,portals = awall(crnp1,crnp2,fh,[])
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
                        wpys,portals = awall(wp1,wp2,wh,hs)
                        for wpy in wpys:m.asurf(wpy,tm)
                        for hps in portals:
                            wpys = aportal(hps,wp1,wp2,eww)
                            for wpy in wpys:m.asurf(wpy,tm)
                        sp1,sp2 = wp1.cp().ztrn(-sh),wp2.cp().ztrn(-sh)
                        cp1,cp2 = wp1.cp().ztrn( wh),wp2.cp().ztrn( wh)
                        wpys,portals = awall(sp1,sp2,sh,[])
                        for wpy in wpys:m.asurf(wpy,tm)
                        wpys,portals = awall(cp1,cp2,ch,[])
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

class AAblggraph:

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
        if ax is None:ax = plot_axes_xy(50)
        for rmv in self.vs:
            if rmv is None:continue
            if rmv[2]['level'] > 0:continue
            bcol = 'b' if rmv[2]['shaft'] else None
            rconbs = 0.5 if pym.bccw(rmv[2]['bound']) else -0.5
            rconb = pym.contract(rmv[2]['bound'],rconbs)
            ax = plot_polygon_xy(rconb,ax,lw = 4,col = bcol)
            rc = vec3(0,0,0).com(rmv[2]['bound'])
            rs = str(rmv[0])+','+str(rmv[1])
            ax = plot_point_xy(rc,plot_point_xy_annotate(rc,ax,rs))
            if rmv[1]:
                for re in rmv[1]:
                    if self.vs[re] is None:continue
                    rt = (rc,vec3(0,0,0).com(self.vs[re][2]['bound']))
                    ax = plot_edges_xy(rt,ax,col = 'g')
            for exit in rmv[2]['exits']:
                if exit is True:
                    exitp = rc
                    ax = plot_point_xy(exitp,
                        plot_point_xy_annotate(exitp,ax,'exit'))
                else:
                    ax = plot_point_xy(exit,
                        plot_point_xy_annotate(exit,ax,'exit'))
        if not fp is None:
            ax = plot_polygon_xy(fp,ax,lw = 2,col = 'r')
        return ax

    def plot(self,fp = None,ax = None):
        if ax is None:ax = plot_axes(50)
        for rmv in self.vs:
            if rmv is None:continue
            ax = plot_polygon(rmv[2]['bound'],ax,lw = 4)
            rc = vec3(0,0,0).com(rmv[2]['bound'])
            rs = str(rmv[0])+','+str(rmv[1])
            ax = plot_point(rc,plot_point_xy_annotate(rc,ax,rs))
            if rmv[1]:
                for re in rmv[1]:
                    if self.vs[re] is None:continue
                    rt = (rc,vec3(0,0,0).com(self.vs[re][2]['bound']))
                    ax = plot_edges(rt,ax,col = 'g')
            for exit in rmv[2]['exits']:
                if exit is True:
                    exitp = rc
                    ax = plot_point(exitp,
                        plot_point_xy_annotate(exitp,ax,'exit'))
                else:
                    ax = plot_point(exit,
                        plot_point_xy_annotate(exit,ax,'exit'))
        if not fp is None:
            ax = plot_polygon(fp,ax,lw = 2,col = 'r')
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


def AAAnew(*ags,**kws):
    blgg = blggraph(**kws)
    blgg.graph()
    #blgg.plotxy()
    #plt.show()
    kws['bgraph'] = blgg
    n = building(*ags,**kws)
    return n
