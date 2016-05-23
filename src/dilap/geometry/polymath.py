from dilap.geometry.vec3 import vec3
import dilap.geometry.tools as gtl

import dilap.core.plotting as dtl
import matplotlib.pyplot as plt

import numpy,math

import pdb



# compute union of two polygons
def polyu(p1,p2):
    raise NotImplemented
    return p1

# compute intersect of two polygons
def polyi(p1,p2):
    raise NotImplemented
    return p1

# compute difference of two polygons
def polyd(p1,p2):
    # if ext p1 contains ext p2, consider holes of p1 and ext of p2
    eb1,eb2 = list(p1[0]),list(p2[0])

    #if pdisjp(eb1,eb2):
    #    print('p2 is nonoverlapping p1')
    #    return p1

    if binbxy(eb1,eb2,ie = False):
        print('eb1 contained by eb2')
        raise ValueError

    elif binbxy(eb2,eb1,ie = False):
        print('eb2 contained by eb1')
        print('if ib1s contact eb2, union, else just add to ib1')
        return (tuple(eb1),(tuple(eb2),))

    elif bintbxy(eb1,eb2,ie = True):
        print('eb1 and eb2 intersect!')
        neb = ebdxy(eb1,eb2)
        return (tuple(neb),p1[1])

    print('off record')
    return p1


    np1,np2 = [],[]
    for x1 in range(len(eb1)):
        ep11,ep12 = eb1[x1-1].cpxy(),eb1[x1].cpxy()

        if ep11.inbxy(eb2):
            print('ep11 inside boundary')
        elif ep11.onbxy(eb2):
            print('ep11 on boundary')
            np1.append(ep11)
        else:
            print('ep11 outside boundary')
            np1.append(ep11)
            if ep12.inbxy(eb2):
                print('find intersection and add that; look for p2 points?')

        #for x2 in range(len(eb2)):
        #    ep21,ep22 = eb2[x2-1].cpxy(),eb2[x2].cpxy()

    for x2 in range(len(eb2)):
        ep21,ep22 = eb2[x2-1].cpxy(),eb2[x2].cpxy()

        if ep21.inbxy(eb1):
            print('ep21 inside boundary')
            np2.append(ep21)
        elif ep21.onbxy(eb1):
            print('ep21 on boundary')
        else:
            print('ep21 outside boundary')

    if not np1:
        print('polyd:p1 contained by p2')
        raise ValueError
        
    #ax = dtl.plot_axes()
    #ax = dtl.plot_polygon(eb2,ax)
    #ax = dtl.plot_point(ep11,ax)
    #plt.show()

    print('compute polygon difference!')
    nibs = []

    if len(np2) == len(eb2):nibs.append(np2)

    return tuple(np1),tuple(tuple(n) for n in nibs)


###############################################################################
### functions of line segments with line segments or boundary polygons
###############################################################################

# does one segment intersect another
# ie  : do endpoint intersections count (end of one segment only)
# ieb : do endpoint intersections count (end of both segments)
# col : do colinear intersections counts
def sintsxy(s11,s12,s21,s22,ie = True,ieb = True,col = True,skew = True):
    # if segments are colinear and overlapping
    # elif segments are colinear and nonoverlapping
    # elif segments are skew and possibly overlapping
    #s11,s12,s21,s22 = s11.cpxy(),s12.cpxy(),s21.cpxy(),s22.cpxy()
    s1tn,s2tn = s11.tov(s12),s21.tov(s22)
    ds = s11.tov(s21)
    s1crss2 = s1tn.crs(s2tn)
    s1s2area = gtl.near(s1crss2.mag(),0)
    dscrss1 = ds.crs(s1tn)
    dss1area = gtl.near(dscrss1.mag(),0)
    if s1s2area == 0 and dss1area == 0:
        if not col:return 0
        # if uncontained overlap
        # elif contained overlap
        # elif perfect overlap
        # elif single point overlap
        # else disjoint
        s1l = s1tn.mag()
        s1dots2 = s1tn.dot(s2tn)
        apara = gtl.near(s1dots2,0) < 0
        t0 = ds.dot(s1tn)/s1l**2
        t1 = t0+s1dots2/s1l**2
        t0 = gtl.near(gtl.near(t0,0),1)
        t1 = gtl.near(gtl.near(t1,0),1)
        if apara:t0,t1 = t1,t0
        if gtl.inrng(t0,0,1) or gtl.inrng(t1,0,1):return 1
        elif gtl.inrng(0,t0,t1) and gtl.inrng(1,t0,t1):return 1
        elif (t0 == 0 and t1 == 1) or (t0 == 1 and t1 == 0):return 1 if ieb else 0
        elif (t0 == 0 and t1 > t0) or (t1 == 0 and t0 > t1):return 1
        elif (t0 < t1 and t1 == 1) or (t1 < t0 and t0 == 1):return 1
        elif (t0 == 1 and t1 > t0) or (t0 < t1 and t1 == 0):return 1 if ie else 0
        else:return 0
    elif s1s2area == 0 and not dss1area == 0:return 0
    elif not s1s2area == 0:
        if not skew:return 0
        # if intersection at endpoints of both segments but not counting
        # elif intersection at endpoint of only one segment but not counting
        # elif no proper intersection at neither segments endpoints
        # else must have intersected by one of three possible cases
        dscrss2 = ds.crs(s2tn)
        u = gtl.near(gtl.near(dscrss1.z/s1crss2.z,0),1)
        t = gtl.near(gtl.near(dscrss2.z/s1crss2.z,0),1)
        if not ieb and ((u == 0 or u == 1) and (t == 0 or t == 1)):return 0
        elif not ie and ((u == 0 or t == 0) or (u == 1 or t == 1)):return 0
        elif not (0 <= u and u <= 1) or not (0 <= t and t <= 1):return 0
        else:return 1

# given 4 points and a tangent vector, 
# extract the two whose projections are in between the extreme projections
def prjmed(p1,p2,p3,p4,tn):
    sp = [p1,p2,p3,p4]
    ss = [p.dot(tn) for p in sp]
    sx = list(range(len(ss)))
    mm = ss.index(min(ss)),ss.index(max(ss))
    ix1,ix2 = [j for j in sx if not j in mm]
    ip1,ip2 = sp[ix1],sp[ix2]
    return ip1,ip2

# does one segment intersect another
# ie  : do endpoint intersections count (end of one segment only)
# ieb : do endpoint intersections count (end of both segments)
# col : do colinear intersections counts
def sintsxyp(s11,s12,s21,s22,ie = True,ieb = True,col = True,skew = True):
    # if segments are colinear and overlapping
    # elif segments are colinear and nonoverlapping
    # elif segments are skew and possibly overlapping
    #s11,s12,s21,s22 = s11.cpxy(),s12.cpxy(),s21.cpxy(),s22.cpxy()
    s1tn,s2tn = s11.tov(s12),s21.tov(s22)
    ds = s11.tov(s21)
    s1crss2 = s1tn.crs(s2tn)
    s1s2area = gtl.near(s1crss2.mag(),0)
    dscrss1 = ds.crs(s1tn)
    dss1area = gtl.near(dscrss1.mag(),0)
    if s1s2area == 0 and dss1area == 0:
        if not col:return None
        # if uncontained overlap
        # elif contained overlap without either endpoint
        # elif perfect overlap
        # elif single endpoint overlap and contained
        # elif single point overlap
        # else disjoint
        s1l,s2l = s1tn.mag(),s2tn.mag()
        s1dots2 = s1tn.dot(s2tn)
        apara = gtl.near(s1dots2,0) < 0
        t0 = ds.dot(s1tn)/s1l**2
        t1 = t0+s1dots2/s1l**2
        t0 = gtl.near(gtl.near(t0,0),1)
        t1 = gtl.near(gtl.near(t1,0),1)
        if apara:t0,t1 = t1,t0
        if gtl.inrng(t0,0,1) or gtl.inrng(t1,0,1):
            return prjmed(s11,s12,s21,s22,s1tn)
        elif gtl.inrng(0,t0,t1) and gtl.inrng(1,t0,t1):
            return prjmed(s11,s12,s21,s22,s1tn)
        elif (t0 == 0 and t1 == 1) or (t0 == 1 and t1 == 0):
            return None if not ieb else (s11.cp(),s12.cp())
        elif (t0 == 0 and t1 > t0) or (t1 == 0 and t0 > t1):
            return prjmed(s11,s12,s21,s22,s1tn)
        elif (t0 < t1 and t1 == 1) or (t1 < t0 and t0 == 1):
            return prjmed(s11,s12,s21,s22,s1tn)
        elif (t0 == 1 and t1 > t0) or (t0 < t1 and t1 == 0):
            if not ie:return None
            else:
                if s11.isnear(s21):return s11.cp()
                elif s11.isnear(s22):return s11.cp()
                else:return s12.cp()
        else:return None
    elif s1s2area == 0 and not dss1area == 0:return None
    elif not s1s2area == 0:
        if not skew:return 0
        # if intersection at endpoints of both segments but not counting
        # elif intersection at endpoint of only one segment but not counting
        # elif no proper intersection at neither segments endpoints
        # else must have intersected by one of three possible cases
        dscrss2 = ds.crs(s2tn)
        u = gtl.near(gtl.near(dscrss1.z/s1crss2.z,0),1)
        t = gtl.near(gtl.near(dscrss2.z/s1crss2.z,0),1)
        if not ieb and ((u == 0 or u == 1) and (t == 0 or t == 1)):return None
        elif not ie and ((u == 0 or t == 0) or (u == 1 or t == 1)):return None
        elif not (0 <= u and u <= 1) or not (0 <= t and t <= 1):return None
        else:return s11.cp().trn(s1tn.cp().uscl(t))

# does a segment intersect a boundary polygon
# ie : if the intersection is strictly along the boundary, does it count
def sintbxy(s1,s2,b,ie = True,ieb = True,col = True):
    for x in range(len(b)):
        b1,b2 = b[x-1],b[x]
        if sintsxy(s1,s2,b1,b2,ie = ie,ieb = ieb,col = col):
            return 1
    return 0

# where does a segment intersect a boundary polygon
# ie : if the intersection is strictly along the boundary, does it count
def sintbxyp(s1,s2,b,ie = True,ieb = True,col = True):
    def isnew(np):
        for p in ips:
            if p.isnear(np):return
        ips.append(np)
    ips = []
    for x in range(len(b)):
        b1,b2 = b[x-1],b[x]
        ip = sintsxyp(s1,s2,b1,b2,ie = ie,ieb = ieb,col = col)
        if not ip is None:
            if type(ip) == type(()):
                isnew(ip[0])
                isnew(ip[1])
            else:isnew(ip)
    return ips

###############################################################################
###############################################################################
###############################################################################

###############################################################################
### functions of boundary polygons
###############################################################################

# is a boundary polygon contained within another polygon
# ie : do edge intersections mean no containment
# NOTE: WILL NOT WORK FOR ALL CONCAVE BOUNDARIES AS IS...
def binbxy(b1,b2,ie = True):
    for p in b1:
        if ie and p.onbxy(b2):return 0
        elif not p.inbxy(b2):return 0
    return 1

# is a boundary polygon intersecting another boundary polygon
# ie : do edge intersections count
def bintbxy(b1,b2,ie = True,ieb = True,col = True):
    for b1x in range(len(b1)):
        b1p1,b1p2 = b1[b1x-1],b1[b1x]
        for b2x in range(len(b2)):
            b2p1,b2p2 = b2[b2x-1],b2[b2x]
            if sintsxy(b1p1,b1p2,b2p1,b2p2,ie = ie,ieb = ieb,col = col):
                return 1
    return 0

# given a sequence of segments with potential gaps, 
# determine where the sg could first be attached
#
def placeseg(segs,s1,s2):
    s1xy,s2xy = s1.cpxy(),s2.cpxy()
    for sx in range(len(segs)):
        os1,os2 = segs[sx]
        os1xy,os2xy = os1.cpxy(),os2.cpxy()
        if os1xy.isnear(s2xy):
            segs.insert(sx,(s1,s2))
            return segs
        elif os1xy.isnear(s1xy):
            segs.insert(sx,(s2,s1))
            return segs

    print('failed to placeseg')
    ax = dtl.plot_axes_xy(60)
    for s in segs:
        ax = dtl.plot_edges_xy(s,ax,lw = 2,col = 'b')
    ax = dtl.plot_edges_xy((s1xy,s2xy),ax,lw = 4,col = 'g')
    plt.show()

    pdb.set_trace()

# segment a boundary polygon based on intersections with a line segment
# return two new boundary polygons
def bsegsxy(b,s1,s2):
    #bes = [(b[x-1].cp(),b[x].cp()) for x in range(len(b))]
    bes = bsegbxy(b,(s1,s2))
    left,right = [],[]
    while bes:
        u1,u2 = bes.pop(0)
        u1o = gtl.orient2d(s1,s2,u1)
        u2o = gtl.orient2d(s1,s2,u2)
        if u1o == 0 and u2o == 0:pass
        elif u1o > 0 or u2o > 0:left.insert(0,(u2,u1))
        elif u1o < 0 or u2o < 0:right.append((u1,u2))
        else:raise ValueError

    if not left or not right:
        print('seg missed bpoly?')
        pdb.set_trace()

    ips = sintbxyp(s1,s2,b,col = False)
    left = placeseg(left,ips[0],ips[1])
    right = placeseg(right,ips[1],ips[0])
    leftb = [left[x][1] for x in range(len(left))][::-1]
    rightb = [right[x][1] for x in range(len(right))]

    '''#
    ax = dtl.plot_axes_xy(50)
    ax = dtl.plot_polygon_xy(b,ax,lw = 2,col = 'g')
    ax = dtl.plot_edges_xy((s1,s2),ax,lw = 2)
    for e in left:
        ax = dtl.plot_edges_xy(e,ax,lw = 2,col = 'b')
    for e in right:
        ax = dtl.plot_edges_xy(e,ax,lw = 2,col = 'r')
    ax = dtl.plot_polygon_xy(contract(leftb,0.25),ax,lw = 2,col = 'g')
    ax = dtl.plot_polygon_xy(contract(rightb,0.25),ax,lw = 2,col = 'g')
    plt.show()
    '''#

    return leftb,rightb

# segment a boundary polygon based on intersections with another
def bsegbxy(b1,b2):
    b1es = []
    unfn = [(b1[x-1],b1[x]) for x in range(len(b1))]
    while unfn:
        u1,u2 = unfn.pop(0)
        ips = sintbxyp(u1,u2,b2,ieb = False)
        ips = [p for p in ips if not p.isnear(u1) and not p.isnear(u2)]
        if len(ips) == 0:
            b1es.append((u1.cp(),u2.cp()))
        else:
            u1u2 = u1.tov(u2)
            ipsprj = [ip.dot(u1u2) for ip in ips]
            lst = u1
            while ips:
                ipx = ipsprj.index(min(ipsprj))
                nxtprj = ipsprj.pop(ipx)
                nxt = ips.pop(ipx)
                b1es.append((lst.cp(),nxt.cp()))
                lst = nxt
            b1es.append((lst.cp(),u2.cp()))
    #nb = tuple(be[1] for be in b1es)
    return b1es

# produce some number of closed loops from a set of edges
def loopseg(es):
    unfn = es[:]
    nxt = unfn.pop(0)
    lst = nxt
    loops = [[nxt[0]]]

    def pl():
        ax = dtl.plot_axes_xy(10)
        for e in unfn:ax = dtl.plot_edges_xy(e,ax,lw = 2,col = 'r')
        for l in loops:
            ax = dtl.plot_polygon_xy(l,ax,lw = 2)
        ax = dtl.plot_edges_xy([nxt[1],u1],ax,col = 'g',lw = 5.0)
        ax = dtl.plot_edges_xy([lst[1],u1],ax,col = 'b',lw = 3.0)
        plt.show()

    lost = False
    while unfn:

        if lost:
            print('im lost!!',len(unfn))
            pl()
            #pdb.set_trace()
            return loops

        lost = True
        nxtpotentials = []
        lstpotentials = []
        for ufx in range(len(unfn)):
            u1,u2 = unfn[ufx]
            if u1.isnear(nxt[1]):
                nxtpotentials.append((u1,ufx,False))
                lost = False
                #pl()
            elif u2.isnear(nxt[1]):
                nxtpotentials.append((u2,ufx,True))
                lost = False
                #pl()
            elif u2.isnear(lst[0]):
                lstpotentials.append((u1,ufx,False))
                lost = False
                #pl()
            elif u1.isnear(lst[0]):
                lstpotentials.append((u2,ufx,True))
                lost = False
                #pl()

        if len(nxtpotentials) == 0:
            if len(lstpotentials) == 1:
                p = lstpotentials[0]
                loops[-1].insert(0,p[0])
                if p[2]:lst = unfn.pop(p[1])[::-1]
                else:lst = unfn.pop(p[1])
            elif len(lstpotentials) > 1:
                pdb.set_trace()
        if len(nxtpotentials) == 1:
            p = nxtpotentials[0]
            loops[-1].append(p[0])
            if p[2]:nxt = unfn.pop(p[1])[::-1]
            else:nxt = unfn.pop(p[1])
            #pl()
        elif len(nxtpotentials) > 1:

            ax = dtl.plot_axes_xy(10)
            for p in nxtpotentials:
                print('nxtpotential',p)
                if p[2]:ax = dtl.plot_edges_xy([nxt[1],p[0]],ax,lw = 6,col = 'g')
                else:ax = dtl.plot_edges_xy([nxt[1],p[0]],ax,lw = 6,col = 'g')
            plt.show()
            #pdb.set_trace()
            
            p = nxtpotentials[0]
            loops[-1].append(p[0])
            if p[2]:nxt = unfn.pop(p[1])[::-1]
            else:nxt = unfn.pop(p[1])
    return loops

# compute difference of boundary polygons
def ebdxy(b1,b2):
    b1segs = bsegbxy(b1,b2)
    b2segs = bsegbxy(b2,b1)
    bo = lambda s1,s2,b : s1.inbxy(b) or s2.inbxy(b) or (s1.onbxy(b) and s2.onbxy(b))
    b1only = [p for p in b1segs if not bo(p[0],p[1],b2)]

    #b2inb1 = [p for p in b2segs if bo(p[0],p[1],b1)]
    #b2inb1 = [p for p in b2segs if p[0].inbxy(b1) or p[1].inbxy(b1) or (p[0].onbxy(b2) and p[1].onbxy(b2) and p[0].onbxy(b1) and p[1].onbxy(b1))]
    #b2inb1 = [p for p in b2segs if p[0].inbxy(b1) or p[1].inbxy(b1) or (p[0].onbxy(b1) and p[1].onbxy(b1))]
    b2inb1 = [p for p in b2segs if p[0].inbxy(b1) or p[1].inbxy(b1)]

    #b2only = [p for p in b2segs if not p[0].inbxy(b1) and not p[1].inbxy(b1)]
    #b1inb2 = [p for p in b1segs if p[0].inbxy(b2) or p[1].inbxy(b2)]

    #ax = dtl.plot_axes_xy(10)
    #for bs in b1segs:
    #    ax = dtl.plot_edges_xy(bs,ax,lw = 1.0,col = 'r')
    #for bs in b2segs:
    #    ax = dtl.plot_edges_xy(bs,ax,lw = 3.0,col = 'r')
    #for bs in b1only:
    #    ax = dtl.plot_edges_xy(bs,ax,lw = 5.0,col = 'b')
    #for bs in b2inb1:
    #    ax = dtl.plot_edges_xy(bs,ax,lw = 5.0,col = 'g')
    #plt.show()

    dfs = loopseg(b1only+b2inb1)

    '''#
    ax = dtl.plot_axes_xy(10)
    ax = dtl.plot_polygon_xy(b1,ax,lw = 1.0,col = None)
    ax = dtl.plot_polygon_xy(b2,ax,lw = 1.0,col = None)
    for df in dfs:
        ax = dtl.plot_polygon_xy(df,ax,lw = 8.0,col = None)
    plt.show()
    '''#

    if len(dfs) == 0:return None
    elif len(dfs) == 1:return dfs[0]
    else:return dfs

###############################################################################
###############################################################################
###############################################################################

# return a vector normal to the boundary polygon b
def bnrm(b):
    pcnt = len(b)
    pn = vec3(0,0,0)
    for px in range(pcnt):
        c1,c2,c3  = b[px-2],b[px-1],b[px]
        c1c2 = c1.tov(c2)
        c2c3 = c2.tov(c3)
        cang = c1c2.ang(c2c3)
        if cang > -numpy.pi+0.001 and cang < numpy.pi-0.001:
            pn.trn(c1c2.crs(c2c3).nrm())
    return pn.nrm()

# make a new boundary polygon including the midpoints of every edge
def bisectb(b):
    nb = []
    for x in range(len(b)):
        nb.append(b[x-1])
        nb.append(b[x-1].mid(b[x]))
    return nb

# THIS SHOULD BE ELSEWHERE
def contract(ps,ds):
    fbnd = []
    pns = []
    for x in range(len(ps)):
        p1,p2,p3 = ps[x-2],ps[x-1],ps[x]
        w1t = p1.tov(p2).nrm()
        w2t = p2.tov(p3).nrm()
        w1n = vec3(0,0,1).crs(w1t)
        w2n = vec3(0,0,1).crs(w2t)
        s11 = p1.cp().trn(w1n).trn(w1t.cp().uscl(-100))
        s12 = p2.cp().trn(w1n).trn(w1t.cp().uscl( 100))
        s21 = p2.cp().trn(w2n).trn(w2t.cp().uscl(-100))
        s22 = p3.cp().trn(w2n).trn(w2t.cp().uscl( 100))
        ip = sintsxyp(s11,s12,s21,s22,col = False)
        if ip is None:pn = p2.cp().trn(w1n)
        else:pn = ip.cp()
        pns.append(pn)
    for x in range(len(ps)):
        p1,p2,p3 = ps[x-2],ps[x-1],ps[x]
        fbnd.append(p2.lerp(pns[x],ds))
    fbnd.append(fbnd.pop(0))
    return fbnd



