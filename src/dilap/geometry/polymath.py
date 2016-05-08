from dilap.geometry.vec3 import vec3
import dilap.geometry.tools as gtl

import dilap.core.plotting as dtl
import matplotlib.pyplot as plt

import numpy

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
def sintsxy(s11,s12,s21,s22,ie = True,ieb = True,col = True):
    # if segments are colinear and overlapping
    # elif segments are colinear and nonoverlapping
    # elif segments are skew and possibly overlapping
    s1tn,s2tn = s11.tov(s12),s21.tov(s22)
    ds = s11.tov(s21)
    s1crss2 = s1tn.crs(s2tn)
    s1s2area = gtl.near(s1crss2.mag(),0)
    dscrss1 = ds.crs(s1tn)
    dss1area = gtl.near(dscrss1.mag(),0)
    if s1s2area == 0 and dss1area == 0:
        # if uncontained overlap
        # elif contained overlap
        # elif perfect overlap
        # elif single point overlap
        # else disjoint
        if not col:return 0
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
        elif (t0 == 1 and t1 > t0) or (t0 < t1 and t1 == 0):return 1 if ie else 0
        else:return 0
    elif s1s2area == 0 and not dss1area == 0:return 0
    elif not s1s2area == 0:
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

# does one segment intersect another
# ie  : do endpoint intersections count (end of one segment only)
# ieb : do endpoint intersections count (end of both segments)
# col : do colinear intersections counts
def sintsxyp(s11,s12,s21,s22,ie = True,ieb = True,col = True):
    # if segments are colinear and overlapping
    # elif segments are colinear and nonoverlapping
    # elif segments are skew and possibly overlapping
    s1tn,s2tn = s11.tov(s12),s21.tov(s22)
    ds = s11.tov(s21)
    s1crss2 = s1tn.crs(s2tn)
    s1s2area = gtl.near(s1crss2.mag(),0)
    dscrss1 = ds.crs(s1tn)
    dss1area = gtl.near(dscrss1.mag(),0)
    if s1s2area == 0 and dss1area == 0:
        # if uncontained overlap
        # elif contained overlap
        # elif perfect overlap
        # elif single point overlap
        # else disjoint
        if not col:return None
        s1l,s2l = s1tn.mag(),s2tn.mag()
        s1dots2 = s1tn.dot(s2tn)
        apara = gtl.near(s1dots2,0) < 0
        t0 = ds.dot(s1tn)/s1l**2
        t1 = t0+s1dots2/s1l**2
        t0 = gtl.near(gtl.near(t0,0),1)
        t1 = gtl.near(gtl.near(t1,0),1)
        if apara:t0,t1 = t1,t0
        if gtl.inrng(t0,0,1) or gtl.inrng(t1,0,1):
            sp = [s11,s12,s21,s22]
            ss = [s11.dot(s1tn),s12.dot(s1tn),s21.dot(s1tn),s22.dot(s1tn)]
            sx = [0,1,2,3]
            mm = ss.index(min(ss)),ss.index(max(ss))
            sx = [j for j in sx if not j in mm]
            ip1,ip2 = sp[sx[0]],sp[sx[1]]
            return ip1,ip2
        elif gtl.inrng(0,t0,t1) and gtl.inrng(1,t0,t1):
            sp = [s11,s12,s21,s22]
            ss = [s11.dot(s1tn),s12.dot(s1tn),s21.dot(s1tn),s22.dot(s1tn)]
            sx = [0,1,2,3]
            mm = ss.index(min(ss)),ss.index(max(ss))
            sx = [j for j in sx if not j in mm]
            ip1,ip2 = sp[sx[0]],sp[sx[1]]
            return ip1,ip2
        elif (t0 == 0 and t1 == 1) or (t0 == 1 and t1 == 0):
            if not ieb:return None
            else:return (s11.cp(),s12.cp())
        elif (t0 == 1 and t1 > t0) or (t0 < t1 and t1 == 0):
            if not ie:return None
            else:
                if s11.isnear(s21):return s11.cp()
                elif s11.isnear(s22):return s11.cp()
                else:return s12.cp()
        else:return None
    elif s1s2area == 0 and not dss1area == 0:return None
    elif not s1s2area == 0:
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
def sintbxy(s1,s2,b,ie = True):
    for x in range(len(b)):
        b1,b2 = b[x-1],b[x]
        if sintsxy(s1,s2,b1,b2,ie = ie,ieb = ie,col = ie):
            return 1
    return 0

# where does a segment intersect a boundary polygon
# ie : if the intersection is strictly along the boundary, does it count
def sintbxyp(s1,s2,b,ie = True):
    for x in range(len(b)):
        b1,b2 = b[x-1],b[x]
        ip = sintsxyp(s1,s2,b1,b2,ie = ie,ieb = ie,col = ie):
        if not ip is None:
            return ip
    return 0

###############################################################################
###############################################################################
###############################################################################

###############################################################################
### functions of boundary polygons
###############################################################################

# is a boundary polygon contained within another polygon
# ie : do edge intersections mean no containment
def binbxy(b1,b2,ie = True):
    for p in b1:
        if ie and p.onbxy(b2):return 0
        elif not p.inbxy(b2):return 0
    return 1

# is a boundary polygon intersecting another boundary polygon
# ie : do edge intersections count
def bintbxy(b1,b2,ie = True):
    for b1x in range(len(b1)):
        b1p1,b1p2 = b1[b1x-1],b1[b1x]
        for b2x in range(len(b2)):
            b2p1,b2p2 = b2[b2x-1],b2[b2x]
            if sintsxy(b1p1,b1p2,b2p1,b2p2,ie = ie,ieb = ie,col = ie):
                return 1
    return 0

# segment a boundary polygon based on intersections with another
def bsegbxy(b1,b2):
    b1es = []
    unfn = [(b1[x-1],b1[x]) for x in range(len(b1))]

    while unfn:
        u1,u2 = unfn.pop(0)
        ip = sintbxyp(u1,u2,b2)
        if ip is None:
            pdb.set_trace()
        else:
            pdb.set_trace()

    #b1es,b2es = [],[]
    #b2es = [(b2[x-1],b2[x]) for x in range(len(b2))]
    for b1x in range(len(b1)):
        b1p1,b1p2 = b1[b1x-1],b1[b1x]
        for b2x in range(len(b2)):
            b2p1,b2p2 = b2[b2x-1],b2[b2x]

            ip = sintsxyp(b1p1,b1p2,b2p1,b2p2,ie = True,ieb = False)
            if ip is None:
                b1es.append((b1p1,b1p2))

# compute difference of boundary polygons
def ebdxy(b1,b2):
    b1es,b2es = [],[]
    #b2es = [(b2[x-1],b2[x]) for x in range(len(b2))]
    for b1x in range(len(b1)):
        b1p1,b1p2 = b1[b1x-1],b1[b1x]
        for b2x in range(len(b2)):
            b2p1,b2p2 = b2[b2x-1],b2[b2x]

            ip = sintsxyp(b1p1,b1p2,b2p1,b2p2,ie = True,ieb = False)
            if ip is None:
                b1es.append((b1p1,b1p2))

            if not ip is None:

                ax = dtl.plot_axes_xy(10)
                ax = dtl.plot_edges_xy([b1p1,b1p2],ax,lw = 2.0,col = 'b')
                ax = dtl.plot_edges_xy([b2p1,b2p2],ax,lw = 2.0,col = 'g')
                plt.show()

                pdb.set_trace()
                
    #ax = dtl.plot_axes_xy(10)
    #for e in b1es:ax = dtl.plot_edges_xy(e,ax,lw = 2.0,col = 'b')
    #for e in b2es:ax = dtl.plot_edges_xy(e,ax,lw = 2.0,col = 'g')
    #plt.show()

    return b1

###############################################################################
###############################################################################
###############################################################################

# is the interior of a boundary polygon disjoint from that of another
def pdisjpxy(b1,b2):
    print('do the interiors of any edge segments intersect?')
    for b1x in range(len(b1)):
        b1p1,b1p2 = b1[b1x-1],b1[b1x]
        for b2x in range(len(b2)):
            b2p1,b2p2 = b2[b2x-1],b2[b2x]
            if gtl.segs_isect_perp(b1p1,b1p2,b2p1,b2p2):
                return False
    return True





