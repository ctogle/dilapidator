from dilap.geometry.vec3 import vec3
import dilap.geometry.tools as gtl

import dilap.topology.planargraph as pgr

import dilap.core.plotting as dtl
import matplotlib.pyplot as plt

import numpy,math

import pdb



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
        if not skew:return None
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
        # NOTE: THE ERROR HERE GETS TOO LARGE WHEN THE GEOMETRIC PROBLEM SCALES UP...
        # NOTE: THE ERROR HERE GETS TOO LARGE WHEN THE GEOMETRIC PROBLEM SCALES UP...
        # NOTE: THE ERROR HERE GETS TOO LARGE WHEN THE GEOMETRIC PROBLEM SCALES UP...
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

# produce some number of closed loops from a set of edges
def sloops(es,epsilon = 0.1):
    #def fp(p):
    #    for j in range(g.vcnt):
    #        #if g.vs[j][1]['p'].isnear(p):
    #        if g.vs[j][1]['p'].d(p) < epsilon:
    #            return j
    #    return g.av(p = p.cp())
    g = pgr.planargraph()
    for e1,e2 in es:
        v1,v2 = g.fp(e1,epsilon),g.fp(e2,epsilon)
        if not v2 in g.rings[v1]:g.ae(v1,v2)
    uls = g.uloops('ccw')
    ols = [[g.vs[j][1]['p'].cp() for j in ul] for ul in uls]
    for ol in ols:
        if bnrm(ol).z < 0:ol.reverse()
    return ols

# segment a boundary polygon based on intersections with a line segment
# return two new boundary polygons
def bsegsxy(b,s1,s2):
    #bes = [(b[x-1].cp(),b[x].cp()) for x in range(len(b))]
    bes = bsegbxy(b,(s1,s2))
    left,right,center = [],[],[]
    while bes:
        u1,u2 = bes.pop(0)
        u1o = gtl.orient2d(s1,s2,u1)
        u2o = gtl.orient2d(s1,s2,u2)
        if u1o == 0 and u2o == 0:center.append((u1,u2))
        elif u1o > 0 or u2o > 0:left.insert(0,(u2,u1))
        elif u1o < 0 or u2o < 0:right.append((u1,u2))
        else:raise ValueError
    if not left or not right:
        print('seg missed bpoly')
        bs = [b]
    else:
        ips = sintbxyp(s1,s2,b,col = False)
        icnt = len(ips)
        iedges = [(ips[x-1].cp(),ips[x].cp()) for x in range(1,icnt) if x % 2 != 0]
        if iedges:bs = sloops(left+iedges)+sloops(right+iedges)
        else:bs = [b]
    return bs

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
    return b1es

# given two boundary polygons, 
# determine if any two edges are colinear and intersecting
def badjbxy(b1,b2,minovlp = 1):
    adjs = []
    for b1x in range(len(b1)):
        b1p1,b1p2 = b1[b1x-1],b1[b1x]
        for b2x in range(len(b2)):
            b2p1,b2p2 = b2[b2x-1],b2[b2x]
            b1p1nr = b1p1.isnear(b2p1) or b1p1.isnear(b2p2)
            b1p2nr = b1p2.isnear(b2p1) or b1p2.isnear(b2p2)
            ips = sintsxyp(b1p1,b1p2,b2p1,b2p2,ie = False,skew = False)
            if type(ips) == type(()):
                if ips[0].d(ips[1]) > minovlp:
                    adjs.append((b1x,b2x))
    return adjs

###############################################################################
###############################################################################

###############################################################################
### csg type functions for boundary polygons
###############################################################################

# compute the union of two boundary polygons
def ebuxy(b1,b2):
    b1segs,b2segs = bsegbxy(b1,b2),bsegbxy(b2,b1)
    bo = lambda s1,s2,b : s1.inbxy(b) or s2.inbxy(b) or (s1.onbxy(b) and s2.onbxy(b))
    b1inb2 = [p for p in b1segs if bo(p[0],p[1],b2)]
    b2inb1 = [p for p in b2segs if bo(p[0],p[1],b1)]
    b1only = [p for p in b1segs if not p in b1inb2]
    b2only = [p for p in b2segs if not p in b2inb1]
    dfs = sloops(b1only+b2only)
    if len(dfs) == 0:return None
    elif len(dfs) == 1:return dfs[0]
    else:return dfs

# compute difference of boundary polygons
def ebdxy(b1,b2):
    b1segs,b2segs = bsegbxy(b1,b2),bsegbxy(b2,b1)
    bo = lambda s1,s2,b : s1.inbxy(b) or s2.inbxy(b) or (s1.onbxy(b) and s2.onbxy(b))
    b1only = [p for p in b1segs if not bo(p[0],p[1],b2)]
    b2inb1 = [p for p in b2segs if p[0].inbxy(b1) or p[1].inbxy(b1)]
    dfs = sloops(b1only+b2inb1)
    if len(dfs) == 0:return None
    elif len(dfs) == 1:return dfs[0]
    else:return dfs

# compute the intersection of two boundary polygons
def ebixy(b1,b2):
    b1segs,b2segs = bsegbxy(b1,b2),bsegbxy(b2,b1)
    bo = lambda s1,s2,b : s1.inbxy(b) or s2.inbxy(b) or (s1.onbxy(b) and s2.onbxy(b))
    b1inb2 = [p for p in b1segs if bo(p[0],p[1],b2)]
    b2inb1 = [p for p in b2segs if bo(p[0],p[1],b1)]
    dfs = sloops(b1inb2+b2inb1)
    if len(dfs) == 0:return None
    elif len(dfs) == 1:return dfs[0]
    else:return dfs

###############################################################################
###############################################################################

# transform a point based on the xy projection of a boundary polygon
def ptob(b,p):
    x,y,z = vec3(1,0,0),vec3(0,1,0),vec3(0,0,1)
    bpx,bpy,bpz = x.prjps(b),y.prjps(b),z.prjps(b)
    bp = vec3(
        bpx[0]+p.x*(bpx[1]-bpx[0]),
        bpy[0]+p.y*(bpy[1]-bpy[0]),
        bpz[0]+p.z*(bpz[1]-bpz[0]))
    return bp

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

# return outward normals for each edge of a boundary polygon
def bnrmsxy(b):
    nms = []
    for x in range(len(b)):
        b1,b2 = b[x-1],b[x]
        nms.append(b1.tov(b2).crs(vec3(0,0,1)).nrm())
    return nms

# compute the distance of a point to the line containing 
# each edge of a boundary polygon
def bdistpxy(b,p):
    bnms = pym.bnrmsxy(b)
    bds = []
    for x in range(len(b)):
        bpp,bpn = b[x-1],bnms[x]
        bds.append(abs(p.dot(bpn)-bpp.dot(bpn)))
    return bds

# find the nearest edge of a boundary polygon to a point
# return the index of the edge and the associated distance
def bnearpxy(b,p):
    mind = 1e10
    minx = None
    for x in range(len(b)):
        b1,b2 = b[x-1],b[x]
        btn = b1.tov(b2)
        b1p = b1.dot(btn)
        b2p = b2.dot(btn)
        ptn = p.dot(btn)
        pd1 = p.d(b1)
        pd2 = p.d(b2)
        if pd1 < mind:mind,minx = pd1,x-1
        if pd2 < mind:mind,minx = pd2,x
        if (b1p <= ptn and ptn <= b2p) or (b1p >= ptn and ptn >= b2p):
            bnm = vec3(0,0,1).crs(btn).nrm()
            bed = abs(p.dot(bnm)-b1.dot(bnm))
            if bed < mind:
                mind = bed
                minx = x
    return minx,mind

# make a new boundary polygon including the midpoints of every edge
def bisectb(b):
    nb = []
    for x in range(len(b)):
        nb.append(b[x-1])
        nb.append(b[x-1].mid(b[x]))
    nb.append(nb.pop(0))
    nb.append(nb.pop(0))
    return nb

# evenly contract a boundary polygon
def contract(b,ds):
    fbnd = []
    pns = []
    for x in range(len(b)):
        p1,p2,p3 = b[x-2],b[x-1],b[x]
        w1t = p1.tov(p2).nrm()
        w2t = p2.tov(p3).nrm()
        w1n = vec3(0,0,1).crs(w1t)
        w2n = vec3(0,0,1).crs(w2t)
        s11 = p1.cp().trn(w1n).trn(w1t.cp().uscl(-1000))
        s12 = p2.cp().trn(w1n).trn(w1t.cp().uscl( 1000))
        s21 = p2.cp().trn(w2n).trn(w2t.cp().uscl(-1000))
        s22 = p3.cp().trn(w2n).trn(w2t.cp().uscl( 1000))
        ip = sintsxyp(s11,s12,s21,s22,ie = False,col = False)
        if ip is None:pn = p2.cp().trn(w1n)
        else:pn = ip.cp()
        pns.append(pn)
    for x in range(len(b)):
        p1,p2,p3 = b[x-2],b[x-1],b[x]
        fbnd.append(p2.lerp(pns[x],ds))
    fbnd.append(fbnd.pop(0))
    return fbnd

# return an xy plane laplacian smoothed version of a boundary polygon
def smoothxy(b,w = 0.1):
    db = []
    for x in range(len(b)):
        b0,b1,b2 = b[x-2],b[x-1],b[x]
        ecom = vec3(0,0,0).com((b0,b2))
        db.append(b1.tov(ecom).uscl(w))
    db.append(db.pop(0))
    sb = [p.cp().trn(ds) for p,ds in zip(b,db)]
    return sb

###############################################################################
###############################################################################





