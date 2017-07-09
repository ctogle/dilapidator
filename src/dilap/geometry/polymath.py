from .vec3 import vec3
from .planargraph import planargraph
from .tools import isnear, near, inrng, orient2d, epsilon
import numpy
import math
import pdb


def sintsxy(s11,s12,s21,s22,ie = True,ieb = True,col = True,skew = True,e = 0.01):
    '''Determine if two line segments intersect in the xy-plane
    ie  : do endpoint intersections count (end of one segment only)
    ieb : do endpoint intersections count (end of both segments)
    col : do colinear intersections counts
    if segments are colinear and overlapping
    elif segments are colinear and nonoverlapping
    elif segments are skew and possibly overlapping'''
    s1tn,s2tn = s11.tov(s12),s21.tov(s22)
    ds = s11.tov(s21)
    s1crss2 = s1tn.crs(s2tn)
    dscrss1 = ds.crs(s1tn)

    #s1s2area = gtl.near(s1crss2.mag(),0)
    #dss1area = gtl.near(dscrss1.mag(),0)

    s1s2area = s1crss2.mag()
    if s1s2area < e:s1s2area = 0
    dss1area = dscrss1.mag()
    if dss1area < e:dss1area = 0

    if s1s2area == 0 and dss1area == 0:
        if not col:return 0
        # if uncontained overlap
        # elif contained overlap
        # elif perfect overlap
        # elif single point overlap
        # else disjoint
        s1l = s1tn.mag()
        s1dots2 = s1tn.dot(s2tn)
        apara = near(s1dots2,0,e) < 0
        t0 = ds.dot(s1tn)/s1l**2
        t1 = t0+s1dots2/s1l**2
        t0 = near(near(t0,0,e),1,e)
        t1 = near(near(t1,0,e),1,e)
        if apara:t0,t1 = t1,t0
        if inrng(t0,0,1,e) or inrng(t1,0,1,e):return 1
        elif inrng(0,t0,t1,e) and inrng(1,t0,t1,e):return 1
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
        u = near(near(dscrss1.z/s1crss2.z,0,e),1,e)
        t = near(near(dscrss2.z/s1crss2.z,0,e),1,e)
        if not ieb and ((u == 0 or u == 1) and (t == 0 or t == 1)):return 0
        elif not ie and ((u == 0 or t == 0) or (u == 1 or t == 1)):return 0
        elif not (0 <= u and u <= 1) or not (0 <= t and t <= 1):return 0
        else:return 1


# given 4 points and a tangent vector, 
# extract the two whose projections are in between the extreme projections
def prjmed(p1,p2,p3,p4,tn):
    ps = [p1,p2,p3,p4]
    ss = [p1.dot(tn),p2.dot(tn),p3.dot(tn),p4.dot(tn)]
    minss,maxss,minfnd,maxfnd = min(ss),max(ss),False,False
    ssn = [0,1,2,3]
    for j in range(4):
        if not minfnd and ss[j] == minss:
            ssn.remove(j)
            minfnd = True
        elif not maxfnd and ss[j] == maxss:
            ssn.remove(j)
            maxfnd = True
    return ps[ssn[0]],ps[ssn[1]]


# does one segment intersect another
# ie  : do endpoint intersections count (end of one segment only)
# ieb : do endpoint intersections count (end of both segments)
# col : do colinear intersections counts
def sintsxyp(s11,s12,s21,s22,ie = True,ieb = True,col = True,skew = True,e = 0.01):
    # if segments are colinear and overlapping
    # elif segments are colinear and nonoverlapping
    # elif segments are skew and possibly overlapping
    #s11,s12,s21,s22 = s11.cpxy(),s12.cpxy(),s21.cpxy(),s22.cpxy()
    s1tn,s2tn = s11.tov(s12),s21.tov(s22)
    ds = s11.tov(s21)
    s1crss2 = s1tn.crs(s2tn)
    #s1s2area = gtl.near(s1crss2.mag(),0)
    s1s2area = s1crss2.mag()
    if s1s2area < e:
        #print('ROUND s1s2area',s1s2area)
        s1s2area = 0
    dscrss1 = ds.crs(s1tn)
    #dss1area = gtl.near(dscrss1.mag(),0)
    dss1area = dscrss1.mag()
    if dss1area < e:
        #print('ROUND dss1area',dss1area)
        dss1area = 0
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
        apara = near(s1dots2,0,e) < 0
        t0 = ds.dot(s1tn)/s1l**2
        t1 = t0+s1dots2/s1l**2
        t0 = near(near(t0,0,e),1,e)
        t1 = near(near(t1,0,e),1,e)
        if apara:t0,t1 = t1,t0
        if inrng(t0,0,1,e) or inrng(t1,0,1,e):
            return prjmed(s11,s12,s21,s22,s1tn)
        elif inrng(0,t0,t1,e) and inrng(1,t0,t1,e):
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
        #u = gtl.near(gtl.near(dscrss1.z/s1crss2.z,0),1)
        #t = gtl.near(gtl.near(dscrss2.z/s1crss2.z,0),1)
        u = dscrss1.z/s1crss2.z
        t = dscrss2.z/s1crss2.z
        if not ieb and ((u == 0 or u == 1) and (t == 0 or t == 1)):return None
        elif not ie and ((u == 0 or t == 0) or (u == 1 or t == 1)):return None
        elif not (0 <= u and u <= 1) or not (0 <= t and t <= 1):return None
        # NOTE: THE ERROR HERE GETS TOO LARGE WHEN THE GEOMETRIC PROBLEM SCALES UP...
        # NOTE: THE ERROR HERE GETS TOO LARGE WHEN THE GEOMETRIC PROBLEM SCALES UP...
        # NOTE: THE ERROR HERE GETS TOO LARGE WHEN THE GEOMETRIC PROBLEM SCALES UP...
        else:return s11.cp().trn(s1tn.cp().uscl(t))


# return acceptable endpoints wrt s1,s2 and s3,s4
def ssegsxy(s1,s2,s3,s4):
    endpoint = lambda p : p.isnear(s1) or p.isnear(s2) or p.isnear(s3) or p.isnear(s4)
    ips = sintsxyp(s1,s2,s3,s4)
    if isinstance(ips,tuple):
        tn = s1.tov(s2).nrm()
        prjs = (s1.dot(tn),s2.dot(tn),s3.dot(tn),s4.dot(tn))

        if not s4.onsxy(s1,s2,0):prjs.pop(3)
        if not s3.onsxy(s1,s2,0):prjs.pop(2)

        prjs = sorted(prjs)
        line = []
        for prj in prjs:
            p = s1.cp().trn(tn.cp().uscl(prj-s1.dot(tn)))
            clear = True
            for l in line:
                if l.isnear(p.cp()):
                    clear = False
                    break
            if clear:
                line.append(p.cp())
        #ax = dtl.plot_axes_xy(200)
        #ax = dtl.plot_points_xy(line,ax,number = True)
        #plt.show()
        return [(line[x-1],line[x]) for x in range(1,len(line))]
    elif isinstance(ips,vec3):
        if endpoint(ips):
            # split one of the segments in two...
            if   s1.onsxy(s3,s4,0):
                #return [(s1.cp(),s2.cp()),(s3.cp(),s1.cp()),(s4.cp(),s1.cp())]
                return [(s1.cp(),s2.cp())]
            elif s2.onsxy(s3,s4,0):
                #return [(s1.cp(),s2.cp()),(s3.cp(),s2.cp()),(s4.cp(),s2.cp())]
                return [(s1.cp(),s2.cp())]
            elif s3.onsxy(s1,s2,0):
                #return [(s3.cp(),s4.cp()),(s3.cp(),s2.cp()),(s3.cp(),s1.cp())]
                return [(s3.cp(),s2.cp()),(s3.cp(),s1.cp())]
            elif s4.onsxy(s1,s2,0):
                #return [(s3.cp(),s4.cp()),(s4.cp(),s2.cp()),(s4.cp(),s1.cp())]
                return [(s4.cp(),s2.cp()),(s4.cp(),s1.cp())]
            else:
                #return [(s1.cp(),s2.cp()),(s3.cp(),s4.cp())]
                return [(s1.cp(),s2.cp())]
        else:
            #return [
            #    (s1.cp(),ips.cp()),(s2.cp(),ips.cp()),
            #    (s3.cp(),ips.cp()),(s4.cp(),ips.cp())]
            return [(s1.cp(),ips.cp()),(s2.cp(),ips.cp())]
    else:
        #return [(s1.cp(),s2.cp()),(s3.cp(),s4.cp())]
        return [(s1.cp(),s2.cp())]


# produce some number of closed loops from a set of edges
def sloops(es,epsilon = 0.1):
    g = planargraph.segstopg(es,epsilon)
    uls = g.uloops('ccw')
    ols = [[g.vs[j][1]['p'].cp() for j in ul] for ul in uls]
    #for ol in ols:
    #    if bnrm(ol).z < 0:
    #        ol.reverse()
    return ols


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
            if type(ip) == type(1):
                pdb.set_trace()
            elif type(ip) == type(()) or type(ip) == type([]):
                isnew(ip[0])
                isnew(ip[1])
            else:isnew(ip)
    return ips


# segment a boundary polygon based on intersections with a line segment
# return some number of new boundary polygons
def bsegsxy(b,s1,s2,epsilon = 0.1):
    bes = bsegbxy(b,(s1,s2))
    left,right = [],[]
    while bes:
        u1,u2 = bes.pop(0)
        u1o = orient2d(s1,s2,u1)
        u2o = orient2d(s1,s2,u2)
        if u1o == 0 and u2o == 0:pass
        elif u1o > 0 or u2o > 0:left.insert(0,(u2,u1))
        elif u1o < 0 or u2o < 0:right.append((u1,u2))
        else:raise ValueError
    if not left or not right:
        print('seg missed bpoly')
        bs = [b]
    else:
        ips = sintbxyp(s1,s2,b,col = False)
        icnt = len(ips)
        ies = [(ips[x-1].cp(),ips[x].cp()) for x in range(1,icnt) if x % 2 != 0]
        if ies:bs = sloops(left+ies,epsilon)+sloops(right+ies,epsilon)
        else:bs = [b]
        if len(bs) > 2:
            print('multi bsegsxy...')

        '''#
        ax = dtl.plot_axes_xy(200)
        ax = dtl.plot_polygon_xy(b,ax,ls = '--',lw = 2,col = 'k')
        ax = dtl.plot_edges_xy((s1,s2),ax,col = 'g')
        for e in left:ax = dtl.plot_edges_xy(e,ax,col = 'r')
        for e in right:ax = dtl.plot_edges_xy(e,ax,col = 'b')
        ax = dtl.plot_points_xy(b,ax,number = True)
        for b in bs:ax = dtl.plot_polygon_xy(b,ax,col = 'g')
        plt.show()
        '''#

    return bs


# determine if a boundary polygon intersects itself
def bintselfxy(b):
    for b1x in range(len(b)):
        b1p1,b1p2 = b[b1x-1],b[b1x]
        for b2x in range(len(b)):
            if abs(b1x-b2x) < 2:continue
            b2p1,b2p2 = b[b2x-1],b[b2x]
            if sintsxy(b1p1,b1p2,b2p1,b2p2,True,True,True):
                return 1
    return 0


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
    ie = ie if col else col
    for b1x in range(len(b1)):
        b1p1,b1p2 = b1[b1x-1],b1[b1x]
        for b2x in range(len(b2)):
            b2p1,b2p2 = b2[b2x-1],b2[b2x]
            if sintsxy(b1p1,b1p2,b2p1,b2p2,ie = ie,ieb = ieb,col = col):
                return 1
    return 0


# segment a boundary polygon based on intersections with another
def bsegbxy(b1,b2):
    b1es = []
    unfn = [(b1[x-1],b1[x]) for x in range(len(b1))]
    while unfn:
        u1,u2 = unfn.pop(0)
        ips = sintbxyp(u1,u2,b2,ieb = False)
        ips = [p for p in ips if not (p.isnear(u1) or p.isnear(u2))]
        if len(ips) == 0:b1es.append((u1.cp(),u2.cp()))
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


# translate and scale (in place) b1 to be roughly inscribed by b2
def bfitbxy(b1,b2):
    def prj(b):
        prjx = vec3(1,0,0).prjps(b)
        prjy = vec3(0,1,0).prjps(b)
        return prjx,prjy
    # scale the polygon to the b
    b1prjx,b1prjy = prj(b1)
    b2prjx,b2prjy = prj(b2)
    cx = sum(b1prjx)-sum(b2prjx)
    cy = sum(b1prjy)-sum(b2prjy)
    recenter = vec3(cx,cy,0).uscl(-0.5)
    for p in b1:
        p.trn(recenter)
    b1prjx,b1prjy = prj(b1)
    #b2prjx,b2prjy = prj(b2)
    d1x = (b1prjx[1]-b1prjx[0])+(b1prjy[1]-b1prjy[0])
    d2x = (b2prjx[1]-b2prjx[0])+(b2prjy[1]-b2prjy[0])
    scale = 1. * d2x / d1x
    for p in b1:
        p.scl(vec3(scale,scale,0))
    return b1


# compute the union of many adjacent boundary polygons
# similar to ebuxy except using a planar graph
def ebuxy_special(bs,epsilon = 5,cellperimlength = 2):

    # do i need this function at all???
    # do i need this function at all???
    # do i need this function at all???

    #return bsuxy(bs,epsilon)

    # bs -> set of edges
    # -> ssegsxy until none overlap
    # -> sstopg using segments
    # -> counterclockwise loop around pg
    new = lambda s1,s2 : not ((s1,s2) in unfin or (s1,s2) in unfin)
    unfin = [(b[x-1],b[x]) for b in bs for x in range(len(b))]
    fin = []
    while unfin:
        u1,u2 = unfin.pop(0)
        clear = True
        for e1,e2 in unfin:
            if e1.onsxy(u1,u2,0):
                if new(u1,e1):unfin.append((u1,e1))
                if new(e1,u2):unfin.append((e1,u2))
                clear = False
                break
            elif e2.onsxy(u1,u2,0):
                if new(u1,e2):unfin.append((u1,e2))
                if new(e2,u2):unfin.append((e2,u2))
                clear = False
                break
            elif sintsxy(u1,u2,e1,e2,ie = False,ieb = False):
                plt.show()
                '''#
                ax = dtl.plot_axes_xy(200)
                for u in unfin:
                    ax = dtl.plot_edges_xy(u,ax,col = 'r',lw = 8)
                for f in fin:
                    ax = dtl.plot_edges_xy(f,ax,col = 'g',lw = 7)
                ax = dtl.plot_edges_xy((u1,u2),ax,col = 'c',lw = 5)
                ax = dtl.plot_edges_xy((e1,e2),ax,col = 'b',lw = 2)
                plt.show()
                '''#
                if sintsxy(u1,u2,e1,e2,ie = False,ieb = False,col = False):
                    ip = sintsxyp(u1,u2,e1,e2)
                    if new(u1,ip):unfin.append((u1,ip))
                    if new(ip,u2):unfin.append((ip,u2))
                    clear = False
                    break
                else:
                    if new(u1,u2):unfin.append((u1,u2))
                    clear = False
        if clear:
            fin.append((u1,u2))
    '''#
    unfin2 = [(b[x-1],b[x]) for b in bs for x in range(len(b))]
    ax = dtl.plot_axes_xy(200)
    for e in unfin2:
        ax = dtl.plot_edges_xy(e,ax,lw = 4,col = 'b')
    for e in fin:
        ax = dtl.plot_edges_xy(e,ax,lw = 2,col = 'g')
    plt.show()
    '''#
    uloops = bsuxy(sloops(fin,epsilon),epsilon)
    '''#
    ax = pg.plotxy(l = 150,aspect = 'equal')
    #for lp in loops:
    colors = ['b','g','r','c','m','y','k','g']
    for j in range(len(uloops)):
        lpp = uloops[j]
        #lpp = [pg.vs[lpx][1]['p'] for lpx in lp]
        lpp = contract(lpp,2)
        ax = dtl.plot_polygon_xy(lpp,ax,lw = (j+1)*2,col = colors[j])
        #ax = dtl.plot_points_xy(lpp,ax,number = True)
    plt.show()
    '''#
    return uloops


# compute the union of two full polygons
# return a list of full polygons
def epuxy(py1,py2,epsilon = 0.1,debug = False):
    b1,h1 = py1
    b2,h2 = py2
    pys = ebuxy(b1,b2,epsilon,debug,holes = True)
    if h1 or h2:
        print('epuxy not implemented warning...')
    #ax = dtl.plot_axes_xy(200)
    #for py in pys:
    #    ax = dtl.plot_polygon_full_xy(py,ax,lw = 2,col = 'b')
    #plt.show()
    return pys


# compute the union of two boundary polygons
def ebuxy(b1,b2,epsilon = 0.1,debug = False,holes = False):
    b1segs,b2segs = bsegbxy(b1,b2),bsegbxy(b2,b1)
    bo = lambda s1,s2,b : s1.mid(s2).inbxy(b)
    def inn(p):
        pn = vec3(0,0,1).crs(p[0].tov(p[1])).nrm().uscl(2*epsilon)
        pm = p[0].mid(p[1])
        ep = pm.cp().trn(pn)
        if not ep.inbxy(b1) and not ep.inbxy(b2):return False
        ep = pm.cp().trn(pn.flp())
        if not ep.inbxy(b1) and not ep.inbxy(b2):return False
        return True
    b1inb2 = [p for p in b1segs if bo(p[0],p[1],b2)]
    b2inb1 = [p for p in b2segs if bo(p[0],p[1],b1)]
    b1only = [p for p in b1segs if not p in b1inb2 and not inn(p)]
    b2only = [p for p in b2segs if not p in b2inb1 and not inn(p)]
    dfs = sloops(b1only+b2only,epsilon)

    if debug:
        ax = dtl.plot_axes_xy(200)
        ax = dtl.plot_polygon_xy(b2,ax,lw = 2,col = 'm')
        for edge in b1segs:
            ax = dtl.plot_edges_xy(edge,ax,lw = 5,col = 'b')
        for edge in b1only:
            ax = dtl.plot_edges_xy(edge,ax,lw = 8,col = 'b')
        for edge in b1inb2:
            ax = dtl.plot_edges_xy(edge,ax,lw = 3,col = 'g')
        #for edge in b2only:
        #    ax = dtl.plot_edges_xy(edge,ax,lw = 5,col = 'g')
        plt.show()

    if len(dfs) > 1:
        edfs = 0
        for jdfs,lp in enumerate(dfs):
            if binbxy(dfs[edfs],lp,ie = 0):
                edfs = jdfs
        if holes:
            print('genus of ebuxy changed!')
            return [(dfs.pop(edfs),dfs)]
        else:
            print('discarding holes resulting from ebuxy!')
            return [dfs[edfs]]
        '''
        print('promisedland?')
        ax = dtl.plot_axes_xy(200)
        ax = dtl.plot_polygon_xy(dfs[edfs],ax,lw = 7,col = 'r')
        ax = dtl.plot_polygon_xy(dfs[0],ax,lw = 5,col = 'b')
        ax = dtl.plot_polygon_xy(dfs[1],ax,lw = 3,col = 'g')
        #ax = dtl.plot_polygon_xy(b1,ax,ls = '-.',lw = 5,col = 'k')
        #ax = dtl.plot_polygon_xy(b2,ax,ls = '--',lw = 5,col = 'k')
        plt.show()
        '''
    if holes:return [(dfs[0],[])]
    else:return dfs


# compute difference of boundary polygons
def ebdxy(b1,b2,epsilon = 0.1):
    b1segs,b2segs = bsegbxy(b1,b2),bsegbxy(b2,b1)
    bo = lambda s1,s2,b : s1.mid(s2).inbxy(b)
    b1only = [p for p in b1segs if not bo(p[0],p[1],b2) and not p in b2segs]
    b2inb1 = [p for p in b2segs if bo(p[0],p[1],b1)]
    dfs = sloops(b1only+b2inb1,epsilon)
    # need to throw away nonsatisfactory loops
    #dfs = [l for l in dfs]
    '''#
    ax = dtl.plot_axes(100)
    ax = dtl.plot_polygon(b1,ax,lw = 2,col = 'b')
    ax = dtl.plot_polygon(b2,ax,lw = 2,col = 'g')
    for x in range(len(dfs)):
        vec3(0,0,10*x+2).trnos(dfs[x])
        ax = dtl.plot_polygon(dfs[x],ax,col = 'r')
    plt.show()
    '''#
    return dfs


# compute the intersection of two boundary polygons
def ebixy(b1,b2,epsilon = 0.1):
    b1segs,b2segs = bsegbxy(b1,b2),bsegbxy(b2,b1)
    #bo = lambda s1,s2,b : s1.inbxy(b) or s2.inbxy(b) or (s1.onbxy(b) and s2.onbxy(b))
    #bo = lambda s1,s2,b : s1.mid(s2).inbxy(b)
    bo = lambda s1,s2,b : s1.mid(s2).inbxy(b) or s1.mid(s2).onbxy(b)
    b1inb2 = [p for p in b1segs if bo(p[0],p[1],b2)]
    b2inb1 = [p for p in b2segs if bo(p[0],p[1],b1)]
    dfs = sloops(b1inb2+b2inb1,epsilon)
    return dfs


# consolidate a list of boundary polygons using ebuxy
def bsuxy(bs,epsilon = 0.1):
    unfn = bs[:]
    nbs = []
    while unfn:
        ub = unfn.pop(0)
        fnd = False
        for nbx in range(len(nbs)):
            b = nbs[nbx]
            if   binbxy(ub,b):pass
            elif binbxy(b,ub):nbs[nbx] = ub
            #elif bintbxy(ub,b):
            elif bintbxy(ub,b,ie = True,ieb = True):
                nbbb = nbs.pop(nbx)
                newunfn = ebuxy(nbbb,ub,epsilon)

                '''#
                if len(newunfn) == 3:
                    print('len((((',len(newunfn))
                    ax = dtl.plot_axes_xy(500)
                    ax = dtl.plot_polygon_xy(contract(ub,2),ax,col = 'r')
                    ax = dtl.plot_polygon_xy(contract(nbbb,2),ax,col = 'b')
                    for u in newunfn:
                        ax = dtl.plot_polygon_xy(contract(u,2),ax)
                    plt.show()
                '''#

                for u in newunfn:unfn.insert(-1,u)
                fnd = True
                break
        if not fnd:nbs.append(ub)

        '''#
        ax = dtl.plot_axes_xy(110)
        for b in bs:ax = dtl.plot_polygon_xy(b,ax,lw = 1,col = 'b')
        for b in unfn:ax = dtl.plot_polygon_xy(b,ax,lw = 2,col = 'r')
        for b in nbs:ax = dtl.plot_polygon_xy(b,ax,lw = 2,col = 'g')
        plt.show()
        '''#

    return nbs


# like vec3.inbxy except uses a ray intersection.. seeems to work better?
def pinb(b,p):
    s2 = p.cp().xtrn(10000)
    ips = sintbxyp(p,s2,b,col = False)
    ins = len(ips) % 2 == 1
    return ins


# transform a point based on the xy projection of a boundary polygon
def ptob(b,p):
    x,y,z = vec3(1,0,0),vec3(0,1,0),vec3(0,0,1)
    bpx,bpy,bpz = x.prjps(b),y.prjps(b),z.prjps(b)
    bp = vec3(
        bpx[0]+p.x*(bpx[1]-bpx[0]),
        bpy[0]+p.y*(bpy[1]-bpy[0]),
        bpz[0]+p.z*(bpz[1]-bpz[0]))
    return bp


# determine the properness of a boundary polygon
def bvalidxy(b,epsilon = 0.1):
    # polygon should have at least 3 points
    if len(b) < 3:
        print('polygon contains fewer than 3 points')
        return 0
    # polygon should be counter clockwise ordered
    if not bccw(b):
        print('polygon is clockwise')
        return -1
    # normal should be pointing in the positive z direction
    #if bnrm(b).z <= 0:
    #    #print('polygon normal is improper')
    #    print('polygon normal is improper... letting it slide...')
    #    #return -1
    # polygon should not have duplicate points
    for x in range(len(b)):
        bp = b[x]
        for y in range(len(b)):
            if x == y:continue
            if b[x].d(b[y]) < epsilon:
                print('polygon points are within epsilon:',epsilon)
                return -2
    # polygon should not have self intersections
    for x in range(len(b)):
        b1,b2 = b[x-1],b[x]
        for y in range(len(b)):
            if x == y:continue
            b3,b4 = b[y-1],b[y]
            if sintsxy(b1,b2,b3,b4,ie = False):
                print('polygon is self-intersecting')
                return -3
    return 1


# fix an invalid polygon such that it bvalidxy(b) returns True
def cleanbxy(b,epsilon = 0.1):
    r = bvalidxy(b)
    while not r > 0:
        if   r == 0:raise ValueError
        elif r == -1:b.reverse()
        elif r == -2:b = aggregate(b,epsilon)
        elif r == -3:b = pinchb(b,epsilon)
        r = bvalidxy(b)
    return b


# compute the area of a boundary polygon in the xy plane
def bareaxy(b,allowzero = False):
    area = 0.0
    for x in range(len(b)):
        b1,b2 = b[x-1],b[x]
        area += (b1.x+b2.x)*(b1.y-b2.y)
    if not allowzero and area == 0:raise ValueError
    return -area/2.0


# return a vector normal to the boundary polygon b
def bnrm(b):
    pcnt = len(b)
    pn = vec3(0,0,0)
    for px in range(pcnt):
        c1,c2,c3  = b[px-2],b[px-1],b[px]
        c1c2 = c1.tov(c2)
        c2c3 = c2.tov(c3)
        #if gtl.isnear(c2c3.mag(),0) or gtl.isnear(c1c2.mag(),0):
        if c2c3.mag() == 0 or c1c2.mag() == 0:
            print('bnrm adjacency error!!',c1c2,c2c3)
            #ax = dtl.plot_axes_xy(200)
            #ax = dtl.plot_polygon_xy(b,ax,lw = 2,col = 'g')
            #ax = dtl.plot_points_xy(b,ax,number = True)
            #plt.show()
            raise ValueError
        cang = c1c2.ang(c2c3)
        if not isnear(cang,0,epsilon) and not isnear(cang,numpy.pi,epsilon):
            if cang > -numpy.pi+0.001 and cang < numpy.pi-0.001:
                pn.trn(c1c2.crs(c2c3).nrm())
    return pn.nrm()


# determine if a polygon is ordered counterclockwise in the xy plane
def bccw(b):
    s = 0
    for x in range(len(b)):
        b1,b2 = b[x-1],b[x]
        s += (b2.x-b1.x)*(b2.y+b1.y)
    return s < 0


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
    bnms = bnrmsxy(b)
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


# make a new boundary polygon including the midpoints of every edge
def splitb(b,l):
    nb = []
    for x in range(len(b)):
        b1,b2 = b[x-1],b[x]
        nb.append(b1.cp())
        el = b1.d(b2)
        if el > l:
            n = 1
            while el/n > l:n += 1
            divpts = b1.pline(b2,n-1)
            divcnt = len(divpts)
            nb.extend([p.cp() for p in divpts])

    '''#
    ax = dtl.plot_axes_xy(500)
    ax = dtl.plot_polygon_xy(nb,ax,lw = 2,col = 'b')
    ax = dtl.plot_points_xy(nb,ax,number = True)
    plt.show()
    '''#

    return nb


# return a copy of a boundary polygon with some points removed
# possibly replace dissolved points with points in rps
def bdissolvep(b,xs,rps):
    nb = []
    for y in range(len(b)):
        if y in xs:
            j = xs.index(y)
            if rps[j] is None:continue
            else:nb.extend(rps[j])
        else:nb.append(b[y])
        #nb.append(b[y])
    return nb


# modify a boundary polygon such that its hmin will be above minhmin
# NOTE: the resulting polygon should be contained by the input polygon
def blimithmin(b,minhmin,maxhmin):
    m = lambda : min([b[x-1].d(b[x]) for x in range(len(b))])*math.sqrt(3)
    els = [b[x-1].d(b[x]) for x in range(len(b))]
    dps,rps = [],[]
    for elx in range(len(els)):
        el = els[elx]
        if el < minhmin:
            rps.append(None);rps.append([b[elx-1].mid(b[elx])])
            dps.append(elx-1);dps.append(elx)
    if dps:b = bdissolvep(b,dps,rps)
    b = splitb(b,maxhmin)

    #b = smoothxy(b,1.0,1.0)
    #b = smoothxy(b,1.0,1.0)
    #b = smoothxy(b,1.0,1.0)
    #b = smoothxy(b,1.0,1.0)
    #b = smoothxy(b,1.0,1.0)

    '''#
    g = sstopg([(b[x-1],b[x]) for x in range(len(b))],1,True)
    uls = g.uloops('ccw')
    #ax = dtl.plot_axes_xy(500)
    ax = g.plotxy(l = 500)
    for lp in uls[1:]:
        lpps = [g.vs[j][1]['p'] for j in lp]
        print('lPPPS',len(lpps))
        #lpps = contract(lpps,0.1)
        ax = dtl.plot_polygon_xy(lpps,ax,lw = 3,col = 'g')
        ax = dtl.plot_points_xy(lpps,ax,number = True)
    plt.show()
    '''#

    return b


# remove additional loops created by intersections 
# in an improper boundary polygon
def pinchb(b,epsilon = 5):
    pde = lambda p,e1,e2 : abs()
    new = lambda s1,s2 : not ((s1,s2) in unfin or (s1,s2) in unfin)
    unfin = [(b[x-1],b[x]) for x in range(len(b))]
    fin = []
    while unfin:
        u1,u2 = unfin.pop(0)
        clear = True
        for e1,e2 in unfin:
            e1u12d = e1.dexy(u1,u2,0)
            e2u12d = e2.dexy(u1,u2,0)
            if False and e1u12d > -1 and e1u12d < epsilon:

                ax = dtl.plot_axes_xy(200)
                ax = dtl.plot_polygon_xy(b,ax,lw = 2,ls = '--',col = 'k')
                ax = dtl.plot_edges_xy((u1,u2),ax,lw = 3,col = 'b')
                ax = dtl.plot_edges_xy((e1,e2),ax,lw = 3,col = 'g')
                plt.show()

                pdb.set_trace()
            elif False and e2u12d > -1 and e2u12d < epsilon:

                ax = dtl.plot_axes_xy(200)
                ax = dtl.plot_polygon_xy(b,ax,lw = 2,ls = '--',col = 'k')
                ax = dtl.plot_edges_xy((u1,u2),ax,lw = 3,col = 'b')
                ax = dtl.plot_edges_xy((e1,e2),ax,lw = 3,col = 'g')
                plt.show()

                pdb.set_trace()
            elif e1.onsxy(u1,u2,0):
                if new(u1,e1):unfin.append((u1,e1))
                if new(e1,u2):unfin.append((e1,u2))
                clear = False
                break
            elif e2.onsxy(u1,u2,0):
                if new(u1,e2):unfin.append((u1,e2))
                if new(e2,u2):unfin.append((e2,u2))
                clear = False
                break
            elif sintsxy(u1,u2,e1,e2,ie = False,ieb = False):
                '''#
                ax = dtl.plot_axes_xy(200)
                for u in unfin:
                    ax = dtl.plot_edges_xy(u,ax,col = 'r',lw = 8)
                for f in fin:
                    ax = dtl.plot_edges_xy(f,ax,col = 'g',lw = 7)
                ax = dtl.plot_edges_xy((u1,u2),ax,col = 'c',lw = 5)
                ax = dtl.plot_edges_xy((e1,e2),ax,col = 'b',lw = 2)
                plt.show()
                '''#
                if sintsxy(u1,u2,e1,e2,ie = False,ieb = False,col = False):
                    ip = sintsxyp(u1,u2,e1,e2)
                    if new(u1,ip):unfin.append((u1,ip))
                    if new(ip,u2):unfin.append((ip,u2))
                    clear = False
                    break
                else:
                    if new(u1,u2):unfin.append((u1,u2))
                    clear = False
        if clear:
            fin.append((u1,u2))

    uloops = sloops(fin,epsilon)
    if len(uloops) > 1:
        #las = [bareaxy(ol) for ol in uloops]
        las = [bareaxy(ol,True) for ol in uloops]
        ulas = [(u,a) for u,a in zip(uloops,las) if a > 10]
        if ulas:
            uloops,las = zip(*ulas)
            uloops,las = list(uloops),list(las)
            uloops.insert(0,uloops.pop(las.index(max(las))))
    #print('PINCHB ACCOMPLISHED')
    return uloops


# evenly contract a boundary polygon
def contract(b,ds,epsilon = 0.1):
    if not type(ds) is type([]):
        dss = ds
        ds = [1 for x in range(len(b))]
    else:dss = 1
    fbnd = []
    pns = []
    for x in range(len(b)):
        p1,p2,p3 = b[x-2],b[x-1],b[x]
        w1t = p1.tov(p2).nrm()
        w2t = p2.tov(p3).nrm()
        w1n = vec3(0,0,1).crs(w1t).nrm().uscl(ds[x-2])
        w2n = vec3(0,0,1).crs(w2t).nrm().uscl(ds[x-1])
        s11 = p1.cp().trn(w1n).trn(w1t.cp().uscl(-10000))
        s12 = p2.cp().trn(w1n).trn(w1t.cp().uscl( 10000))
        s21 = p2.cp().trn(w2n).trn(w2t.cp().uscl(-10000))
        s22 = p3.cp().trn(w2n).trn(w2t.cp().uscl( 10000))

        ip = sintsxyp(s11,s12,s21,s22,ie = False,col = False)

        #if ip is None:
        if ip is None or isnear(abs(w1t.dot(w2t)),1,epsilon):
            pn = p2.cp().trn(w1n)

        else:pn = ip.cp()
        pns.append(pn)
    for x in range(len(b)):
        p1,p2,p3 = b[x-2],b[x-1],b[x]
        fbnd.append(p2.lerp(pns[x],dss))
        #fbnd.append(p2.lerp(pns[x],random.uniform(0.5*ds,ds)))

        '''#
        print('amen5')
        ax = dtl.plot_axes_xy(500)
        ax = dtl.plot_polygon_xy(fbnd,ax,lw = 2)
        ax = dtl.plot_point_xy_annotate(p1,ax,'p1')
        ax = dtl.plot_point_xy(p1,ax)
        ax = dtl.plot_point_xy_annotate(p2,ax,'p2')
        ax = dtl.plot_point_xy(p2,ax)
        ax = dtl.plot_point_xy_annotate(p3,ax,'p3')
        ax = dtl.plot_point_xy(p3,ax)
        plt.show()
        print('amen6')
        '''#

    fbnd.append(fbnd.pop(0))

    '''#
    print('amen3')
    ax = dtl.plot_axes_xy(500)
    ax = dtl.plot_polygon_xy(fbnd,ax,lw = 2)
    plt.show()
    print('amen4')
    '''#

    #if len(fbnd) > 3:fbnd = pinchb(fbnd,epsilon)
    return fbnd


# return an xy plane laplacian smoothed version of a boundary polygon
def smoothxy(b,w = 0.1,epsilon = 0.1,constraint = 0):
    db = []
    for x in range(len(b)):
        b0,b1,b2 = b[x-2],b[x-1],b[x]
        ecom = vec3(0,0,0).com((b0,b2))
        if constraint:
            closer = b1.cp().trn(ecom.cp().uscl(2*epsilon)).inbxy(b)
            if constraint == -1 and closer:db.append(vec3(0,0,0))
            elif constraint == 1 and not closer:db.append(vec3(0,0,0))
            else:db.append(b1.tov(ecom).uscl(w))
        else:db.append(b1.tov(ecom).uscl(w))
    db.append(db.pop(0))
    sb = [p.cp().trn(ds) for p,ds in zip(b,db)]
    if len(sb) > 3:sb = pinchb(sb,epsilon)[0]
    return sb


def smoothxyi(b,w = 0.1,epsilon = 0.1,i = 10,constraint = 0):
    for j in range(i):
        b = smoothxy(b,w,epsilon,constraint)
        print('smoothxyi',i)
    return b


# return a version of a boundary polygon where 
# points within ds of one another are merged
def aggregate(b,ds,da = numpy.pi/4):
    if len(b) <= 3:
        #print('polygon could not be further aggregated')
        return 
        #return b[:]

    edists = [b[x-1].d(b[x]) for x in range(len(b))]
    edas = [b[x-1].tov(b[x-2]).angxy(b[x-1].tov(b[x])) for x in range(len(b))]
    edas.append(edas.pop(0))

    dmerges = [1 if d < ds else 0 for d in edists]
    amerges = [1 if a < da else 0 for a in edas]

    nb = []
    piece = []
    for x in range(len(b)):
        if not dmerges[x] == 0:
            piece.append(b[x])
        if not amerges[x] == 0:
            u,v,w = (x-1,x,x+1) if x < len(b)-1 else (x-1,x,0)
            piece.extend([b[u],b[v],b[w]])
        if piece:
            #nb.append(vec3(0,0,0).com(piece))
            #nb.append(piece[int(len(piece)/2.0)])
            piece = []
        else:nb.append(b[x])
    return nb


# generate a full polygon from a planar graph
##### why is it sometimes useful to flip z???
def pgtopy(pg,r,epsilon = 0.1,z = vec3(0,0,1),findeseam = False):
    loops,seams,eseam = pg.uloops('ccw'),[],None
    for lp in loops:
        seam = []
        seams.append(seam)

        # IF PG IS NONPLANAR, PY SHOULD BE NONPLANAR
        # NONPLANAR PG CAUSES ISSUES IN THIS COMPUTATION CURRENTLY
        # MODIFY TO USE XY PROJECTION OF PG INSTEAD!!!

        for lpx in range(len(lp)):
            lp0,lp1,lp2 = lp[lpx-2],lp[lpx-1],lp[lpx]
            lpp0 = pg.vs[lp0][1]['p']
            lpp1 = pg.vs[lp1][1]['p']
            lpp2 = pg.vs[lp2][1]['p']
            lptn1,lptn2 = lpp0.tov(lpp1).nrm(),lpp1.tov(lpp2).nrm()
            lpnm1,lpnm2 = z.crs(lptn1).uscl(r),z.crs(lptn2).uscl(r)
            if lp0 == lp2:
                stemoffset = lptn1.cp().uscl(r)
                seam.append(lpp1.cp().trn(lpnm1+stemoffset))
                seam.append(lpp1.cp().trn(lpnm2+stemoffset))
            else:
                if lpnm1.isnear(lpnm2):cnm = lpnm1
                else:
                    s1 = lpp0.cp().trn(lpnm1).trn(lptn1.cp().uscl(-1000))
                    s2 = lpp1.cp().trn(lpnm1).trn(lptn1.cp().uscl( 1000))
                    s3 = lpp1.cp().trn(lpnm2).trn(lptn2.cp().uscl(-1000))
                    s4 = lpp2.cp().trn(lpnm2).trn(lptn2.cp().uscl( 1000))
                    ip = sintsxyp(s1,s2,s3,s4,col = 0)
                    cnm = lpp1.tov(ip)
                seam.append(lpp1.cp().trn(cnm))
        if not bccw(seam):seam.reverse()
        if eseam is None:eseam = 0
        else:
            if binbxy(seams[eseam],seam):
                eseam = len(seams)-1
    if eseam is None:return None
    if findeseam:return eseam,seams,loops
    py = [seams.pop(eseam),seams]
    return py


# generate a full polygon from a planar graph (different from pgtopy)
def pgbleed(pg,r,epsilon = 0.1):

    # get a copy of the points of the graph without 
    # duplicates within neighborhood of epsilon
    ps = []
    for vx in range(pg.vcnt):
        ps.append(pg.fp(pg.vs[vx][1]['p'],epsilon))

    pdb.set_trace()


    '''#
                    ax = dtl.plot_axes_xy(200)
                    ax = dtl.plot_edges_xy(seam,ax,lw = 2,col = 'k')
                    ax = dtl.plot_edges_xy((s1,s2),ax,col = 'b')
                    ax = dtl.plot_edges_xy((s3,s4),ax,col = 'g')
                    ax = dtl.plot_point_xy(ip,ax,mk = 's')
                    plt.show()
    '''#

    py = [eb,ibs]

    ax = dtl.plot_axes_xy(200)
    ax = dtl.plot_polygon_full_xy(py,ax,lw = 2,col = 'k')
    plt.show()

    return py

