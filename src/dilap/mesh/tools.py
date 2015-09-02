import dilap.core.tools as dpr
import dilap.core.vector as dpv
import dilap.core.quaternion as dpq
import dilap.core.bbox as dpb
import dilap.core.profiler as dprf
import dilap.mesh.triangulate as dtg

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import math,numpy

import pdb

###############################################################################

def plot_axes_xy():
    ax = plt.figure().add_subplot(111)
    ax.set_aspect('equal')
    return ax

def plot_axes(x = 5):
    ax = plt.figure().add_subplot(111,projection = '3d')
    ax.set_xlim([-x,x])
    ax.set_ylim([-x,x])
    ax.set_zlim([-(9.0/16.0)*x,(9.0/16.0)*x])
    return ax

def plot_point_xy(pt,ax,marker = 'o'):
    ax.plot([pt.x],[pt.y],marker = marker)
    return ax

def plot_point_xy_annotate(pt,ax,text):
    ax.annotate(text,xy = (pt.x,pt.y),xytext = (-20, 20),
        textcoords = 'offset points',ha = 'right',va = 'bottom',
        arrowprops = dict(arrowstyle = '->',
        connectionstyle = 'arc3,rad=0'))
    return ax

def plot_point(pt,ax,marker = 'o'):
    ax.plot([pt.x],[pt.y],zs = [pt.z],marker = marker)
    return ax

def plot_points_xy(points,ax = None,ms = None,number = False):
    if ax is None:ax = plot_axes_xy()
    if ms is None:ms = ['o']*len(points)
    for pdx in range(len(points)):
        plot_point_xy(points[pdx],ax,ms[pdx])
        if number:plot_point_xy_annotate(points[pdx],ax,str(pdx+1))
    return ax

def plot_points(points,ax = None,ms = None):
    if ax is None:ax = plot_axes()
    if ms is None:ms = ['o']*len(points)
    for pdx in range(len(points)):plot_point(points[pdx],ax,ms[pdx])  
    return ax

def plot_edges_xy(points,ax = None,mk = None,lw = 1.0,center = False):
    if ax is None:ax = plot_axes_xy()
    if mk is None:mk = '+'
    pts = [p.to_tuple() for p in points]
    xs,ys,zs = zip(*pts)
    ax.plot(xs,ys,marker = mk,lw = lw)
    if center:
        centers = [dpv.midpoint(points[x-1],points[x]) 
                        for x in range(1,len(points))]
        plot_points_xy(centers,ax)
    return ax

def plot_edges(points,ax = None,lw = 1.0,center = False):
    if ax is None:ax = plot_axes()
    pts = [p.to_tuple() for p in points]
    xs,ys,zs = zip(*pts)
    ax.plot(xs,ys,zs,marker = '+',lw = lw)
    if center:
        centers = [dpv.midpoint(points[x-1],points[x]) 
                        for x in range(1,len(points))]
        plot_points(centers,ax)
    return ax

def plot_polygon_xy(points,ax = None,center = False,lw = 1.0):
    epts = points[:]
    epts.append(points[0])
    ax = plot_edges_xy(epts,ax,lw = lw)
    if center:plot_point_xy(dpv.center_of_mass(points),ax,marker = 's')
    return ax

def plot_polygon(points,ax = None,center = False,lw = 1.0):
    epts = points[:]
    epts.append(points[0])
    ax = plot_edges(epts,ax,lw = lw)
    if center:plot_point(dpv.center_of_mass(points),ax,marker = 's')
    return ax

def plot_polygon_full_xy(poly,ax = None,center = False,lw = 1.0):
    if ax is None:ax = plot_axes_xy()
    ebnd,ibnds = poly
    plot_polygon_xy(list(ebnd),ax,center = True,lw = lw)
    for ib in ibnds:plot_polygon_xy(list(ib),ax,center = True,lw = lw)
    return ax

def plot_polygon_full(poly,ax = None,center = False,lw = 1.0):
    if ax is None:ax = plot_axes()
    ebnd,ibnds = poly
    plot_polygon(list(ebnd),ax,center = True,lw = lw)
    for ib in ibnds:plot_polygon(list(ib),ax,center = True,lw = lw)
    return ax

def plot_tetrahedron(points,ax = None):
    raise NotImplemented

def plot_circle_xy(c,r,ax = None,center = False):
    circ = dpv.translate_coords(dpr.point_ring(r,32),c)
    ax = plot_polygon_xy(circ,ax,center)
    return ax

def plot_circle(c,r,ax = None,center = False):
    circ = dpv.translate_coords(dpr.point_ring(r,32),c)
    ax = plot_polygon(circ,ax,center)
    return ax

###############################################################################

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
# given poly, a sequence of points, return the lengths of 
# the edges of the polygon that poly describes
def edge_lengths(*poly):
    elengs = []
    for x in range(len(poly)):
        p1,p2 = poly[x-1],poly[x]
        elengs.append(dpv.distance(p1,p2))
    return elengs

# given poly, a sequence of points, 
# return the length of the shortest edge
def shortest_edge_tri(*tri):
    p1,p2,p3 = tri
    e1 = dpv.distance(p1,p2)
    e2 = dpv.distance(p2,p3)
    e3 = dpv.distance(p3,p1)
    return min((e1,e2,e3))

# given poly, a sequence of points, 
# return the length of the shortest edge
def shortest_edge(*poly):
    elengs = edge_lengths(*poly)
    return min(elengs)

# given line segment s1, line segment s2
# a segment is a tuple of two points
# return the point of intersection if there is one
# otherwise return None
#   note: currently assumes s1,s2 in xy-plane
def segments_intersect_at(s11,s12,s21,s22):
    p,q = s11,s21
    r = dpv.v1_v2(s11,s12)
    s = dpv.v1_v2(s21,s22)
    qmp = q-p
    rcs = r.cross(s)
    rcsmag = rcs.magnitude()
    qmpcr = qmp.cross(r)
    qmpcrmag = qmpcr.magnitude()
    rmag2 = r.magnitude2()
    if dpr.isnear(rcsmag,0) and dpr.isnear(qmpcrmag,0):

        if dpr.isnear(rmag2,0):
          print('wtffffffffffffffffff',s11,s12,s21,s22,rmag2)
          ax = plot_axes()
          plot_edges([s11,s12],ax,center = True,lw = 5.0)
          plot_edges([s21,s22],ax,center = True)
          plt.show()
          pdb.set_trace()

        t0 = dpr.near(dpr.near(dpv.dot(qmp,r)/rmag2,0),1)
        t1 = dpr.near(dpr.near( t0 + s.dot(r)/rmag2,0),1)
        if dpv.dot(s,r) < 0.0:t0,t1 = t1,t0
        t0pt = p + r.copy().scale_u(t0)
        t1pt = p + r.copy().scale_u(t1)
        if dpr.inrange(t0,0,1) or dpr.inrange(t1,0,1):return t0pt,t1pt
        elif dpr.inrange(0,t0,t1) and dpr.inrange(1,t0,t1):return s11,s12
        elif (t0 == 0 and t1 == 1) or (t0 == 1 and t1 == 0):return t0pt,t1pt
    elif dpr.isnear(rcsmag,0) and not dpr.isnear(qmpcrmag,0):return None
    #elif not dpr.isnear(rcsmag,0) and not dpr.isnear(qmpcrmag,0):
    #elif not dpr.isnear(rcsmag,0):
    elif not dpr.isnear(rcs.z,0):

        '''#
        if dpr.isnear(rcs.z,0):
          ax = plot_axes()
          plot_edges([s11,s12],ax,center = True,lw = 5.0)
          plot_edges([s21,s22],ax,center = True)
          plt.show()
          pdb.set_trace()
        '''#

        u = dpr.near(dpr.near(       qmpcr.z/rcs.z,0),1)
        t = dpr.near(dpr.near(qmp.cross(s).z/rcs.z,0),1)
        if (u == 0 or u == 1) and (t == 0 or t == 1):return
        #print('uvt',u,t)
        if not dpr.inrange(u,0,1) or not dpr.inrange(t,0,1):return None
        else:return q + s.scale_u(u)

# generate a plc for an icosphere with faces split n times
def icosphere(r = 1,n = 1):
    import dilap.mesh.piecewisecomplex as pwc
    plc = pwc.piecewise_linear_complex()

    # create 12 vertices of a icosahedron
    t = (1.0+math.sqrt(5.0))/2.0
    pts = []
    pts.append(dpv.vector(-1, t, 0))
    pts.append(dpv.vector( 1, t, 0))
    pts.append(dpv.vector(-1,-t, 0))
    pts.append(dpv.vector( 1,-t, 0))
    pts.append(dpv.vector( 0,-1, t))
    pts.append(dpv.vector( 0, 1, t))
    pts.append(dpv.vector( 0,-1,-t))
    pts.append(dpv.vector( 0, 1,-t))
    pts.append(dpv.vector( t, 0,-1))
    pts.append(dpv.vector( t, 0, 1))
    pts.append(dpv.vector(-t, 0,-1))
    pts.append(dpv.vector(-t, 0, 1))
    dpv.translate_coords(pts,dpv.center_of_mass(pts).flip())

    triangles = []
    # 5 faces around point 0
    triangles.append((0,11,5))
    triangles.append((0,5,1))
    triangles.append((0,1,7))
    triangles.append((0,7,10))
    triangles.append((0,10,11))
    # 5 adjacent faces
    triangles.append((1,5,9))
    triangles.append((5,11,4))
    triangles.append((11,10,2))
    triangles.append((10,7,6))
    triangles.append((7,1,8))
    # 5 faces around point 3
    triangles.append((3,9,4))
    triangles.append((3,4,2))
    triangles.append((3,2,6))
    triangles.append((3,6,8))
    triangles.append((3,8,9))
    # 5 adjacent faces
    triangles.append((4,9,5))
    triangles.append((2,4,11))
    triangles.append((6,2,10))
    triangles.append((8,6,7))
    triangles.append((9,8,1))

    def esplit(lk,u,v):
        k = (u,v)
        if not k in lk:
            pts.append(dpv.midpoint(pts[u],pts[v]))
            lk[k] = len(pts)-1
        return lk[k]
    def tsplit(lk,u,v,w):
        m1 = esplit(lk,u,v)
        m2 = esplit(lk,v,w)
        m3 = esplit(lk,w,u)
        return (m1,m3,u),(m2,m1,v),(m3,m2,w),(m1,m2,m3)
    def split(itris):
        otris = []
        nwpts = {}
        for t in itris:otris.extend(tsplit(nwpts,*t))
        return otris
    for nx in range(n):triangles = split(triangles)
    for p in pts:p.normalize().scale_u(r)

    polygons = tuple((tuple(pts[x].copy() for x in t),()) for t in triangles)
    plc.add_polygons(*polygons)

    #dprf.profile_function(plc.triangulate)
    plc.simplices = [(pts[u],pts[v],pts[w]) for u,v,w in triangles]
    plc.ghostbnds = []
    ax = plc.plot()
    plt.show()
    return plc

# generate a plc for a triangulated cube
def box(l = 1,w = 1,h = 2):
    import dilap.mesh.piecewisecomplex as pwc
    plc = pwc.piecewise_linear_complex()

    pts = []
    pts.append(dpv.vector( 0, 0, 0))
    pts.append(dpv.vector( l, 0, 0))
    pts.append(dpv.vector( l, w, 0))
    pts.append(dpv.vector( 0, w, 0))
    pts.append(dpv.vector( 0, 0, h))
    pts.append(dpv.vector( l, 0, h))
    pts.append(dpv.vector( l, w, h))
    pts.append(dpv.vector( 0, w, h))
    dpv.translate_coords(pts,dpv.vector(-0.5*l,-0.5*w,0.0))

    polygons = ((3,2,1,0),(4,5,6,7),(0,1,5,4),(1,2,6,5),(2,3,7,6),(3,0,4,7))
    polygons = tuple((tuple(pts[x].copy() for x in p),()) for p in polygons)

    plc.add_polygons(*polygons)
    #dprf.profile_function(plc.triangulate)
    #ax = plc.plot()
    #plt.show()
    return plc

def roof():
    import dilap.mesh.piecewisecomplex as pwc
    plc = pwc.piecewise_linear_complex()

    h = 5
    ebnd = (
        dpv.vector(0,0,0),dpv.vector(10,0,0),
        dpv.vector(10,10,0),dpv.vector(20,10,0),
        dpv.vector(20,20,0),dpv.vector(0,20,0))
    ibnd = (
        dpv.vector(2,2,h),dpv.vector(8,2,h),
        dpv.vector(8,12,h),dpv.vector(18,12,h),
        dpv.vector(18,18,h),dpv.vector(2,18,h))
    polygons = ((ebnd,(ibnd,)),)

    plc.add_polygons(*polygons)
    dprf.profile_function(plc.triangulate)
    ax = plc.plot()
    plt.show()
    return plc

def facade():
    import dilap.mesh.piecewisecomplex as pwc
    rim = [
        dpv.vector(0,0,0),dpv.vector(15,0,0),
        dpv.vector(15,0,4),dpv.vector(0,0,4)]

    door = [
        dpv.vector(-1,0,0.2),dpv.vector(1,0,0.2),
        dpv.vector(1,0,3),dpv.vector(-1,0,3)]
    dpv.translate_coords(door,dpv.vector(10,0,0))

    wspline = dpv.vector_spline(
        dpv.vector(2,0,3),dpv.vector(1.8,0,3.2),
        dpv.vector(-1.8,0,3.2),dpv.vector(-2,0,3),10)
    window = [
        dpv.vector(-2,0,0.5),dpv.vector(2,0,0.5),
        dpv.vector(2,0,3)]+wspline+[dpv.vector(-2,0,3)]
    dpv.translate_coords(window,dpv.vector(5,0,0))

    beam = [
        dpv.vector(1,0,1),dpv.vector(2,0,1),
        dpv.vector(2,0,9),dpv.vector(1,0,9)]

    #fac = (rim,(door,window,beam)),(beam,())
    fac = (rim,(door,window,beam)),
    
    plc = pwc.piecewise_linear_complex()
    plc.add_polygons(*fac)

    #plc.translate_polygon(2,dpv.vector(0,0,10))
    #plc.extrude_polygon(1,dpv.vector(0,1,0))

    plc.triangulate()

    #ax = plc.plot_xy()
    #ax = plc.plot()

    #plt.show()
    return plc

# profile the test function
def profile_triangulation():
    import dilap.mesh.piecewisecomplex as pwc

    ps1 = tuple(dpr.point_ring(100,32))
    ps2 = tuple(dpr.point_ring(20,8))
    ps3 = tuple(dpr.point_ring(100,64))
    ps4 = tuple([p.copy() for p in ps2])
    q = dpq.q_from_av(numpy.pi/2.0,dpv.xhat)
    for p in ps1:p.rotate(q)
    for p in ps2:p.rotate(q)
    dpv.translate_coords(list(ps1),dpv.vector(0,100,0))
    dpv.translate_coords(list(ps2),dpv.vector(0,100,0))
    dpv.translate_coords(list(ps3),dpv.vector(0,0,-100))
    dpv.translate_coords(list(ps4),dpv.vector(0,0,-100))

    polygons = ((ps1,(ps2,)),(ps3,()))
    #polygons = ((ps3,(ps4,)),)
    #polygons = ((ps3,()),)

    plc = pwc.piecewise_linear_complex()
    plc.add_polygons(*polygons)
    
    dprf.profile_function(plc.triangulate)
    #dprf.profile_function(plc.triangulate_xy)

    ax = plc.plot()
    plt.show()

###############################################################################

# does the line segment interior given by s1,s2 intersect the polygon py
def intersect_segment_polygon(s1,s2,py):
    if dpr.inconcave_xy(s1,py) and dpr.inconcave_xy(s2,py):return True
    for ex in range(len(py)):
        ep1,ep2 = py[ex-1],py[ex]
        isect = segments_intersect_at(s1,s2,ep1,ep2)
        if not isect is None:
            if type(isect) == type(()):return True
            #    if dpv.distance(*isect) > 0.001:return True
    return False

# given a line segment and a polygon, produce
# all points of intersection for edge segmentation
def segsects(s1,s2,py):
    unn = False
    isects = []
    for pyx in range(len(py)):
        ep1,ep2 = py[pyx-1],py[pyx]
        isect = segments_intersect_at(s1,s2,ep1,ep2)
        if not isect is None:
            if type(isect) == type(()):
                isects.append(isect[0])
                isects.append(isect[1])
                unn = True
            elif isect.near(s1) or isect.near(s2):pass
            elif isect.near(ep1) or isect.near(ep2):isects.append(isect)
            else:
                isects.append(isect)
                unn = True
    return isects,unn

# given polygon 1 and 2, return segmented edges of polygon 1
# segmented based on intersection with py2
def segpoly(py1,py2):
    unfinished,unionize = [],False
    for e1x in range(len(py1)):
        ep11,ep12 = py1[e1x-1],py1[e1x]
        isects,union = segsects(ep11,ep12,py2)
        if union:unionize = True
        isects = segclean(ep11,ep12,isects)
        for x in range(1,len(isects)):
            i1,i2 = (isects[x-1],isects[x])
            #if dpr.inconcave_xy(i1,py2) and dpr.inconcave_xy(i2,py2):pass
            #elif not intersect_segment_polygon(i1,i2,py2):
            #if not intersect_segment_polygon(i1,i2,py2):
            if True:
                unfinished.append((i1,i2))
        #else:unfinished.append((ep11,ep12))
    return unfinished,unionize

# return an ordered nonduplicative list of points from s1 to s2
def segclean(s1,s2,isects):
    if   len(isects) == 0:return [s1,s2]
    elif len(isects) == 1:return [s1,isects[0],s2]
    clean = [s1]
    iordr = list(dpv.proximity_order(s1,isects))
    while iordr:
        nxtx = iordr.pop(0)
        nxti = isects[nxtx]
        if not clean[-1].near(nxti):clean.append(nxti)
    if not clean[-1].near(s2):clean.append(s2)
    return clean

# given two concave polygons with holes
# return one polygon if p1 and p2 can be merged
# or return None
#def merge_two_polygons(p1,p2):
def union_polygon(p1,p2):
    eb1,eb2 = p1[0],p2[0]
    ebn1 = dpr.polygon_normal(eb1)
    pj1 = dpv.project_coords(list(eb1),ebn1)
    pj2 = dpv.project_coords(list(eb2),ebn1)
    nr = dpr.isnear
    if not (nr(pj1.x,pj1.y) and nr(pj2.x,pj2.y) and nr(pj1.x,pj2.x)):return
    if       ebn1.near(dpv.nz()):prot = dpq.q_from_av(dpr.PI,dpv.x())
    elif not ebn1.near(dpv.z() ):prot = dpq.q_from_uu(ebn1,dpv.z())
    else:                        prot = dpq.zero()
    dpr.rotate_polygon(p1,prot)
    dpr.rotate_polygon(p2,prot)
    unionize = False
    unfinished = []
    eb1seg,unionize1 = segpoly(eb1,eb2)
    eb2seg,unionize2 = segpoly(eb2,eb1)
    unfinished.extend(eb1seg)
    unfinished.extend(eb2seg)


    ax = plot_axes_xy()
    #plot_polygon_full_xy(p1,ax)
    #plot_polygon_full_xy(p2,ax)
    for u in unfinished:plot_edges_xy(u,ax,lw = 4.0)
    plt.show()


    unionize = unionize1 or unionize2
    if not unionize:
        prot.flip()
        dpr.rotate_polygon(p1,prot)
        dpr.rotate_polygon(p2,prot)
        return
    finished = []
    last = unfinished[0][0]
    while unfinished:
        for ufx in range(len(unfinished)):
            unfin = unfinished[ufx]
            if   last.near(unfin[0]):break
            elif last.near(unfin[1]):break
        unfin = unfinished.pop(ufx)
        which = unfin[1] if last.near(unfin[0]) else unfin[0]
        finished.append(which)
        last = which
    merged = (tuple(finished),tuple(p1[1]+p2[1]))
    return dpr.rotate_polygon(merged,prot.flip())

# given a collection of concave polygons with holes
# return a new set of polygons where holes are maintained
# but exerior bounds which intersect in a tractable way are
# merged into one polygon
def merge_polygons(polys):
    islands = [p for p in polys]
    checked = False
    while not checked:
        icnt = len(islands)
        if icnt == 1:checked = True
        else:
            for x in range(icnt):
                i1 = islands[x]
                for y in range(icnt):
                    if x == y:continue
                    i2 = islands[y]
                    #mp = merge_two_polygons(i1,i2)
                    mp = polygon_union(i1,i2)
                    if not mp is None:
                        islands.remove(i1)
                        islands.remove(i2)
                        islands.append(mp)
                        break
                if not mp is None:break
                if x == icnt-1:checked = True
    return islands

# i need the ability to do csg with polygons in the xy plane
# all 3 operations begin the same
# verify the two polygons are in the same plane
#   if they are not, then cgs between them is not valid
# find all points of intersection between the two polygons
#   should use sweepline algorithm?
# segment all segments which contain intersections at each intersection
#   just detect if a point is on the interior of any segment and return two new segments
# determine if each new segment should be in the final polygon
#   this is the more complicated step thats specific to the operation
# use nearness of segments to decide new connectivity of resulting polygon
#   this is where degeneracy/ambiguity may affect things

# given two polygons, determine if they occupy the same plane
# if they do, return the quaternion to rotate them to the xy plane
# if they do not, return None
def valid_pair(py1,py2):
    eb1,eb2 = py1[0],py2[0]
    ebn1 = dpr.polygon_normal(eb1)
    pj1 = dpv.project_coords(list(eb1),ebn1)
    pj2 = dpv.project_coords(list(eb2),ebn1)
    if not (dpr.isnear(pj2.x,pj2.y) and dpr.isnear(pj1.x,pj2.x)):return
    if       ebn1.near(dpv.nz()):prot = dpq.q_from_av(dpr.PI,dpv.x())
    elif not ebn1.near(dpv.z() ):prot = dpq.q_from_uu(ebn1,dpv.z())
    else:                        prot = dpq.zero()
    return prot

# given two polygons, return 0 if they are disjoint
# return 1 if py1 is inside py2, or -1 if py2 is inside py1
def containment(py1,py2):
    raise NotImplemented

# given a set of edge segments, 
# return all points of interior intersection
# NOTE: this is brute-force now.... consider sweepline in the future
def intersections(segments):
    def newpoint(ip):
        for p in ipts:
            if p.near(ip):return
        ipts.append(ip)
    ipts = []
    scnt = len(segments)
    for x in range(scnt):
        x1,x2 = segments[x]
        for y in range(scnt):
            if x == y:continue
            y1,y2 = segments[y]
            isect = segments_intersect_at(x1,x2,y1,y2)
            if isect is None:continue
            elif type(isect) == type(()):
                i1,i2 = isect
                newpoint(i1)
                newpoint(i2)
            else:newpoint(isect)
    return ipts

# given a set of edge segments and a set of points, return a new set
# of edge segments where all edges whose interiors intersect a point 
# in points is broken at the point of intersection
def break_segments(segments,points):
    broken = []
    for seg in segments:
        s1,s2 = seg
        relev = []
        for ipt in points:
            if dpr.insegment(ipt,s1,s2):
                relev.append(ipt)
        clean = segclean(s1,s2,relev)
        for x in range(1,len(clean)):
            i1,i2 = (clean[x-1],clean[x])
            broken.append((i1,i2))
    return broken

# given a set of edge segments, 
# return a boundary represented by the segments
def construct_loop(segments):
    srange = [x for x in range(len(segments))]
    loop = []
    last = segments[0][0]
    while srange:
        for x in srange:
            s1,s2 = segments[x]
            if   last.near(s1):found = s2;break
            elif last.near(s2):found = s1;break
        srange.remove(x)
        loop.append(found)
        last = found
    return tuple(loop)

# given a concave boundary, return the set of edge segments
# which generate the boundary
def boundary_segments(bnd):
    segments = []
    for x in range(len(bnd)):segments.append((bnd[x-1],bnd[x]))
    return segments

# given a concave polygon with holes, return the set of 
# edge segments which generate the polygon
def polygon_segments(py):
    eb,ibs = py
    segments = []
    segments.extend(boundary_segments(eb))
    #for ib in ibs:
    #    print('considering polygon segments from hole!!!')
    #    segments.extend(boundary_segments(ib))
    return segments

# given a set of segments, 
# return the subset which intersects a polygon
def segments_inpolygon(segments,py):
    eb,ibs = py
    inpoly = []
    for x in range(len(segments)):
        i1,i2 = segments[x]
        if intersect_segment_polygon(i1,i2,eb):
            inpoly.append((i1,i2))
    return inpoly

# given a set of segments, 
# return the subset which does not intersect a polygon
def segments_outpolygon(segments,py):
    eb,ibs = py
    inpoly = []
    for x in range(len(segments)):
        i1,i2 = segments[x]
        if not intersect_segment_polygon(i1,i2,eb):
            inpoly.append((i1,i2))
    return inpoly

# given two polygons, return the union of their edge segments
# broken upon every point of intersection between polygons
# or return None if there are no intersections
def break_polygons(p1,p2):
    p1segs = polygon_segments(p1)
    p2segs = polygon_segments(p2)
    isects = intersections(p1segs+p2segs)
    if len(isects) == 0:return p1segs,p2segs
    p1segs = break_segments(p1segs,isects)
    p2segs = break_segments(p2segs,isects)
    return p1segs,p2segs

def polygon_union(p1,p2):
    prot = valid_pair(p1,p2)
    if prot is None:print('invalid polygon csg operation...');return
    dpr.rotate_polygon(p1,prot)
    dpr.rotate_polygon(p2,prot)
    broken = break_polygons(p1,p2)
    p1segs,p2segs = break_polygons(p1,p2)
    p1inp2 = segments_outpolygon(p1segs,p2)
    p2inp1 = segments_outpolygon(p2segs,p1)
    #union = (construct_loop(p1inp2+p2inp1),())
    union = (construct_loop(p1inp2+p2inp1),(p1[1]+p2[1]))
    dpr.rotate_polygon(union,prot.flip())
    return union

def polygon_intersection(p1,p2):
    prot = valid_pair(p1,p2)
    if prot is None:print('invalid polygon csg operation...');return
    dpr.rotate_polygon(p1,prot)
    dpr.rotate_polygon(p2,prot)
    p1segs,p2segs = break_polygons(p1,p2)
    p1inp2 = segments_inpolygon(p1segs,p2)
    p2inp1 = segments_inpolygon(p2segs,p1)
    inter = (construct_loop(p1inp2+p2inp1),())
    dpr.rotate_polygon(inter,prot.flip())
    return inter

def polygon_difference(p1,p2):
    prot = valid_pair(p1,p2)
    if prot is None:print('invalid polygon csg operation...');return
    dpr.rotate_polygon(p1,prot)
    dpr.rotate_polygon(p2,prot)
    p1segs,p2segs = break_polygons(p1,p2)
    p1inp2 = segments_outpolygon(p1segs,p2)
    p2inp1 = segments_inpolygon(p2segs,p1)
    diffr = (construct_loop(p1inp2+p2inp1),())
    dpr.rotate_polygon(diffr,prot.flip())
    return diffr

def cgstest():
    p1 = (tuple(dpr.square(5,5)),())
    p2 = (tuple(dpr.square(5,4,dpv.vector(3,0,0))),())

    #union = polygon_union(p1,p2)
    #inter = polygon_intersection(p1,p2)
    diffr = polygon_difference(p1,p2)




    ax = plot_axes_xy()
    #for seg in p1segs:ax = plot_edges_xy(seg,ax)
    #for seg in p2segs:ax = plot_edges_xy(seg,ax)
    plot_polygon_full_xy(p1,ax)
    plot_polygon_full_xy(p2,ax)
    #plot_polygon_full_xy(union,ax,lw = 4.0)
    #plot_polygon_full_xy(inter,ax,lw = 4.0)
    plot_polygon_full_xy(diffr,ax,lw = 4.0)
    plt.show()

    '''#
    prot = valid_pair(p1,p2)
    if prot is None:print('invalid polygon csg operation...');return
    dpr.rotate_polygon(p1,prot)
    dpr.rotate_polygon(p2,prot)
  
    p1segs,p2segs = break_polygons(p1,p2)

    # union
    p1inp2 = segments_outpolygon(p1segs,p2)
    p2inp1 = segments_outpolygon(p2segs,p1)
    union = construct_loop(p1inp2+p2inp1)

    # intersect
    p1inp2 = segments_inpolygon(p1segs,p2)
    p2inp1 = segments_inpolygon(p2segs,p1)
    inter = construct_loop(p1inp2+p2inp1)

    # difference
    p1inp2 = segments_outpolygon(p1segs,p2)
    p2inp1 = segments_inpolygon(p2segs,p1)
    diffr = construct_loop(p1inp2+p2inp1)

    ax = plot_axes_xy()
    for seg in p1segs:ax = plot_edges_xy(seg,ax)
    for seg in p2segs:ax = plot_edges_xy(seg,ax)
    #plot_polygon_full_xy(p1,ax)
    #plot_polygon_full_xy(p2,ax)

    #plot_polygon_xy(union,ax,lw = 4.0)
    #plot_polygon_xy(inter,ax,lw = 4.0)
    plot_polygon_xy(diffr,ax,lw = 4.0)
    plt.show()

    pdb.set_trace()
    '''#





