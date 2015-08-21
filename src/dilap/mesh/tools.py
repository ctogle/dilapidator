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

def plot_edges_xy(points,ax = None,mk = None,lw = 1.0):
    if ax is None:ax = plot_axes_xy()
    if mk is None:mk = '+'
    pts = [p.to_tuple() for p in points]
    xs,ys,zs = zip(*pts)
    ax.plot(xs,ys,marker = mk,lw = lw)
    return ax

def plot_edges(points,ax = None,lw = 1.0):
    if ax is None:ax = plot_axes()
    pts = [p.to_tuple() for p in points]
    xs,ys,zs = zip(*pts)
    ax.plot(xs,ys,zs,marker = '+',lw = lw)
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
# does s1 overlap the interior or s2?
# a segment is a tuple of two points
# return the point of intersection if there is one
# otherwise return None
#   note: currently assumes s1,s2 in xy-plane
def segments_intersect_at(s1,s2):
    p,q = s1[0],s2[0]
    r = dpv.v1_v2(*s1)
    s = dpv.v1_v2(*s2)
    qmp = q-p
    rcs = r.cross(s)
    rcsmag = rcs.magnitude()
    qmpcr = qmp.cross(r)
    qmpcrmag = qmpcr.magnitude()
    rmag2 = r.magnitude2()

    if isnear(rcsmag,0) and isnear(qmpcrmag,0):
        t0 = qmp.dot(r)/rmag2
        t1 = t0 + s.dot(r)/rmag2
        if s.dot(r) < 0.0:t0,t1 = t1,t0
        if dpb.overlap(dpv.vector2d(0,1),dpv.vector2d(t0,t1)):
            if t0 == 1 or t1 == 0:return None
            t0pt = p + r.copy().scale_u(t0)
            t1pt = p + r.copy().scale_u(t1)
            print('colinear and overlapping!')
            return t0pt,t1pt
        else:return None
    elif isnear(rcsmag,0) and not isnear(qmpcrmag,0):return None
    elif not isnear(rcsmag,0) and not isnear(qmpcrmag,0):
        u = qmpcr.z/rcs.z
        t = qmp.cross(s).z/rcs.z
        if (u == 0 or u == 1) and (t == 0 or t == 1):return None
        if dpb.p_in_rng(u,0,1) and dpb.p_in_rng(t,0,1):return q + s.scale_u(u)

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
    dprf.profile_function(plc.triangulate)
    ax = plc.plot()
    plt.show()
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






