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

def plot_point(pt,ax,marker = 'o'):
    ax.plot([pt.x],[pt.y],zs = [pt.z],marker = marker)
    return ax

def plot_points_xy(points,ax = None,ms = None):
    if ax is None:ax = plot_axes_xy()
    if ms is None:ms = ['o']*len(points)
    for pdx in range(len(points)):plot_point_xy(points[pdx],ax,ms[pdx])
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

'''#
# if a is within c of b, return True
# else return False
def isnear(a,b,c = 0.000001):
    if abs(a-b) < c:return True
    else:return False

# if a is within c of b, return b
# else return a
def near(a,b,c = 0.000001):
    if abs(a-b) < c:return b
    else:return a

def orient2d(a,b,c):
    m11,m12 = a.x-c.x,a.y-c.y
    m21,m22 = b.x-c.x,b.y-c.y
    det = m11*m22-m12*m21
    return near(det,0)

# return the signed volume of the parallelpiped created
# by the vectors a-d,b-d,c-d
# return 0 if a,b,c,d are coplanar
def orient3d(a,b,c,d):
    m11 = a.x-d.x;m12 = a.y-d.y;m13 = a.z-d.z
    m21 = b.x-d.x;m22 = b.y-d.y;m23 = b.z-d.z
    m31 = c.x-d.x;m32 = c.y-d.y;m33 = c.z-d.z
    det = m11*(m22*m33-m23*m32)-m12*(m21*m33-m23*m31)+m13*(m21*m32-m22*m31)
    return near(det,0)

def incircle(a,b,c,d):
    m11,m12 = a.x-d.x,a.y-d.y
    m13 = m11*m11 + m12*m12
    m21,m22 = b.x-d.x,b.y-d.y
    m23 = m21*m21 + m22*m22
    m31,m32 = c.x-d.x,c.y-d.y
    m33 = m31*m31 + m32*m32
    det1 = m11*(m22*m33-m23*m32)
    det2 = m12*(m21*m33-m23*m31)
    det3 = m13*(m21*m32-m22*m31)
    incirc = near(det1 - det2 + det3,0)
    #return incirc*orient2d(a,b,c)
    return incirc

# let a,b,c,d be such that orient3d(a,b,c,d) is nonnegative
# return > 0 if e is inside sphere passing through a,b,c,d
# return < 0 if e is outside sphere passing through a,b,c,d
# return 0 if e is on the surface of the sphere passing through a,b,c,d
# return 0 if all five points are coplanar
def insphere(a,b,c,d,e):
    m11 = a.x-e.x;m12 = a.y-e.y;m13 = a.z-e.z
    m14 = m11*m11 + m12*m12 + m13*m13
    m21 = b.x-e.x;m22 = b.y-e.y;m23 = b.z-e.z
    m24 = m21*m21 + m22*m22 + m23*m23
    m31 = c.x-e.x;m32 = c.y-e.y;m33 = c.z-e.z
    m34 = m31*m31 + m32*m32 + m33*m33
    m41 = d.x-e.x;m42 = d.y-e.y;m43 = d.z-e.z
    m44 = m41*m41 + m42*m42 + m43*m43
    det1 = m11*(m22*(m33*m44-m34*m43)-m23*(m32*m44-m34*m42)+m24*(m32*m43-m33*m42))
    det2 = m12*(m21*(m33*m44-m34*m43)-m23*(m31*m44-m34*m41)+m24*(m31*m43-m33*m41))
    det3 = m13*(m21*(m32*m44-m34*m42)-m22*(m31*m44-m34*m41)+m24*(m31*m42-m32*m41))
    det4 = m14*(m21*(m32*m43-m33*m42)-m22*(m31*m43-m33*m41)+m23*(m31*m42-m32*m41))
    insphr = near(det1 - det2 + det3 - det4,0)
    #return insphr*orient3d(a,b,c,d)
    return insphr

def inconcave(pt,poly):
    angle = 0.0
    for x in range(len(poly)):
        p1 = poly[x-1]
        p2 = poly[x]
        e1 = p1 - pt
        e2 = p2 - pt
        angle += dpv.signed_angle_between_xy(e1,e2)
    if abs(angle) < numpy.pi:return False
    else:return True
'''#

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


'''#
/*
 *    Return the angle between two vectors on a plane
 *       The angle is from vector 1 to vector 2, positive anticlockwise
 *          The result is between -pi -> pi
 *          */
double Angle2D(double x1, double y1, double x2, double y2)
{
    double dtheta,theta1,theta2;

    theta1 = atan2(y1,x1);
    theta2 = atan2(y2,x2);
    dtheta = theta2 - theta1;
    while (dtheta > PI)
        dtheta -= TWOPI;
    while (dtheta < -PI)
        dtheta += TWOPI;

    return(dtheta);
}
'''#

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

'''#
# given line segment s1, line segment s2
# does s1 overlap the interior or s2?
# a segment is a tuple of two points
def segments_intersect(s1,s2,err = 0.001):
    pi2 = numpy.pi/2.0
    n1 = dpv.v1_v2(*s1).rotate_z(pi2)
    proj1 = dpv.project_coords([s1[0],s1[1]],n1)
    proj2 = dpv.project_coords([s2[0],s2[1]],n1)
    if proj1.x - proj2.x > err and proj2.y - proj1.x > err:
        n2 = dpv.v1_v2(*s2).rotate_z(pi2)
        proj1 = dpv.project_coords([s1[0],s1[1]],n2)
        proj2 = dpv.project_coords([s2[0],s2[1]],n2)
        if proj2.x - proj1.x > err and proj1.y - proj2.x > err:
            #ax = plot_edges_xy(s1)
            #plot_edges_xy(s2,ax)
            #plt.show()
            return 1
    #if s1[0].near_xy(s2[0]) and s1[1].near_xy(s2[1]):return 1
    #elif s1[0].near_xy(s2[1]) and s1[1].near_xy(s2[0]):return 1
    return 0

# given concave polygon p1, concave polygon p2
# does p1 overlap the interior p2?
# a polygon is a tuple of points
def concaves_intersect(p1,p2):
    isegsectfound = False
    for px in range(len(p1)):
        if isegsectfound:break
        e1 = (p1[px-1],p1[px])
        for py in range(len(p2)):
            e2 = (p2[py-1],p2[py])
            if dpr.segments_intersect(e1,e2):
                isegsectfound = True
                break
    i2 = dpv.center_of_mass(p2)
    if inconcave(i2,p1) or isegsectfound:ins = True
    else:ins = False
    return ins

# given concave polygon p1, concave polygon p2
# does p1 overlap the interior p2?
# a polygon is a tuple of points
def concaves_contains(p1,p2):
    isegsectfound = False
    for px in range(len(p1)):
        if isegsectfound:break
        #e1 = (p1[px-1],p1[px])
        for py in range(len(p2)):
            #e2 = (p2[py-1],p2[py])
            #if dpr.segments_intersect(e1,e2):
            if dpr.segments_intersect(p1[px-1],p1[px],p2[py-1],p2[py]):
                isegsectfound = True
                break
    i2 = dpv.center_of_mass(p2)
    if inconcave(i2,p1) and not isegsectfound:ins = True
    else:ins = False
    return ins
'''#

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

    '''#
    for nx in range(1):
        for egx in range(plc.edgecount):
            eg = plc.edges[egx]
            if eg is None:continue
            plc.split_edge(*eg)
    '''#

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






