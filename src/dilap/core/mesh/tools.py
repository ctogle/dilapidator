import dilap.core.tools as dpr

import dp_vector as dpv

import matplotlib.pyplot as plt
import numpy

import pdb

def plot_point_xy(pt,ax,marker = 'o'):
    ax.plot([pt.x],[pt.y],marker = marker)

def plot_point(pt,ax,marker = 'o'):
    ax.plot([pt.x],[pt.y],zs = [pt.z],marker = marker)

def plot_points_xy(points,ax = None,ms = None):
    if ax is None:ax = plt.figure().add_subplot(111)
    if ms is None:ms = ['o']*len(points)
    for pdx in range(len(points)):plot_point_xy(points[pdx],ax,ms[pdx])
    return ax

def plot_points(points,ax = None,ms = None):
    if ax is None:ax = plt.figure().add_subplot(111,projection = '3d')
    if ms is None:ms = ['o']*len(points)
    for pdx in range(len(points)):plot_point(points[pdx],ax,ms[pdx])  
    return ax

def plot_edges_xy(points,ax = None):
    if ax is None:ax = plt.figure().add_subplot(111)
    pts = [p.to_tuple() for p in points]
    xs,ys,zs = zip(*pts)
    ax.plot(xs,ys,marker = '+')
    return ax

def plot_edges(points,ax = None):
    if ax is None:ax = plt.figure().add_subplot(111,projection = '3d')
    pts = [p.to_tuple() for p in points]
    xs,ys,zs = zip(*pts)
    ax.plot(xs,ys,zs,marker = '+')
    return ax

def plot_polygon_xy(points,ax = None,center = False):
    epts = points[:]
    epts.append(points[0])
    ax = plot_edges_xy(epts,ax)
    if center:plot_point_xy(dpv.center_of_mass(points),ax,marker = 's')
    return ax

def plot_polygon(points,ax = None,center = False):
    epts = points[:]
    epts.append(points[0])
    ax = plot_edges(epts,ax)
    if center:plot_point(dpv.center_of_mass(points),ax,marker = 's')
    return ax

def plot_tetrahedron(points,ax = None):
    raise NotImplemented

def plot_circle_xy(c,r,ax = None,center = False):
    circ = dpv.translate_coords(dpr.point_ring(r,32),c)
    plot_polygon_xy(circ,ax,center)

def plot_circle(c,r,ax = None,center = False):
    circ = dpv.translate_coords(dpr.point_ring(r,32),c)
    plot_polygon(circ,ax,center)

###############################################################################

# if a is within c of b, return b
# else return a
def near(a,b,c = 0.0000000000001):
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

    print('inconcave call!!!!')
    ax = plot_polygon_xy(poly)
    plot_point_xy(pt,ax)
    plt.show()

    pdb.set_trace()

'''#
    int InsidePolygon(Point *polygon,int n,Point p)
    {
           int i;
              double angle=0;
                 Point p1,p2;

                    for (i=0;i<n;i++) {
                            p1.h = polygon[i].h - p.h;
                                  p1.v = polygon[i].v - p.v;
                                        p2.h = polygon[(i+1)%n].h - p.h;
                                              p2.v = polygon[(i+1)%n].v - p.v;
                                                    angle += Angle2D(p1.h,p1.v,p2.h,p2.v);
                                                       }

                       if (ABS(angle) < PI)
                             return(FALSE);
                                else
                                      return(TRUE);
                                      }

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
def segments_intersect(s1,s2):
    pi2 = numpy.pi/2.0
    n1 = dpv.v1_v2(*s1).rotate_z(pi2)
    proj1 = dpv.project_coords([s1[0],s1[1]],n1)
    proj2 = dpv.project_coords([s2[0],s2[1]],n1)
    if proj1.x > proj2.x and proj1.x < proj2.y:
        n2 = dpv.v1_v2(*s2).rotate_z(pi2)
        proj1 = dpv.project_coords([s1[0],s1[1]],n2)
        proj2 = dpv.project_coords([s2[0],s2[1]],n2)
        if proj1.x < proj2.x and proj1.y > proj2.x:
            print('these intersect!!!!')
            ax = plot_edges_xy(s1)
            plot_edges_xy(s2,ax)
            plt.show()
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
            if segments_intersect(e1,e2):
                isegsectfound = True
                break
    i2 = dpv.center_of_mass(p2)
    if inconcave(i2,p1) or isegsectfound:ins = True
    else:ins = False
    print('INZZZZ',ins)


    ax = plot_polygon_xy(p1)
    plot_polygon_xy(p2,ax)
    plt.show()


    pdb.set_trace()










