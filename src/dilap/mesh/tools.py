import dilap.core.tools as dpr
import dilap.core.vector as dpv
import dilap.core.quaternion as dpq
import dilap.core.bbox as dpb
import dilap.core.profiler as dprf

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

def plot_points(points,ax = None,ms = None,marker = None):
    if ax is None:ax = plot_axes()
    if marker is None:marker = 'o'
    if ms is None:ms = [marker]*len(points)
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
    epts = list(points[:])
    epts.append(points[0])
    ax = plot_edges_xy(epts,ax,lw = lw)
    if center:plot_point_xy(dpv.center_of_mass(points),ax,marker = 's')
    return ax

def plot_polygon(points,ax = None,center = False,lw = 1.0):
    epts = list(points[:])
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

def plot_line_xy(l1,l2,r = 25,ax = None,center = False,lw = 1.0):
    ltan = dpv.v1_v2_xy(l1,l2)
    l1far = l1.copy().translate(ltan.copy().scale_u(-r))
    l2far = l2.copy().translate(ltan.copy().scale_u( r))
    ax = plot_edges_xy([l1far,l2far],ax = ax,center = center,lw = lw)
    return ax

def plot_line(l1,l2,r = 25,ax = None,center = False,lw = 1.0):
    ltan = dpv.v1_v2(l1,l2)
    l1far = l1.copy().translate(ltan.copy().scale_u(-r))
    l2far = l2.copy().translate(ltan.copy().scale_u( r))
    ax = plot_edges([l1far,l2far],ax = ax,center = center,lw = lw)
    return ax

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

# given a polygon, find a quaternion rotating it to the xy plane
def prot_to_xy(py):
    eb,ibs = py
    ebn = dpr.polygon_normal(eb)
    if       ebn.near(dpv.nz()):prot = dpq.q_from_av(dpr.PI,dpv.x())
    elif not ebn.near(dpv.z() ):prot = dpq.q_from_uu(ebn,dpv.z())
    else:                       prot = dpq.zero()
    return prot

# given two points on a line and a polygon, 
# return the set of intersections points between the two
#def line_intersects_polygon_at(l1,l2,py):
def line_intersects_polygon_at(l1,l2,py,intins):
    eb,ibs = py

    ltn = dpv.v1_v2(l1,l2).normalize()
    pyprj = dpv.project_coords(list(eb),ltn)
    midpt = dpv.midpoint(l1,l2)
    mpprj = dpv.dot(midpt,ltn)
    offst = mpprj-((pyprj.x+pyprj.y)/2.0)
    l1del = ltn.copy().scale_u(offst+100)
    l2del = ltn.copy().scale_u(offst-100)
    l1far = midpt.copy().translate(l1del)
    l2far = midpt.copy().translate(l2del)
    pysegs = polygon_segments(py)
    lnsegs = [(l1far,l2far)]

    prot = prot_to_xy(py)

    '''#
    print('before you rotate!!',prot.__str__())
    ax = plot_axes()
    for p in pysegs:plot_edges(p,ax)
    plot_edges(lnsegs[-1],ax,lw = 5.0)
    #for p in isects:plot_point(p,ax)
    #for p in isects2:plot_point(p,ax,marker = 's')
    plt.show()
    print('AMEN')
    '''#

    dpr.rotate_segments(pysegs,prot)
    dpr.rotate_coords(list(lnsegs[0]),prot)
    isects = intersections(pysegs+lnsegs)
    #isects = [x for x in isects if dpr.orient3d(eb[0],eb[1],eb[2],x) == 0]

    if len(isects) == 1:
      print('while you were rotated!!')
      ax = plot_axes()
      for p in pysegs:plot_edges(p,ax)
      plot_edges(lnsegs[-1],ax,lw = 5.0)
      for p in isects:plot_point(p,ax)
      #for p in isects2:plot_point(p,ax,marker = 's')
      plt.show()
      print('AMEN.5',len(isects))

    if len(isects) == 0:
        prot.flip()
        dpr.rotate_segments(pysegs,prot)
        dpr.rotate_coords(list(lnsegs[0]),prot)
        dpv.rotate_coords(isects,prot)
        return
    else:
        #print('actuallyinsssssisects',isects)
        actuallyintins = dpr.inconcave_xy(dpv.midpoint(*isects),eb)

        '''#
        print('wattttch',actuallyintins)
        ax = plot_axes()
        plot_polygon(list(eb),ax)
        plot_edges(isects,ax,lw = 4.0)
        plt.show()
        '''#

        prot.flip()
        dpr.rotate_segments(pysegs,prot)
        dpr.rotate_coords(list(lnsegs[0]),prot)
        dpv.rotate_coords(isects,prot)
        if intins and not actuallyintins:return

        if len(isects) == 1:
          print('after you rotate!!')
          ax = plot_axes()
          for p in pysegs:plot_edges(p,ax)
          plot_edges(lnsegs[-1],ax,lw = 5.0)
          for p in isects:plot_point(p,ax)
          #for p in isects2:plot_point(p,ax,marker = 's')
          plt.show()
          print('AMEN2',len(isects))

        #ebn = dpr.polygon_normal(eb)
        #pyj = dpv.project_coords(list(eb),ebn)
        #isects2 = [x for x in isects if dpr.isnear(dpv.dot(x,ebn),pyj.x)]
        #isects2 = [x for x in isects if dpr.orient3d(eb[0],eb[1],eb[2],x) == 0]

        #if intins and not dpr.inconcave_xy(dpv.midpoint(*isects),eb):return
        #intins = dpr.inconcave_xy(dpv.midpoint(*isects),eb)
        #if not intins:return

        '''#
        print('breakin some shit!!',intins)
        ax = plot_axes()
        for p in pysegs:plot_edges(p,ax)
        plot_edges(lnsegs[-1],ax,lw = 5.0)
        for p in isects:plot_point(p,ax)
        #for p in isects2:plot_point(p,ax,marker = 's')
        plt.show()
        '''#

        return isects

# given line segment s1, line segment s2
# a segment is a tuple of two points
# return the point of intersection if there is one
# otherwise return None
#   note: currently assumes s1,s2 in xy-plane
def segments_intersect_at(s11,s12,s21,s22,include_endpoints = False):
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
        elif dpr.inrange(0,t0,t1) and dpr.inrange(1,t0,t1):
            return s11.copy(),s12.copy()
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
        if (u == 0 or u == 1) and (t == 0 or t == 1):
            #print('ENDPOINT')
            if include_endpoints:
                return q + s.scale_u(u)
            return
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
    #plc.simplices = [(pts[u],pts[v],pts[w]) for u,v,w in triangles]
    #plc.ghostbnds = []
    #ax = plc.plot()
    #plt.show()
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

# find the line of intersection between two planes
def planes_intersection(pn1,p01,pn2,p02):
    u = pn1.cross(pn2)
    ax = u.x if u.x >= 0 else -u.x
    ay = u.y if u.y >= 0 else -u.y
    az = u.z if u.z >= 0 else -u.z
    # pn1 and pn2 are near parallel
    if ((ax+ay+az) < 0.0001):
        v = p02 - p01
        if dpr.isnear(dpv.dot(pn1,v),0):return 
        else:return
    # pn1 and pn2 intersect in a line
    # first determine max abs coordinate of cross product
    if (ax > ay):
        if (ax > az):maxc =  1
        else:maxc = 3
    else:
        if (ay > az):maxc =  2
        else:maxc = 3

    # next, to get a point on the intersect line
    # zero the max coord, and solve for the other two
    d1 = -dpv.dot(pn1,p01)
    d2 = -dpv.dot(pn2,p02)
    # select max coordinate
    if maxc == 1: # intersect with x=0 
        y = (d2*pn1.z - d1*pn2.z) / u.x
        z = (d1*pn2.y - d2*pn1.y) / u.x
        ip = dpv.vector(0,y,z)
    elif maxc == 2: # intersect with y=0 
        x = (d1*pn2.z - d2*pn1.z) / u.y
        z = (d2*pn1.x - d1*pn2.x) / u.y
        ip = dpv.vector(x,0,z)
    elif maxc == 3: # intersect with z=0 
        x = (d2*pn1.y - d1*pn2.y) / u.z
        y = (d1*pn2.x - d2*pn1.x) / u.z
        ip = dpv.vector(x,y,0)
    return ip,ip + u

# given a polygon and a coplanar line segment, return 
# a new set of polygons which is split upon intersection with 
# the line or None if no intersection is found
def segment_split_polygon(s1,s2,py,plot = False):
    eb,ibs = py

    ltn = dpv.v1_v2(s1,s2).normalize()
    pyprj = dpv.project_coords(list(eb),ltn)
    midpt = dpv.midpoint(s1,s2)
    mpprj = dpv.dot(midpt,ltn)
    offst = mpprj-((pyprj.x+pyprj.y)/2.0)
    l1del = ltn.copy().scale_u(offst+100)
    l2del = ltn.copy().scale_u(offst-100)
    l1far = midpt.copy().translate(l1del)
    l2far = midpt.copy().translate(l2del)
    pysegs = polygon_segments(py)
    lnsegs = [(l1far,l2far)]
    #pysegs = polygon_segments(py)+[(s1,s2)]
        
    ebn = dpr.polygon_normal(eb)
    prot = dpr.q_to_xy(ebn)
    dpr.rotate_segments(pysegs,prot)
    dpr.rotate_coords(list(lnsegs[0]),prot)

    isects = intersections(pysegs+lnsegs)

    '''#
    print('truly time to line split!!!',len(isects))
    ax = plot_axes()
    for p in pysegs:plot_edges(p,ax)
    plot_edges([l1far,l2far],ax,lw = 4.0)
    for p in isects:plot_point(p,ax)
    plt.show()
    '''#

    if len(isects) == 0:
        prot.flip()
        dpr.rotate_segments(pysegs,prot)
        #dpr.rotate_coords(list(lnsegs[0]),prot)
        return

    pysegs = break_segments(pysegs,isects)
    lnsegs = break_segments(lnsegs,isects)

    # THIS MISSES TANGENTIAL INTERSECTION WITH THE LINE
    # THIS MISSES TANGENTIAL INTERSECTION WITH THE LINE
    # THIS MISSES TANGENTIAL INTERSECTION WITH THE LINE
    # THIS MISSES TANGENTIAL INTERSECTION WITH THE LINE
    # THIS MISSES TANGENTIAL INTERSECTION WITH THE LINE

    lnsegs  = segments_inpolygon(lnsegs,py)
    leftof  = segments_leftofline( pysegs,l1far,l2far)
    rightof = segments_rightofline(pysegs,l1far,l2far)

    if len(leftof) == 0 or len(rightof) == 0:
        #print('maybe??')
        prot.flip()
        dpr.rotate_segments(pysegs,prot)
        #dpr.rotate_coords(list(lnsegs[0]),prot)
        return

    if len(leftof) < 3 or len(rightof) < 3:
        print('HERE IT IS',len(pysegs),s1,s2)
        ax = plot_axes()
        #for p in pysegs:plot_edges(p,ax)
        #for s in leftof:plot_edges(s,ax,lw = 3.0)
        #for s in rightof:plot_edges(s,ax,lw = 6.0)
        for s in lnsegs:plot_edges(s,ax,lw = 6.0)
        #plot_edges([l1far,l2far],ax,lw = 4.0)
        for p in isects:plot_point(p,ax)
        plt.show()
        pdb.set_trace()

    leftloops  = construct_loops( leftof+lnsegs,plot)
    rightloops = construct_loops(rightof+lnsegs,plot)

    newpolys = []
    prot.flip()
    for seg in  leftloops:
        newpolys.append(dpr.rotate_polygon((seg,()),prot))
    for seg in rightloops:
        newpolys.append(dpr.rotate_polygon((seg,()),prot))
    #dpr.rotate_coords(list(lnsegs[0]),prot)

    '''#
    print('MADE NEW POLYSSSS',len(newpolys))
    ax = plot_axes()
    for ll in newpolys:plot_polygon_full(ll,ax)
    plt.show()
    '''#

    return newpolys

#######################################
#######################################


'''#
# does the line segment interior given by s1,s2 intersect the polygon py
### DOES NOT WORK .... and is unused!!!!!!!!!!!!
def intersect_segment_polygon(s1,s2,py):

    if dpr.inconcave_xy(dpv.midpoint(s1,s2),py):return True

    if dpr.inconcave_xy(s1,py) and dpr.inconcave_xy(s2,py):return True
    for ex in range(len(py)):
        ep1,ep2 = py[ex-1],py[ex]
        isect = segments_intersect_at(s1,s2,ep1,ep2)

        #if s1.near(ep1) and s2.near(ep2):
        #    print('here it is?')

        if not isect is None:
            if type(isect) == type(()):return True
            #    if dpv.distance(*isect) > 0.001:return True
    return False
'''#

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

                    '''#
                    print('INPUT')
                    ax = plot_axes(x = 40)
                    plot_polygon_full(i1,ax)
                    plot_polygon_full(i2,ax)
                    plt.show()
                    '''#

                    mp = polygon_union(i1,i2)

                    '''#
                    print('OUTPUT')
                    if not mp is None:
                        ax = plot_axes(x = 40)
                        plot_polygon_full(mp,ax)
                        plt.show()
                    '''#

                    '''#
                    for x in range(len(i2[0])):
                        one,two = i2[0][x-1],i2[0][x]
                        if one.near(two):
                            print('wtwfadfafdf:w',mp)
                            pdb.set_trace()
                    if not mp is None:
                        for x in range(len(mp[0])):
                            one,two = mp[0][x-1],mp[0][x]
                            if one.near(two):
                                print('mpmpmpmp',mp)
                                pdb.set_trace()
                    '''#

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

    for x in range(len(eb2)):
        one,two = eb2[x-1],eb2[x]
        if one.near(two):
            print('invalid!!!')
            pdb.set_trace()

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
    if dpr.inconcave_xy(py1[0],py2):return  1
    if dpr.inconcave_xy(py2[0],py1):return -1
    print('CONTAINMENT IS SKETCHY')
    return 0

    raise NotImplemented

def newpoint(p,ps):
    for q in ps:
        if p.near(q):return
    ps.append(p)

# given a set of edge segments, 
# return all points of interior intersection
# NOTE: this is brute-force now.... consider sweepline in the future
def intersections(segments):
    ipts = []
    scnt = len(segments)
    for x in range(scnt):
        x1,x2 = segments[x]
        for y in range(scnt):
            if x == y:continue
            y1,y2 = segments[y]
            #isect = segments_intersect_at(x1,x2,y1,y2,include_endpoints)
            #isect = segments_intersect_at(x1,x2,y1,y2,True)
            isect = segments_intersect_at(x1,x2,y1,y2)
            if isect is None:continue
            elif type(isect) == type(()):
                i1,i2 = isect
                newpoint(i1,ipts)
                newpoint(i2,ipts)
            else:newpoint(isect,ipts)
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
            if dpr.insegment_xy(ipt,s1,s2):
                relev.append(ipt)
        clean = segclean(s1,s2,relev)
        for x in range(1,len(clean)):
            i1,i2 = (clean[x-1],clean[x])
            broken.append((i1,i2))
    return broken

# remove points that do not form an angle of less than 180 degrees
# with the two edges associated with each point
def clean_loop(loop,plot = False):

    if plot and len(loop) < 3:
      print('ENTER LOOP',len(loop))
      ax = plot_axes()
      plot_polygon(loop,ax)
      plt.show()
      pdb.set_trace()

    def eang(x,y,z):
        e1 = dpv.v1_v2(loop[x],loop[y])
        e2 = dpv.v1_v2(loop[x],loop[z])
        a = dpr.angle_between(e1,e2)
        return a
    lcnt = len(loop)
    for x in range(lcnt):
        y = x-1
        z = x+1 if x < lcnt-2 else 0
        a = eang(x,y,z)
        print('x',x,a)
        if dpr.isnear(a,dpr.PI):loop[x] = None

        #if dpr.isnear(a,dpr.PI):marked.append(x-1)
    #for x in marked:loop.pop(x)
    #print('looptyloop',marked)
    for p in loop:print(p)

    loop = [l for l in loop if not l is None]

    if plot and len(loop) < 3:
      print('EXIT LOOP',len(loop))
      ax = plot_axes()
      plot_polygon(loop,ax)
      plt.show()

    return loop
    #clean = [loop[x-1] for x in range(lcnt) if not x-1 in marked]
    #return clean

# given a set of edge segments, 
# return a boundary represented by the segments
def construct_loops(segments,plot = False):
    srange = [x for x in range(len(segments))]
    loops = [[]]
    last = segments[0][0].copy()
    while srange:
        foundone = False
        for x in srange:
            s1,s2 = segments[x]
            if   last.near(s1):
                found = s2
                foundone = True
                break
            elif last.near(s2):
                found = s1
                foundone = True
                break
        #print('srange',srange)
        if foundone:
            #print('foundone',foundone)
            srange.remove(x)
            newpoint(found.copy(),loops[-1])
            last = found
            '''#
            ax = plot_axes()
            ax = plot_edges(segments[-1],ax,lw = 5.0)
            ax = plot_edges(loops[-1],ax)

            plt.show()
            '''#
        else:
            #print('foundone',foundone)
            last = segments[x][0].copy()
            loops.append([])

    if min([len(l) for l in loops]) < 3:
      print('HALT\n'*5,len(segments))
      ax = plot_axes()
      for l in loops:plot_edges(l,ax)
      plt.show()
      pdb.set_trace()

    loops = [tuple(clean_loop(l,plot)) for l in loops]

    if plot and False:
      print('got loops?!')
      ax = plot_axes()
      for lp in loops:
          ax = plot_polygon(list(lp),ax)
      plt.show()

    return loops

# given a set of edge segments, 
# return a boundary represented by the segments
def construct_loop(segments):
    srange = [x for x in range(len(segments))]
    loop = []
    last = segments[0][0].copy()
    while srange:
        for x in srange:
            s1,s2 = segments[x]
            if   last.near(s1):found = s2;break
            elif last.near(s2):found = s1;break
        srange.remove(x)
        newpoint(found.copy(),loop)
        last = found
    loop = clean_loop(loop)
    return tuple(loop)

# given a concave boundary, return the set of edge segments
# which generate the boundary
def boundary_segments(bnd):
    segments = []
    for x in range(len(bnd)):
        segments.append((bnd[x-1],bnd[x]))
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

# given two polygons, return the union of their edge segments
# broken upon every point of intersection between polygons
# or return None if there are no intersections
def break_polygons(p1,p2):
    p1segs = polygon_segments(p1)
    p2segs = polygon_segments(p2)
    isects = intersections(p1segs+p2segs)
    if len(isects) == 0:
        tpt1,tpt2 = p1segs[0][0],p2segs[0][0]
        p1inp2 = dpr.inconcave_xy(tpt1,p2[0])
        p2inp1 = dpr.inconcave_xy(tpt2,p1[0])
        if not (p1inp2 or p2inp1):
            return
    else:
        p1segs = break_segments(p1segs,isects)
        p2segs = break_segments(p2segs,isects)
    '''#
    ax = plot_axes_xy()
    plot_polygon_full_xy(p1,ax)
    plot_polygon_full_xy(p2,ax)
    for seg in p1segs:ax = plot_edges_xy(seg,ax,lw = 4.0)
    for seg in p2segs:ax = plot_edges_xy(seg,ax,lw = 4.0)
    plt.show()
    '''#
    return p1segs,p2segs

# given a set of segments, 
# return the subset on one side of a line
def segments_leftofline(segments,l1,l2):
    leftof = []
    tangs = 0
    for x in range(len(segments)):
        s1,s2 = segments[x]
        ori1 = dpr.orient2d(s1,l1,l2)
        ori2 = dpr.orient2d(s2,l1,l2)

        #print('oriririe left',ori1,ori2)
        if (ori1 > 0 or ori2 > 0):
            leftof.append((s1,s2))
        elif (ori1 == 0 and ori2 == 0):
            leftof.append((s1,s2))
            tangs += 1
    if tangs == len(leftof):return []
    '''#
    print('ithinkimleft')
    ax = plot_axes()
    for s in leftof:plot_edges(s,ax)
    plot_edges([l1,l2],ax,lw = 5.0)
    plt.show()
    '''#
    return leftof

# given a set of segments, 
# return the subset on one side of a line
def segments_rightofline(segments,l1,l2):
    rightof = []
    tangs = 0
    for x in range(len(segments)):
        s1,s2 = segments[x]
        ori1 = dpr.orient2d(s1,l1,l2)
        ori2 = dpr.orient2d(s2,l1,l2)
        #print('orieeright',ori1,ori2)
        #if (ori1 < 0 or ori2 < 0) or (ori1 == 0 and ori2 == 0):
        #    rightof.append((s1,s2))

        if (ori1 < 0 or ori2 < 0):
            rightof.append((s1,s2))
        elif (ori1 == 0 and ori2 == 0):
            rightof.append((s1,s2))
            tangs += 1
    if tangs == len(rightof):return []

    '''#
    print('ithinkimright')
    ax = plot_axes()
    for s in segments:plot_edges(s,ax,lw = 4.0)
    for s in rightof:plot_edges(s,ax)
    plot_edges([l1,l2],ax,lw = 8.0)
    plt.show()
    '''#
    return rightof

# given a set of segments, 
# return the subset which intersects a polygon
def segments_inpolygon(segments,py):
    eb,ibs = py
    inpoly = []
    for x in range(len(segments)):
        i1,i2 = segments[x]
        eisect = False
        for ex in range(len(eb)):
            ep1,ep2 = eb[ex-1],eb[ex]
            isect = segments_intersect_at(i1,i2,ep1,ep2)
            if not isect is None:eisect = True
        if eisect or dpr.inconcave_xy(dpv.midpoint(i1,i2),eb):
            inpoly.append((i1,i2))
    return inpoly

# given a set of segments, 
# return the subset which does not intersect a polygon
def segments_outpolygon(segments,py):
    eb,ibs = py
    outpoly = []
    for x in range(len(segments)):
        i1,i2 = segments[x]
        eisect = False
        for ex in range(len(eb)):
            ep1,ep2 = eb[ex-1],eb[ex]
            isect = segments_intersect_at(i1,i2,ep1,ep2)
            if not isect is None:eisect = True 
        if not eisect and not dpr.inconcave_xy(dpv.midpoint(i1,i2),eb):
            outpoly.append((i1,i2))
    return outpoly

def polygon_union(p1,p2):
    prot = valid_pair(p1,p2)
    if prot is None:return
    dpr.rotate_polygon(p1,prot)
    dpr.rotate_polygon(p2,prot)
    broken = break_polygons(p1,p2)
    if broken is None:
        dpr.rotate_polygon(p1,prot.flip())
        dpr.rotate_polygon(p2,prot)
        return
    else:p1segs,p2segs = broken
    p1inp2 = segments_outpolygon(p1segs,p2)
    p2inp1 = segments_outpolygon(p2segs,p1)
    union = (construct_loop(p1inp2+p2inp1),(p1[1]+p2[1]))
    dpr.rotate_polygon(union,prot.flip())
    dpr.rotate_polygon(p1,prot)
    dpr.rotate_polygon(p2,prot)
    return union

def polygon_intersection(p1,p2):
    prot = valid_pair(p1,p2)
    if prot is None:return
    dpr.rotate_polygon(p1,prot)
    dpr.rotate_polygon(p2,prot)
    broken = break_polygons(p1,p2)
    if broken is None:
        dpr.rotate_polygon(p1,prot.flip())
        dpr.rotate_polygon(p2,prot)
        return
    else:p1segs,p2segs = broken
    p1inp2 = segments_inpolygon(p1segs,p2)
    p2inp1 = segments_inpolygon(p2segs,p1)
    inter = (construct_loop(p1inp2+p2inp1),(p1[1]+p2[1]))
    dpr.rotate_polygon(inter,prot.flip())
    dpr.rotate_polygon(p1,prot)
    dpr.rotate_polygon(p2,prot)
    return inter

def polygon_difference(p1,p2):
    prot = valid_pair(p1,p2)
    if prot is None:return
    dpr.rotate_polygon(p1,prot)
    dpr.rotate_polygon(p2,prot)
    broken = break_polygons(p1,p2)
    if broken is None:
        dpr.rotate_polygon(p1,prot.flip())
        dpr.rotate_polygon(p2,prot)
        return
    else:p1segs,p2segs = broken
    p1inp2 = segments_outpolygon(p1segs,p2)
    p2inp1 = segments_inpolygon(p2segs,p1)

    print('THE CONTAINMENT PROBLEM ARISES!')
    print('THE CONTAINMENT PROBLEM ARISES!')
    print('THE CONTAINMENT PROBLEM ARISES!',p1inp2,p2inp1)
    if not p1inp2 and not p2inp1:
        dpr.rotate_polygon(p1,prot.flip())
        dpr.rotate_polygon(p2,prot)
        return

    #diffr = (construct_loop(p1inp2+p2inp1),(p1[1]+p2[1]))
    loops = construct_loops(p1inp2+p2inp1)
    if len(loops) == 0:
        print('FAILED TO MAKE LOOPS?!')
        pdb.set_trace()
    elif len(loops) == 1:loop,holes = tuple(loops[0]),(p1[1]+p2[1])
    elif len(loops) > 1:
        ct = containment(loops[0],loops[1])
        if   ct == -1:loop,holes = tuple(loops[0]),(tuple(loops[1]),)
        elif ct ==  1:loop,holes = tuple(loops[1]),(tuple(loops[0]),)
        else:loop,holes = tuple(loops[0]),(p1[1]+p2[1])

    #diffr = (construct_loops(p1inp2+p2inp1),(p1[1]+p2[1]))
    diffr = (loop,holes)
    dpr.rotate_polygon(diffr,prot.flip())
    dpr.rotate_polygon(p1,prot)
    dpr.rotate_polygon(p2,prot)
    return diffr

def polygon_xor(p1,p2):
    raise NotImplemented

def cgstest():
    p1 = (tuple(dpr.square(5,5)),())
    p2 = (tuple(dpr.square(4,4,dpv.vector(0,0,0))),())
    #p2 = (tuple(dpr.square(5,4,dpv.vector(4,0,0))),())
    #p3 = (tuple(dpr.square(5,3,dpv.vector(8,0,0))),())
    #dpr.rotate_polygon(p1,dpq.q_from_av(dpr.PI/2.0,dpv.x()))
    #dpr.rotate_polygon(p2,dpq.q_from_av(dpr.PI/2.0,dpv.x()))
    #dpr.rotate_polygon(p3,dpq.q_from_av(dpr.PI/2.0,dpv.x()))

    ptest = polygon_union(p1,p2)
    #ptest = polygon_union(ptest,p3)

    #ptest = polygon_intersection(p1,p2)

    #ptest = polygon_difference(p1,p2)


    #ptest = merge_polygons([p1,p2,p3])

    ax = plot_axes()
    #for seg in p1segs:ax = plot_edges_xy(seg,ax)
    #for seg in p2segs:ax = plot_edges_xy(seg,ax)
    #plot_polygon_full(p1,ax)
    #plot_polygon_full(p2,ax)
    plot_polygon_full(ptest,ax)
    #for pt in ptest:plot_polygon_full(pt,ax,lw = 4.0)
    plt.show()

    if ptest is None:print('unacceptable union')
    else:
        print('soo whats happeningjj?')
        ax = plot_axes()
        #for seg in p1segs:ax = plot_edges_xy(seg,ax)
        #for seg in p2segs:ax = plot_edges_xy(seg,ax)
        plot_polygon_full(p1,ax)
        plot_polygon_full(p2,ax)
        plot_polygon_full(ptest,ax,lw = 4.0)
        #for pt in ptest:plot_polygon_full(pt,ax,lw = 4.0)
        plt.show()





