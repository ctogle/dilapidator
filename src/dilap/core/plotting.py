from dilap.geometry.vec3 import vec3

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

import pdb

###############################################################################

# create a mpl 2d axes object
def plot_axes_xy(x = 5,f = None):
    if f is None:ax = plt.figure().add_subplot(111)
    else:ax = f.add_subplot(111)
    ax.set_xlim([-x,x])
    ax.set_ylim([-x,x])
    ax.set_aspect('equal')
    return ax

# create a mpl 3d axes object
def plot_axes(x = 5,f = None):
    if f is None:ax = plt.figure().add_subplot(111,projection = '3d')
    else:ax = f.add_subplot(111,projection = '3d')
    ax.set_xlim([-x,x])
    ax.set_ylim([-x,x])
    ax.set_zlim([-(9.0/16.0)*x,(9.0/16.0)*x])
    return ax

def plot_point_xy(pt,ax,mk = 'o',col = None):
    if col is None:col = 'black'
    ax.plot([pt.x],[pt.y],marker = mk,color = col)
    return ax

def plot_vector_xy(pt,tn,ax,mk = 'o',lw = 1.0,col = None):
    tip = pt.cp().trn(tn)
    ax = plot_point_xy(pt,ax,mk = mk,col = col)
    ax = plot_edges_xy([pt,tip],ax,mk = mk,lw = lw,col = col)
    ax = plot_point_xy(tip,ax,mk = 'd',col = col)
    return ax

def plot_point_xy_annotate(pt,ax,text):
    ax.annotate(text,xy = (pt.x,pt.y),xytext = (-20, 20),
        textcoords = 'offset points',ha = 'right',va = 'bottom',
        arrowprops = dict(arrowstyle = '->',
        connectionstyle = 'arc3,rad=0'))
    return ax

def plot_point(pt,ax,mk = 'o',col = None):
    if col is None:col = 'black'
    ax.plot([pt.x],[pt.y],zs = [pt.z],marker = mk,color = col)
    return ax

def plot_points_xy(points,ax = None,ms = None,cs = None,number = False):
    if ax is None:ax = plot_axes_xy()
    if ms is None:ms = ['o']*len(points)
    if cs is None:cs = [None]*len(points)
    for pdx in range(len(points)):
        plot_point_xy(points[pdx],ax,ms[pdx],cs[pdx])
        if number:plot_point_xy_annotate(points[pdx],ax,str(pdx+1))
    return ax

def plot_points(points,ax = None,ms = None,cs = None,marker = None):
    if ax is None:ax = plot_axes()
    if marker is None:marker = 'o'
    if ms is None:ms = [marker]*len(points)
    if cs is None:cs = [None]*len(points)
    for pdx in range(len(points)):plot_point(points[pdx],ax,ms[pdx],cs[pdx])  
    return ax

def plot_edges_xy(points,ax = None,mk = None,lw = 1.0,center = False,col = None):
    if ax is None:ax = plot_axes_xy()
    if mk is None:mk = '+'
    if col is None:col = 'black'
    pts = [p.__iter__() for p in points]
    xs,ys,zs = zip(*pts)
    ax.plot(xs,ys,marker = mk,lw = lw,color = col)
    if center:
        centers = [points[x-1].mid(points[x]) for x in range(1,len(points))]
        plot_points_xy(centers,ax)
    return ax

def plot_edges(points,ax = None,mk = None,lw = 1.0,center = False,col = None):
    if ax is None:ax = plot_axes()
    if mk is None:mk = '+'
    if col is None:col = 'black'
    pts = [p.__iter__() for p in points]
    xs,ys,zs = zip(*pts)
    ax.plot(xs,ys,zs,marker = mk,lw = lw,color = col)
    if center:
        centers = [points[x-1].mid(points[x]) for x in range(1,len(points))]
        plot_points(centers,ax)
    return ax

def plot_polygon_xy(points,ax = None,center = False,mk = None,lw = 1.0,col = None):
    epts = list(points[:])
    epts.append(points[0])
    ax = plot_edges_xy(epts,ax,mk = mk,lw = lw,col = col)
    if center:plot_point_xy(vec3(0,0,0).com(points),ax,mk = 's',col = col)
    return ax

def plot_polygon(points,ax = None,center = False,mk = None,lw = 1.0,col = None):
    epts = list(points[:])
    epts.append(points[0])
    ax = plot_edges(epts,ax,mk = mk,lw = lw,col = col)
    if center:plot_point(vec3(0,0,0).com(points),ax,mk = 's',col = col)
    return ax

def plot_polygon_full_xy(poly,ax = None,center = False,lw = 1.0,col = None):
    if ax is None:ax = plot_axes_xy()
    ebnd,ibnds = poly
    plot_polygon_xy(list(ebnd),ax,center = True,lw = lw,col = col)
    for ib in ibnds:plot_polygon_xy(list(ib),ax,center = True,lw = lw,col = col)
    return ax

def plot_polygon_full(poly,ax = None,center = False,lw = 1.0,col = None):
    if ax is None:ax = plot_axes()
    ebnd,ibnds = poly
    plot_polygon(list(ebnd),ax,center = True,lw = lw,col = col)
    for ib in ibnds:plot_polygon(list(ib),ax,center = True,lw = lw,col = col)
    return ax

def plot_tetrahedron(points,ax = None):
    raise NotImplemented

def plot_line_xy(l1,l2,r = 25,ax = None,center = False,lw = 1.0,col = None):
    ltan = l1.tov(l2).cpxy()
    l1far = l1.cp().trn(ltan.cp().uscl(-r))
    l2far = l2.cp().trn(ltan.cp().uscl( r))
    ax = plot_edges_xy([l1far,l2far],ax = ax,center = center,lw = lw,col = col)
    return ax

def plot_line(l1,l2,r = 25,ax = None,center = False,lw = 1.0,col = None):
    ltan = l1.tov(l2)
    l1far = l1.cp().trn(ltan.cp().uscl(-r))
    l2far = l2.cp().trn(ltan.cp().uscl( r))
    ax = plot_edges([l1far,l2far],ax = ax,center = center,lw = lw,col = col)
    return ax

def plot_circle_xy(c,r,ax = None,center = False,lw = 1.0,col = None):
    circ = c.pring(r,32)
    ax = plot_polygon_xy(circ,ax,center,lw,col)
    return ax

def plot_circle(c,r,ax = None,center = False,lw = 1.0,col = None):
    circ = c.pring(r,32)
    ax = plot_polygon(circ,ax,center,lw,col)
    return ax

def plot_ray_xy(r,ax = None,col = None):
    ax = plot_point_xy(r.o,ax,mk = 's')
    ax = plot_edges_xy((r.o,r.o.cp().trn(r.d.cp().scl(100))),ax,lw = 2.0,col = col)
    return ax

def plot_ray(r,ax = None,col = None):
    ax = plot_point(r.o,ax,mk = 's')
    ax = plot_edges((r.o,r.o.cp().trn(r.d.cp().scl(100))),ax,lw = 2.0,col = col)
    return ax

###############################################################################





