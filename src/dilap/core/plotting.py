#import dilap.core.tools as dpr
import dilap.geometry.tools as dpr
#import dilap.core.vector as dpv
from dilap.geometry.vec3 import vec3
#import dilap.core.quaternion as dpq

#import dilap.core.bbox as dpb

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import math,numpy

import pdb

###############################################################################

# create a mpl 2d axes object
def plot_axes_xy():
    ax = plt.figure().add_subplot(111)
    ax.set_aspect('equal')
    return ax

# create a mpl 3d axes object
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
    #pts = [p.to_tuple() for p in points]
    pts = [p.__iter__() for p in points]
    xs,ys,zs = zip(*pts)
    ax.plot(xs,ys,marker = mk,lw = lw)
    if center:
        centers = [dpv.midpoint(points[x-1],points[x]) 
                        for x in range(1,len(points))]
        plot_points_xy(centers,ax)
    return ax

def plot_edges(points,ax = None,lw = 1.0,center = False):
    if ax is None:ax = plot_axes()
    #pts = [p.to_tuple() for p in points]
    pts = [p.__iter__() for p in points]
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

def plot_ray_xy(r,ax = None):
    ax = plot_point_xy(r.o,ax,marker = 's')
    ax = plot_edges_xy((r.o,r.o.cp().trn(r.d.cp().scl(100))),ax,lw = 2.0)
    return ax

def plot_ray(r,ax = None):
    ax = plot_point(r.o,ax,marker = 's')
    ax = plot_edges((r.o,r.o.cp().trn(r.d.cp().scl(100))),ax,lw = 2.0)
    return ax

###############################################################################





