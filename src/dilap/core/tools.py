import dilap.core.base as db

import dp_vector as dpv
import dp_quaternion as dpq
import dp_ray as dr

import numpy,os,appdirs,pdb

PI = numpy.pi

# reduce a list of models to a single model
def combine(models):
    final = models[0]
    if len(models) > 1:
        for m in models[1:]:
            final._consume(m)
    return final

# return a polygon of radius 1 with n sides
# similar to point_ring
def polygon(n):
    angle = 360.0/n
    turns = [x*angle for x in range(n)]
    poly = [dpv.zero()]
    current_angle = 0.0
    for si in range(n):
        l,t = 1.0,turns[si]
        current_angle = t
        dx = l*numpy.cos(rad(current_angle))
        dy = l*numpy.sin(rad(current_angle))
        new = poly[-1].copy().translate_x(dx).translate_y(dy)
        poly.append(new)
    poly.pop()
    dpv.translate_coords(poly,dpv.center_of_mass(poly).flip())
    return poly

# given start point s, end point e, and n segments, 
# return a colinear set of points equally spaced between s and e
def point_line(s,e,n):
    line = [s.copy()]
    tn = dpv.v1_v2(s,e)
    l = tn.magnitude()
    tn.normalize()
    tn.scale_u(float(l)/n)
    for x in range(n):
        line.append(line[-1].copy().translate(tn))
    return line

# return a ring of points of radius r with n corners
def point_ring(r,n):
    st = dpv.zero().translate_x(r)
    alpha = PI*(2.0/n)
    points = []
    for x in range(n):
        points.append(st.copy().rotate_z(x*alpha))
    return points

def extrude_edge(c1,c2,length,direction):
    c1c2n = direction.copy().normalize().scale_u(length)
    c3 = c2.copy().translate(c1c2n)
    c4 = c1.copy().translate(c1c2n)
    return c3,c4

def extrude_edge_normal(c1,c2,length):
    c1c2 = cv.v1_v2(c1,c2)
    c1c2n = cv.cross(cv.zhat,c1c2)
    return extrude_edge(c1,c2,length,c1c2n)

# return a square of length,width l,w at position p and zrot phi
def corners(l,w,p = None,phi = None):
    l2,w2 = l/2.0,w/2.0
    cs = [
        dpv.vector(-l2,-w2,0),dpv.vector( l2,-w2,0), 
        dpv.vector( l2, w2,0),dpv.vector(-l2, w2,0)]
    if not phi is None:dpv.rotate_coords_z(cs,phi)
    if not p is None:dpv.translate_coords(cs,p)
    return cs

def project_coords_plane_nearest(coords,p0,n):
    projected = []
    for c in coords:
        v = c - p0
        d = v.dot(n)
        nd = n.copy().scale_u(d)
        pj = c - nd
        projected.append(pj)
    return projected

def project_coords_plane_along(coords,p0,n,direc):
    projected = []
    for c in coords:
        r = dr.ray(c,direc)
        iplane = r.intersect_plane(p0,n)
        if iplane:pj = r.cast.copy()
        else:raise ValueError
        projected.append(pj)
    return projected

# origin is the point where the loop should be anchored to
# loop starts in the space of the origin, which is usually world space
# control is the point in the loops local space which is the anchor
# targetnormal is the vector which the normal of loop must coincide with
#   the loop must reside in exactly one plane!
#
# orient the loop by rotating around its control point
# then translate by a vector pointing from control to origin
def orient_loop(loop,targetnormal,control = None):
    if control is None:control = dpv.center_of_mass(loop)
    n = normal(*loop[:3])
    if n == targetnormal.copy().flip():
        #print('HACK?')
        #qrot = dpq.q_from_av(numpy.pi,dpv.zhat)
        qrot = dpq.q_from_av(numpy.pi,dpv.yhat)
    else:qrot = dpq.q_from_uu(n,targetnormal)
    looprot = dpq.rotate_coords(loop,control,qrot)
    return looprot

# return a vector normal to the plane containing c1,c2,c3
def normal(c1,c2,c3):
    c1c2 = dpv.v1_v2(c1,c2).normalize()
    c2c3 = dpv.v1_v2(c2,c3).normalize()
    cn = c1c2.cross(c2c3).normalize()
    return cn

# return a vector tanget to the plane containing c1,c2,c3
def tangent(c1,c2,c3):
    tn = dpv.v1_v2(c1,c2).normalize()
    return tn

# add index offset to a list of faces
def offset_faces(faces,offset):
    for fdx in range(len(faces)):
        fa = faces[fdx]
        tfcnt = len(fa)
        for tfdx in range(tfcnt):
            fa[tfdx] += offset
    return faces

#
def clamp(v,f,c):
    if v < f: return f
    elif v > c: return c
    else: return v

def rad(deg):return numpy.pi*deg/180.0
def deg(rad):return 180.0*rad/numpy.pi

# return the path to a safe resource directory, 
# or a full path to a file therein
res_path = os.path.join(appdirs.user_data_dir(),'dilap_resources')
def resource_path(res = None):
    if res is None:rpath = res_path[:]
    else:rpath = os.path.join(res_path,res)
    return rpath


