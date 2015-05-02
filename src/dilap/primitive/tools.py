import dilap.core.base as db
import dp_vector as dpv

import numpy

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
        dx = l*numpy.cos(db.rad(current_angle))
        dy = l*numpy.sin(db.rad(current_angle))
        new = poly[-1].copy().translate_x(dx).translate_y(dy)
        poly.append(new)
    poly.pop()
    dpv.translate_coords(poly,dpv.center_of_mass(poly).flip())
    return poly

# return a ring of points of radius r with n corners
def point_ring(r,n):
    st = dpv.zero().translate_x(r)
    alpha = PI*(2.0/n)
    points = []
    for x in range(n):
        points.append(st.copy().rotate_z(x*alpha))
    return points

# return a square of length,width l,w at position p and zrot phi
def corners(l,w,p = None,phi = None):
    l2,w2 = l/2.0,w/2.0
    cs = [
        dpv.vector(-l2,-w2,0),dpv.vector( l2,-w2,0), 
        dpv.vector( l2, w2,0),dpv.vector(-l2, w2,0)]
    if not phi is None:dpv.rotate_coords_z(cs,phi)
    if not p is None:dpv.translate_coords(cs,p)
    return cs

# return a vector normal to the plane containing c1,c2,c3
def normal(c1,c2,c3):
    c1c2 = dpv.v1_v2(c1,c2).normalize()
    c2c3 = dpv.v1_v2(c2,c3).normalize()
    cn = c1c2.cross(c2c3).normalize()
    return cn

# add index offset to a list of faces
def offset_faces(faces,offset):
    for fdx in range(len(faces)):
        fa = faces[fdx]
        tfcnt = len(fa)
        for tfdx in range(tfcnt):
            fa[tfdx] += offset
    return faces


