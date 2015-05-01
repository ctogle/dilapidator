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

# return a ring of points of radius r with n corners
def point_ring(r,n):
    st = dpv.zero().translate_x(r)
    alpha = PI*(2.0/n)
    points = []
    for x in range(n):
        points.append(st.copy().rotate_z(x*alpha))
    return points

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


