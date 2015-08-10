import dilap.core.base as db

import dp_vector as dpv
import dp_quaternion as dpq
import dp_ray as dr

import math,numpy,os,appdirs,pdb
import matplotlib.pyplot as plt

PI = numpy.pi

def flatten(unflat):
    return [item for sublist in unflat for item in sublist]

# return a copy of seq without duplicates; preserves order
def uniqfy(seq):
    result = []
    for item in seq:
        if item in result:continue
        result.append(item)
    return result

# is seq1 a cyclic permutation of seq2?
def cyclic_permutation(seq1,seq2):
    s1cnt,s2cnt = len(seq1),len(seq2)
    if not s1cnt == s2cnt:return False
    for pmdx in range(s1cnt):
        perm = seq1[pmdx:] + seq1[:pmdx]
        if perm == seq2:return True
    return False

def plot_points(pts,proj = '3d',edges = True,mark = 'o',postpone = False):
    xs = [v.x for v in pts]
    ys = [v.y for v in pts]
    zs = [v.z for v in pts]
    fig = plt.figure()
    ax = fig.add_subplot(111,projection = proj)
    if edges:ax.plot(xs,ys,zs = zs,marker = mark,zdir = 'z')
    else:
        for x,y,z in zip(xs,ys,zs):
            ax.plot([x],[y],[z],marker = mark,zdir = 'z')
        #[ax.plot([x],[y],[z],marker = mark,zdir = 'z') for x,y,z in zip(xs,ys,zs)]
    #if not postpone:plt.show()
    plt.show()

# reduce a list of models to a single model
def combine(models):
    final = models[0]
    if len(models) > 1:
        for m in models[1:]:
            final._consume(m)
    return final

# given 3 points and a plane, determine the center and radius
# of the circumcenter found in the plane plane by projecting p1,p2,p3
# in the plane
def circumscribe_tri(p1,p2,p3,plane = None):
    if not plane is None:
        p1 = p1.project_plane(*plane)
        p2 = p2.project_plane(*plane)
        p3 = p3.project_plane(*plane)
    e1 = p1 - p3
    e2 = p2 - p3
    th = dpv.angle_between(e1,e2)
    th = numpy.pi if th < 0.0001 else th
    cr = dpv.distance(p1,p2)/(2*numpy.sin(th))
    cp = e2.copy().scale_u(e1.magnitude2())-e1.copy().scale_u(e2.magnitude2())
    cp = cp.cross(e1.cross(e2)).scale_u(1.0/(2.0*(e1.cross(e2).magnitude2())))
    return cp+p3,cr

def inside_circle(pt,center,radius,plane):
    pt = pt.project_plane(*plane)
    center = center.project_plane(*plane)
    ins = not dpv.distance(pt,center) > radius
    return ins

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

def dice_edges(verts,dices = 3):
    for di in range(dices):
        newpts = []
        vcnt = len(verts)
        for tdx in range(vcnt):
            p1 = verts[tdx-1]
            p2 = verts[tdx]
            mpt = dpv.midpoint(p1,p2)
            newpts.append(p1)
            newpts.append(mpt)
        newpts.append(newpts.pop(0))
        newpts.append(newpts.pop(0))
        verts = newpts
    return verts

def inflate(convex,radius):
    enorms = dpv.edge_normals_xy(convex)
    for cdx in range(len(convex)):
        lead = enorms[cdx]
        rear = enorms[cdx-1]
        norm = dpv.midpoint(lead,rear).normalize()
        convex[cdx].translate(norm.scale_u(radius))
    return convex

def pts_to_convex_xy(pts):
    # return the corners of the polygon, a subset of pts
    # it could be that pts is all one point or is colinear
    new = find_extremes_x(pts)[1]
    tang = None
    shape = []
    while not new in shape:
        shape.append(new)
        if len(shape) > 1:
            tang = dpv.v1_v2(shape[-2],shape[-1])
        new = sweep_search(pts,new,tang)
    return shape

def find_extremes_y(pts):
    lo = pts[0]
    hi = pts[0]
    for pt in pts[1:]:
        if pt.y < lo.y:lo = pt
        if pt.y > hi.y:hi = pt
    return lo,hi

def find_extremes_x(pts):
    lo = pts[0]
    hi = pts[0]
    for pt in pts[1:]:
        if pt.x < lo.x:lo = pt
        if pt.x > hi.x:hi = pt
    return lo,hi

def sweep_search(pts,center,tangent = None):
    # this will behave oddly when `which` would
    #  be a colinear set
    offset = center.copy().flip()
    dpv.translate_coords(pts,offset)
    if not tangent is None:
        tangent_rot = dpv.angle_from_xaxis_xy(tangent)
        dpv.rotate_z_coords(pts,-tangent_rot)
    which = center
    pang = 2*PI
    pcnt = len(pts)
    for adx in range(pcnt):
        pt = pts[adx]
        if pt is center: continue
        tpang = dpv.angle_from_xaxis_xy(pt)
        if tpang < pang:
            pang = tpang
            which = pt
    if not tangent is None:
        dpv.rotate_z_coords(pts,tangent_rot)
    dpv.translate_coords(pts,offset.flip())
    return which

def triangle_cover(boundary,side_length):
    side_length = float(side_length)
    perp_side_length = math.sqrt(3)*side_length/2.0
    bproj_x = dpv.project_coords(boundary,dpv.xhat)
    bproj_y = dpv.project_coords(boundary,dpv.yhat)
    xrng = bproj_x.y - bproj_x.x
    yrng = bproj_y.y - bproj_y.x
    total_offset = dpv.vector(
        bproj_x.x - side_length,
        bproj_y.x - perp_side_length,0)

    xtcnt = int(xrng/side_length+0.5) + 3
    ytcnt = int(yrng/perp_side_length+0.5) + 3

    broffset = dpv.vector(side_length/2.0,perp_side_length,0)
    corners = []
    xs = [x*side_length for x in range(xtcnt)]

    pts = []
    seedrow = [dpv.vector(x,0,0).translate(total_offset) for x in xs]
    rowleng = len(seedrow)
    pts.extend(seedrow)
    bottomrow = [x for x in range(rowleng)]
    pts.extend([b.copy().translate(broffset) for b in seedrow])
    toprow = [x+rowleng for x in range(rowleng)]
    
    def next_rows(bottom,top,even):
        sign = -1 if even else 1
        broffset.translate_x(sign*side_length)
        newpts = [pts[t].copy().translate(broffset) for t in top]
        pcnt = len(pts)
        pts.extend(newpts)
        newrng = range(pcnt,len(pts))
        newtop = [x for x in newrng]
        newbottom = [x for x in top]
        return newbottom,newtop

    for rdx in range(ytcnt):
        even = rdx % 2 == 0
        for vdx in range(1,len(bottomrow)):
            c1 = bottomrow[vdx-1]
            c2 = bottomrow[vdx]
            c3 = toprow[vdx-1]
            c4 = toprow[vdx]
            if even:
                corners.append((c1,c2,c3))
                corners.append((c3,c2,c4))
            else:
                corners.append((c1,c2,c4))
                corners.append((c3,c1,c4))
        if not rdx == ytcnt - 1:
            bottomrow,toprow = next_rows(bottomrow,toprow,even)

    extras = []
    keeppts = []
    for cnx in range(len(corners)):
        ptri = [pts[x] for x in corners[cnx]]
        ins = [dpv.inside(c,boundary) for c in ptri]
        if not ins.count(True) > 0:extras.append(cnx)
        else:
            for x in corners[cnx]:
                if not x in keeppts:
                    keeppts.append(x)

    extras = sorted(extras)
    extras.reverse()
    for x in extras:corners.pop(x)
    extrapts = [x for x in range(len(pts)) if not x in keeppts]
    extrapts = sorted(extrapts)
    extrapts.reverse()
    # THIS MESSES UP THE INDICES OF CORNERS!!!
    #for x in extrapts:pts.pop(x)
    #plot_points(pts+boundary)
    return pts,corners

def extrude_edge(c1,c2,length,direction):
    c1c2n = direction.copy().normalize().scale_u(length)
    c3 = c2.copy().translate(c1c2n)
    c4 = c1.copy().translate(c1c2n)
    return c3,c4

def extrude_edge_normal(c1,c2,length):
    c1c2 = dpv.v1_v2(c1,c2)
    c1c2n = dpv.cross(dpv.zhat,c1c2)
    return extrude_edge(c1,c2,length,c1c2n)

# return a square of length,width l,w at position p and zrot phi
def corners(l,w,p = None,phi = None):
    l2,w2 = l/2.0,w/2.0
    cs = [
        dpv.vector(-l2,-w2,0),dpv.vector( l2,-w2,0), 
        dpv.vector( l2, w2,0),dpv.vector(-l2, w2,0)]
    if not phi is None:dpv.rotate_z_coords(cs,phi)
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
        qrot = dpq.q_from_av(PI,dpv.yhat)
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

# keep the value val bounded by f and c by flooring
def clamp(v,f,c):
    if v < f: return f
    elif v > c: return c
    else: return v

# keep the value val bounded by f and c by wrapping around
def clamp_periodic(val,f,c):
    period = c - f
    while val < f: val += period
    while val > c: val -= period
    else: return val

def rad(deg):return PI*deg/180.0
def deg(rad):return 180.0*rad/PI

# return the path to a safe resource directory, 
# or a full path to a file therein
res_path = os.path.join(appdirs.user_data_dir(),'dilap_resources')
def resource_path(res = None):
    if res is None:rpath = res_path[:]
    else:rpath = os.path.join(res_path,res)
    return rpath


