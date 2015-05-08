#imports
# cython: profile=True
#cimport cython
cimport dp_quaternion as dpq

import matplotlib.pyplot as plt

from libc.math cimport sqrt
from libc.math cimport cos
from libc.math cimport sin
from libc.math cimport tan
from libc.math cimport hypot
import numpy as np
 
#import make_places.core.support.mp_utils as mpu



stuff = 'hi'

# 2d classes/functions

cdef class vector2d:

    def __cinit__(self,float x,float y):
        self.x = x
        self.y = y

    def __str__(self):
        strr = 'vector2d:' + str((self.x,self.y))
        return strr

    cpdef vector2d copy(self):
        cdef vector2d new = vector2d(self.x,self.y)
        return new

    cpdef list to_list(self):
        cdef list new = [self.x,self.y]
        return new

    cpdef tuple to_tuple(self):
        cdef tuple new = (self.x,self.y)
        return new

    cpdef vector2d flip(self):
        self.x *= -1.0
        self.y *= -1.0
        return self

    cpdef vector2d reciprocate(self):
        self.x = 1.0/self.x
        self.y = 1.0/self.y
        return self

    cpdef vector2d translate(self, vector2d tv):
        self.x += tv.x
        self.y += tv.y
        return self

    cpdef vector2d translate_x(self, float tx):
        self.x += tx
        return self

    cpdef vector2d translate_y(self, float ty):
        self.y += ty
        return self

    cpdef vector2d scale(self, vector2d sv):
        self.x *= sv.x
        self.y *= sv.y
        return self

    cpdef vector2d scale_x(self, float sx):
        self.x *= sx
        return self

    cpdef vector2d scale_y(self, float sy):
        self.y *= sy
        return self

cdef vector2d midpoint2d_c(vector2d v1, vector2d v2):
    cdef float x = (v1.x + v2.x)/2.0
    cdef float y = (v1.y + v2.y)/2.0
    cdef vector2d new = vector2d(x,y)
    return new

cpdef vector2d midpoint2d(vector2d v1, vector2d v2):
    cdef vector2d pt = midpoint2d_c(v1, v2)
    return pt

cdef void translate_coords2d_c(list coords, vector2d t):
    cdef int ccnt = len(coords)
    cdef int cdx
    cdef vector2d coo
    for cdx in range(ccnt):
        coo = <vector2d>coords[cdx]
        coo.translate(t)

cpdef translate_coords2d(list coords, vector2d t):
    translate_coords2d_c(coords,t)

cdef void scale_coords2d_c(list coords, vector2d s):
    cdef int ccnt = len(coords)
    cdef int cdx
    cdef vector2d coo
    for cdx in range(ccnt):
        coo = <vector2d>coords[cdx]
        coo.scale(s)

cpdef scale_coords2d(list coords, vector2d s):
    scale_coords2d_c(coords,s)

# 3d classes/functions

cdef class vector:

    def plot_xy(self):
        return [(self.x,self.y)]

    def __cinit__(self,float x,float y,float z):
        self.x = x
        self.y = y
        self.z = z

    def __str__(self):
        strr = 'vector:' + str((self.x,self.y,self.z))
        return strr

    def __add__(self,other):
        new = vector(self.x+other.x,self.y+other.y,self.z+other.z)
        return new

    def __sub__(self,other):
        new = vector(self.x-other.x,self.y-other.y,self.z-other.z)
        return new

    def __mul__(self,other):
        new = vector(self.x*other.x,self.y*other.y,self.z*other.z)
        return new

    def __richcmp__(self, other, comparator):
        if self.x == other.x:
            if self.y == other.y:
                if self.z == other.z: return True
        return False

    cpdef vector2d xy2d(self):
        cdef vector2d new = vector2d(self.x,self.y)
        return new

    cpdef vector2d xz2d(self):
        cdef vector2d new = vector2d(self.x,self.z)
        return new

    cpdef vector2d yz2d(self):
        cdef vector2d new = vector2d(self.y,self.z)
        return new

    cpdef vector xy(self):
        cdef vector new = vector(self.x,self.y,0.0)
        return new

    cpdef vector xz(self):
        cdef vector new = vector(self.x,0.0,self.z)
        return new

    cpdef vector yz(self):
        cdef vector new = vector(0.0,self.y,self.z)
        return new

    cpdef vector copy(self):
        cdef vector new = vector(self.x,self.y,self.z)
        return new

    def __iter__(self):
        yield self.x
        yield self.y
        yield self.z

    cpdef list to_list(self):
        cdef list new = [self.x,self.y,self.z]
        return new

    cpdef tuple to_tuple(self):
        cdef tuple new = (self.x,self.y,self.z)
        return new

    # return the magintude of self squared
    cpdef float magnitude2(self):
        cdef float xx = self.x**2
        cdef float yy = self.y**2
        cdef float zz = self.z**2
        cdef float ss = xx + yy + zz
        return ss

    cpdef float magnitude(self):
        cdef float xx = self.x**2
        cdef float yy = self.y**2
        cdef float zz = self.z**2
        cdef float ss = sqrt(xx + yy + zz)
        return ss

    cpdef vector normalize(self):
        cdef float mag = self.magnitude()
        if not mag == 0.0:
            self.x /= mag
            self.y /= mag
            self.z /= mag
        return self

    cpdef vector flip(self):
        self.x *= -1.0
        self.y *= -1.0
        self.z *= -1.0
        return self

    cpdef vector reciprocate(self):
        self.x = 1.0/self.x
        self.y = 1.0/self.y
        self.z = 1.0/self.z
        return self

    # linearly interpolate between self and other proportionally to delta
    cpdef vector linterpolate(self,vector other,float delta):
        cdef float dx = self.x + (other.x - self.x)*delta
        cdef float dy = self.y + (other.y - self.y)*delta
        cdef float dz = self.z + (other.z - self.z)*delta
        cdef vector new = vector(dx,dy,dz)
        return new

    # return a copy of self in basis of [b1,b2,b3]
    cpdef vector in_basis(self,vector b1,vector b2,vector b3):
        cdef float bx = self.dot(b1)/b1.magnitude2()
        cdef float by = self.dot(b2)/b2.magnitude2()
        cdef float bz = self.dot(b3)/b3.magnitude2()
        cdef vector new = vector(bx,by,bz)
        return new 

    # rotate self in place by q, without modifying q
    cpdef vector rotate(self, dpq.quaternion q):
        cdef float row1x = q.w**2 + q.x**2 - q.y**2 - q.z**2
        cdef float row1y = 2*(q.x*q.y - q.w*q.z)
        cdef float row1z = 2*(q.x*q.z + q.w*q.y)
        cdef vector row1 = vector(row1x,row1y,row1z)
        cdef float row2x = 2*(q.x*q.y + q.w*q.z)
        cdef float row2y = q.w**2 - q.x**2 + q.y**2 - q.z**2
        cdef float row2z = 2*(q.y*q.z - q.w*q.x)
        cdef vector row2 = vector(row2x,row2y,row2z)
        cdef float row3x = 2*(q.x*q.z - q.w*q.y)
        cdef float row3y = 2*(q.y*q.z + q.w*q.x)
        cdef float row3z = q.w**2 - q.x**2 - q.y**2 + q.z**2
        cdef vector row3 = vector(row3x,row3y,row3z)
        cdef float rotx = row1.dot(self)
        cdef float roty = row2.dot(self)
        cdef float rotz = row3.dot(self)
        self.x = rotx
        self.y = roty
        self.z = rotz
        return self

    cpdef vector rotate_x(self, float zang):
        cdef float cosz = cos(zang)
        cdef float sinz = sin(zang)
        cdef float newy = cosz*self.y - sinz*self.z
        cdef float newz = sinz*self.y + cosz*self.z
        self.y = newy
        self.z = newz
        return self

    cpdef vector rotate_y(self, float zang):
        cdef float cosz = cos(zang)
        cdef float sinz = sin(zang)
        cdef float newx = cosz*self.x - sinz*self.z
        cdef float newz = sinz*self.x + cosz*self.z
        self.x = newx
        self.z = newz
        return self

    cpdef vector rotate_z(self, float zang):
        cdef float cosz = cos(zang)
        cdef float sinz = sin(zang)
        cdef float newx = cosz*self.x - sinz*self.y
        cdef float newy = sinz*self.x + cosz*self.y
        self.x = newx
        self.y = newy
        return self

    cpdef vector translate_x(self, float tx):
        self.x += tx
        return self

    cpdef vector translate_y(self, float ty):
        self.y += ty
        return self

    cpdef vector translate_z(self, float tz):
        self.z += tz
        return self

    cpdef vector translate(self, vector tv):
        self.x += tv.x
        self.y += tv.y
        self.z += tv.z
        return self

    cpdef vector scale(self, vector sv):
        self.x *= sv.x
        self.y *= sv.y
        self.z *= sv.z
        return self

    cpdef vector scale_u(self, float s):
        self.x *= s
        self.y *= s
        self.z *= s
        return self

    cpdef vector cross(self, vector v):
        return cross_c(self,v)

    cpdef float dot(self, vector v):
        return dot_c(self,v)

cdef vector zero_c():
    cdef vector new = vector(0,0,0)
    return new

cpdef vector zero():
    return zero_c()

cdef vector2d zero2d_c():
    cdef vector2d new = vector2d(0,0)
    return new

cpdef vector2d zero2d():
    return zero2d_c()

cdef vector one_c():
    cdef vector new = vector(1,1,1)
    return new

cpdef vector one():
    return one_c()

cdef vector2d one2d_c():
    cdef vector2d new = vector2d(1,1)
    return new

cpdef vector2d one2d():
    return one2d_c()

cdef vector flip_c(vector f):
    cdef vector new = vector(-1.0*f.x,-1.0*f.y,-1.0*f.z)
    return new

cpdef vector flip(vector f):
    return flip_c(f)

cdef vector2d flip2d_c(vector2d f):
    cdef vector2d new = vector2d(-1.0*f.x,-1.0*f.y)
    return new

cpdef vector2d flip2d(vector2d f):
    return flip2d_c(f)

cdef vector normalize_c(vector v):
    cdef vector new = v.copy().normalize()
    return new

cpdef vector normalize(vector v):
    return normalize_c(v)

cpdef float magnitude(vector v):
    cdef float xx = v.x**2
    cdef float yy = v.y**2
    cdef float zz = v.z**2
    cdef float ss = sqrt(xx + yy + zz)
    return ss

cdef float dot_c(vector v1, vector v2):
    return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z

cpdef float dot(vector v1, vector v2):
    return dot_c(v1,v2)

cdef vector cross_c(vector v1, vector v2):
    cdef float cx = v1.y*v2.z-v1.z*v2.y
    cdef float cy = v1.z*v2.x-v1.x*v2.z
    cdef float cz = v1.x*v2.y-v1.y*v2.x
    cdef vector res = vector(cx,cy,cz)
    return res

cpdef vector cross(vector v1, vector v2):
    return cross_c(v1,v2)

cdef vector v1_v2_c(vector v1, vector v2):
    cdef float dx = v2.x - v1.x
    cdef float dy = v2.y - v1.y
    cdef float dz = v2.z - v1.z
    cdef vector new = vector(dx,dy,dz)
    return new

cpdef vector v1_v2(vector v1, vector v2):
    cdef vector pt = v1_v2_c(v1, v2)
    return pt

cdef vector vzip_c(vector v1, vector v2):
    cdef float x = v1.x*v2.x
    cdef float y = v1.y*v2.y
    cdef float z = v1.z*v2.z
    return vector(x,y,z)

cpdef vector vzip(vector v1, vector v2):
    cdef vector pt = vzip_c(v1, v2)
    return pt

cdef vector midpoint_c(vector v1, vector v2):
    cdef float x = (v1.x + v2.x)/2.0
    cdef float y = (v1.y + v2.y)/2.0
    cdef float z = (v1.z + v2.z)/2.0
    cdef vector new = vector(x,y,z)
    return new

cpdef vector midpoint(vector v1, vector v2):
    cdef vector pt = midpoint_c(v1, v2)
    return pt

cdef vector barymetric_to_world_c(float u,float v,vector v0,vector v1,vector v2):
    cdef vector b0 = v0.copy().scale_u(1-u-v)
    cdef vector b1 = v1.copy().scale_u(u)
    cdef vector b2 = v2.copy().scale_u(v)
    cdef vector iw = b0.translate(b1).translate(b2)
    return iw

cpdef vector barymetric_to_world(float u,float v,vector v0,vector v1,vector v2):
    return barymetric_to_world_c(u,v,v0,v1,v2)

cdef vector2d project_coords_c(list coords, vector axis):
    cdef ccnt = len(coords)
    cdef int cdx
    cdef float proj = coords[0].dot(axis)
    cdef vector2d projv = vector2d(proj,proj)
    for cdx in range(1,ccnt):
        proj = coords[cdx].dot(axis)
        if proj < projv.x:projv.x = proj
        if proj > projv.y:projv.y = proj
    return projv

cpdef vector2d project_coords(list coords, vector axis):
    return project_coords_c(coords,axis)

cdef float distance_to_edge_c(vector pt,vector e1,vector e2,vector nm):
    print('not implemented?')
    eproj = project_coords_c([e1,e2],nm)
    pproj = project_coords_c([pt],nm)
    return abs(eproj.x - pproj.x)

cpdef float distance_to_edge(vector pt,vector e1,vector e2,vector nm):
    return distance_to_edge_c(pt,e1,e2,nm)

cdef list edge_tangents_c(list verts):
    cdef list tangs = []
    cdef int vcnt = len(verts)
    cdef int vdx
    cdef vector v1
    cdef vector v2
    cdef float dx
    cdef float dy
    cdef float dv
    cdef vector tang
    for vdx in range(1,vcnt):
        v1,v2 = verts[vdx-1],verts[vdx]
        tang = v1_v2(v1,v2).normalize()
        tangs.append(tang)
    return tangs

cpdef list edge_tangents(list verts):
    return edge_tangents_c(verts)

cdef list edge_normals_xy_c(list verts):
    cdef list norms = []
    cdef int vcnt = len(verts)
    cdef int vdx
    cdef vector v1
    cdef vector v2
    cdef float dx
    cdef float dy
    cdef float dv
    cdef vector norm
    for vdx in range(vcnt):
        v1,v2 = verts[vdx-1],verts[vdx]
        dx = v2.x - v1.x
        dy = v2.y - v1.y
        dv = sqrt(dx**2 + dy**2)
        norm = vector(dy/dv,-dx/dv,0)
        norms.append(norm)
    norms.append(norms.pop(0))
    return norms

cpdef list edge_normals_xy(list verts):
    return edge_normals_xy_c(verts)

cdef float distance_to_border_xy_c(vector pt,list border):
    edgenorms = edge_normals_xy_c(border)
    dists = []
    for edx in range(len(border)):
        e1 = border[edx-1]
        e2 = border[edx]
        norm = edgenorms[edx-1]
        dists.append(distance_to_edge(pt,e1,e2,norm))
    dists.append(dists.pop(0))
    distance = min(dists)
    return distance

cpdef float distance_to_border_xy(vector pt,list border):
    return distance_to_border_xy_c(pt,border)

#########################################################################
### spline interpolation business
#########################################################################

cpdef list parameterize_time(list points,list time,float alpha):
    cdef float total = 0.0
    cdef int idx
    cdef list v1v2
    for idx in range(1,4):
        total += v1_v2_c(points[idx-1],points[idx]).magnitude()**(2.0*alpha)
        time[idx] = total

cpdef list catmull_rom(list P,list T,int tcnt):
    cdef int j
    cdef int t
    cdef float tt
    cdef float p
    cdef list spl = P[:1]
    for j in range(1, len(P)-2):  # skip the ends
        for t in range(tcnt):  # t: 0 .1 .2 .. .9
            tt = float(t)/tcnt
            tt = T[1] + tt*(T[2]-T[1])
            p = spline([P[j-1], P[j], P[j+1], P[j+2]],
                    [T[j-1], T[j], T[j+1], T[j+2]],tt)
            spl.append(p)
    spl.extend(P[-2:])
    return spl

cpdef float spline(list p,list time,float t):
    L01 = p[0] * (time[1] - t) / (time[1] - time[0]) + p[1] * (t - time[0]) / (time[1] - time[0])
    L12 = p[1] * (time[2] - t) / (time[2] - time[1]) + p[2] * (t - time[1]) / (time[2] - time[1])
    L23 = p[2] * (time[3] - t) / (time[3] - time[2]) + p[3] * (t - time[2]) / (time[3] - time[2])
    L012 = L01 * (time[2] - t) / (time[2] - time[0]) + L12 * (t - time[0]) / (time[2] - time[0])
    L123 = L12 * (time[3] - t) / (time[3] - time[1]) + L23 * (t - time[1]) / (time[3] - time[1])
    C12 = L012 * (time[2] - t) / (time[2] - time[1]) + L123 * (t - time[1]) / (time[2] - time[1])
    return C12

cpdef list vector_spline(vector c1,vector c2,vector c3,vector c4,int scnt):
    cox = [c1.x,c2.x,c3.x,c4.x]
    coy = [c1.y,c2.y,c3.y,c4.y]
    coz = [c1.z,c2.z,c3.z,c4.z]
    tim = [0.0,1.0,2.0,3.0]
    alpha = 1.0/2.0
    parameterize_time([c1,c2,c3,c4],tim,alpha)
    cox = catmull_rom(cox,tim,scnt)[1:-1]
    coy = catmull_rom(coy,tim,scnt)[1:-1] 
    coz = catmull_rom(coz,tim,scnt)[1:-1] 
    filled = [vector(*i) for i in zip(cox,coy,coz)]
    return filled

###########################################################################

cdef void rotate_x_coords_c(list coords, float ang):
    cdef int ccnt = len(coords)
    cdef int cdx
    cdef vector coo
    for cdx in range(ccnt):
        coo = <vector>coords[cdx]
        coo.rotate_x(ang)

cpdef rotate_x_coords(list coords, float ang):
    rotate_x_coords_c(coords,ang)

cdef void rotate_y_coords_c(list coords, float ang):
    cdef int ccnt = len(coords)
    cdef int cdx
    cdef vector coo
    for cdx in range(ccnt):
        coo = <vector>coords[cdx]
        coo.rotate_y(ang)

cpdef rotate_y_coords(list coords, float ang):
    rotate_y_coords_c(coords,ang)

cdef void rotate_z_coords_c(list coords, float ang):
    cdef int ccnt = len(coords)
    cdef int cdx
    cdef vector coo
    for cdx in range(ccnt):
        coo = <vector>coords[cdx]
        coo.rotate_z(ang)

cpdef rotate_z_coords(list coords, float ang):
    rotate_z_coords_c(coords,ang)

cdef void rotate_z_coords_about_c(list coords, vector about, float ang):
    cdef int ccnt = len(coords)
    cdef int cdx
    cdef vector coo
    for cdx in range(ccnt):
        coo = <vector>coords[cdx]
        coo.translate(about.flip()).rotate_z(ang).translate(about.flip())

cpdef rotate_z_coords_about(list coords, vector about, float ang):
    rotate_z_coords_about_c(coords,about,ang)

cdef void translate_coords_x_c(list coords, float tv):
    cdef int ccnt = len(coords)
    cdef int cdx
    cdef vector coo
    for cdx in range(ccnt):
        coo = <vector>coords[cdx]
        coo.translate_x(tv)

cpdef translate_coords_x(list coords, float tv):
    translate_coords_x_c(coords,tv)

cdef void translate_coords_y_c(list coords, float tv):
    cdef int ccnt = len(coords)
    cdef int cdx
    cdef vector coo
    for cdx in range(ccnt):
        coo = <vector>coords[cdx]
        coo.translate_y(tv)

cpdef translate_coords_y(list coords, float tv):
    translate_coords_y_c(coords,tv)

cdef void translate_coords_z_c(list coords, float tv):
    cdef int ccnt = len(coords)
    cdef int cdx
    cdef vector coo
    for cdx in range(ccnt):
        coo = <vector>coords[cdx]
        coo.translate_z(tv)

cpdef translate_coords_z(list coords, float tv):
    translate_coords_z_c(coords,tv)

cdef void translate_coords_c(list coords, vector t):
    cdef int ccnt = len(coords)
    cdef int cdx
    cdef vector coo
    for cdx in range(ccnt):
        coo = <vector>coords[cdx]
        coo.translate(t)

cpdef translate_coords(list coords, vector t):
    translate_coords_c(coords,t)

cdef void scale_coords_x_c(list coords, float s):
    cdef int ccnt = len(coords)
    cdef int cdx
    cdef vector coo
    for cdx in range(ccnt):
        coo = <vector>coords[cdx]
        coo.x *= s

cpdef scale_coords_x(list coords, float s):
    scale_coords_x_c(coords,s)

cdef void scale_coords_y_c(list coords, float s):
    cdef int ccnt = len(coords)
    cdef int cdx
    cdef vector coo
    for cdx in range(ccnt):
        coo = <vector>coords[cdx]
        coo.y *= s

cpdef scale_coords_y(list coords, float s):
    scale_coords_y_c(coords,s)

cdef void scale_coords_z_c(list coords, float s):
    cdef int ccnt = len(coords)
    cdef int cdx
    cdef vector coo
    for cdx in range(ccnt):
        coo = <vector>coords[cdx]
        coo.z *= s

cpdef scale_coords_z(list coords, float s):
    scale_coords_z_c(coords,s)
    
cdef void scale_coords_c(list coords, vector t):
    cdef int ccnt = len(coords)
    cdef int cdx
    cdef vector coo
    for cdx in range(ccnt):
        coo = <vector>coords[cdx]
        coo.scale(t)

cpdef scale_coords(list coords, vector t):
    scale_coords_c(coords,t)

cdef vector com(list coords):
    cdef int ccnt = len(coords)
    cdef float ccntf = float(ccnt)
    cdef int cdx = 0
    cdef float x = 0.0
    cdef float y = 0.0
    cdef float z = 0.0
    cdef vector coo
    cdef vector new
    for cdx in range(ccnt):
        coo = <vector>coords[cdx]
        x += coo.x
        y += coo.y
        z += coo.z
    new = vector(x/ccntf,y/ccntf,z/ccntf)
    return new

cpdef vector center_of_mass(list coords):
    return com(coords)

cdef float distance_xy_c(vector v1, vector v2):
    cdef float dx = v2.x - v1.x
    cdef float dy = v2.y - v1.y
    cdef float ds = sqrt(dx**2 + dy**2)
    return ds

cpdef float distance_xy(vector v1, vector v2):
    return distance_xy_c(v1, v2)

cdef float distance_c(vector v1, vector v2):
    cdef vector v1v2 = v1_v2(v1,v2)
    cdef float mag = v1v2.magnitude()
    return mag

cpdef float distance(vector v1, vector v2):
    return distance_c(v1, v2)

cpdef bint near(vector v1, vector v2):
    cdef bint isnear = distance_c(v1,v2) < 0.01
    return isnear

cpdef bint near_xy(vector v1, vector v2):
    cdef bint isnear = distance_xy_c(v1,v2) < 0.01
    return isnear

cdef int find_closest_c(vector one,list bunch,int bcnt,float close_enough):
    cdef float nearest = 100000000.0
    cdef float ds = nearest
    cdef int bdx
    cdef int ndx = 0
    cdef vector which
    for bdx in range(bcnt):
        which = <vector>bunch[bdx]
        ds = distance_c(one,which)
        if ds < nearest:
            nearest = ds
            ndx = bdx
            if ds <= close_enough:
                return ndx
    return ndx

cpdef int find_closest(vector one,list bunch,int bcnt,float close_enough):
    return find_closest_c(one,bunch,bcnt,close_enough)

cdef int find_closest_xy_c(vector one,list bunch,int bcnt,float close_enough):
    cdef float nearest = 100000000.0
    cdef float ds = nearest
    cdef int bdx
    cdef int ndx = 0
    cdef vector which
    for bdx in range(bcnt):
        which = <vector>bunch[bdx]
        ds = distance_xy_c(one,which)
        if ds < nearest:
            nearest = ds
            ndx = bdx
            if ds <= close_enough:
                return ndx
    return ndx

cpdef int find_closest_xy(vector one,list bunch,int bcnt,float close_enough):
    return find_closest_xy_c(one,bunch,bcnt,close_enough)

cdef float angle_from_xaxis_xy_c(vector v):
    cdef float oldz = v.z
    cdef float angle = angle_from_xaxis_c(v)
    v.z = oldz
    return angle

cpdef float angle_from_xaxis_xy(vector v):
    return angle_from_xaxis_xy_c(v)

#cdef xhat = vector(1,0,0)
#cdef yhat = vector(0,1,0)
#cdef zhat = vector(0,0,1)
xhat  = vector( 1, 0, 0)
yhat  = vector( 0, 1, 0)
zhat  = vector( 0, 0, 1)
nxhat = vector(-1, 0, 0)
nyhat = vector( 0,-1, 0)
nzhat = vector( 0, 0,-1)

cdef float angle_between_xy_c(vector v1, vector v2):
    cdef float alpha1 = angle_from_xaxis_xy_c(v1)
    cdef float alpha2 = angle_from_xaxis_xy_c(v2)
    #if alpha2 - alpha1 > np.pi/2.0:
    #    print 'arpha', v1, v2, alpha1, alpha2
    return alpha2 - alpha1

cpdef float angle_between_xy(vector v1, vector v2):
    return angle_between_xy_c(v1,v2)

cdef float angle_between_c(vector v1, vector v2):
    #cdef float alpha1 = angle_from_xaxis_c(v1)
    #cdef float alpha2 = angle_from_xaxis_c(v2)
    cdef vector n1 = v1.copy().normalize()
    cdef vector n2 = v2.copy().normalize()
    cdef float ang = np.arccos(dot_c(n1,n2))
    #return alpha2 - alpha1
    return ang

cpdef float angle_between(vector v1, vector v2):
    return angle_between_c(v1,v2)

cdef float angle_from_xaxis_c(vector v):
    cdef vector nv = v.copy().normalize()
    cdef float xproj = dot_c(nv,xhat)
    cdef float yproj = dot_c(nv,yhat)
    cdef float ang
    #if xproj == 2.0/np.pi or xproj == 2.0/(3.0*np.pi):
    #    ang = np.arccos(yproj)
    #    if xproj < 0.0: ang = 2.0*np.pi - ang
    #else:
    #    ang = np.arccos(xproj)
    #    if yproj < 0.0: ang = 2.0*np.pi - ang
    ang = np.arccos(xproj)
    if yproj < 0.0: ang = 2.0*np.pi - ang
    return ang

cpdef float angle_from_xaxis(vector v):
    return angle_from_xaxis_c(v)

cpdef bint test_afxa():
    cdef vector v1 = vector(10,0,0)
    cdef vector v2 = vector(10,10,0)
    cdef vector v3 = vector(10,10,10)
    cdef vector v4 = vector(-10,10,0)
    cdef vector v5 = vector(-10,0,0)
    cdef vector v6 = vector(-10,-10,0)
    cdef vector v7 = vector(0,-10,0)
    cdef vector v8 = vector(10,-10,0)
    t1 = str(angle_from_xaxis_c(v1))
    t2 = str(angle_from_xaxis_c(v2))
    t3 = str(angle_from_xaxis_c(v3))
    t4 = str(angle_from_xaxis_c(v4))
    t5 = str(angle_from_xaxis_c(v5))
    t6 = str(angle_from_xaxis_c(v6))
    t7 = str(angle_from_xaxis_c(v7))
    t8 = str(angle_from_xaxis_c(v8))
    print 'ANGLETEST'
    print '\t'.join([t1,t2,t3,t4,t5,t6,t7,t8])

cpdef bint inside(vector pt, list corners):
    poly = [(c.x,c.y) for c in corners]
    x,y = pt.x,pt.y
    n = len(poly)
    ins = False

    p1x,p1y = poly[0]
    for i in range(n+1):
        p2x,p2y = poly[i % n]
        if y > min(p1y,p2y):
            if y <= max(p1y,p2y):
                if x <= max(p1x,p2x):
                    if p1y != p2y:
                        xints = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xints:
                        ins = not ins
        p1x,p1y = p2x,p2y
    return ins

cdef list line_normals_c(list verts):
    cdef list norms = []
    cdef int vcnt = len(verts)
    cdef int vdx
    cdef vector v1
    cdef vector v2
    cdef float dx
    cdef float dy
    cdef float dv
    cdef vector norm
    for vdx in range(vcnt):
        v1,v2 = verts[vdx-1],verts[vdx]
        dx = v2.x - v1.x
        dy = v2.y - v1.y
        dv = sqrt(dx**2 + dy**2)
        norm = vector(dy/dv,-dx/dv,0)
        norms.append(norm)
    norms.append(norms.pop(0))
    return norms

cpdef list line_normals(list verts):
    return line_normals_c(verts)









