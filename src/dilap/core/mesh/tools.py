import matplotlib.pyplot as plt

def plot_points_xy(points,ax = None):
    if ax is None:ax = plt.figure().add_subplot(111)
    for pdx in range(len(points)):
        p = points[pdx]
        ax.plot([p.x],[p.y],marker = 'o')
    return ax

def plot_points(points,ax = None):
    if ax is None:ax = plt.figure().add_subplot(111,projection = '3d')
    for pdx in range(len(points)):
        p = points[pdx]
        ax.plot([p.x],[p.y],zs = [p.z],marker = 'o')
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

def plot_polygon_xy(points,ax = None):
    epts = points[:]
    epts.append(points[0])
    plot_edges_xy(epts,ax)

def plot_polygon(points,ax = None):
    epts = points[:]
    epts.append(points[0])
    plot_edges(epts,ax)

def plot_tetrahedron(points,ax = None):
    raise NotImplemented

# if a is within c of b, return b
# else return a
def near(a,b,c = 0.000001):
    if abs(a-b) < c:return b
    else:return a

def orient2d(a,b,c):
    m11,m12 = a.x-c.x,a.y-c.y
    m21,m22 = b.x-b.x,b.y-b.y
    return near(m11*m22-m12*m21,0)

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
    return near(det1 - det2 + det3,0)

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
    return near(det1 - det2 + det3 - det4,0)


