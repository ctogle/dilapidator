
def combine(models):
    final = models[0]
    if len(models) > 1:
        for m in models[1:]:
            final._consume(m)
    return final

def point_ring(r,n):
    st = cv.zero().translate_x(r)
    alpha = fu.PI*(2.0/n)
    points = []
    for x in range(n):
        points.append(st.copy().rotate_z(x*alpha))
    return points

def normal(c1,c2,c3):
    c1c2 = cv.v1_v2(c1,c2).normalize()
    c2c3 = cv.v1_v2(c2,c3).normalize()
    cn = c1c2.cross(c2c3).normalize()
    return cn

