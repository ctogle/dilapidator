import dilap.core.tools as dpr
import dilap.core.vector as dpv
import dilap.core.quaternion as dpq

import dilap.mesh.tools as dtl
import matplotlib.pyplot as plt
import pdb



def doorhole(l,wlw,dw,dh,lp):
    dhole = [
        dpv.vector(-dw/2.0, 0,0),dpv.vector(-dw/2.0,dh,0),
        dpv.vector( dw/2.0,dh,0),dpv.vector( dw/2.0, 0,0)]
    dpv.translate_coords_x(dhole,l*lp)
    return dhole

def doorframe(l,wlw,dw,dh,lp):
    ps = doorhole(l,wlw,dw,dh,lp)
    ps.extend([p.copy() for p in ps])
    dpv.translate_coords(ps[:4],dpv.z().scale_u(-wlw/2.0))
    dpv.translate_coords(ps[4:],dpv.z().scale_u( wlw/2.0))
    p0,p1,p2,p3,p4,p5,p6,p7 = ps
    fpolys = (
        ((p4.copy(),p0.copy(),p1.copy(),p5.copy()),()),
        ((p5.copy(),p1.copy(),p2.copy(),p6.copy()),()),
        ((p6.copy(),p2.copy(),p3.copy(),p7.copy()),()),
        ((p7.copy(),p3.copy(),p0.copy(),p4.copy()),()),
            )
    return fpolys

def windowhole(l,wlw,ww,wh,wz,lp):
    whole = [
        dpv.vector(-ww/2.0,   wz,0),dpv.vector(-ww/2.0,wh+wz,0),
        dpv.vector( ww/2.0,wh+wz,0),dpv.vector( ww/2.0,   wz,0)]
    dpv.translate_coords_x(whole,l*lp)
    return whole

def windowframe(l,wlw,ww,wh,wz,lp):
    ps = windowhole(l,wlw,ww,wh,wz,lp)
    ps.extend([p.copy() for p in ps])
    dpv.translate_coords(ps[:4],dpv.z().scale_u(-wlw/2.0))
    dpv.translate_coords(ps[4:],dpv.z().scale_u( wlw/2.0))
    p0,p1,p2,p3,p4,p5,p6,p7 = ps
    fpolys = (
        ((p4.copy(),p0.copy(),p1.copy(),p5.copy()),()),
        ((p5.copy(),p1.copy(),p2.copy(),p6.copy()),()),
        ((p6.copy(),p2.copy(),p3.copy(),p7.copy()),()),
        ((p7.copy(),p3.copy(),p0.copy(),p4.copy()),()),
            )
    return fpolys

# create polygons for a plc representing a wall from p1 to p2
def wall(p1,p2,wh1,wh2,ww,doors = (),windows = ()):
    l = dpv.distance(p1,p2)
    a = dpr.angle_from_xaxis_xy(dpv.v1_v2(p1,p2))
    ww2 = ww/2.0
    bnd = [
        dpv.vector(  ww2,  0,0),dpv.vector(l-ww2, 0,0),
        dpv.vector(l-ww2,wh2,0),dpv.vector( ww2,wh1,0)]
    for d in doors:
        dhp = doorhole(l,*d)
        for x in range(len(dhp)-1,-1,-1):
            bnd.insert(1,dhp[x])
    wholes = [tuple(windowhole(l,*w)) for w in windows]
    extras = []
    for d in doors:
        dpolys = doorframe(l,*d)
        for dp in dpolys:extras.append(dp)
    for w in windows:
        wpolys = windowframe(l,*w)
        for wp in wpolys:extras.append(wp)
    x = dpv.z().scale_u(ww2)
    main0 = (tuple(bnd),tuple(wholes))
    main1 = dpr.translate_polygon(dpr.copy_polygon(main0),x)
    dpr.translate_polygon(main0,x.flip())
    wallgons = (main0,main1)+tuple(extras)
    for wgon in wallgons:
        dpr.rotate_x_polygon(wgon,dpr.PI/2.0)
        dpr.rotate_z_polygon(wgon,a)
        dpr.translate_polygon(wgon,p1)
    return wallgons

def post(p,s,w,h):
    #base = dpr.point_ring(w/2.0,s)
    base = dpr.square(w,w)
    dpv.translate_coords(base,p)
    top = [b.copy().translate_z(h) for b in base]
    wls = (
        ((base[0].copy(),base[1].copy(),top[1].copy(),top[0].copy()),()),
        ((base[1].copy(),base[2].copy(),top[2].copy(),top[1].copy()),()),
        ((base[2].copy(),base[3].copy(),top[3].copy(),top[2].copy()),()),
        ((base[3].copy(),base[0].copy(),top[0].copy(),top[3].copy()),()),
            )
    return wls

# id like to answer several questions
# can these be done for full 3d or xy only?
# is a point on a line segment, including endpoints
# is a point inside a polygon, including edges
# is a point on a line segment, not including enpoints
# is a point inside a polygon, not including edges

# does the line segment interior given by s1,s2 intersect the polygon py
def intersect_segment_polygon(s1,s2,py):
    for ex in range(len(py)):
        ep1,ep2 = py[ex-1],py[ex]
        isect = dtl.segments_intersect_at(s1,s2,ep1,ep2)
        if not isect is None:
            if type(isect) == type(()):return True
            #    if dpv.distance(*isect) > 0.001:return True
    return False

# given two concave polygons with holes
# return one polygon if p1 and p2 can be merged
# or return None
def merge_two_polygons(p1,p2):

    import matplotlib.pyplot as plt
    import dilap.mesh.tools as dtl

    # return an ordered nonduplicative list of points from s1 to s2
    def segclean(s1,s2,isects):
        clean = [s1]
        if len(isects) < 2:
            clean.append(s2)
            return clean
        maxd = dpv.distance(s1,s2)
        iordr = list(dpv.proximity_order(s1,isects))
        while iordr:
            nxtx = iordr.pop(0)
            nxti = isects[nxtx]
            if not clean[-1].near(nxti):clean.append(nxti)
        if not clean[-1].near(s2):clean.append(s2)
        return clean

    # given a line segment and a polygon, produce
    # all points of intersection for edge segmentation
    def segsects(s1,s2,py):
        unn = False
        isects = []
        for pyx in range(len(py)):
            ep1,ep2 = py[pyx-1],py[pyx]
            isect = dtl.segments_intersect_at(s1,s2,ep1,ep2)
            if not isect is None:
                if type(isect) == type(()):
                    isects.append(isect[0])
                    isects.append(isect[1])
                    unn = True
                elif isect.near(s1) or isect.near(s2):pass
                elif isect.near(ep1) or isect.near(ep2):isects.append(isect)
                else:
                    isects.append(isect)
                    unn = True
        return isects,unn

    unionize = False
    eb1,eb2 = p1[0],p2[0]
    ebn1 = dpr.polygon_normal(eb1)
    pj1 = dpv.project_coords(list(eb1),ebn1)
    pj2 = dpv.project_coords(list(eb2),ebn1)
    nr = dpr.isnear
    if not (nr(pj1.x,pj1.y) and nr(pj2.x,pj2.y) and nr(pj1.x,pj2.x)):return

    if       ebn1.near(dpv.nz()):prot = dpq.q_from_av(dpr.PI,dpv.x())
    elif not ebn1.near(dpv.z() ):prot = dpq.q_from_uu(ebn1,dpv.z())
    else:                        prot = dpq.zero()
    dpr.rotate_polygon(p1,prot)
    dpr.rotate_polygon(p2,prot)
    unfinished = []
    for e1x in range(len(eb1)):
        ep11,ep12 = eb1[e1x-1],eb1[e1x]
        isects,union = segsects(ep11,ep12,eb2)
        if union:unionize = True
        isects = segclean(ep11,ep12,isects)
        for x in range(1,len(isects)):
            iseg = (isects[x-1],isects[x])
            if not intersect_segment_polygon(iseg[0],iseg[1],eb2):
                unfinished.append(iseg)
    for e2x in range(len(eb2)):
        ep21,ep22 = eb2[e2x-1],eb2[e2x]
        isects,union = segsects(ep21,ep22,eb1)
        if union:unionize = True
        isects = segclean(ep21,ep22,isects)
        for x in range(1,len(isects)):
            iseg = (isects[x-1],isects[x])
            if not intersect_segment_polygon(iseg[0],iseg[1],eb1):
                unfinished.append(iseg)
    if not unionize:
        prot.flip()
        dpr.rotate_polygon(p1,prot)
        dpr.rotate_polygon(p2,prot)
        return

    finished = []
    last = unfinished[0][0]
    while unfinished:
        for ufx in range(len(unfinished)):
            unfin = unfinished[ufx]
            if   last.near(unfin[0]):break
            elif last.near(unfin[1]):break
        unfin = unfinished.pop(ufx)
        which = unfin[1] if last.near(unfin[0]) else unfin[0]
        finished.append(which)
        last = which
        
    ibnds = p1[1]+p2[1]
    merged = (tuple(finished),tuple(ibnds))

    prot.flip()
    dpr.rotate_polygon(merged,prot)
    return merged

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
                    mp = merge_two_polygons(i1,i2)
                    if not mp is None:
                        islands.remove(i1)
                        islands.remove(i2)
                        islands.append(mp)
                        break
                if not mp is None:break
                if x == icnt-1:checked = True
    return islands

def mergetest():
    p1 = tuple(dpr.square(5,5))
    p2 = tuple(dpr.square(5,5,dpv.vector(4,0,0)))

    isect = dpr.concaves_intersect(p1,p2)

    print('isect',isect)
    ax = dtl.plot_axes_xy()
    dtl.plot_polygon_xy(list(p1),ax)
    dtl.plot_polygon_xy(list(p2),ax)
    plt.show()

    pdb.set_trace()



