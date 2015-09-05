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
    #wallgons = (main1,)+tuple(extras)
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





