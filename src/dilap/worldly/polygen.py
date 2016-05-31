import dilap.core.base as db
import dilap.geometry.tools as gtl
from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat
import dilap.geometry.polymath as pym

import pdb



###############################################################################
### functions which create geometry in the form of polygons
###############################################################################

# return geometry describing a wall from p1 to p2 of height wh with holes hs
#   return polygons constituting a trimesh for the wall
#   return boundary polygons where portals will connect to the wall
def awall(p1,p2,wh,hs,hp1 = None,hp2 = None):
    if hp1 is None:hp1 = p1
    if hp2 is None:hp2 = p2
    wn = p1.tov(p2).nrm().crs(vec3(0,0,1))
    hp1.prj(p1,wn);hp2.prj(p2,wn)
    polys,portals = [],[]

    eb = [p1.cp(),p2.cp(),p2.cp().ztrn(wh),p1.cp().ztrn(wh)]
    ibs = []

    for hx in range(len(hs)):
        h = hs[hx]
        dx,dp,dw,dh,dz = h
        ddp = dw/(2.0*hp1.d(hp2))
        d1,d2 = hp1.lerp(hp2,dp-ddp),hp1.lerp(hp2,dp+ddp)
        if dh+dz < wh:
            h = [d2.cp().ztrn(dz),d2.cp().ztrn(dz+dh),
                d1.cp().ztrn(dz+dh),d1.cp().ztrn(dz)]
            portals.append(h)
            if dz == 0:
                for hp in h:eb.insert(1,hp)
            else:ibs.append(h)
        else:
            tdz,tdh = 0,wh
            h = [d2.cp().ztrn(tdz),d2.cp().ztrn(tdz+tdh),
                d1.cp().ztrn(tdz+tdh),d1.cp().ztrn(tdz)]
            portals.append(h)
            ebextra = [d2.cp()]+eb[1:-3]+\
                [eb[-3].cp(),eb[-3].cp().ztrn(wh),d2.cp().ztrn(wh)]
            wpy = (tuple(ebextra),())
            polys.append(wpy)
            eb = [p1.cp(),d1.cp(),d1.cp().ztrn(wh),p1.cp().ztrn(wh)]
    if eb:
        wpy = (tuple(eb),tuple(tuple(x) for x in ibs))
        polys.append(wpy)
    return polys,portals

# return polygons constituting a trimesh of a portal 
def aportal(oloop,p1,p2,ww):
    polys = []
    pn = vec3(0,0,1).crs(p1.tov(p2).nrm()).uscl(ww)
    iloop = [p.cp().trn(pn) for p in oloop]
    for lx in range(len(oloop)):
        polys.append(((iloop[lx-1],oloop[lx-1],oloop[lx],iloop[lx]),()))
    return polys

###############################################################################
###############################################################################
###############################################################################



