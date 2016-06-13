from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat
import dilap.geometry.tools as gtl
import dilap.geometry.polymath as pym

import dilap.core.plotting as dtl
import matplotlib.pyplot as plt

import math,numpy,random,pdb



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

# produce a highly irregular boundary polygon using 
#   pym.ebdxy (e.g. for landmasses)
# the result must fit within b without intersections
# the result must meet criterion: pym.bvalidxy(r) > 1
def ajagged(b,epsilon):
    j = random.randint(0,len(b)-1)
    l = b[j-1].d(b[j])
    t = random.uniform(0,2.0*numpy.pi)
    t = random.uniform(0,2.0*numpy.pi)
    n = random.randint(3,8)

    stamp = b[j-1].mid(b[j]).pring(l/2.0,n)
    q = quat(1,0,0,0).av(t,vec3(0,0,1))
    q.rotps(stamp)

    #xpj = vec3(1,0,0).prjps(b)
    #xsc = xpj[1]-xpj[0]
    #r = random.uniform(xsc/8.0-100,xsc/8.0+100)
    #r = xsc/8.0
    #stamp = vec3(xsc/2.0*math.cos(t),xsc/2.0*math.sin(t),0).pring(r,n)
    
    if pym.bintbxy(b,stamp,col = False,ie = False):
        nbs = pym.ebdxy(b,stamp,epsilon)
        nbas = [pym.bareaxy(nb) for nb in nbs]
        nb = nbs[nbas.index(max(nbas))]
        nb = pym.aggregate(nb,1)
        if pym.bvalidxy(nb) > 0:b = nb

        print('ISVALID',pym.bvalidxy(b),nbas,max(nbas))
        '''#
        ax = dtl.plot_axes_xy(600)
        for nb in nbs:
            ax = dtl.plot_polygon_xy(nb,ax,lw = 1,col = 'k')
        ax = dtl.plot_polygon_xy(b,ax,lw = 4,col = 'b')
        ax = dtl.plot_polygon_xy(stamp,ax,lw = 2,col = 'r')
        plt.show()
        '''#

    bval = pym.bvalidxy(b)
    if bval == -1:b.reverse()

    if not pym.bvalidxy(b) > 0:
        ax = dtl.plot_axes_xy(700)
        ax = dtl.plot_polygon_xy(b,ax,lw = 4,col = 'b')
        ax = dtl.plot_points_xy(b,ax,number = True)
        ax = dtl.plot_polygon_xy(stamp,ax,lw = 2,col = 'r')
        plt.show()
        #pdb.set_trace()
        raise ValueError
    
    return b

###############################################################################





