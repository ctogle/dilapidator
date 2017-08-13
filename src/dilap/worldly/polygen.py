from dilap.core import *
from dilap.geometry import *
import dilap.geometry.tools as gtl
import dilap.geometry.polymath as pym
import math
import numpy
import random
import pdb


# return geometry describing a wall from p1 to p2 of height wh with holes hs
#   return polygons constituting a trimesh for the wall
#   return boundary polygons where portals will connect to the wall
def awall(p1,p2,wh,hs,hp1 = None,hp2 = None):
    if hp1 is None:hp1 = p1
    if hp2 is None:hp2 = p2
    wn = p1.tov(p2).nrm().crs(vec3(0,0,1))
    hp1.prj(p1,wn);hp2.prj(p2,wn)

    if isinstance(hp1,vec3) and hp1.isnear(hp2):
        # the prj is the problem!!
        pdb.set_trace()
        raise ValueError('go mow the yard sucka')

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
    n = random.randint(3,8)

    stamp = b[j-1].mid(b[j]).pring(l/2.0,n)
    q = quat(1,0,0,0).av(t,vec3(0,0,1))
    q.rotps(stamp)
    
    if pym.bintbxy(b,stamp,col = False,ie = False):
        nbs = pym.ebdxy(b,stamp,epsilon)
        nbas = [pym.bareaxy(nb) for nb in nbs]
        nb = nbs[nbas.index(max(nbas))]
        nb = pym.aggregate(nb,1)
        if pym.bvalidxy(nb) > 0:b = nb

    bval = pym.bvalidxy(b)
    if bval == -1:b.reverse()
    if not pym.bvalidxy(b) > 0:
        ax = plot_axes_xy(700)
        ax = plot_polygon_xy(b,ax,lw = 4,col = 'b')
        ax = plot_points_xy(b,ax,number = True)
        ax = plot_polygon_xy(stamp,ax,lw = 2,col = 'r')
        plt.show()
        #pdb.set_trace()
        raise ValueError
    
    return b


def chunk(b,epsilon,lscl = 1.0,j1 = None,j2 = None,edge = False):
    if j1 is None:j1 = random.randint(0,len(b)-1)
    if j2 is None:j2 = j1-1
    l = b[j1].d(b[j2])*lscl
    t = random.uniform(0,2.0*numpy.pi)
    n = random.randint(3,8)

    if edge:
        stamp = b[j1].mid(b[j2]).pring(l/2.0,n)
    else:
        stamp = vec3(
            random.randint(int(l/4),int(l)),
            random.randint(int(l/4),int(l)),
            0).com(b).pring(l/2.0,n)
    q = quat(1,0,0,0).av(t,vec3(0,0,1))
    q.rotps(stamp)
    
    if pym.binbxy(stamp,b):
        return stamp
    elif pym.bintbxy(b,stamp,col = False,ie = False):
        nbs = pym.ebixy(b,stamp,epsilon)
        nbas = [pym.bareaxy(nb) for nb in nbs]
        nb = nbs[nbas.index(max(nbas))]
        nb = pym.aggregate(nb,1)
        if pym.bvalidxy(nb) > 0:
            return nb
    else:
        ax = plot_axes_xy(200)
        ax = plot_polygon_xy(b,ax,col = 'r',lw = 2)
        ax = plot_polygon_xy(stamp,ax,col = 'b',lw = 2)
        plt.show()
        raise ValueError


def splotch(vb,easement = 10):
    # WIP
    vbes = [(vb[0][x-1],vb[0][x]) for x in range(len(vb[0]))]
    vbels = [vb[0][x-1].d(vb[0][x]) for x in range(len(vb[0]))]
    vbe_primary = vbels.index(max(vbels))
    vbe_primary_tn = vb[0][vbe_primary-1].tov(vb[0][vbe_primary]).nrm()
    bq = quat(0,0,0,1).uu(vec3(1,0,0),vbe_primary_tn)
    if not bq:bq = quat(1,0,0,0)

    #vbcon = pym.aggregate(pym.contract(vb[0],easement),5)
    vbcon = pym.aggregate(pym.contract(vb[0],easement))

    vbxpj = vbe_primary_tn.prjps(vb[0])
    vbypj = vbe_primary_tn.cp().zrot(gtl.PI2).prjps(vb[0])
    vbxl,vbyl = vbxpj[1]-vbxpj[0],vbypj[1]-vbypj[0]

    dx,dy = 20,20
    xn,yn = int(1.5*vbxl/dx),int(1.5*vbyl/dy)
    print('xn,yn',xn,yn)
    o = vec3(0,0,0).com(vbcon)
    vgrid = [o.cp().xtrn(dx*(x-(xn/2.0))).ytrn(dy*(y-(yn/2.0)))
        for x in range(xn) for y in range(yn)]
    #vgrid = [p for p in vgrid if abs(p.x-10) > 10]
    boxes = [p.cp().trn(o.flp()).rot(bq).trn(o.flp()).sq(dx,dy) for p in vgrid]

    for b in boxes:
        bcom = vec3(0,0,0).com(b).flp()
        bcom.trnps(b)
        bq.rotps(b)
        bcom.flp().trnps(b)
    boxes = [b for b in boxes if pym.binbxy(b,vbcon)]
    for ib in vb[1]:
        boxes = [b for b in boxes if not pym.bintbxy(b,ib) 
            and not b[0].inbxy(ib) and not ib[0].inbxy(b[0])]

    blgfps = pym.ebuxy_special(boxes,5,4)

    ps = [vec3(0,0,0).com(blgfp) for blgfp in blgfps]
    qs = [bq.cp() for blgfp in blgfps]
    ss = [vec3(1,1,1) for blgfp in blgfps]

    for p,q,s,blgfp in zip(ps,qs,ss,blgfps):
        p.cpxy().flp().trnps(blgfp)
        q.cpf().rotps(blgfp)

    return zip(ps,qs,ss,blgfps)


def lsystopg(b,lsys,e = 2):
    ''''''
    # lsystem -> planargraph
    pg = planargraph()
    leaves = []
    for piece in lsys:
        if isinstance(piece,tuple):
            p1,p2 = piece
            #v1,v2 = pg.fp(p1,10),pg.fp(p2,10)
            v1,v2 = pg.fp(p1,e),pg.fp(p2,e)
            e12 = pg.fe(v1,v2)
        elif isinstance(piece,vec3):
            leaves.append(piece)
            pass
    return pg, leaves


def lsystopy(b,lsys,e = 2):
    '''Return an outline of a fractal scaled to boundary polygon'''
    # lsystem -> planargraph
    pg = planargraph()
    for piece in lsys:
        if isinstance(piece,tuple):
            p1,p2 = piece
            v1,v2 = pg.fp(p1,10),pg.fp(p2,10)
            e12 = pg.fe(v1,v2)
        elif isinstance(piece,vec3):
            pass
    # planargraph -> polygon -> smooth -> pinch -> fit
    py = pym.pgtopy(pg,5)[0]
    py = pym.smoothxy(py,0.5,2,0)
    py = pym.pinchb(py,5)[0]
    py = pym.bfitbxy(py,b)
    # return the intersection of the outline and the b
    #py = pym.ebixy(b,py)[0]
    #py = pym.pinchb(py,5)[0]
    return py


def blgfpFFFF(b, e, s=4):
    '''Create a full polygon representing the footprint of a building'''
    entrance = b[0].mid(b[1])
    enttn = b[0].tov(b[1])
    entnm = enttn.crs(vec3(0,0,1)).nrm()
    entln = enttn.mag()
    seedlo = entrance.cp().trn(enttn.cp().nrm().uscl( s))
    seedro = entrance.cp().trn(enttn.cp().nrm().uscl(-s))
    seedli = seedlo.cp().trn(entnm.cp().uscl(-s))
    seedri = seedro.cp().trn(entnm.cp().uscl(-s))
    seed = (seedlo, seedro, seedri, seedli)

    pg = planargraph()
    fairgame = []
    for x in range(len(seed)):
        s1 = pg.fp(seed[x-1],1,1)
        s2 = pg.fp(seed[x],1,1)
        e12 = pg.fe(s1,s2)
        fairgame.append((s1,s2))

    while fairgame:
        which = random.choice(fairgame)
        w1,w2 = which
        p1,p2 = pg.vs[w1][1]['p'],pg.vs[w2][1]['p']
        whichtn = p1.tov(p2)
        whichnm = whichtn.crs(vec3(0,0,1)).nrm().uscl(-s)
        p3,p4 = p2.cp().trn(whichnm),p1.cp().trn(whichnm)
        
        btest = pym.bintbxy((p1, p2, p3, p4), b)

        if btest:
            fairgame.remove(which)
        else:
            n1 = pg.fp(p3,1,1)
            n2 = pg.fp(p4,1,1)
            e23 = pg.fe(w2,n1)
            e34 = pg.fe(n1,n2)
            e41 = pg.fe(n2,w1)
            fairgame.append((w2,n1))
            fairgame.append((n1,n2))
            fairgame.append((n2,w1))

        print('btest',btest)

        #edge = (fpb[which-1],fpb[which])

        # can i operate on which
        # if not; remove from fairgame
        # else; which operations can i perform
        #   is there space to grow the edge?
        #   modify the edge
        #     translate, extend, trim
        #   grow graph from edge

        ax = pg.plotxy(l = 50)
        ax = plot_polygon_xy(b, ax, col = 'r')
        ax = plot_polygon_xy((p1, p2, p3, p4), ax, col = 'g')
        ax = plot_polygon_xy((seedlo, seedro), ax, col = 'b', lw = 4)
        plt.show()

        pdb.set_trace()





    def clean(py):
        es = ((py[x-1],py[x]) for x in range(len(py)))
        py = pym.pgtopy(planargraph.segstopg(es,1),0)[0]
        return py


    fpb,fphs = [s.cp() for s in seed],[]
    fairgame = list(range(1,len(fpb)))
    while True:
        which = random.choice(fairgame)
        edge = (fpb[which-1],fpb[which])
        edgenm = edge[1].tov(edge[0]).crs(vec3(0,0,1)).nrm()
        tip = (edge[0].cp().trn(edgenm.cp().uscl(s)),
               edge[1].cp().trn(edgenm.cp().uscl(s)), )
        fpb.insert(which,tip[1])
        fpb.insert(which,tip[0])

        fpb = clean(fpb)
        if not pym.bnrm(fpb).z > 0:
            fpb.reverse()

        ax = plot_axes_xy(100)
        ax = plot_polygon_xy(b, ax, col = 'r')
        ax = plot_polygon_xy(fpb, ax, col = 'g', lw = 2)
        for h in fphs:
            ax = plot_polygon_xy(h, ax, col = 'm', lw = 2)
        ax = plot_polygon_xy(seed, ax, col = 'b', lw = 4)
        ax = plot_edges_xy(edge, ax, col = 'c', lw = 4)
        ax = plot_edges_xy(tip, ax, col = 'y', lw = 4)
        plt.show()

    
    return (fpb,fphs)


def blgfp(b, e, s=4):
    '''Create a full polygon representing the footprint of a building'''

    i,p,d = 3,vec3(0,0,0),vec3(0,1,0)
    #axiom,rules = 'X',dict([('X','{[[X}{]XF}]X'),('F','FA'),('A','F')])
    axiom,rules = 'X',dict([('X','{[[X}{]F}]X'),('F','FA'),('A','F')])
    params = dict(dazimuthal = gtl.rad(90),drho = 16)
    lsys = lsystem(i, p, d, axiom, rules, **params)


    stems, leaves = lsystopg(b, lsys, 2)
    py = pym.pgtopy(stems,7)
    #py[0] = pym.pinchb(py[0],5)[0]


    ax = plot_axes_xy(100)
    ax = plot_polygon_xy(py[0], ax, col = 'b')
    ax = plot_points_xy(py[0], ax, number = True)
    for h in py[1]:
        ax = plot_polygon_xy(h, ax, col = 'r')
    ax = stems.plotxy(ax)
    for leaf in leaves:
        ax = plot_point_xy(leaf, ax, col = 'b')
    ax = plot_polygon_xy(b, ax, col = 'r', ls = '--')
    plt.show()

    eds = pym.bdistpxy(b,p)

    #feat = lsystopy(b, lsys, e)
    #vec3(1,1,0).sclps(feat)


    return (py[0],())


