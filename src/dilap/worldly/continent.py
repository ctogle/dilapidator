from .polygen import lsystopy
from .terrain import terrain
from .infrastructure import roadmap
from .partitiongraph import world
from dilap.geometry import *
import dilap.geometry.polymath as pym
import dilap.geometry.planarmesh as pme
import dilap.geometry.tools as gtl

from dilap.core.scenegraph import scenegraph
from dilap.worldly.polygen import boundaries

from dilap.core import *
import numpy
import random
import math
import pdb


def land_model(b, hs=[], z=(lambda x, y: p.z)):
    m = model()
    tm = m.agfxmesh()
    ngvs = m.asurf((b, hs), tm, fm='grass2',
        ref=True, hmin=8, hmax=8, minhmin=0.1,
        zfunc=z, rv=pym.bnrm(b).z < 0, uvstacked=None, autoconnect=True, e=0.001)
    lockf = lambda p: p.onpxy((b, hs))
    m.subdiv(tm, False, True, lockf)
    m.uvs(tm)
    return m

def ocean_model(b, hs=[], z=10):
    m = model()
    tm = m.agfxmesh()
    ngvs = m.asurf((b, hs), tm, fm='generic',
                   ref=False, hmin=100, zfunc=(lambda x, y: z),
                   rv=pym.bnrm(b).z < 0)
    return m

def road_models(b, hs, z):
    m = model()
    tm = m.agfxmesh()
    ngvs = m.asurf((b, hs), tm, fm='concrete1', ref=True, hmin=8, hmax=8, minhmin=0.1,
                   zfunc=z, rv=pym.bnrm(b).z < 0, uvstacked=None, autoconnect=True, e=0.001)
    lockf = lambda p: p.onpxy((b, hs))
    m.subdiv(tm, False, True, lockf)
    m.uvs(tm)
    return [m]

def land_and_sea(t_b, fp, z_f, sealevel=1, sg=None):
    fp_b, fp_hs = fp
    lands = [land_model(t_b, [fp_b], z_f)]
    lands.extend([land_model(ifp, [], z_f) for ifp in fp_hs])
    waters = [ocean_model(t_b, [], sealevel)]
    roads = road_models(fp_b, fp_hs, z_f)
    sg = sg if sg else scenegraph()
    sgv = sg.avert(None, None, None, models=lands, parent=sg.root)
    sgv = sg.avert(None, None, None, models=waters, parent=sg.root)
    sgv = sg.avert(None, None, None, models=roads, parent=sg.root)
    return sg





def terraininput(b,e):


    def regionize(boundary, max_unsplittable_size, show=False):
        # create set of convex non-overlapping polygons touching b
        wide_boundary = pym.contract(boundary, -e*100)
        tb = pym.splitb(wide_boundary, 20)

        i = 3
        p,d = vec3(0,0,0),vec3(0,1,0)
        axiom = 'X'
        rules = dict([('X','F]{{X}[X}[F{[FX}]X'),('F','FF')])
        params = dict(dazimuthal=gtl.rad(25), drho=32)
        lsys = lsystem(i, p, d, axiom, rules, **params)
        pg = planargraph()
        for piece in lsys:
            if isinstance(piece, tuple):
                p1, p2 = piece
                v1 = pg.fp(p1, 16, w=0.5)
                v2 = pg.fp(p2, 16, w=0.5)
                e12 = pg.fe(v1, v2)
            elif isinstance(piece, vec3):
                pass
        pg.smooth_sticks()
        pg.fitbxy(wide_boundary, 1.1)
        py = pym.pgtopy(pg, 10, epsilon=1, z=vec3(0,0,1))
        smooth = lambda b: pym.smoothxy(pym.splitb(b, 20), 0.5)
        py = [pym.smoothxy(smooth(py[0]), 0.5), [smooth(h) for h in py[1]]]

        waters = [py]
        lands = pym.ebdxy(tb, py[0])
        islands = py[1][:]

        highs, lows = [], []
        highs.extend([l for l in lands if pym.bareaxy(l) > ((e * 100) ** 2)])
        #highs.extend([i for i in islands if pym.bareaxy(i) > ((e * 100) ** 2)])
        lows.append(py[0])

        print('regions', len(lands), len(islands), len(waters))
        if show:
            ax = plot_axes_xy(300)                 
            ax = plot_polygon_xy(wide_boundary, ax, lw=2, col='m', mk='o')
            ax = plot_polygon_xy(boundary, ax, lw=2, col='r', mk='o')
            #for land in lands:
            #    ax = plot_polygon_xy(land, ax, lw=3, col='g', mk='o')
            #for water in waters:
            #    ax = plot_polygon_full_xy(water, ax, lw=2, col='b', mk='o')
            #for island in islands:
            #    ax = plot_polygon_xy(island, ax, lw=3, col='g', mk='o')
            for l in lows:
                ax = plot_polygon_xy(l, ax, lw=3, col='b', mk='o')
            for h in highs:
                ax = plot_polygon_xy(h, ax, lw=3, col='g', mk='o')
            plt.show()

            for land in lands:
                ax = plot_axes_xy(300)                 
                ax = plot_polygon_xy(wide_boundary, ax, lw=2, col='m')
                ax = plot_polygon_xy(boundary, ax, lw=2, col='r')
                for island in islands:
                    ax = plot_polygon_xy(island, ax, lw=3, col='y')
                for water in waters:
                    ax = plot_polygon_full_xy(water, ax, lw=2, col='b')
                ax = plot_polygon_xy(land, ax, lw=5, col='g')
                plt.show()

        return [(pym.bareaxy(h), h) for h in highs], [(pym.bareaxy(l), l) for l in lows]

    def map_regions(boundary, max_unsplittable_size, show=False):
        # create set of convex non-overlapping polygons touching b
        wide_boundary = pym.contract(boundary, -e*100)

        #pm = pme.planarmesh()
        #pm.al(wide_boundary)
        #pm.pattern_split()
        #pm.plotxy(l=400)
        #plt.show()

        unfinished = [tuple(p.cp() for p in wide_boundary)]
        unfinished = [(pym.bareaxy(py), py) for py in unfinished]
        finished = []
        while unfinished:
            a, py = unfinished.pop(0)
            if a < max_unsplittable_size:
                finished.append((a, py))
            else:
                v = vec3(random.uniform(-1, 1), random.uniform(-1, 1), 0).nrm()
                l, r = pym.vsplitb(v, py)
                if pym.bnrm(l).z < 0:
                    l.reverse()
                if pym.bnrm(r).z < 0:
                    r.reverse()
                unfinished.append((pym.bareaxy(l), l))
                unfinished.append((pym.bareaxy(r), r))
        dl = 5
        #finish = lambda a, b: pym.smart_contract(pym.splitb(b, dl * 3), dl)
        #finished = [(a, finish(a, b)) for a, b in finished]
        print('regions')
        if show:
            ax = plot_axes(300)                 
            ax = plot_polygon_xy(wide_boundary, ax, lw=2, col='m')
            ax = plot_polygon_xy(boundary, ax, lw=2, col='b')
            for a, fpy in finished:
                ax = plot_polygon_xy(fpy, ax, lw=2, col='g')
            for a, upy in unfinished:
                ax = plot_polygon_xy(upy, ax, lw=2, col='r')
            plt.show()
        return finished

    def designate_regions(regions, feature_scale, show=True):
        i = 3
        p,d = vec3(0,0,0),vec3(0,1,0)
        axiom = 'X'
        rules = dict([('X','F]{{X}[X}[F{[FX}]X'),('F','FF')])
        params = dict(dazimuthal=gtl.rad(25), drho=20)
        lsys = lsystem(i, p, d, axiom, rules, **params)

        pg = planargraph()
        for piece in lsys:
            if isinstance(piece, tuple):
                p1, p2 = piece
                v1 = pg.fp(p1, 16, w=0.5)
                v2 = pg.fp(p2, 16, w=0.5)
                e12 = pg.fe(v1, v2)
            elif isinstance(piece, vec3):
                pass
        pg.smooth_sticks()

        highs = regions[:]
        lows = []
        # find regions which touch the pg!
        for v in pg.vs:
            p = v[1]['p']
            while True:
                for j, (a, region) in enumerate(highs):
                    if p.inbxy(region) or p.onbxy(region):
                        lows.append(highs.pop(j))
                        break
                else:
                    break
        py = pym.pgtopy(pg, 10, epsilon=1, z=vec3(0,0,1))

        newlows = []
        for a, l in lows:
            newlows.extend(pym.ebdxy(l, py[0]))

        if show:
            ax = plot_axes(300)
            ax = pg.plotxy(ax)
            for a, h in highs:
                ax = plot_polygon_xy(h, ax, lw=2, col='g')
            for l in newlows:
                ax = plot_polygon_xy(l, ax, lw=4, col='r')
            for a, l in lows:
                ax = plot_polygon_xy(l, ax, lw=2, col='r')
            ax = plot_polygon_full_xy(py, ax, lw=3, col='b')
            ax = plot_polygon_xy(py[0], ax, lw=5, col='m')
            plt.show()

        return highs, [(pym.bareaxy(b), b) for b in newlows]

    # do all this stuff with a planarmesh instance...?

    feature_scale = 0.02
    wide_boundary = pym.contract(b, -e*100)
    max_unsplittable_size = feature_scale*pym.bareaxy(wide_boundary)

    #regions = map_regions(b, max_unsplittable_size, True)
    #highs, lows = designate_regions(regions, max_unsplittable_size)
    highs, lows = regionize(b, max_unsplittable_size, False)

    ax = plot_axes(300)
    ax = plot_polygon_xy(wide_boundary,ax,lw = 2,col = 'm')
    ax = plot_polygon_xy(b,ax,lw = 2,col = 'r')
    for a, l in lows:
        ax = plot_polygon(l,ax,col = 'b')
    for a, h in highs:
        ax = plot_polygon(h,ax,col = 'g')
    plt.show()

    print('figure out highs and lows')

    scale = 10.0
    steps = [
        ([h[1] for h in highs],  scale, ((0.1 , 2), (0.5 , 3), (0.1, None))), 
        #(lows , -scale, ((0.08, 5), (0.02, 3), (0.1, 5), (0.2, None))), ]
        ([l[1] for l in lows], -scale, ((0.2, None), )), ]
    sealevel = 9
    return (wide_boundary, e, steps, sealevel)




    def feature(b):
        i,p,d = 6,vec3(0,0,0),vec3(0,1,0)
        axiom,rules = 'X',dict([('X','{[[X}{]X}F]X'),('F','FA'),('A','F')])
        params = dict(dazimuthal = gtl.rad(25.7),drho = 20)
        lsys = lsystem(i, p, d, axiom, rules, **params)
        feat = lsystopy(b, lsys, e)
        vec3(1,1,0).sclps(feat)
        return feat

    def scuttle(f,b):
        fnd = True
        while fnd:
            fnd = False
            for j in range(len(f)):
                p1,p2,p3 = f[j-2],f[j-1],f[j]
                if not (p2.inbxy(b) or p2.onbxy(b)):
                    if not (p1.inbxy(b) or p3.inbxy(b)):
                        if not pym.sintbxy(p1,p3,b):
                            fnd = True
                            break
            if fnd:
                f.pop(j-1)
        es = [(f[k-1],f[k]) for k in range(len(f))]
        return pym.sloops(es,2)


    #l,r = pym.vsplitb(vec3(0,1,0),b)
    #tl,bl = pym.vsplitb(vec3(1,0,0),l)
    #tr,br = pym.vsplitb(vec3(1,0,0),r)

    #trnsclf = lambda b : vec3(0,0,0).com(b).uscl(0.1).trnps(vec3(10,10,0).sclps(feature(b)))

    #highs = [feature(tb)]

    #highs = []
    #for b in (tl,bl,tr,br):
    #    highs.append(feature(b,highs))

    #highs = [pym.bfitbxy(feature(tl),tb)]
    #highs = [pym.bisectb(feature(tb))]
    highs = [feature(tb)]
    #highs = [x for y in highs for x in y]
    #highs = [trnsclf(feature(tl))]
    #highs = [feature(tl),feature(bl),feature(tr),feature(br)]
    #highs = [feature(tl),feature(br)]

    #highs = [scuttle(trnsclf(b),tb) for b in (tl,bl,tr,br)]
    #highs = [trnsclf(b) for b in (tl,bl,tr,br)]
    #highs = [trnsclf(b) for b in (tl,br)]

    #highs = [vec3(0,0,0).com(tl).pring(150,8),vec3(0,0,0).com(br).pring(150,8)]
    #highs = pym.bsuxy(highs,e)

    # HIGHS SHOULD BE A LIST OF DISJOINT POLYGONS 

    #lows = [tb]
    #for h in highs:
    #    newlows = []
    #    for l in lows:
    #        newlows.extend(pym.ebdxy(tb, h))
    #    lows = newlows

    #lows = pym.ebdxy(tb,pym.bfitbxy(highs[0],tb))
    #valleys = []
    #lows.extend(valleys)

    lows  = [pym.ebdxy(tb, h) for h in highs]
    lows  = [x for y in lows for x in y]

    ax = plot_axes(300)
    ax = plot_polygon_xy(tb,ax,lw = 2,col = 'm')
    ax = plot_polygon_xy(b,ax,lw = 2,col = 'r')
    for l in lows:
        ax = plot_polygon(l,ax,col = 'b')
    for h in highs:
        ax = plot_polygon(h,ax,col = 'g')
    plt.show()

    scale = 10.0
    steps = [
        (highs,  scale, ((0.1 , 2), (0.5 , 3), (1.0, None))), 
        #(lows , -scale, ((0.08, 5), (0.02, 3), (0.1, 5), (0.2, None))), ]
        (lows , -scale, ((0.2, None), )), ]
    sealevel = 9
    return (tb, e, steps, sealevel)


# fix blender plugin!!!
# requirements.txt
# fix blender plugin!!!
# requirements.txt
# fix blender plugin!!!
# requirements.txt
# fix blender plugin!!!
# requirements.txt
# fix blender plugin!!!
# requirements.txt
# fix blender plugin!!!
# requirements.txt


def continent(b,e, 
              includeocean=True, 
              includeland=True, 
              includeroads=False, 
              includebuildings=False, **kws):
    '''Create a partition of the world and generate a scenegraph
    Given a boundary and/or epsilon (can determine epsilon from boundary), 
    create a terrain object; from that create a roadmap; from those create
    a partition; vertices in the partition can be deterministically 
    represented by models in a scenegraph (self)'''
    if includeland:
        print('generating terrain')
        #targs = terraininput(b, e)
        #t = terrain(*targs,**kws)
        t = terrain.from_boundary(b, e)
    else:
        print('skipping terrain')
        t = None

    r = roadmap(t,e,**kws) if includeroads else None
    a = world(b,t,r,e,**kws)
    sg = scenegraph()

    sealevel = -10
    if includeland:
        for v,depth in t.enum():
            if hasattr(v,'style') and 'sealevel' in v.style:
                sealevel = v.loop[0].z
                break
    if includeocean:
        ocean(sg,b,sealevel)

    '''#
    for v,depth in a.enum():
        if 'world' in v.style:pass
        elif 'ocean' in v.style:
            #ocean(sg,b,sealevel)
            pass
        elif 'land' in v.style:
            pass
            #natural(w,v,a)
        elif 'bounded' in v.style:
            pass
            #developed(w,v,a)
        elif 'infrastructure' in v.style:
            pass
            #infrastructure(sg,v,t,a)
    '''#

    if includeland:
        for v,depth in a.enum():
            if 'land' == v.style or 'bounded' == v.style:
                print('generate land for style: \'%s\'' % v.style)
                land(sg,v,t.interpolate,a)

    return sg


def deform(b, e, dr_0=50):
    sb = pym.splitb(b, dr_0)
    nb = []
    for p in sb:
        n = p.cpxy().nrm()
        theta = 2 * numpy.arctan(n.y / n.x)
        dr = p.mag() + ((dr_0 / 2.0) * numpy.sin(theta))
        nb.append(n.uscl(dr))
    return nb


def ocean(sg,loop,sealevel,depth = 2.0,bleed = 5.0):
    ### create an unrefined flat surface for the bottom of the ocean
    ###  and an unrefined flat surface with contracted holes for the 
    ###  surface of the ocean
    m = model()
    sgv = sg.avert(None,None,None,models = [m],parent = sg.root)
    tm = m.agfxmesh()
    #holes = [pym.contract([p.cp() for p in ib],bleed) for ib in holes]
    holes = []
    gb = [[p.cp().ztrn(sealevel) for p in loop],holes]
    wb = [[p.cp().ztrn(sealevel+depth) for p in gb[0]],
         [[p.cp().ztrn(sealevel+depth) for p in lh] for lh in holes]]
    ngvs = m.asurf(gb,tm,fm = 'generic',ref = False,hmin = 100,zfunc = None)
    ngvs = m.asurf(wb,tm,fm = 'generic',ref = False,hmin = 100,zfunc = None)


def infrastructure(sg,v,t,a):
    print('infrastructure vertex',v.style,len(v.loop),len(v.holes))
    m = model()
    sgv = sg.avert(None,None,None,models = [m],parent = sg.root)
    tm = m.agfxmesh()

    rh = 0.0
    vb = [b.cp().ztrn(rh) for b in v.loop]
    vibs = [[b.cp().ztrn(rh) for b in h] for h in v.holes]
    vb = vstitch(vb,v,a)
    vibs = [vstitch(vib,v,a) for vib in v.holes]
    v.loop = vb
    v.holes = vibs

    ngvs = m.asurf((vb,vibs),tm,
        fm = 'concrete1',ref = True,hmin = 100,zfunc = t,
        uvstacked = None,autoconnect = True)
    lockf = lambda p : p.onpxy((v.loop,v.holes)) 
    m.subdiv(tm,False,True,lockf)

    return m

    rh = 0.0
    vb = [b.cp().ztrn(rh) for b in v[1]['b'][0]]
    vibs = [[b.cp().ztrn(rh) for b in v[1]['b'][1][x]] 
        for x in range(len(v[1]['b'][1]))]
    vb = vstitch(vb,v,a)
    vibs = [vstitch(vib,v,a) for vib in v[1]['b'][1]]
    v[1]['b'] = vb,vibs

    ngvs = m.asurf((vb,vibs),tm,
        fm = 'concrete1',ref = True,hmin = 100,zfunc = v[1]['tmesh'],
        uvstacked = None,autoconnect = True)
    lockf = lambda p : p.onpxy(v[1]['b']) 
    m.subdiv(tm,False,True,lockf)
    #m.subdiv(tm,False,True,lockf)
    #m.subdiv(tm,False,True,lockf)
    #m.uvs(tm)
    print('generated infrastructure')
    return m


def land(sg,v,z,a):
    m = model()
    sgv = sg.avert(None,None,None,models = [m],parent = sg.root)
    tm = m.agfxmesh()
    #vstitch isnt being used!!!
    # need to enforce that refinement is no problem for stitching!!!
    #  need to stitch from holes to interior regions!!!

    vb = [b.cp() for b in v.loop]
    #vibs = [[b.cp() for b in h] for h in v.holes]
    vibs = []

    #vb = stitchareas(vb,v,a)
    #vibs = [stitchareas(h,v,a) for h in v.holes]
    #v.loop,v.holes = vb,vibs

    # WHY DO I GET COLINEAR TRIANGLES WITH HMIN == 4 (did i fix this?)
    print('generating land')
    ngvs = m.asurf((vb,vibs),tm,
        fm = 'grass2',ref = True,hmin = 16,hmax = 32,minhmin = 0.1,
        zfunc = z,rv = pym.bnrm(vb).z < 0,
        uvstacked = None,autoconnect = True)
    lockf = lambda p : p.onpxy((vb,vibs)) 
    m.subdiv(tm,False,True,lockf)
    m.uvs(tm)
    print('generated land')
    return m


def stitchareas(vb,v,a,l = 10):
    #vb = [b.cp() for b in v[1]['b'][0]]
    if not pym.bccw(vb):vb.reverse()
    fnd = True
    while fnd:
        fnd = False
        for o,depth in a.enum():
            if o.ix == v.ix:continue
            #if ox == v.ix:continue
            ob = o.loop
            adjs = pym.badjbxy(vb,ob,0.1)
            for adj in adjs:
                vbx,obx = adj
                vb1,vb2,ob1,ob2 = vb[vbx-1],vb[vbx],ob[obx-1],ob[obx]
                ips = pym.sintsxyp(vb1,vb2,ob1,ob2,ieb = 0,skew = 0,ie = 0)
                if ips is None:continue
                ip1,ip2 = ips
                if not ip1 in vb or not ip2 in vb:
                    if ip1.onsxy(vb1,vb2,0):vb.insert(vbx,ip1);fnd = True
                    if ip2.onsxy(vb1,vb2,0):vb.insert(vbx,ip2);fnd = True
            if fnd:break
    vbx = 0
    while vbx < len(vb):
        vb1,vb2 = vb[vbx-1],vb[vbx]
        if vb1.d(vb2) < l:vbx += 1
        else:vb.insert(vbx,vb1.mid(vb2))
    return vb
