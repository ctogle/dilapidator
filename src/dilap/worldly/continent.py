from .polygen import lsystopy
from .terrain import terrain
from .infrastructure import roadmap
from .partitiongraph import world
from dilap.geometry import *
import dilap.geometry.polymath as pym
import dilap.geometry.tools as gtl
from dilap.core import *
import pdb


def terraininput(b,e):

    def feature(b):
        i,p,d = 5,vec3(0,0,0),vec3(0,1,0)
        axiom,rules = 'X',dict([('X','{[[X}{]X}F]X'),('F','FA'),('A','F')])
        params = dict(dazimuthal = gtl.rad(25.7),drho = 20)
        lsys = lsystem(i, p, d, axiom, rules, **params)
        feat = lsystopy(b, lsys, e)
        vec3(1,1,0).sclps(feat)
        #feat = lsystopy(b, lsys, e)
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

    tb = pym.contract(b, -e*100)
    l,r = pym.vsplitb(vec3(0,1,0),b)
    tl,bl = pym.vsplitb(vec3(1,0,0),l)
    tr,br = pym.vsplitb(vec3(1,0,0),r)

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
    for l in lows:
        ax = plot_polygon(l,ax,col = 'b')
    for h in highs:
        ax = plot_polygon(h,ax,col = 'g')
    ax = plot_polygon_xy(tb,ax,lw = 2,col = 'm')
    ax = plot_polygon_xy(b,ax,lw = 2,col = 'r')
    plt.show()

    scale = 1.0
    steps = [
        (highs,  scale, ((0.1 , 2), (0.5 , 3), (1.0, None))), 
        #(lows , -scale, ((0.08, 5), (0.02, 3), (0.1, 5), (0.2, None))), ]
        (lows , -scale, ((0.2, None), )), ]
    sealevel = 8
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
        targs = terraininput(b,e)
        t = terrain(*targs,**kws)
    else:
        t = None

    r = roadmap(t,e,**kws) if includeroads else None
    a = world(b,t,r,e,**kws)
    sg = scenegraph()

    sealevel = 0
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
        fm = 'grass2',ref = True,hmin = 8,hmax = 16,minhmin = 0.1,
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
