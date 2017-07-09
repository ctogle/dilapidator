from .polygen import lsystopy
from .terrain import terrain
from .infrastructure import roadmap
from .partitiongraph import world
from dilap.geometry import *
import dilap.geometry.polymath as pym
import dilap.geometry.tools as gtl
from dilap.core import *
import pdb


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

        tb = pym.contract(b, -e*100)
        i,p,d = 6,vec3(0,0,0),vec3(0,1,0)
        axiom,rules = 'X',dict([('X','{[[X}{]X}F]X'),('F','FA'),('A','F')])
        params = dict(dazimuthal = gtl.rad(25.7),drho = 20)
        lsys = lsystem(i, p, d, axiom, rules, **params)
        highs = [pym.bisectb(lsystopy(tb, lsys, e))]
        lows  = [pym.ebdxy(tb, h) for h in highs]
        lows  = [x for y in lows for x in y]
        scale = 1.0
        steps = [
            (highs,  scale, ((0.1 , 2), (0.5 , 3), (1.0, None))), 
            (lows , -scale, ((0.08, 5), (0.02, 3), (0.1, 5), (0.2, None))), ]
        sealevel = 8

        t = terrain(tb,e,steps,sealevel,**kws)
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
    #holes = [pym.contract([p.cp() for p in ib],bleed,5) for ib in holes]
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

    print('generating land')
    ngvs = m.asurf((vb,vibs),tm,
        fm = 'grass2',ref = True,hmin = 16,hmax = 16,minhmin = 0.1,
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
