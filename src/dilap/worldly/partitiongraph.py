import dilap.geometry.polymath as pym
from dilap.geometry import *
from dilap.core.plotting import *
import pdb


class partition(topography):


    def plot(self,ax = None):
        if ax is None:
            ax = plot_axes(300)
        for v,depth in self.enum(False):
            col = 'g' if self.below(v) else 'b'
            vloop = [p.ztrn(depth*10) for p in pym.contract(v.loop,2)]
            #vloop = [p.cp().ztrn(depth*10) for p in v.loop]
            vanchor = vloop[0]
            if depth > 0:
                panchor = self.above(v).loop[0].cp().ztrn((depth-1)*10)
            else:panchor = vec3(0,0,0)
            ax = plot_polygon(vloop,ax,lw = 3,col = col)
            ax = plot_edges((vanchor,panchor),ax,col = col,lw = 5)
            for hole in v.holes:
                hole = [p.ztrn(depth*10) for p in pym.contract(hole,-2)]
                #hole = [p.cp().ztrn(depth*10) for p in hole]
                ax = plot_polygon(hole,ax,lw = 2,col = 'r')
        return ax


    def split(self,sv,nvloop,svloop = None,**kws):
        '''Add two vertices with parent vertex sv; ensure resulting loops do not overlap'''
        nsvs = []
        nv1holes,nv2holes = sv.holes[:],[]
        if svloop is None:
            svloop = sv.loop[:]
            if pym.binbxy(nvloop,svloop):
                nv1holes.append(nvloop)
                nsvs.append(topography.avert(self,sv,svloop,holes = nv1holes))
            elif pym.binbxy(svloop,nvloop):
                nv2holes.append(svloop)
                nsvs.append(topography.avert(self,sv,svloop,holes = nv1holes))
            elif pym.bintbxy(nvloop,svloop,ie = False):
                svloop = pym.ebdxy(svloop,nvloop)
                if len(svloop) > 1:
                    for svl in svloop[:-1]:
                        nsvs.append(topography.avert(self,sv,svl,holes = nv1holes))
                svloop = svloop[-1]
                nsvs.append(topography.avert(self,sv,svloop,holes = nv1holes))
        nv = topography.avert(self,sv,nvloop,holes = nv2holes)
        for k in kws:
            for nsv in nsvs:
                nsv.__setattr__(k,sv.__getattribute__(k))
            nv.__setattr__(k,kws[k])
        return nsvs,nv


def world(b,t,r,e,s = 'land'):
    '''Create a partition representing a world from a terrain, roadmap, and epsilon'''
    p = partition(loop = b,holes = [],style = s)

    if r is None:
        blgfps,rpy = [],None
    else:
        blgfps,rpy = r.footprints,pym.pgtopy(r.roads,5)

        # could calculate landmass boundary from sealevel style t.verts
        lmpv = p.root
        lmpvs,rdpv = p.split(lmpv,rpy[0],style = 'infrastructure')

        for fp in blgfps[0]:
            fnd = False
            for lmpv in lmpvs:
                if pym.binbxy(fp[0],lmpv.loop):
                    #blgpvs,blgpv = p.split(lmpv,fp[0],style = 'structure')
                    #found = True
                    #break
                    pass
            if fnd:break

        if len(rpy[1]) > 0:
            irdpvs,irdpv = p.split(rdpv,rpy[1][0],style = 'bounded')
            if len(rpy[1]) > 1:
                for irpy in rpy[1][1:]:
                    irdpvs,irdpv = p.split(irdpvs[0],irpy,style = 'bounded')

        # here is where loop correction is applied?
        rtv = t.avert(t.locloop(rpy[0], 1)[0], loop=rpy[0])
        for dpy in rpy[1]:
            rtv = t.avert(t.locloop(dpy, 1)[0], loop=dpy)

    '''#
    for lmtv in t.below(t.root):

        rootpvs,lmpv = p.split(p.root,lmtv.loop,None,style = 'land')

        if not rpy is None:
            lmpvs,rdpv = p.split(lmpv,rpy[0],style = 'infrastructure')

            for fp in blgfps[0]:
                fnd = False
                for lmpv in lmpvs:
                    if pym.binbxy(fp[0],lmpv.loop):
                        #blgpvs,blgpv = p.split(lmpv,fp[0],style = 'structure')
                        #found = True
                        #break
                        pass
                if fnd:break

            if len(rpy[1]) > 0:
                irdpvs,irdpv = p.split(rdpv,rpy[1][0],style = 'bounded')
                if len(rpy[1]) > 1:
                    for irpy in rpy[1][1:]:
                        irdpvs,irdpv = p.split(irdpvs[0],irpy,style = 'bounded')

            # here is where loop correction is applied?
            rtv = t.avert(t.locloop(rpy[0], 1)[0], loop=rpy[0])
            for dpy in rpy[1]:
                rtv = t.avert(t.locloop(dpy, 1)[0], loop=dpy)
    '''#
    return p
