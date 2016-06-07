import dilap.core.base as db

from dilap.geometry.vec3 import vec3
from dilap.geometry.quat import quat
import dilap.geometry.polymath as pym

import dilap.worldly.blockletters as dbl

import dilap.core.plotting as dtl
import matplotlib.pyplot as plt

import pdb





###############################################################################
### utility functions
###############################################################################

def checkseq(fp,h,seq,show = False,**kws):
    print('check-pseq:',seq)
    kws['footprint'] = fp
    kws['height'] = h
    kws['sequence'] = seq
    p = pgraph(**kws).graph()
    if show:
        #p.plot()
        p.plotxy()
        plt.show()
    return p

###############################################################################
### fundamental partition graph operations
###############################################################################

def split(g,subseq):
    srx,spx,spy,sdx,sdy,sframe = [v for v in subseq.split(',')]
    spz,sdz = g.vs[int(srx)][1]['fp'][0].z,0
    lp,ld = vec3(float(spx),float(spy),spz),vec3(float(sdx),float(sdy),sdz)
    fp = g.vs[int(srx)][1]['fp']
    if   sframe == 'root':  lp = g.tfp(lp,g.footprint)
    elif sframe == 'vertex':lp = g.tfp(lp,fp)
    sp1 = lp.cp().trn(ld.cp().uscl( 10000))
    sp2 = lp.cp().trn(ld.cp().uscl(-10000))
    nfps = pym.bsegsxy(fp,sp1,sp2)
    l = nfps.pop(0)
    #for r in nfps:new = g.sv(int(srx),l,r,False,False)
    for r in nfps:new = g.sv(int(srx),l,r,True,True)
    return g

# split a vertex by nesting a union of polygons as a void of the vertex
def splotch(g,subseq):
    irx,easement,splotchseq = subseq.split(',')
    irx,easement = int(irx),float(easement)
    print('SPLOTCH!',subseq,splotchseq)
    iv = g.vs[irx]
    fp = iv[1]['fp']

    xprj = vec3(1,0,0).prjps(fp)
    yprj = vec3(0,1,0).prjps(fp)
    xlen = xprj[1]-xprj[0]
    ylen = yprj[1]-yprj[0]
    rad = max(xlen,ylen)/2.0
    
    #r1 = vec3(0,10,0).com(fp).sq(30,40)
    #r2 = vec3(-10,-20,0).com(fp).sq(20,40)
    #r3 = vec3(0,-10,0).com(fp).sq(50,30)
    #rs = [r1,r2,r3]

    #r = rs.pop(0)
    #while rs:r = pym.ebuxy(r,rs.pop(0))

    r = dbl.block('H',rad/3.0,rad,rad)
    r = r[0]

    for j in range(3):r = pym.smoothxy(r,0.1)

    #r = pym.ebixy(fp,r)
    l = fp

    nv = g.sv(irx,l,r)

    ax = dtl.plot_axes_xy(2*rad)
    ax = dtl.plot_polygon_xy(fp,ax,lw = 6,col = None)
    ax = dtl.plot_polygon_xy(r,ax,lw = 4,col = 'b')
    ax = dtl.plot_polygon_xy(l,ax,lw = 2,col = 'g')
    plt.show()

    return g

# contract the voids of a vertex
def bleed(g,subseq):
    irx,rad = subseq.split(',')
    irx,rad = int(irx),float(rad)
    voids = g.vs[irx][1]['voids']
    for vx in range(len(voids)):
        voids[vx] = pym.contract(voids[vx],rad)
    return g

def merge(g,subseq):

    print('merge verts',subseq)
    
    pdb.set_trace()
    return g

def insert(g,subseq):
    irx = int(subseq[:subseq.find(',')])
    subgseq = subseq[subseq.find(',')+1:]
    iv = g.vs[irx][1]
    fp,h = [p.cp() for p in iv['fp']],iv['height']
    subg = checkseq(fp,h,subgseq,False,grammer = g.grammer)
    g.rvag(irx,subg,True)
    return g

def addexit(g,subseq):
    vx,ex,ey,ez,frame = subseq.split(',')
    vx = int(vx)
    exloc = vec3(float(ex),float(ey),float(ez))
    if   frame == 'root':  fp = g.footprint
    elif frame == 'vertex':fp = g.vs[vx][1]['fp']
    exloc = g.tfp(exloc,fp)
    g.vs[int(vx)][1]['exits'].append(exloc)
    return g

def addedge(g,subseq):
    ss = subseq.split(',')
    g.ae(int(ss[0]),int(ss[1]))
    return g

def settype(g,subseq):
    ss = subseq.split(',')
    g.vs[int(ss[0])][1]['type'] = ss[1:]
    return g



###############################################################################
###
###############################################################################



# pgraph represents a partitioning of 3d space
class pgraph(db.base):



    def plotxy(self,ax = None):
        scale = vec3(1,0,0).prjps(self.footprint)
        scale = scale[1]-scale[0]
        if ax is None:ax = dtl.plot_axes_xy(scale)
        for vx in range(self.vcnt):
            v = self.vs[vx]
            if v is None:continue
            vb = v[1]['fp']
            vrtls = '-' if 'developed' in v[1]['type'] else '--'
            #vrtls = '--' if 'natural' in v[1]['type'] else '-'
            ax = dtl.plot_polygon_xy(
                pym.contract(vb,0.01*scale),ax,lw = 2,ls = vrtls,col = None)
            vbc = vec3(0,0,0).com(vb)
            rs = str(vx)+','+str(v[1]['edges'])+',\n'+str(v[1]['type'])
            ax = dtl.plot_point_xy_annotate(vbc,ax,rs)
            for ve in v[1]['edges']:
                if self.vs[ve] is None:continue
                ovbc = vec3(0,0,0).com(self.vs[ve][1]['fp'])
                ax = dtl.plot_edges_xy((vbc,ovbc),ax,lw = 3,col = 'c')
            for ep in v[1]['exits']:
                ax = dtl.plot_point_xy(ep,ax,mk = 's',col = 'r')
        ax = dtl.plot_polygon_xy(self.footprint,ax,lw = 1,col = 'g')
        return ax



    def plot(self,ax = None):
        scale = vec3(1,0,0).prjps(self.footprint)
        scale = scale[1]-scale[0]
        if ax is None:ax = dtl.plot_axes(scale)
        for vx in range(self.vcnt):
            v = self.vs[vx]
            if v is None:continue
            vb = v[1]['fp']
            vrtls = '--' if 'natural' in v[1]['type'] else '-'
            ax = dtl.plot_polygon(
                pym.contract(vb,0.01*scale),ax,lw = 2,ls = vrtls,col = None)
            vbc = vec3(0,0,0).com(vb)
            for ve in v[1]['edges']:
                if self.vs[ve] is None:continue
                ovbc = vec3(0,0,0).com(self.vs[ve][1]['fp'])
                ax = dtl.plot_edges((vbc,ovbc),ax,lw = 3,col = 'c')
            for ep in v[1]['exits']:
                ax = dtl.plot_point(ep,ax,mk = 's',col = 'r')
        top = [p.cp().ztrn(self.height) for p in self.footprint]
        ax = dtl.plot_polygon(self.footprint,ax,lw = 1,col = 'g')
        ax = dtl.plot_polygon(top,ax,lw = 1,col = 'b')
        return ax



    # add a vertex to the graph
    def av(self,vkws,null = False):
        vx = self.vcnt
        self.vcnt += 1
        if null:self.vs.append(None)
        else:self.vs.append([vx,vkws])
        return vx



    # merge two vertices
    def mv(self,u,v):
        if u == v:return
        uv,vv = self.vs[u][1],self.vs[v][1]

        print('merge vertices!',u,v)

        raise NotImplemented



    # replace a vertex with the vertices of another graph
    def rvag(self,vx,og,connect_all = False):
        vc = self.vcnt
        v = self.vs[vx]
        nvs = []
        for ox in range(og.vcnt):
            ov = og.vs[ox]
            if ov is None:
                self.av(None,null = True)
                continue
            okws = ov[1]
            newring = [ovx+vc for ovx in okws['edges']]
            okws['edges'] = newring
            if 'terrainmesh' in v[1]['info']:
                okws['info']['terrainmesh'] = v[1]['info']['terrainmesh']
            new = self.av(okws)
            nvs.append(new)
        for res in v[1]['edges']:
            ringv = self.vs[res]
            if ringv is None:continue
            for nvx in nvs:
                nv = self.vs[nvx]
                for ep in nv[1]['exits']:
                    if ep.onbxy(ringv[1]['fp']):
                        self.ae(nvx,res)
                if connect_all:
                    aws = pym.badjbxy(nv[1]['fp'],ringv[1]['fp'],0)
                    if aws:self.ae(nvx,res)
        self.vs[vx] = None
        for vx in range(self.vcnt):self.vve(vx)
        return nvs



    # verify the geometric viability of each edge of a vertex
    # when an exit resides geometrically between two vertices, add an edge
    def vve(self,vx):
        v = self.vs[vx]
        if v is None:return
        vb = v[1]['fp']
        rem = []
        for adj in v[1]['edges']:
            ov = self.vs[adj]
            if ov is None:
                rem.append(adj)
                continue
            ob = ov[1]['fp']
            aws = pym.badjbxy(vb,ob)
            if not aws:self.re(adj,vx)
        for rm in rem:v[1]['edges'].remove(rm)
        for vexit in v[1]['exits']:
            for ox in range(self.vcnt):
                ov = self.vs[ox]
                if ov is None:continue
                if vexit.onbxy(ov[1]['fp']):
                    self.ae(ox,vx)



    # add an edge to the graph
    def ae(self,u,v):
        if u == v:return
        ur,vr = self.vs[u][1]['edges'],self.vs[v][1]['edges']
        if not v in ur:ur.append(v)
        if not u in vr:vr.append(u)



    # remove the edge between two vertices if it exists
    def re(self,u,v):
        ur,vr = self.vs[u][1]['edges'],self.vs[v][1]['edges']
        if v in ur:ur.remove(v)
        if u in vr:vr.remove(u)



    # split a vertex based on intersection with a polyline (adds a vertex)
    def sv(self,vx,l,r,connect_split = True,connect_ring = True):
        sv = self.vs[vx]
        if pym.binbxy(r,l):sv[1]['voids'].append(r)
        sv[1]['fp'] = l
        newes = []
        if connect_ring:newes.extend(sv[1]['edges'])
        if connect_split:newes.append(vx)
        vinfo = {}
        if 'terrainmesh' in sv[1]['info']:
            vinfo['terrainmesh'] = sv[1]['info']['terrainmesh']
        newkws = {
            'fp':r,'voids':[],
            'height':sv[1]['height'],'type':sv[1]['type'][:],
            'edges':newes,'exits':[],'info':vinfo,
                }
        new = self.av(newkws)
        if pym.binbxy(l,r):self.vs[new][1]['voids'].append(l)
        if connect_ring:
            for ringvx in sv[1]['edges']:
                self.ae(new,ringvx)
        if connect_split:self.ae(new,vx)
        return new



    # transform a relative location to an absolute location 
    #   based on the xy projection of the footprint
    def tfp(self,p,fp):
        bprjx = vec3(1,0,0).prjps(fp)
        bprjy = vec3(0,1,0).prjps(fp)
        bprjz = vec3(0,0,1).prjps(fp)
        bprjz = bprjz[0],bprjz[0]+self.height
        lp = vec3(
            bprjx[0]+p.x*(bprjx[1]-bprjx[0]),
            bprjy[0]+p.y*(bprjy[1]-bprjy[0]),
            bprjz[0]+p.z*(bprjz[1]-bprjz[0]))
        return lp



    # return a set of vertices constituting a topological 
    #   path between two distinct vertices in the graph
    def path(self,u,v):
        p = []
        uv,vv = self.vs[u],self.vs[v]
        if u == v:return p

        raise NotImplemented



    def __init__(self,*ags,**kws):
        self.vs = []
        self.vcnt = 0

        self._def('footprint',None,**kws)
        self._def('height',10,**kws)
        self._def('stitchsize',2,**kws)
        self._def('sequence','',**kws)
        
        self.grammer = {
            'S':split,'J':splotch,'M':merge,'I':insert,
            'X':addexit,'E':addedge,'R':settype,'V':bleed,
                }
        if 'grammer' in kws:
            for eg in kws['grammer']:
                self.grammer[eg] = kws['grammer'][eg]

        rootkws = {
            'fp':[p.cp() for p in self.footprint],
            'height':self.height,'edges':[],'exits':[],
            'type':['ocean'],'info':{},'voids':[],
                }
        self.av(rootkws)



    def graph(self):
        seq = self.sequence[:]
        scnt = len(seq)
        sx = 0
        while sx < scnt:
            c = seq[sx]
            if c in self.grammer:
                ex = db.seqread(seq,sx)
                self.grammer[c](self,seq[sx+2:ex])
            else:ex = sx+1
            sx = ex
        for vx in range(self.vcnt):self.vve(vx)
        return self



###############################################################################
###############################################################################



