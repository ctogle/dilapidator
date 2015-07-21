import dilap.core.base as db
import dilap.core.uinfo as di
import dilap.core.mesh.tools as dtl
import dilap.core.tools as dpr

import dp_vector as dpv

import xml.etree.cElementTree
import matplotlib.pyplot as plt
import math

import pdb

# p1,2 are latitudes, l1,2 are longitudes
def haverside(p1,p2,l1,l2,earthradius = 6371000):
    p1,p2,l1,l2 = dpr.rad(p1),dpr.rad(p2),dpr.rad(l1),dpr.rad(l2) 
    dp = abs(p1-p2)
    dl = abs(l1-l2)
    p = math.sin(dp/2)**2
    l = math.cos(p1)*math.cos(p2)*math.sin(dl/2)**2
    h = 2*earthradius*math.asin(math.sqrt(p+l))
    return h

class osm_world(db.base):

    def plot_boundary(self,ax):
        l1,l2,p1,p2,lonext,latext,xext,yext = self.extents
        b1x,b1y = self._transform(l1,p1)
        b2x,b2y = self._transform(l2,p1)
        b3x,b3y = self._transform(l2,p2)
        b4x,b4y = self._transform(l1,p2)
        b1 = dpv.vector(b1x,b1y,0)
        b2 = dpv.vector(b2x,b2y,0)
        b3 = dpv.vector(b3x,b3y,0)
        b4 = dpv.vector(b4x,b4y,0)
        dtl.plot_polygon([b1,b2,b3,b4],ax)
        return ax

    def plot_way(self,which,ax):
        way = self.ways[which]
        nds = [c.get('ref') for c in way if c.tag == 'nd']
        nds = [self.nodes[self.nodes_lookup[n]] for n in nds] 
        npos = []
        for nd in nds:
            nlon,nlat = nd.get('lon'),nd.get('lat')
            nx,ny = self._transform(float(nlon),float(nlat))
            np = dpv.vector(nx,ny,0.0)
            npos.append(np)
        dtl.plot_edges(npos,ax)
        return ax

    def plot_all_routes(self,ax):
        rnames = self.find_routes()
        for name in rnames:ax = self.plot_route_by_name(name,ax)
        return ax

    def plot_route_by_name(self,name,ax):
        which = self.find_route(name)
        if which is None:
            print('failed to plot route:',name)
            return ax
        whichms = which.findall('member')
        whichws = []
        for wm in whichms:
            mref = wm.get('ref')
            if mref in self.ways_lookup:
                whichws.append(self.ways_lookup[mref])
            print('failed to find route member:',mref)
        for w in whichws:ax = self.plot_way(w,ax)
        return ax

    def find_routes(self):
        routes = {}
        nonroutes = {}
        for relx in range(self.relationcount):
            rtags = self.relations[relx].findall('tag')
            rname = 'None'
            rtype = 'None'
            for t in rtags:
                if t.get('k') == 'name':rname = t.get('v')
                if t.get('k') == 'type':rtype = t.get('v')
                if not rname == 'None' and not rtype == 'None':break
            if True:routes[rname] = relx
            #if rtype == 'route':routes[rname] = relx
            else:nonroutes[rname] = relx
        self.routes = routes
        self.nonroutes = nonroutes
        return routes

    def find_route(self,name):
        routes = self.routes
        candidates = [r for r in routes 
            if r.count(name) and not r.count('(super)')]
        candcnt = len(candidates)
        if candcnt == 0:return
        elif candcnt == 1:cx = 0
        else:
            print('choose from:'+str(candidates)+':')
            #cx = int(input('\t\t::'))
            cx = 0
        rrel = self.relations[routes[candidates[cx]]]
        return rrel

    def __init__(self,bounds,nodes,ways,relations):
        self._set_transform(bounds)
        self.nodes = []
        self.ways = []
        self.relations = []
        self.nodes_lookup = {}
        self.ways_lookup = {}
        self.relations_lookup = {}
        self.nodecount = 0
        self.waycount = 0
        self.relationcount = 0
        self._add_nodes(*nodes)
        self._add_ways(*ways)
        self._add_relations(*relations)
        self.find_routes()

    def _set_transform(self,bnds):
        self.bounds = bnds
        p1,p2 = float(bnds.get('minlat')),float(bnds.get('maxlat'))
        l1,l2 = float(bnds.get('minlon')),float(bnds.get('maxlon'))
        lavg = (l1+l2)/2.0
        pavg = (p1+p2)/2.0
        xextent = haverside(pavg,pavg,l1,l2)
        yextent = haverside(p1,p2,lavg,lavg)
        lonextent = l2-l1
        latextent = p2-p1
        self.extents = (l1,l2,p1,p2,lonextent,latextent,xextent,yextent)

    def _transform(self,lon,lat):
        x = self.extents[6]*(lon-self.extents[0])/self.extents[4]
        y = self.extents[7]*(lat-self.extents[2])/self.extents[5]
        return x,y

    def _add_nodes(self,*nds):
        for nd in nds:
            self.nodes.append(nd)
            self.nodes_lookup[nd.get('id')] = self.nodecount
            self.nodecount += 1

    def _add_ways(self,*wys):
        for wy in wys:
            self.ways.append(wy)
            self.ways_lookup[wy.get('id')] = self.waycount
            self.waycount += 1

    def _add_relations(self,*rls):
        for rl in rls:
            self.relations.append(rl)
            self.relations_lookup[rl.get('id')] = self.relationcount
            self.relationcount += 1

def read_xml(xmlfile):
    tree = xml.etree.cElementTree.parse(xmlfile)
    root = tree.getroot()
    obnds = root.find('bounds')
    onods = root.findall('node')
    oways = root.findall('way')
    orels = root.findall('relation')
    oworld = osm_world(obnds,onods,oways,orels)
    return oworld

def test():
    user_info = di.fetch_info()
    osmdata = user_info['osmdata']
    oworld = read_xml(osmdata)

    ax = dtl.plot_axes()
    ax = oworld.plot_boundary(ax)
    ax = oworld.plot_all_routes(ax)
    #ax = oworld.plot_way_by_name('I 95',ax)
    #ax = oworld.plot_way_by_name('I 64',ax)
    #ax = oworld.plot_way_by_name('360',ax)
    plt.show()



    pdb.set_trace()

    

test()



