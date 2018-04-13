import pdb

'''
wire graph class (purely topological)

wiregraph is a simple data structure for vertices and edges
connecting them it also provides a means for ordering edges
'''

class wiregraph(object):

    ###################################
    ### fundamental topological methods
    ###################################

    def orphan(self, vx):
        return len(self.orings[vx]) == 0 if vx in self.orings else False

    # add a new intersection to the graph
    def av(self,**vkws):
        j = self.vcnt
        i = (j,vkws)
        self.vs.append(i)
        self.rings[j] = {}
        self.orings[j] = []
        self.vcnt += 1
        return j

    # remove vertex u
    def rv(self,u):
        i = self.vs[u]
        if i is None:
            print('cannot remove missing vertex', u)
        else:
            ur = self.rings[u]
            for v in set(ur.keys()):
                if u == v:
                    print('u == v')
                else:
                    self.re(u,v)
            self.vs[u] = None
            del self.rings[u]
            del self.orings[u]
            return i

    # and a new road to the graph
    def ae(self,u,v,**ekws):
        ur,vr = self.rings[u],self.rings[v]
        uor,vor = self.orings[u],self.orings[v]
        m = self.ecnt
        r = (m,ekws)
        self.es.append(r)
        self.elook[(u,v)] = m
        self.elook[(v,u)] = m
        if not v in ur:ur[v] = r
        if not u in vr:vr[u] = r

        urcnt = len(uor)
        uor.append(v)
        '''#
        if urcnt < 2:uor.append(v)
        else:
            w = uor[0]
            nea = self.ea(u,w,v)
            f = False
            for eax in range(1,urcnt):
                if nea < self.ea(u,w,uor[eax]):
                    uor.insert(eax,v)
                    f = True
                    break
            if not f:uor.append(v)
        '''#

        vrcnt = len(vor)
        vor.append(u)
        '''#
        if vrcnt < 2:vor.append(u)
        else:
            w = vor[0]
            nea = self.ea(v,w,u)
            f = False
            for eax in range(1,vrcnt):
                if nea < vor[eax]:
                    vor.insert(eax,u)
                    f = True
                    break
            if not f:vor.append(u)
        '''#

        self.ecnt += 1
        return m

    # remove an edge between u and v
    def re(self,u,v):
        r = self.rings[u][v]
        if r is None:return r
        self.es[r[0]] = None
        del self.elook[(u,v)]
        del self.elook[(v,u)]
        del self.rings[u][v]
        del self.rings[v][u]
        if v in self.orings[u]:self.orings[u].remove(v)
        if u in self.orings[v]:self.orings[v].remove(u)
        return r

    ###################################
    ### additional topological methods
    ###################################

    # add a new vertex and edge to another vertex
    def mev(self,ov,vkws,ekws):
        nv = self.av(**vkws)
        ne = self.ae(ov,nv,**ekws)
        return nv,ne

    # split an edge/road into two edges/roads
    def se(self,u,v,w):
        ruv = self.re(u,v)
        ruw = self.ae(u,w,**ruv[1])
        rwv = self.ae(w,v,**ruv[1])
        return ruw,rwv

    ###################################
    ### edge ordering mechanisms
    ###################################

    # provide an ordering mechanism for edges in the graph
    #  the default behavoir is the order of edge addition
    def ea(self,u,w,v):
        return 0

    def cw(self,u,v):
        #uor = self.orings[u]
        vor = self.orings[v]
        uori = vor.index(u)
        ror = vor[uori+1:]+vor[:uori]
        if ror:tip = ror[0]
        else:tip = u
        return tip

    def ccw(self,u,v):
        #uor = self.orings[u]
        vor = self.orings[v]
        uori = vor.index(u)
        ror = vor[uori+1:]+vor[:uori]
        if ror:tip = ror[-1]
        else:tip = u
        return tip

    # return a list of vertex indices which form a loop
    #   the first edge will be from u to v, turns of direction 
    #   d (clockwise or counterclockwise) form the loop
    def loop(self,u,v,d = 'cw',usematch = False):
        if not v in self.rings[u]:raise ValueError
        lp = [u,v]

        c = 0

        while True:

            c += 1
            if c > self.vcnt*5:
                #print('LOOPWARNING',d,u,v,len(lp))
                return self.loop(u,v,d,True)

            if   d == 'cw': tip = self.cw( lp[-2],lp[-1])
            elif d == 'ccw':tip = self.ccw(lp[-2],lp[-1])
            else:raise ValueError
            lp.append(tip)
            if lp[-1] == lp[1] and lp[-2] == lp[0]:
                lp.pop(-1)
                lp.pop(-1)
                return lp
            if usematch:
                lseqmatch = seqmatch(lp)
                if lseqmatch[0]:
                    #print('LOOPWARNINGRECONCILED')
                    if lseqmatch[0]:
                        lp.pop(-1)
                        lp.pop(-1)
                        return lp

    # return a list of all unique loops of the graph
    def uloops(self,d = 'cw'):
        loops = {}
        unfn = [x for x in range(self.ecnt)]
        for vx, v in self:
            for ox in self.orings[vx]:
                r = self.rings[vx][ox]
                rx,rkws = r
                if rx in unfn:
                    #lp = self.loop(vx,ox,'ccw')
                    lp = self.loop(vx,ox,d)
                    #lpk = tuple(set(lp))
                    lpk = set(lp)
                    #if not lpk in loops:
                    if newloopkey(lpk,loops):
                        loops[tuple(lpk)] = lp
                    unfn.remove(rx)
                else:continue
            if not unfn:break

        for key in loops:
            if 149 in key and not 150 in key and 121 in key:
                print('key', key)

        lps = [loops[lpk] for lpk in loops]
        return lps

    ###################################

    def __init__(self):
        self.vs = []
        self.vcnt = 0
        self.es = []
        self.ecnt = 0
        self.elook = {}
        self.rings = {}
        self.orings = {}

    def __iter__(self):
        for vx in range(self.vcnt):
            v = self.vs[vx]
            if not v is None:
                yield vx, v

    ###################################

def newloopkey(key,loops):
    for loop in loops:
        if set(loop) == set(key):
            return False
    return True

def seqmatch(l = list(range(10))+list(range(5))):
    #print('seqmatch: %s' % str(l))
    longest,longestseq = None,None
    gperm = lambda lx : l[lx-n:]+l[:lx] if lx-n < 0 else l[lx-n:lx]
    for n in range(2,int(len(l)+1/2)):
        fnd = False
        perms = [gperm(lx) for lx in range(len(l))]
        uniq = []
        for p in perms:
            if not p in uniq:
                uniq.append(p)
            else:
                longest = n
                fnd = True
        plen,ulen = len(perms),len(uniq)
        if plen == ulen:pass
        elif plen > ulen:
            for p in uniq:
                if perms.count(p) > 1:
                    longestseq = p
                    #print('found repeated permutation: %s' % p)
        else:pdb.set_trace()
        if not fnd:break
    #print('seqmatch (of length %d) had longest match: %s' % (len(l),str(longest)))
    #print(n)
    return longest,longestseq
 
