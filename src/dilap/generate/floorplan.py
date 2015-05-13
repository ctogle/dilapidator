import dilap.core.base as db
import dilap.core.tools as dpr
import dilap.core.sgraph as dsg
import dilap.core.context as dgc
import dilap.io.io as dio
import dilap.primitive.cube as dcu
import dilap.primitive.cone as dco
import dilap.primitive.cylinder as dcyl
import dilap.primitive.wall as dw
import dilap.primitive.floor as df
import dilap.generate.room as dgr
import dilap.generate.shaft as dsh
import dilap.generate.roof as drf

import dp_vector as dpv
import dp_bbox as dbb

import random as rm

import pdb

class floorplan(dgc.context):

    def __init__(self,building,*args,**kwargs):
        dgc.context.__init__(self,*args,**kwargs)
        self.bldg = building
        self._def('max_rooms',2,**kwargs)

    def should_shaft(self,rmplan):
    #def should_shaft(self,newpos,newl,neww):
        #sdists = [dpv.distance(newpos,sc.position) for sc in self.sectors]
        #if not sdists:sdmin = 100
        #else:sdmin = min(sdists)

        #if sdmin > 20 and newl >= 20 and neww >= 24:
        splan = None
        newl,neww = rmplan[1]['l'],rmplan[1]['w']
        if newl >= 20 and neww >= 24:
            shx,shy = rmplan[1]['x'],rmplan[1]['y']
            shl = 8
            shw = 10
            sps = dpv.vector(shx,shy,0)
            gap = (dpv.zero(),shl,shw)
            rmplan[1]['shafted'] = True
            rmplan[1]['fgap'] = gap
            rmplan[1]['cgap'] = gap
            splan = ((),{'p':sps,'l':gap[1],'w':gap[2]})
        return splan

    # given the plan for a room, return verts for its corners
    def wall_verts(self,rmargs):
        rmkws = rmargs[1]
        cs = dpr.corners(rmkws['l'],rmkws['w'],
            dpv.vector(rmkws['x'],rmkws['y'],0.0))
        pairs = []
        for cdx in range(1,len(cs)):
            c1,c2 = cs[cdx-1].copy(),cs[cdx].copy()
            pairs.append((c1,c2))
        c1,c2 = cs[-1].copy(),cs[0].copy()
        pairs.append((c1,c2))
        return pairs

    # initializes the set of plans with a single room 
    # this room is always considered the main entry point
    def entry(self):
        l,w = self.bldg.l,self.bldg.w
        subl = dpr.clamp(rm.choice([0.25*l,0.5*l,0.75*l]),20,l)
        subw = dpr.clamp(rm.choice([0.25*w,0.5*w,0.75*w]),24,w)
        mx = 0.0
        my = - w/2.0 + subw/2.0 + 10.0
        margs = ((),{'x':mx,'y':my,'l':subl,'w':subw,'shafted':False})
        splan = self.should_shaft(margs)
        cpairs = self.wall_verts(margs)
        ewargs = [[cp,
            {'h':4.0,'w':0.5,'walltype':'exterior','room':margs}] 
                for cp in cpairs]
        # 0 may not be the front wall, i did not check
        ewargs[0][1]['walltype'] = 'entryway'
        #ewargs[0][1]['doorgaps'] = [(0.5,0.25)]

        rmplans = [margs]
        ewplans = ewargs[:]
        iwplans = []
        shplans = [splan] if splan else []
        return rmplans,ewplans,iwplans,shplans

    # choose a side if necessary
    def grow_side(self,plans,side):
        if side is None:
            options = plans[1][:]
            while options:
                which = rm.choice(options)
                if not which[1]['walltype'] == 'entryway':
                    side = which
                    break
                options.remove(which)
        return side

    # choose a grow length if necessary
    def grow_length(self,plans,length,side):
        if length is None:

            #bdist = side._distance_to_border(self.corners)
            #if bdist < 8 and not force:
            #    #print 'too close to a border to grow'
            #    return False
            #elif gleng > bdist: gleng = bdist

            length = rm.choice([8,12,16,20,24,28,32])
        return length

    def face_away(self,side):
        # THIS IS WHERE IT GOES OFF THE RAILS
        rp = side[1]['room']

        intpt = dpv.vector(rp[1]['x'],rp[1]['y'],0.0)
        midpt = dpv.midpoint(*side[0])
        tangt = dpv.v1_v2(*side[0]).normalize()
        norml = tangt.copy().rotate_z(dpr.rad(90)).normalize()
        tstpt = midpt.copy().translate(norml)

        side[1]['normal'] = norml
        if dpv.distance(intpt,midpt) > dpv.distance(intpt,tstpt):
            side[1]['normal'].flip()

    # plans is a tuple of 4 lists of plans
    #def grow(self,length = None,side = None,force = False):
    def grow(self,plans,length = None,side = None,force = False):
        side = self.grow_side(plans,side)
        if side is None:return False
        gleng = self.grow_length(plans,length,side)
        if gleng is None:return False

        self.face_away(side)
        v1,v2 = side[0]
        c1 = v2.copy()
        c2 = v1.copy()
        c3,c4 = dpr.extrude_edge(c1,c2,gleng,side[1]['normal'])

        newcorners = [c1,c2,c3,c4]
        x,y,z = dpv.center_of_mass(newcorners)
        #l,w = c2.x-c1.x,c4.y-c1.y
        xpj = dpv.project_coords(newcorners,dpv.xhat)
        ypj = dpv.project_coords(newcorners,dpv.yhat)
        l,w = xpj.y-xpj.x,ypj.y-ypj.x

        margs = ((),{'x':x,'y':y,'l':l,'w':w,'shafted':False})
        #gaps = self.should_shaft(sect.position,sect.length,sect.width)
        #margs['fgaps'] = gaps[:]
        #margs['cgaps'] = gaps[:]
        #splan = self.should_shaft(margs)
        #if gaps:sect.shafted = True

        #for esect in self.sectors:
        #    ebb = esect.get_bboxes()
        #    nbb =  sect.get_bboxes()
        #    if mpbb.intersects(ebb,nbb):
        #        if gaps: self.shaft_kwargs.pop(-1)
        #        #print 'new sect intersected!'
        #        return False

        cpairs = [(c2,c3),(c3,c4),(c4,c1)]
        ewargs = [[cp,
            {'h':4.0,'w':0.5,'walltype':'exterior','room':margs}] 
                for cp in cpairs]
        iwargs = []

        #if self.resolve_walls(extwalls,intwalls,sect):
        if True:
            #self.switch_wall_sort(side)

            rps,eps,ips,sps = plans
            rps.append(margs)
            eps.extend(ewargs)
            ips.extend(iwargs)
            return True
        else:return False

    def plan_specific(self,specific):
        if specific > 0:self.allplans[1][0][1]['walltype'] = 'exterior'
        if specific == self.bldg.stories-1:
            for sp in self.allplans[0]:
                sp[1]['cgap'] = None

    def plan(self):
        plans = self.entry()
        
        #self.allplans = plans
        #return

        for sect in range(self.max_rooms):
            grew = self.grow(plans)
            if grew is False:
                print('growth rejected')
                pass
            elif grew is None:
                print('division aborted')
                break
        self.allplans = plans

    # generate and return shafts
    def generate_shafts(self,worn = 0):
        rooms,ewalls,iwalls,shafts = [],[],[],[]
        rps,eps,ips,sps = self.allplans
        for sp in sps:shafts.append(dsh.shaft(self.bldg,*sp[0],**sp[1]))
        for s in shafts:s.generate(worn)
        return shafts

    # generate and return roof pieces
    def generate_roof(self,worn = 0):
        roof = drf.roof(self)
        roof.generate(worn)
        return roof

    # generates a floor without associated shafts
    # other contexts are expected to use this contexts nodes directly
    def generate(self,wallheight = 4.0,worn = 0):
        rooms,ewalls,iwalls,shafts = [],[],[],[]
        rps,eps,ips,sps = self.allplans
        for rp in rps:rooms.append(dgr.room(*rp[0],**rp[1]))
        for ep in eps:ewalls.append(dw.wall(*ep[0],**ep[1]))
        for ip in ips:iwalls.append(dw.wall(*ip[0],**ip[1]))
        for r in rooms:
            r.generate(wallheight,worn)
            for n in r.sgraph.nodes:n._assign_material('concrete2')
            self._consume(r)
        for ew in ewalls:ew._assign_material('brick2')
        for iw in iwalls:iw._assign_material('concrete2')
        self._nodes_to_graph(
            self._node_wrap(*ewalls),
            self._node_wrap(*iwalls))


