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

    def _terrain_points(self):
        def addtpt(tp):
            for etp in tpts:
                if etp.near(tp):
                    return
            tpts.append(tp)
        ewplans = self.allplans[1]
        tpts = []
        for ew in ewplans:
            addtpt(ew[0][0])
            addtpt(ew[0][1])
        return tpts

    def __init__(self,building,*args,**kwargs):
        dgc.context.__init__(self,*args,**kwargs)
        self.bldg = building
        self._def('max_rooms',0,**kwargs)
        self._def('min_shaft_distance',50,**kwargs)

    # return a plan if a shaft would be good, otherwise None
    def should_shaft(self,plans,rmplan):
        rps,eps,ips,sps = plans
        sdist = 1000.0
        rpos = dpv.vector(rmplan[1]['x'],rmplan[1]['y'],0)
        for sp in sps:
            spos = sp[1]['p']
            sd = dpv.distance(spos,rpos)
            if sd < sdist:sdist = sd
        splan = None
        newl,neww = rmplan[1]['l'],rmplan[1]['w']
        if sdist > self.min_shaft_distance and newl >= 16 and neww >= 18:
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
        self.create_bbox(margs)
        splan = self.should_shaft([[],[],[],[]],margs)
        cpairs = self.wall_verts(margs)
        ewargs = [[cp,
            {'h':4.0,'w':0.5,'walltype':'exterior','room':margs}] 
                for cp in cpairs]
        ewargs[0][1]['walltype'] = 'entryway'

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

    # add a bbox to a rooms plan based on its plan
    def create_bbox(self,roomplan):
        x,y = roomplan[1]['x'],roomplan[1]['y']
        l,w = roomplan[1]['l']-0.01,roomplan[1]['w']-0.01
        cns = dpr.corners(l,w,dpv.vector(x,y,0))
        xpj = dpv.project_coords(cns,dpv.xhat)
        ypj = dpv.project_coords(cns,dpv.yhat)
        zpj = dpv.project_coords(cns,dpv.zhat)
        bb = dbb.bbox(xpj,ypj,zpj)
        roomplan[1]['bbox'] = bb 

    # choose a grow length if necessary
    def grow_length(self,plans,length,side):
        if length is None:

            l,w = self.bldg.l,self.bldg.w
            corners = dpr.corners(l,w)

            sidept = dpv.midpoint(*side[0])
            bdist = dpv.distance_to_border_xy(sidept,corners)
            length = rm.choice([8,12,16,20,24,28,32])
            if bdist < 8:return
            elif length > bdist:length = bdist
        return length

    def verify_growth(self,plans,margs):
        mbb = margs[1]['bbox']
        for rp in plans[0]:
            bb = rp[1]['bbox']
            print('bb',bb)
            #if not dbb.separating_axis(mbb,bb):
            if not mbb.separating_axis(bb):
                print('mbb',mbb)
                print('bb',bb)
                print('bbs intersected')
                return False
        return True

    # based on the position of the room associated with this wall
    # recalculate its normal so that growths work intuitively
    def face_away(self,side):
        rp = side[1]['room']
        intpt = dpv.vector(rp[1]['x'],rp[1]['y'],0.0)
        midpt = dpv.midpoint(*side[0])
        tangt = dpv.v1_v2(*side[0]).normalize()
        norml = tangt.copy().rotate_z(dpr.rad(90)).normalize()
        tstpt = midpt.copy().translate(norml)
        side[1]['normal'] = norml
        if dpv.distance(intpt,midpt) > dpv.distance(intpt,tstpt):
            side[1]['normal'].flip()

    # if side is an exterior wall, make it an interior one, and visa versa
    def switch_walltype(self,plans,side):
        wtype = side[1]['walltype']
        if   wtype == 'exterior':
            wtype = 'interior'
            plans[1].remove(side)
            plans[2].append(side)
        elif wtype == 'interior':
            wtype = 'exterior'
            plans[1].append(side)
            plans[2].remove(side)
        side[1]['walltype'] = wtype
        return side

    # compare the new walls (ewargs,iwargs), modify them as possible
    # to resolve issues with existing walls in plans
    def resolve_walls(self,plans,ewargs,iwargs,margs):
        rps,eps,ips,sps = plans

#####
        return True

        # ewalls and iwalls are new but will be added soon
        eewalls = self.exterior_walls
        eiwalls = self.interior_walls
        
        def resolve_verts(shared,one,two):
            done = cv.distance(shared,one)
            dtwo = cv.distance(shared,two)
            if done < dtwo:
                ipair = [shared,one]
                epair = [one,two]
            elif done > dtwo:
                ipair = [shared,two]
                epair = [two,one]
            return ipair, epair

        newalls = []
        niwalls = []

        etoiwalls = []
        extrawalls = []

        for ew in ewalls: # for each new potential wall
            ewvs = [ew.v1,ew.v2] # positions of ew endpoints
            ewn = ew.normal.copy()
            #ewt = ew.v1v2.copy().normalize()
            ewt = ew.tangent.copy().normalize()
            ewewnproj = mpbb.project(ewvs,ewn) # projection of ew onto ewn; always a dot

            for eew in eewalls: # for each existing exterior wall
                eewvs = [eew.v1,eew.v2] # positions of eew endpoints
                eewn = eew.normal.copy()
                eewt = eew.tangent.copy().normalize()
                #eewt = eew.v1v2.copy().normalize()
                eeweewnproj = mpbb.project(eewvs,eewn) # projection of eew onto eewn; always a dot

                eewewnproj = mpbb.project(eewvs,ewn) # projection of eew onto ew
                eweewnproj = mpbb.project(ewvs,eewn) # projection of ew onto eew

                if mpbb.overlap(eewewnproj,ewewnproj) or\
                    mpbb.overlap(eweewnproj,eeweewnproj):
                    # if each wall is a dot in the other's normal projection, the overlap is parallel
                    if eewewnproj.x == eewewnproj.y and eweewnproj.x == eweewnproj.y:
                        # note that eewt and ewt are antiparallel typically
                        ewewtproj = mpbb.project(ewvs,ewt) # projection of ew onto ew tangent space
                        eewewtproj = mpbb.project(eewvs,ewt) # projection of eew onto ew tangent space
                        eweewtproj = mpbb.project(ewvs,eewt) # projection of ew onto eew tangent space
                        eeweewtproj = mpbb.project(eewvs,eewt) # projection of eew onto eew tangent space

                        # be sure they overlap in one another's tangent space
                        if mpbb.overlap(ewewtproj,eewewtproj) and\
                                mpbb.overlap(eweewtproj,eeweewtproj):

                            # make sure they are not just sharing an endpoint
                            if not ewewtproj.x == eewewtproj.y and\
                                    not ewewtproj.y == eewewtproj.x:

                                evals = ewvs + eewvs
            
                                # if theres a complete overlap
                                if ewewtproj.x == eewewtproj.x and\
                                        ewewtproj.y == eewewtproj.y:

                                    extrawalls.append(ew)
                                    etoiwalls.append(eew)

                                # if they share a single endpoint
                                elif ewewtproj.x == eewewtproj.x:
                                    if ewt == eewt:
                                        ipair,epair = resolve_verts(evals[0],evals[1],evals[3])
                                    elif ewt == cv.flip(eewt):
                                        ipair,epair = resolve_verts(evals[0],evals[1],evals[2])

                                    extrawalls.append(ew)
                                    extrawalls.append(eew)

                                    if ew.length > eew.length:
                                        epair.reverse()

                                    wmat = self.wall_material
                                    niw = wa.newwall(v1 = ipair[0],v2 = ipair[1],
                                            sector = sect,sort = 'interior',m = wmat,
                                            h = self.wall_height,fh = self.floor_height,
                                            w = self.wall_width)
                                    new = wa.newwall(v1 = epair[0],v2 = epair[1],
                                            sector = sect,sort = 'exterior',m = wmat,
                                            h = self.wall_height,fh = self.floor_height,
                                            w = self.wall_width)
                                    niwalls.append(niw)
                                    newalls.append(new)

                                # or perhaps the other endpoint
                                elif ewewtproj.y == eewewtproj.y:

                                    if ewt == eewt:
                                        ipair,epair = resolve_verts(evals[1],evals[0],evals[2])
                                    elif ewt == cv.flip(eewt):
                                        ipair,epair = resolve_verts(evals[1],evals[0],evals[3])

                                    extrawalls.append(ew)
                                    extrawalls.append(eew)

                                    if ew.length > eew.length:
                                        epair.reverse()

                                    wmat = self.wall_material
                                    niw = wa.newwall(v1 = ipair[0],v2 = ipair[1],
                                            sector = sect,sort = 'interior',m = wmat,
                                            h = self.wall_height,fh = self.floor_height,
                                            w = self.wall_width)
                                    new = wa.newwall(v1 = epair[0],v2 = epair[1],
                                            sector = sect,sort = 'exterior',m = wmat,
                                            h = self.wall_height,fh = self.floor_height,
                                            w = self.wall_width)
                                    niwalls.append(niw)
                                    newalls.append(new)

                                else:
                                    #print 'the complicated case'
                                    return False

                    # overlap is nonparallel
                    else:
                        ewewtproj = mpbb.project(ewvs,ewt) # projection of ew onto ew tangent space
                        eewewtproj = mpbb.project(eewvs,ewt) # projection of eew onto ew tangent space
                        eweewtproj = mpbb.project(ewvs,eewt) # projection of ew onto eew tangent space
                        eeweewtproj = mpbb.project(eewvs,eewt) # projection of eew onto eew tangent space

                        # be sure they overlap in one another's tangent space
                        if mpbb.overlap(ewewtproj,eewewtproj) and\
                                mpbb.overlap(eweewtproj,eeweewtproj):

                            if eweewtproj.x == eweewtproj.y and\
                                    (eweewtproj.x == eeweewtproj.x or\
                                        eweewtproj.x == eeweewtproj.y):
                                #print 'perp endpoint intersection...'
                                pass

                            elif eewewtproj.x == eewewtproj.y and\
                                    (eewewtproj.x == ewewtproj.x or\
                                        eewewtproj.x == ewewtproj.y):
                                #print 'perp endpoint intersection...'
                                pass

                            else:
                                #print 'perp overlap case'
                                return False

        for ew in extrawalls:
            if ew in ewalls: ewalls.remove(ew)
            if ew in iwalls: iwalls.remove(ew)
            if ew in eewalls: eewalls.remove(ew)
            if ew in eiwalls: eiwalls.remove(ew)
        [ewalls.append(new) for new in newalls]
        [iwalls.append(niw) for niw in niwalls]
        for ew in etoiwalls:self.switch_wall_sort(ew)
        return True

    #####



    # plans is a tuple of 4 lists of plans
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
        xpj = dpv.project_coords(newcorners,dpv.xhat)
        ypj = dpv.project_coords(newcorners,dpv.yhat)
        l,w = xpj.y-xpj.x,ypj.y-ypj.x
        margs = ((),{'x':x,'y':y,'l':l,'w':w,'shafted':False})
        self.create_bbox(margs)
        if not self.verify_growth(plans,margs):return False
        splan = self.should_shaft(plans,margs)
        cpairs = [(c2,c3),(c3,c4),(c4,c1)]
        iwargs = []
        ewargs = [[cp,
            {'h':4.0,'w':0.5,'walltype':'exterior','room':margs}] 
                for cp in cpairs]
        self.switch_walltype(plans,side)
        if self.resolve_walls(plans,ewargs,iwargs,margs):
            rps,eps,ips,sps = plans
            rps.append(margs)
            eps.extend(ewargs)
            ips.extend(iwargs)
            if splan:sps.append(splan)
            return True
        else:return False

    # alter plans on a floor specific basis before building them
    def plan_specific(self,specific):
        if specific > 0:self.allplans[1][0][1]['walltype'] = 'exterior'
        if specific == self.bldg.stories-1:
            for sp in self.allplans[0]:
                sp[1]['cgap'] = None

    # generate all the plans necessary to build a floor
    def plan(self):
        plans = self.entry()
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
            #for n in r.sgraph.nodes:n._assign_material('concrete2')
            self._consume(r)
        #for ew in ewalls:ew._assign_material('brick2')
        #for iw in iwalls:iw._assign_material('concrete2')
        self._nodes_to_graph(
            self._node_wrap(*ewalls),
            self._node_wrap(*iwalls))


