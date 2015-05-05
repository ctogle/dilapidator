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

import dp_vector as dpv

import random as rm

import pdb

class floorplan(dgc.context):

    def __init__(self,building,*args,**kwargs):
        dgc.context.__init__(self,*args,**kwargs)
        self.bldg = building
        self._def('max_rooms',10,**kwargs)

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
            shl = 10
            shw = 16
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
        ewargs = [[cp,{'h':4.0,'w':0.5,'switchable':True}] for cp in cpairs]
        # 0 may not be the front wall, i did not check
        ewargs[0][1]['switchable'] = False

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
                if not which[1]['switchable']:
                    side = which
                    break
                options.remove(which)
        return side

    # choose a grow length if necessary
    def grow_length(self,plans,length):
        if length is None:return rm.choice([8,12,16,20,24,28,32])
        else:return length

    # plans is a tuple of 4 lists of plans
    #def grow(self,length = None,side = None,force = False):
    def grow(self,plans,length = None,side = None,force = False):
        side = self.grow_side(plans,side)
        if side is None:return False
        gleng = self.grow_length(plans,length)

        #bdist = side._distance_to_border(self.corners)
        #if bdist < 8 and not force:
        #    #print 'too close to a border to grow'
        #    return False
        #elif gleng > bdist: gleng = bdist

        pdb.set_trace()
        # THIS IS WHERE IT GOES OFF THE RAILS


        side._face_away()
        c1 = side.v2.copy()
        c2 = side.v1.copy()
        cn = side.normal.copy()
        c3,c4 = mbp.extrude_edge(c1,c2,gleng,cn)
        newcorners = [c1,c2,c3,c4]
        sect = fl.floor_sector(corners = newcorners, 
                floor_height = self.floor_height, 
                ceiling_height = self.ceiling_height, 
                wall_height = self.wall_height)
        gaps = self.should_shaft(sect.position,sect.length,sect.width)
        sect.fgaps = gaps[:]
        sect.cgaps = gaps[:]
        if gaps: sect.shafted = True
        for esect in self.sectors:
            ebb = esect.get_bboxes()
            nbb =  sect.get_bboxes()
            if mpbb.intersects(ebb,nbb):
                if gaps: self.shaft_kwargs.pop(-1)
                #print 'new sect intersected!'
                return False

        wmat = self.wall_material
        cpairs = [(c2,c3),(c3,c4),(c4,c1)]
        extwalls = [wa.newwall(v1 = cp[0],v2 = cp[1],sort = 'exterior', 
                    m = wmat,sector = sect,h = self.wall_height,
                    fh = self.floor_height,w = self.wall_width) 
                        for cp in cpairs]
        intwalls = []

        #cpairs = [(c2,c3),(c3,c4),(c4,c1)]
        #extwalls = [wa.wall_plan(*cp, 
        #    sector = sect, 
        #    sort = 'exterior', 
        #    wall_height = self.wall_height) 
        #        for cp in cpairs] 
        #intwalls = []
        #if self.resolve_walls(extwalls,intwalls,sect):
        if True:
            self.switch_wall_sort(side)

            self.sectors.append(sect)
            self.exterior_walls.extend(extwalls)
            self.interior_walls.extend(intwalls)
            return True

        else: return False


    #####



    def plan(self):
        plans = self.entry()
        
        self.allplans = plans
        return

        for sect in range(self.max_rooms):
            grew = self.grow(plans)
            if grew is false:
                print('growth rejected')
                pass
            elif grew is none:
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
            self._consume(r)
        self._nodes_to_graph(
            self._node_wrap(*ewalls),
            self._node_wrap(*iwalls))


