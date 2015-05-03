import dilap.core.base as db
import dilap.core.model as dmo
import dilap.core.tools as dpr
import dilap.core.context as dgc
import dilap.primitive.cube as dcu
import dilap.primitive.stairs as ds

import dp_vector as dpv

class shaft(dgc.context):

    def __init__(self,*args,**kwargs):
        dgc.context.__init__(self,*args,**kwargs)
        self._def('style','uturn',**kwargs)
        self._def('p',dpv.zero(),**kwargs)
        self._def('l',12,**kwargs)
        self._def('w',16,**kwargs)
        if args:
            b = args[0]
            s = b.stories
            fheights,cheights,wheights = b.fheights,b.cheights,b.wheights
        else:
            s = 3
            fheights,cheights,wheights = [0.25]*s,[0.25]*s,[4.0]*s
        self._def('stories',s,**kwargs)
        self._def('fheights',fheights,**kwargs)
        self._def('cheights',cheights,**kwargs)
        self._def('wheights',wheights,**kwargs)

    # helper function to get story-dep. info
    def _params(self,level):
        fh = self.fheights[level]
        wh = self.wheights[level]
        ch = self.cheights[level]
        l,w = self.l,self.w
        return l,w,fh,wh,ch

    # generate and return a node representing stairs of one story
    def generate_uturn(self,level):
        l,w,fh,wh,ch = self._params(level)

        p1h = fh;p1l = l;p1w = 3.0
        p1x = 0.0;p1y = (p1w - w)/2.0;p1z = 0
        pform1 = dcu.cube().translate_z(0.5)
        pform1.scale(dpv.vector(p1l,p1w,p1h))
        pform1.translate(dpv.vector(p1x,p1y,p1z))

        p2h = ch;p2l = l;p2w = 3.0
        p2x = 0.0;p2y = (p2w - w)/2.0;p2z = wh - ch
        pform2 = dcu.cube().translate_z(0.5)
        pform2.scale(dpv.vector(p2l,p2w,p2h))
        pform2.translate(dpv.vector(p2x,p2y,p2z))

        gap = l/5.0
        s = int(wh)
        rw = 2.0*gap
        rl = w - 2.0*p1w
        diff = wh/2.0
        rwoffx = l/2.0 - rw/2.0
        rwoffy = rl/2.0
        sheight = diff/s

        p3h = 2.0*sheight;p3l = l;p3w = 3.0
        p3x = 0.0;p3y = (w - p3w)/2.0;p3z = diff-sheight
        pform3 = dcu.cube().translate_z(0.5)
        pform3.scale(dpv.vector(p3l,p3w,p3h))
        pform3.translate(dpv.vector(p3x,p3y,p3z))

        extra = dcu.cube().translate_z(0.5)
        extra.scale(dpv.vector(gap,rl,wh))
        extra.translate(dpv.vector(0,0,0))

        sopts1 = {'steps':s,'l':rl,'w':rw,'h':diff}
        sopts2 = {'steps':s,'l':rl,'w':rw,'h':diff}
        lside = ds.stairs(**sopts1)
        rside = ds.stairs(**sopts2)
        #lside.rotate_z(dpr.PI).translate_y(rl).translate_z(diff)
        lside.translate_y(rl).translate_z(diff)
        lside.translate_x(-rwoffx).translate_y(-rwoffy).translate_z(fh-sheight)
        rside.translate_x( rwoffx).translate_y(-rwoffy).translate_z(fh-sheight)
        final = dpr.combine([pform1,pform2,pform3,lside,rside,extra])
        return self._node_wrap(final)

    # generate and return a node representing ramps of one story
    def generate_switchback(self,level):
        floor = dmo.model()
        pdb.set_trace()
        return self._node_wrap(floor)

    ##########
    # generate and return a node representing walls of one story
    def generate_walls(self,fldex):
        l,w,fh,wh,ch = self._params(fldex)
        cs = mpu.make_corners(cv.zero(),l,w,0)
        wargs = [
            {'v1':cs[1],'v2':cs[2], 
            'solid':True,'m':'brick2',
            'h':wh,'fh':fh,'w':0.25}, 
            {'v1':cs[2],'v2':cs[3], 
            'solid':True,'m':'brick2',
            'h':wh,'fh':fh,'w':0.25}, 
            {'v1':cs[3],'v2':cs[0], 
            'solid':True,'m':'brick2',
            'h':wh,'fh':fh,'w':0.25},
            {'v1':cs[0],'v2':cs[1],
            'sort':'interior','solid':False,'m':'brick2',
            'h':wh,'fh':fh,'w':0.25}]
        wals = [wa.newwall(**w) for w in wargs]
        [w._face_away(cv.zero()) for w in wals]
        [w._build() for w in wals]
        return pr.sum_primitives([w._primitives() for w in wals])
    ##########

    # generate and return a node representing one story
    def generate_story(self,level):
        if self.style == 'uturn':
            floor = self.generate_uturn(level)
        elif self.style == 'switchback':
            floor = self.generate_switchback(level)
        #floor._add_child(self.generate_walls(level))
        return floor

    def generate(self,worn = 0):
        p = self.p
        zoff = 0.0
        for fx in range(self.stories):
            np = self.generate_story(fx)
            np.translate(p).translate_z(zoff)
            self._nodes_to_graph(np)
            zoff += self.wheights[fx]


