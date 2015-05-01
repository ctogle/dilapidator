import dilap.core.base as db
import dilap.core.tform as dtf

#import make_places.core.primitives as pr

import dp_vector as dpv

###############################################################################
### sgraph represents a scenegraph and contains references to all
### top level nodes in the graph
###############################################################################

class sgraph(db.base):

    def __str__(self):
        strr = '\n' + self.__str__() + '\t'
        strr += '\tscenegraph:\n'
        if not self.nodes:strr += '\tEMPTY\n'
        else:
            strr += '\t\ttop-level-nodes:'
            strr += ','.join([n.name for n in self.nodes])
            strr += '\n'.join([n.__str__() for n in self.nodes])
        return strr

    def __init__(self,*args,**kwargs):
        self._def('nodes',[],**kwargs)

    def make_scene(self, scenetype, center = False):
        for nd in self.nodes:
            nd.make(scenetype = scenetype,center = center)

###############################################################################
###############################################################################

###############################################################################
### dpnode is the basic scenegraph node of dilap
### it references a transform that contains references to children and parent
###############################################################################

unused_dpnode_id = 0
class dpnode(db.base):

    def _dpid(self):
        global unused_dpnode_id
        self.dpid = unused_dpnode_id
        unused_dpnode_id += 1
        return self.dpid
    
    def _name(self):
        self.name = 'node.'+str(self.dpid)
        return self.name

    def _def_tform(self,*args,**kwargs):
        if hasattr(self,'tform'):return
        kweys = kwargs.keys()
        pos = kwargs[ke] if 'pos' in kweys else dpv.zero()
        rot = kwargs[ke] if 'rot' in kweys else dpv.zero()
        scl = kwargs[ke] if 'scl' in kweys else dpv.one()
        tpar = kwargs[ke] if 'parent' in kweys else None
        tchi = kwargs[ke] if 'children' in kweys else []

        ntf = dtf.tform(self,parent = tpar,
            pos = pos,rot = rot,scl = scl,
            children = [ch.tform for ch in tchi])
        self.tform = ntf

    def _def_uv_tform(self,*args,**kwargs):
        if hasattr(self,'uv_tform'):return
        kweys = kwargs.keys()
        pos = kwargs[ke] if 'uv_pos' in kweys else dpv.zero()
        rot = kwargs[ke] if 'uv_rot' in kweys else dpv.zero()
        scl = kwargs[ke] if 'uv_scl' in kweys else dpv.one()
        tpar = kwargs[ke] if 'uv_parent' in kweys else None
        tchi = kwargs[ke] if 'uv_children' in kweys else []

        ntf = dtf.tform(self,parent = tpar,
            pos = pos,rot = rot,scl = scl,
            children = [ch.uv_tform for ch in tchi])
        self.uv_tform = ntf

    def __str__(self,b = -1):
        b += 1
        tf = self.tform
        pcnt = str(len(self.models))
        ccnt = str(len(tf.children))
        strr = '\t'*b + 'node:' + str(self.name) + '\n'
        strr += tf.__str__(b) + '\n'
        strr += '\t'*b + 'with ' + pcnt + ' primitives '
        strr += 'and ' + ccnt + ' children'
        childstr = [c.owner.__str__(b) for c in tf.children]
        strr += '\n' + '\n'.join(childstr)
        return strr

    def __init__(self,*args,**kwargs):
        self._dpid()
        self._name()
        self._def('models',[],**kwargs)
        self._def('lod_models',[],**kwargs)
        self._def_tform(*args,**kwargs)
        self._def_uv_tform(*args,**kwargs)

        self._def('consumption',False,**kwargs)

    def _add_child(self,*chil):
        for ch in chil:
            ch.tform.parent = self.tform
            if not ch in self.tform.children:
                self.tform.children.append(ch.tform)
    
    def _remove_child(self,*chil):
        for ch in chil:
            ch.tform.parent = None
            if ch.tform in self.tform.children:
                self.tform.children.remove(ch.tform)

    def _set_parent(self,parent):
        self.tform.parent = parent.tform
        if not self.tform in parent.tform.children:
            parent.tform.children.append(self.tform)










    def assign_material(self, mat, propagate = True):
        if propagate:
            for ch in self.tform.children:
                ch.owner.assign_material(mat)
        for p in self.primitives: p.assign_material(mat)
        for p in self.lod_primitives: p.assign_material(mat)

    def worldly_primitive(self, prim, ttf, uv_ttf, **kwargs):
        tpm = prim
        tpm.scale(ttf.scales)
        tpm.worldly_uvs(uv_ttf)
        tpm.rotate_z(ttf.rotation.z)
        tpm.translate(ttf.position)
        kwargs['name'] = self.name
        kwargs['rdist'] = self.grit_renderingdistance
        kwargs['lodrdist'] = self.grit_lod_renderingdistance
        return tpm, kwargs

    def worldly_children(self, **kwargs):
        ttf = self.tform.true()
        uv_ttf = self.uv_tform.true()
        newpms = []
        for pm in self.models:
            tpm,kwargs = self.worldly_primitive(pm,ttf,uv_ttf,**kwargs)
            newpms.append(tpm)
        for ch in self.tform.children:
            chpms = ch.owner.worldly_children()
            newpms.extend(chpms)
        return newpms

    def lod_worldly_children(self, **kwargs):
        ttf = self.tform.true()
        uv_ttf = self.uv_tform.true()
        newpms = []
        for pm in self.lod_primitives:
            tpm, kwargs = self.worldly_primitive(pm,ttf,uv_ttf,**kwargs)
            newpms.append(tpm)
        for ch in self.tform.children:
            chpms = ch.owner.lod_worldly_children()
            newpms.extend(chpms)
        return newpms

    # consume geometry of children into this nodes models/lods
    def _consume(self):
        chps = []
        chlps = []
        for ch in self.tform.children:
            ch.owner.consume()
            ch.parent = None
            chps.extend(ch.owner.worldly_children())
            chlps.extend(ch.owner.lod_worldly_children())
        self.tform.children = []
        self.models.extend(chps)
        self.lod_models.extend(chlps)
        if self.models:
            final_prim = pr.sum_primitives(self.models)
            self.models = [final_prim]
        if self.lod_primitives:
            final_lod_prim = pr.sum_primitives(self.lod_primitives)
            self.lod_models = [final_lod_prim]

            if self.models:
                self.models[0].has_lod = True
            self.lod_models[0].is_lod = True

    # return models and lods in 1-1 with Nones filled in
    # also return total model count
    def _model_listing(self):
        pcnt = len(self.models)
        lcnt = len(self.lod_models)
        ms = self.models[:]
        ls = self.lod_models[:]
        if pcnt < lcnt:
            ms.extend([None]*(lcnt-pcnt))
            apcnt = lcnt
        elif pcnt > lcnt:
            ls.extend([None]*(pcnt-lcnt))
            apcnt = pcnt
        else:apcnt = pcnt
        return ms,ls,apcnt

    def _modelize(self,iotype,**kwargs):
        ttf = self.tform.true()
        uv_ttf = self.uv_tform.true()

        models,lods,mcnt = self._model_listing()
        for pmdx in range(mcnt):
            pm = models[pmdx]
            lpm = lods[pmdx]

            if not pm is None:
                tpm,kwargs = self.worldly_primitive(pm,ttf,uv_ttf)
                
                iotype.create_primitive(tpm,**kwargs)

            if not lpm is None:
                tpm,kwargs = self.worldly_primitive(lpm,ttf,uv_ttf)
                tpm.is_lod = True

                iotype.create_primitive(tpm,**kwargs)

    # given an io module, produce model outputs 
    # for this node and its children
    def _realize(self,iotype):
        if self.consumption:self._consume()
        else:
            for ch in self.tform.children:
                ch.owner._realize(iotype)
        self._modelize(iotype)


