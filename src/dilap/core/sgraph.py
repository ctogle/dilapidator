import dilap.core.base as db
import dilap.core.tform as dtf
import dilap.core.tools as dpr

import dp_vector as dpv
import dp_quaternion as dpq

import pdb

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

    def to_world(self):
        for nd in self.nodes:
            for ch in nd.tform.children:
                ch.owner._models_to_world()
            nd._models_to_world()

    def to_local(self):
        for nd in self.nodes:
            for ch in nd.tform.children:
                ch.owner._models_to_local()
            nd._models_to_local()

    def graph(self,iotype):
        iotype.write_materials()
        for nd in self.nodes:
            if nd.consumption:nd._consume()
            else:
                for ch in nd.tform.children:
                    ch.owner._realize(iotype)
            nd._realize(iotype)

###############################################################################
###############################################################################

###############################################################################
### node is the basic scenegraph node of dilap
### it references a transform that contains references to children and parent
###############################################################################

unused_dpnode_id = 0
class node(db.base):

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
        pos = kwargs['pos'] if 'pos' in kweys else dpv.zero()
        #rot = kwargs['rot'] if 'rot' in kweys else dpv.zero()
        rot = kwargs['rot'] if 'rot' in kweys else dpq.zero()
        scl = kwargs['scl'] if 'scl' in kweys else dpv.one()
        tpar = kwargs['parent'] if 'parent' in kweys else None
        tchi = kwargs['children'] if 'children' in kweys else []

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

    def translate(self,v):
        self.tform.pos.translate(v)
        return self

    def translate_x(self,dx):
        self.tform.pos.translate_x(dx)
        return self

    def translate_y(self,dy):
        self.tform.pos.translate_y(dy)
        return self

    def translate_z(self,dz):
        self.tform.pos.translate_z(dz)
        return self

    def scale(self,s):
        self.tform.scl.scale(s)
        return self

    def scale_x(self,sx):
        self.tform.scl.x *= sx
        return self

    def scale_y(self,sy):
        self.tform.scl.y *= sy
        return self

    def scale_z(self,sz):
        self.tform.scl.z *= sz
        return self

    def rotate(self,q):
        print('its time to try the quaternion')
        pdb.set_trace()

        self.tform.pos.translate(v)
        return self

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

    # this changes the model to local space given a world space tform
    def _transform_to_local(self,mod,ttf,uv_ttf):
        mod.translate(ttf.pos.copy().flip())
        #tpm.rotate_z(ttf.rot.copy().flip().z)
        mod.scale(ttf.scl.copy().reciprocate())
        mod._uvs_to_local(uv_ttf)

    # this changes the model to world space given a world space tform
    def _transform_to_world(self,mod,ttf,uv_ttf):
        mod.scale(ttf.scl)
        mod._uvs_to_world(uv_ttf)
        #mod.rotate_z(ttf.rot.z)
        mod.translate(ttf.pos)

    # transform models and childrens models to world space
    def _transform_to_world_walk(self,lod):
        ttf = self.tform.true()
        uv_ttf = self.uv_tform.true()
        newpms = []
        if lod:which = self.lod_models
        else:which = self.models
        for pm in which:
            self._transform_to_world(pm,ttf,uv_ttf)
            newpms.append(pm)
        for ch in self.tform.children:
            chpms = ch.owner._transform_to_world_walk(lod)
            newpms.extend(chpms)
        return newpms

    # assign a material mat to all faces of all models/lods
    # if propagate, also assign the material recursively to children
    def _assign_material(self,mat,propagate = True):
        if propagate:
            for ch in self.tform.children:
                ch.owner._assign_material(mat)
        for p in self.models: p._assign_material(mat)
        for p in self.lod_models: p._assign_material(mat)

    # consume geometry of children into this nodes models/lods
    # combine all models/lods into one model/lod
    # this is irreverisble so long as _transform_to_world_walk is
    def _consume(self):
        chps = []
        chlps = []
        for ch in self.tform.children:
            ch.owner._consume()
            ch.parent = None
            chps.extend(ch.owner._transform_to_world_walk(False))
            chlps.extend(ch.owner._transform_to_world_walk(True))
        self.tform.children = []
        self.models.extend(chps)
        self.lod_models.extend(chlps)
        if self.models:
            final_prim = dpr.combine(self.models)
            self.models = [final_prim]
        if self.lod_models:
            final_lod_prim = dpr.combine(self.lod_models)
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

    # more this nodes models to world space
    def _models_to_world(self):
        ttf = self.tform.true()
        uv_ttf = self.uv_tform.true()
        models,lods,mcnt = self._model_listing()
        for pmdx in range(mcnt):
            pm = models[pmdx]
            lpm = lods[pmdx]
            if not pm is None:self._transform_to_world(pm,ttf,uv_ttf)
            if not lpm is None:self._transform_to_world(lpm,ttf,uv_ttf)

    # more this nodes models to local space
    def _models_to_local(self):
        ttf = self.tform.true()
        uv_ttf = self.uv_tform.true()
        models,lods,mcnt = self._model_listing()
        for pmdx in range(mcnt):
            pm = models[pmdx]
            lpm = lods[pmdx]
            if not pm is None:self._transform_to_local(pm,ttf,uv_ttf)
            if not lpm is None:self._transform_to_local(lpm,ttf,uv_ttf)
    
    # transform models/lods to world space
    # use iotype to build the models
    # transform models/lods back to local spaces
    def _build_models(self,iotype,**kwargs):
        models,lods,mcnt = self._model_listing()
        for pmdx in range(mcnt):
            pm = models[pmdx]
            lpm = lods[pmdx]
            if not pm is None:iotype.build_model(pm,**kwargs)
            if not lpm is None:iotype.build_model(lpm,**kwargs)

    # transform models/lods to world space
    # use iotype to build the models
    # transform models/lods back to local spaces
    def _modelize(self,iotype,**kwargs):
        ttf = self.tform.true()
        uv_ttf = self.uv_tform.true()
        models,lods,mcnt = self._model_listing()
        for pmdx in range(mcnt):
            pm = models[pmdx]
            lpm = lods[pmdx]
            if not pm is None:
                self._transform_to_world(pm,ttf,uv_ttf)
                iotype.build_model(pm,**kwargs)
                self._transform_to_local(pm,ttf,uv_ttf)
            if not lpm is None:
                self._transform_to_world(lpm,ttf,uv_ttf)
                #lpm.is_lod = True
                iotype.build_model(lpm,**kwargs)
                self._transform_to_local(lpm,ttf,uv_ttf)

    # given an io module, produce model outputs 
    # for this node and its children
    def _realize(self,iotype):
        #if self.consumption:self._consume()
        #else:
        #    for ch in self.tform.children:
        #        ch.owner._realize(iotype)
        self._modelize(iotype)


