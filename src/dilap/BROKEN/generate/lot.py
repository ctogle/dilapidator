import dilap.core.tools as dpr
import dilap.core.vector as dpv
import dilap.core.quaternion as dpq
import dilap.core.context as dgc
import dilap.structures.house as dlh

class lot(dgc.context):

    def __init__(self,*args,**kwargs):
        dgc.context.__init__(self,*args,**kwargs)
        self._def('boundary',None,**kwargs)

    def _transform(self,t,q,s):
        q = dpq.q_from_av(q,dpv.z())
        dgc.context._transform(self,t,q,s)
        dpv.scale_coords(self.terrain_points,s)
        dpv.rotate_coords(self.terrain_points,q)
        dpv.translate_coords(self.terrain_points,t)

    def _terrain_points(self):
        tpts = []
        [tpts.extend(s._terrain_points()) for s in self.structures]
        self.terrain_points = tpts

    def generate(self,worn = 0):
        e1 = dpv.v1_v2(self.boundary[0],self.boundary[1])
        e2 = dpv.v1_v2(self.boundary[1],self.boundary[2])
        a = dpr.angle_from_xaxis_xy(e1)
        t = dpv.center_of_mass(list(self.boundary))
        r = dpq.q_from_av(a,dpv.z())
        tf,rf = t.copy().flip(),r.copy().flip()
        hl = e1.magnitude()
        hw = e2.magnitude()

        hbnd = [b.copy().translate(tf).rotate(rf) for b in self.boundary]
        house = dlh.house(hbnd,l = hl,w = hw)
        self.structures = [house]

        for s in self.structures:
            node = self._node_wrap(s.model(t,r))
            self._nodes_to_graph(node)
        return self

   

