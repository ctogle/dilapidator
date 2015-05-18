import dilap.core.model as dmo
import dilap.core.tools as dpr

import dp_vector as dpv

import pdb

###############################################################################
### terrain represents a piece of deterministic terrain
### terrain is a triangle when xy-projected
###############################################################################

class terrain(dmo.model):
    def __init__(self,*args,**kwargs):
        dmo.model.__init__(self,*args,**kwargs)
        self._def('mesh',None,**kwargs)
        self._geo()

    def _geo(self):
        pxy = self.mesh._modeldata()
        self.pcoords = pxy['pcoords']
        self.ncoords = pxy['ncoords']
        self.ucoords = pxy['ucoords']
        self.faces = pxy['faces']
        self.face_mats = pxy['face_mats']

def region_pts_to_boundary(rpts,radius = 100):
    convex = pts_to_convex_xy(rpts)
    inflate(convex,radius)
    return convex

def some_input():
    ls = [5, 10, 15]
    ws = [2, 4, 6, 8]
    base = dpv.vector(50,50,0)
    fixed_pts = []
    for sq in range(5):
        pt = base.copy().translate_x(sq*25).translate_y(sq*10)
        l,w = rm.choice(ls),rm.choice(ws)
        corners = mpu.make_corners(pt,l,w,0)
        fixed_pts.extend(corners)

    hole_corners = [pt.copy() for pt in fixed_pts]

    region_pts = []
    for lpt in range(5):
        pt = dpv.zero().translate_x(lpt*100).translate_y((4*lpt-lpt**2)*100)
        region_pts.append(pt)
    region_pts.append(dpv.vector(200,-100,0))
    for rpt in region_pts:print(rpt)

    target_polygon_edge_length = 10
    target_primitive_edge_length = 200

    theinput = {
        'fixed_pts':fixed_pts, 
        'hole_pts':hole_corners, 
        'region_pts':region_pts, 
        'polygon_edge_length':target_polygon_edge_length, 
        'primitive_edge_length':target_primitive_edge_length, 
            }
    return theinput

def test():
    someinput = some_input()
    boundary = region_pts_to_boundary(someinput['region_pts'])
    fixed = someinput['fixed_pts']

    mesh_corners = triangle_cover(boundary,someinput['primitive_edge_length'])
    mesh_corners = intersects_xy(mesh_corners,boundary)

    map_curves = [(mc,'green',None,0.75) for mc in mesh_corners]
    map_curves.append((boundary,'black','o',1.5))
    make_map(map_curves)



