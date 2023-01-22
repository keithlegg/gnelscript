

from gnelscript.pygfx.render import *
from gnelscript.pygfx.point_ops import *
from gnelscript.pygfx.math_ops import  *
from gnelscript.pygfx.obj3d import  *







pt_op = point_operator()


"""

pt_op.tuple_pop(self, listTuples, tup):


def add_margin_bbox(self, bbox, size):

def extents_fr_bbox(self, bbox, offset=None):

def closest_to_axis(self, points, val, axis):

def calc_square_diag(self, tl, br):

def calc_square_pt(self, size, origin=None ):

def calc_bbox_pt(self, size, origin=None ):

def calc_circle(self, pos=(0,0,0), rot=(0,0,0), dia=1, axis='z', periodic=True, spokes=23):

def triangle_mean_z(self, triangle):

def cvt_2d_to_3d(self, points):

def locate_pt_along3d(self, x1, y1, z1, x2, y2, z2, num):


"""


############################################################################################################################
############################################################################################################################



poly_op = polygon_operator() 

"""


def z_sort(self, reverse=False):

def flush(self):

def inspect_geom(self, geom):

def verify_geom(self, geom):

def calc_bbox(self, prgrp=None, facgrp=None ):

def _reindex_ply(self, f_idx, offset):

def centroid(self, pts):

def calc_tripoly_normal(self, three_pts, unitlen):

def three_vec3_to_normal(self, v1, v2, v3, unitlen=False):

def any_pt_is_near(self, pt_list, pt2, dist ):

def pt_is_near(self, pt1, pt2, dist ):

def select_by_location(self, select_type, pt_two, dist):

def indexer(self, slice=None, ids=None):

def geom_to_ptgrp(self, geom):

def sub_select_geom(self, slice=None, ids=None, reindex=False):

def get_pt_ids(self, fids=None):

def get_face_pts(self, fid):

def insert_pt_grp(self, ptgrp):

def get_pt_grp(self, slice=None, ids=None):

def get_face_group(self, slice=None, ids=None):

def get_geom_edges(self, geom ):

def get_face_edges(self, fid, reindex=False, geom=None):

def get_face_geom(self, fid,  reindex=False, geom=None):

def get_face_normal(self, fid=None, unitlen=False ):

def get_face_centroid(self, fid):

def apply_matrix_pts(self, pts, m33=None, m44=None):

def apply_matrix_ptgrp(self, ptgrp, m33=None, m44=None):

def scale_pts(self, amt, pts=None, ptgrp=None ):

def rotate_pts(self, rot, pts=None, ptgrp=None):

def insert_polygons(self, plyids, points, asnew_shell=True, geom=None):

def extrude_face(self, f_id, distance):

def copy_sop(self, slice=None, ids=None, reindex=False, offset=(0,1,0), rot=(0,0,0), num=2, distance=2):

def xform_pts(self, pos, pts=None, ptgrp=None):

def radial_triangulate_face(self, fid, offset=None, as_new_obj=False ):

def radial_triangulate_obj(self, as_new_obj=False, offset=None ):

def triangulate(self, force=False, offset=(0,0,0)):

def poly_loft(self, obj2, as_new_obj=True):

def repair(self):

def load(self, filename, doflush=True):

def save(self, filename, as_lines=False):
    
"""



############################################################################################################################
############################################################################################################################


obj = object3d() 

""" 
obj.reset() 

obj.copy()

obj.insert( obj, replace=False)

obj.append(otherobj)

obj.show()

obj.numply

obj.numpts

obj.numfacnrml

obj.convert_pts_vec3

obj.calc_face_normals

obj.one_vec_to_arrow(r3, pos=None):

obj.one_vec_to_obj(r3, pos=None, arrowhead=False):

obj.vectorlist_to_obj( vecs, pos=None):


"""

###############################

"""
obj.prim_line(axis, pos, rot, size=1)

obj.prim_triangle(axis, pos, rot, size=1)

obj.prim_quad(axis, pos, rot, size=1)

def prim_circle(self, axis, pos=(0,0,0), rot=(0,0,0), dia=1, spokes=9):

def prim_cone(self, axis, pos=(0,0,0), rot=(0,0,0), dia=1, spokes=8):

def prim_sphere(self, pos, rot, size=1 ):

def prim_locator(self, pos=(0,0,0), rot=(0,0,0), size=1):

def prim_locator_xyz(self, pos, rot, size=1):

def prim_arrow(self, axis='y', pos=(0,0,0), rot=(0,0,0), vec3=None, size=1, pivot='obj'): 

def prim_cube(self, linecolor=(0,128,0), pos=(0,0,0), rot=(0,0,0), size=1, pivot='obj'):

"""


############################################################################################################################
############################################################################################################################






