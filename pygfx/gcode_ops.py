

import math 
import re 
import os 

from gnelscript import SHAPELY_IS_LOADED

from gnelscript.pygfx.obj3d import object3d
from gnelscript.pygfx.obj2d import object2d
from gnelscript.pygfx.point_ops import *
from gnelscript.pygfx.math_ops import *
from gnelscript.pygfx.vector_ops import *
from gnelscript.pygfx.grid_ops import tessellator

#from pygfx.gcode.bridgeport import  parser_commands
from gnelscript.pygfx.gcode.linuxcnc import  parser_commands


##-----------------------------## 

#I resisted using external libraries too long - shapely and trimesh look awesome 
import networkx as nx

if GEOJSON_IS_LOADED:
    from geojson import dump  
    from geojson import Point as gjpt
    from geojson import Polygon as gjply
    from geojson import Feature as gjftr
    from geojson import LineString as gjln
    from geojson import FeatureCollection as gjfc

if SHAPELY_IS_LOADED:
    from shapely import buffer, BufferCapStyle, BufferJoinStyle
    from shapely import Point as shp_pt
    from shapely import MultiPolygon as shp_mply
    from shapely import Polygon as shp_ply
    from shapely import LineString as shp_ln
    #from shapely import Feature as shp_ftr
    #from shapely import FeatureCollection as shp_fc

import trimesh
import numpy as np




