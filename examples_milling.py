
"""
from pygfx.render import *
from pygfx.point_ops import *
from pygfx.math_ops import  *
from pygfx.obj3d import  *
"""

""" 
    A nifty free tool to view gcode:

    https://ncviewer.com/

"""



#from pygfx.obj3d import  *

from pygfx.milling_ops import * 


gc_poly = gcode_to_polyline()


#gc_poly.load_gcode_textfile('gcode/1001.txt')

gc_poly.load_gcode('gcode/king.txt')

#gc_poly.save_gcode('gcode/fubar.txt')

# pts_to_linesegment

#gc_poly.show_data()
gc_poly.save_3d_object('gcode_3d.obj')





