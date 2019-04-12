
"""
from pygfx.render import *
from pygfx.point_ops import *
from pygfx.math_ops import  *
from pygfx.obj3d import  *
"""

#from pygfx.obj3d import  *

from pygfx.milling_ops import * 


gc_asm = gcode_assembly()


#gc_asm.load_gcode_textfile('gcode/1001.txt')

gc_asm.load_gcode_textfile('gcode/engrave.txt')

#gc_asm.show_data()
gc_asm.save_3d_object('gcode_3d.obj')





