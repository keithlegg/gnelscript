
 
#from gnelscript.pygfx.render import *
from gnelscript.pygfx.point_ops import *
from gnelscript.pygfx.math_ops import  *
from gnelscript.pygfx.obj3d import  *
 

""" 
    A nifty free tool to view gcode:
    https://ncviewer.com/

    //gear generator  
    https://woodgears.ca/gear_cutting/index.html

"""

 

from gnelscript.pygfx.kicad_ops import * 
from gnelscript.pygfx.milling_ops import * 

#gc_poly = gcode()

#gc_poly.load_gcode('gcode/ngc/3dtest.ngc')
#gc_poly.load_gcode('gcode/ngc/arcspiral.ngc')
#gc_poly.load_gcode('gcode/ngc/calibrate_5X1.ngc')
#gc_poly.load_gcode('gcode/ngc/3D_Chips.ngc')
#print( gc_poly.commented )
#print( gc_poly.param_names, gc_poly.param_values )
#gc_poly.save_gcode('gcode/fubar.txt')
#gc_poly.show_data()
#gc_poly.save_3d_object('gcode_3d.obj')


class gear_generator(object3d):
    def __init__(self):
        super().__init__()         
        #self.contact_angle
        self.tooth_spacing = 20 
        self.num_teeth = 5
        self.shaft_hole_dia = 4


    def build(self, shaftdia, dia, teeth_height, numteeth ):
        outer = self.calc_circle(pos=(0,0,0), rot=(0,0,0), dia=(dia+teeth_height), axis='z', periodic=True, spokes=36)
        inner = self.calc_circle(pos=(0,0,0), rot=(0,0,0), dia=dia               , axis='z', periodic=True, spokes=36)
        shaft = self.calc_circle(pos=(0,0,0), rot=(0,0,0), dia=shaftdia          , axis='z', periodic=True, spokes=36)

        outer.extend(inner)
        outer.extend(shaft)
        
        return outer 


#-------------------------------

class gear_mesh_generator(object3d):
    def __init__(self):
        self.gear1 = gear_generator()
        self.gear2 = gear_generator()



#-------------------------------

 
#kicad = pcbfile()

#kicad.load_gcode('gcode/ngc/3D_Chips.ngc')
#kicad.load_kicadpcb('gcode/kicad/sample1.kicad_pcb')
#kicad.load_kicadpcb('gcode/kicad/zipper.kicad_pcb')
#kicad.load_kicadpcb('gcode/kicad/simple.kicad_pcb')
#kicad.load_kicadpcb('gcode/kicad/cnc1.kicad_pcb')
# kicad.show_geom()
#kicad.save_3d_obj('kicad.obj') 
#kicad.show_modules()
#-------------------------------



