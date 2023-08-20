
 
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



def dagoptimize(filename):
    kip = vectorflow()
    
    #DEBUG - EXPONENT COORDS ARE GETTING THROUGH AND MESSING UP THINGS 

    #pop = point_operator()
    #ids = pop.indexer(ids, span, unique, nth)
    #kip.load_geojson('%s/images/out/%s'%(GLOBAL_PROJ,'%s.json'%filename), 0, getfids=[0], getids=[2])
    


    ##!! OLD ARGS - FIX THIS!!
    #kip.load_geojson('%s/images/out/%s'%(GLOBAL_PROJ,'%s.json'%filename), 0, getfids=None, getids=None)



    kip._omit_ids(ids=[32])

    kip.gl_move_center()    
    kip.gl_scale( (.3,.3,.3) )
    
    #kip.gl_move( (1,0,0) )
    kip.gl_rotate( (0,0,90) )

    #kip.make_grid(2, '%s/images/out'%(GLOBAL_PROJ),6, 4)
    #kip.export_grid_gfx('gridgfx', '%s/images/out'%(GLOBAL_PROJ) )

    ##kip.tesl.save_graph_file('%s/test.klsys'%GLOBAL_PROJ)

    #kip.save_line_obj('3d_obj/foo.obj')

    #kip.show_buffers()    
    kip.show_setup()    

    # export_ngc(self, rh, ch, cdpi, cmax, filename, do3d=False):
    #kip.export_extents_ngc(.1, None, .1, 2, '%s/images/out/%s.ngc'%(GLOBAL_PROJ, 'extents'), do3d=True )
    kip.export_ngc(.2, None, .1, 2, '%s/images/out/%s.ngc'%(GLOBAL_PROJ, filename), do3d=True )







    





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



