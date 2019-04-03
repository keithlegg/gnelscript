#!/usr/local/bin/python3



import math 

#from pygfx.point_ops import polygon_operator
#from pygfx.math_ops import vec3 

import os 

from pygfx.obj3d import object3d




##------------------------------------------

## GCODE KNOWN DIALECTS 

bridgeport_gcode =  {
  "G00": "??",
  "G01": "??",
  "G80": "drill_cycle",
  "G170": "mill_cycle",
  "G90": "abs_pos",
  "E": "end",
  "M9": "???FOO",

}


##------------------------------------------

 
class gcode_assembly(object3d):

    def __init__(self):
        super().__init__()  
        self.segments  = [] # [index, command] 
        self.commented = [] # [index, string] 

        self.dialect = bridgeport_gcode
        # self.linear_units = 'in'
        # self.z_axis  ?? 

    def load_gcode_textfile(self, filename):

        #if os.path.lexists(filename) == 0:
        #    self.scribe("%s DOES NOT EXIST !! "%filename )
        #    #raise
            
        if os.path.lexists(filename):
            f = open( filename,"r", encoding='utf-8')
            contents = f.readlines()

            for lin in contents:
                n_idx = lin[0:3]
                comm = lin[3:]
                
                if comm:
                   
                    if lin[0]==';' or lin[0]=='\'':  
                        # print('LINE IS COMMENTED OUT ', lin )
                        self.commented.append(lin)
                    else:    
                        if comm in self.dialect:
                            print("COMMAND FOUND!")
                        
                        print(" COM IS ", comm  )

                        #print("index is %s command is %s " % (n_idx, comm ) )




