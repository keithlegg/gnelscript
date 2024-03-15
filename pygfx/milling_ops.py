"""
PYTHON MILLING TO DO 

COMMANDS 
    5. M62 - M65 Digital Output Control
        M62 P- - turn on digital output synchronized with motion. The P- word specifies the digital output number.
        M63 P- - turn off digital output synchronized with motion. The P- word specifies the digital output number.
        M64 P- - turn on digital output immediately. The P- word specifies the digital output number.
        M65 P- - turn off digital output immediately. The P- word specifies the digital output number.

    The P-word ranges from 0 to a default value of 3. If needed the the number of I/O can be increased by using the num_dio parameter when loading the motion controller. See the Motion Section for more information.
    The M62 & M63 commands will be queued. Subsequent commands referring to the same output number will overwrite the older settings. More than one output change can be specified by issuing more than one M62/M63 command.
    The actual change of the specified outputs will happen at the beginning of the next motion command. If there is no subsequent motion command, the queued output changes won’t happen. It’s best to always program a motion G code (G0, G1, etc) right after the M62/63.
    M64 & M65 happen immediately as they are received by the motion controller. They are not synchronized with movement, and they will break blending.


#GENERAL TODO 
   pendant
   calibrate machine 
   learn change tools mid cuts 
   cut to same location 
   flip over, etc 


#PYTHON TODO 

tool definition class 
    diameter 
    shape
    height 
    flutes 
    speed 

cut circles with tool 

drill with tool 
   
gear generator 

gerbox generator 

flip paths (cut on backside)


SPLINES / ARCS/ ETC 

G2, G3 Arc Move
G2 or G3 axes offsets (center format)
G2 or G3 axes R- (radius format)
G2 or G3 offsets|R- <P-> (full circles)

G5 Cubic Spline
G5 X- Y- <I- J-> P- Q-
    I - X incremental offset from start point to first control point
    J - Y incremental offset from start point to first control point
    P - X incremental offset from end point to second control point
    Q - Y incremental offset from end point to second control point

G5.1 Quadratic Spline
G5.1 X- Y- I- J-
    I - X incremental offset from start point to control point
    J - Y incremental offset from start point to control point

ARCS 
https://www.instructables.com/How-to-program-arcs-and-linear-movement-in-G-Code-/

1) Xs=Xc+(R*cos(Theta1))
2) Ys=Yc+(R*sin(Theta1))
3) Xe=Xc+(R*cos(Theta2))
4) Ye=Yc+(R*sin(Theta2))
5) I=(Xc-(R*cos(Theta1)))-Xc
6) J=(Yc-(R*sin(Theta1)))-Yc 


R=Radius of arc
Theta1=angle of the position of the start point relative to the X axis
Theta2=angle of the position of the end point relative to the X axis
Xc=X coordinate of arc center
Yc=Y coordinate of the arc center
Xs=X coordinate of arc start point
Ys=Y coordinate of arc start point
Xe=X coordinate of arc end point
Ye=Y coordinate of arc end point
I=Incremental X coordinate of start point
J=Incremental Y coordinate of start point



"""


import math 
import re 

import os 



from gnelscript import SHAPELY_IS_LOADED

from gnelscript.pygfx.obj3d import object3d
from gnelscript.pygfx.obj2d import object2d
from gnelscript.pygfx.point_ops import *
from gnelscript.pygfx.math_ops import *
from gnelscript.pygfx.gis_vector_ops import *
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
    from shapely import Polygon as shp_ply
    from shapely import LineString as shp_ln
    #from shapely import Feature as shp_ftr
    #from shapely import FeatureCollection as shp_fc

import trimesh
import numpy as np


"""
A  A axis of machine
B  B axis of machine
C  C axis of machine
D  Tool radius compensation number
F  Feed rate
G  General function (See table Modal Groups)
H  Tool length offset index
I  X offset for arcs and G87 canned cycles
J  Y offset for arcs and G87 canned cycles
K  Z offset for arcs and G87 canned cycles.
   Spindle-Motion Ratio for G33 synchronized movements.
L  generic parameter word for G10, M66 and others
M  Miscellaneous function (See table Modal Groups)
N  Line number
P  Dwell time in canned cycles and with G4.
   Key used with G10.
Q  Feed increment in G73, G83 canned cycles
R  Arc radius or canned cycle plane
S  Spindle speed
T  Tool selection
U  U axis of machine
V  V axis of machine
W  W axis of machine
X  X axis of machine
Y  Y axis of machine
Z  Z axis of machine

"""

"""
    FROM http://linuxcnc.org/docs/html/gcode.html


  G0    Rapid Move
  G1    Linear Move
  G2, G3  I J K or R, P Arc Move
  G4  P Dwell
  G5  I J P Q Cubic Spline
  G5.1  I J Quadratic Spline
  G5.2  P L NURBS
  G38.2 - G38.5   Straight Probe
  G33 K ($) Spindle Synchronized Motion
  G33.1 K ($) Rigid Tapping
  G80     Cancel Canned Cycle


  G81 R L (P) Drilling Cycle
  G82 R L (P) Drilling Cycle, Dwell
  G83 R L Q Drilling Cycle, Peck
  G84 R L (P) ($) Right-hand Tapping Cycle, Dwell
  G73 R L Q Drilling Cycle, Chip Breaking
  G74 R L (P) ($) Left-hand Tapping Cycle, Dwell
  G85 R L (P) Boring Cycle, Feed Out
  G89 R L (P) Boring Cycle, Dwell, Feed Out
  G76 P Z I J R K Q H L E ($) Threading Cycle


  G90, G91    Distance Mode
  G90.1, G91.1    Arc Distance Mode
  G7    Lathe Diameter Mode
  G8    Lathe Radius Mode

  M3, M4, M5  S ($)   Spindle Control
  M19 R Q (P) ($) Orient Spindle
  G96, G97    S D ($) Spindle Control Mode

"""

# COMMENT = ';' # bridgeport uses this 

COMMENT = '('    # linuxcnc uses this 
PARAM   = '#<'   # parameter (variable) , followed with brackets 

#<xscale> = 1.0
#<yscale> = 1.0
#<zscale> = 1.0
#<fscale> = 10000.0
#<toolno> = 1
#<rpm>    = 1600

##------------------------------------------

## GCODE KNOWN DIALECTS 



#  X-Axis (table)  
#  Y-Axis (saddle)  
#  Z-Axis (quill)  


##------------------------------------------



#SHAPELY_IS_LOADED
#https://shapely.readthedocs.io/en/stable/manual.html



# Delaunay triangulation
# line = LineString([(0, 0), (1, 1), (0, 2), (2, 2), (3, 1), (1, 0)])
# dilated = line.buffer(0.5)
# eroded = dilated.buffer(-0.3)



"""
anycubic_fmt = {
    "brand": "Ford",
    "model": "Mustang",
    "year": 1964
}


linuxcnc_fmt = {
    "brand": "Ford",
    "model": "Mustang",
    "year": 1964
}
"""




class gcode_object(object3d):
    ## first stab of gcode translator 
    ## kicad ops has the first working simple gcode exporter
    
    def __init__(self):
        super().__init__()  

        self.DEBUG_MODE = True

        #comments get dumped into this array (full or partial line)
        self.commented = [] # [index, string] 
        
        self.param_names  = [] # [name, value ]
        self.param_values = [] # [name, value ]
        # #swappable dialects of gcode commands 
        # self.dialect = parser_commands

        self.coord_words = ['N','G', 'X','Y','Z','U','V','W','I','J','K','R','P','F'] # ,'A','B','C','D']

        self.segments  = [] # [index, command, xyz_pos] 

    ##-------------------------------##
    def _scrub(self, inp):
        """ clean up parsed characters from (kicad) file """

        out = inp.lower()
        out = out.strip()
        out = out.replace('(','')
        out = out.replace(')','')
        out = out.replace(' ','')
        out = out.replace('\"','')
        out = out.replace('\'','')
        return out

    ##-------------------------------##
    def _clean(self, inp):
        out = inp.strip()
        out = out.replace('(','')
        out = out.replace(')','')
        out = out.replace(' ','')
        out = out.replace('\"','')
        out = out.replace('\'','')
        return out

    ##-------------------------------##
    def _seek_coords(self, tokens):
        out = [0,0,0]

        for token in tokens:
            print(token)

    ##-------------------------------##
    def show_data(self):
        for s in self.segments:
            print(s[1])

    ##-------------------------------##
    def contains_coord_words(self, string): 
        # check for any known coordinate words
        has_coords = False 
        for cw in self.coord_words:
            if cw in string:
                has_coords = True
        return  has_coords       

    ##-------------------------------##
    def which_coord_words(self, string):
        words = []
        for cw2 in self.coord_words:
            if cw2 in string:
                words.append(cw2)
        return words

    ##-------------------------------##
    def parse_params(self, string):
        ## discover and retain any parameters in the file 
        if PARAM in string:
            tmp = self.between_token_list(string, [PARAM, '>']) 

            if len(tmp)==2:
                if tmp[0] not in self.param_names:
                    if '=' in tmp[1]:
                        tmp2 = tmp[1].split('=')
                        self.param_names.append(tmp[0])
                        self.param_values.append(tmp2[1].replace(' ',''))

    ##-------------------------------##
    def between_token_list(self, string, tokens, return_match=False):
        ## give a list of tokens, return a list of betweens 
        ## string = 'y987a123b541c307d999'
        ## tokens = ['a','b','c','d']
        ## out = self.between_token_list( string, tokens)
        ## out will be:
        ##      ['123', '541', '307']
        ## unless return_match is True,
        ##     then out will be:
        ##         [ ['a',123'], ['b',541'], ['c',307'] ]

        

        #return re.split( tokens_str ,string)
        match = ''
        if len(tokens)==1:
            match = tokens[0] + '|'
        if len(tokens)>1:    
            for i,t in enumerate(tokens):
                if i<len(tokens)-1:
                    match = match + t + '|'
                if i==len(tokens)-1:
                    match = match + t 
        
        if return_match is False:
            tmp = re.split( match, string)
        else:
            match = '('+match+')'
            tmp = re.split( match, string)

        out = []
        for t in tmp:
            tmp2 = t.replace('\n','')
            tmp2 = tmp2.replace(';','')
            out.append(tmp2)

        if return_match is False:
            return out[1:]  
        else:
            out2 = [] 
            for i,tok in enumerate(out[1:]):
                if i%2==0:
                    out2.append( [tok, out[1:][i+1]] )
            return out2


    ##-------------------------------##
    def save_gcode(self, filename):
        f = open( filename,"w", encoding='utf-8')
        
        for seg in self.segments:
            #f.write(seg[])
            print( seg[0], seg[2])


    ##-------------------------------##
    def calc_arc(self, R, Theta1, Theta2, Xc, Yc):
        """
            appears to work but not fully tested 


            R       - Radius of arc                            
            Theta1  - angle of the position of the start point relative to the X axis 
            Theta2  - angle of the position of the end point relative to the X axis 
            Xc      - X coordinate of arc center 
            Yc      - Y coordinate of the arc center 

            Xs      - X coordinate of arc start point
            Ys      - Y coordinate of arc start point
            Xe      - X coordinate of arc end point
            Ye      - Y coordinate of arc end point
            I       - Incremental X coordinate of start point
            J       - Incremental Y coordinate of start point


            USAGE:
                gcc = cnc_op()
                gcc.calc_arc( .5, 90, 0, 2.25, 3) 

        """

        pts = [] 

        #R=float(R)
        #Theta1=float(Theta1)
        #Theta2=float(Theta2)
        #Xc=float(Xc)
        #Yc=float(Yc)

        Xs = Xc +(R* math.cos(Theta1))       
        Ys = Yc +(R* math.sin(Theta1))       
        Xe = Xc +(R* math.cos(Theta2))       
        Ye = Yc +(R* math.sin(Theta2))       
        I  = (Xc-(R* math.cos(Theta1)))-Xc   
        J  = (Yc-(R* math.sin(Theta1)))-Yc   
        
        print(Xs)
        print(Ys)
        print(Xe)
        print(Ye)
        print(I)
        print(J)


        return pts


    ##-------------------------------##
    def load_gcode(self, filename):

        if os.path.lexists(filename) == 0:
            self._scribe("%s DOES NOT EXIST !! "%filename )
            #raise
            
        if os.path.lexists(filename):
            f = open( filename,"r", encoding='utf-8')
            contents = f.readlines()
          
            # #scan entire file for parametrs first 
            # for lin in contents:
            #     self.parse_params(lin)
            # print("#### param_names  ", self.param_names )
            # print("#### param_values ", self.param_values )

            cleaned_contents = []

            lastx = 0  
            lasty = 0  
            lastz = 0  
            tmp_poly = []

            for lin in contents:

                #if self.DEBUG_MODE:
                #print("#### LINE IS ", lin.replace("\n","") )

                # ignore commented out lines 
                if lin[0]==COMMENT or lin[0]=='\'':  
                    self.commented.append(lin)

                else:
                    # seperate the line index out and split from rest of command
                    lindex = self.between_token_list(lin, self.coord_words, return_match=True)
                    if lindex:
                        #print( lindex)

                        # execute the parameters(variables) if they exist 
                        for i,tok in enumerate(lindex):
                            xval=None;yval=None;zval=None

                            if self._clean(tok[0])=='G':
                                if self._clean(tok[1])=='0':
                                    if len(tmp_poly):
                                        self.segments.append(tmp_poly)
                                        tmp_poly = []

                                if self._clean(tok[1])=='1':
                                        self.segments.append(tmp_poly)
                                        tmp_poly = []

                                if self._clean(tok[1])=='2':
                                        #G02- Clockwise circular interpolation

                                        print('###### G2 FOUND ')
                                        print(lin)

                                if self._clean(tok[1])=='3':
                                        #G03- Counter Clockwise circular interpolation

                                        print('###### G3 FOUND ')
                                        print(lin)

                            if self._clean(tok[0])=='N':
                                if 'G43' not in lin: 
                                    for t in lindex:
                                        if self._clean(t[0])=='X':
                                            xval = float(t[1])
                                            lastx = xval  
                                        if self._clean(t[0])=='Y':
                                            yval = float(t[1])
                                            lasty = yval                             
                                        if self._clean(t[0])=='Z':
                                            zval = float(t[1])
                                            lastz = zval 


                                    tmp_poly.append( (xval if xval else lastx, 
                                                      yval if yval else lasty,
                                                      zval if zval else lastz) )

 
                            


                            #for j,p in enumerate(self.param_names):


                        #print( lindex )


                        ## #check for comments on the line but not at start   
                        ## if COMMENT in lin:
                        ##     tmp = lin.split(COMMENT)
                        ##     lin = tmp[0]
                        ##     if tmp[1] is not '\n':
                        ##         self.commented.append(tmp[1])
                        ## # check for known commands 
                        ## for key in self.dialect :
                        ##     if key in comm:
                        ##         print(" COM FOUND at INDEX %s ! "%n_idx , key, '---', self.dialect [key] )
                        ## # check for any known coordinate words
                        ## has_coords = self.contains_coord_words(comm) 
                        ## # if has coordinate words, determine which ones it has 
                        ## if has_coords is True:             
                        ##     line_contains = ( self.which_coord_words(comm) )
                        ##     output        = self.between_token_list(comm, line_contains)
                        ##     #print( line_contains, output )
                        ##     for i,token in enumerate(line_contains):
                        ##         if token=='X':
                        ##             self.POSX = float(output[i])
                        ##         if token=='Y':
                        ##             self.POSY = float(output[i])
                        ##         if token=='Z':
                        ##             self.POSZ = float(output[i])

                        #self.segments.append( [n_idx, [self.POSX, self.POSY, self.POSZ], comm ] )



##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##

class cnc_op( vectorflow ):
    """ 
     gcode parsing was complex enough to get its own class (gcode_object)  
     vectorflow is the experimental tool to seamlessly integreate OBJ, JSON, and NGC 

     This is the main place for all general CNC stuff 
     
     - hole drilling 
     - dialate/erode 
     - etc 

     the complex cutting edge stuff will go in cam_op 

    """


    def __init__(self):
        super().__init__()  

        # simulated position of the cutting head
        self.HEAD_POSX = 0
        self.HEAD_POSY = 0
        self.HEAD_POSZ = 0

    ##-------------------------------## 
    ##-------------------------------## 
    if SHAPELY_IS_LOADED:
        def test_shapely(self):
            line = shp_ln([(2, 0), (2, 4), (3, 4)])
            print(line)

        ##-------------------------------## 
        def cvt_grsort_shapely(self):
            #               [[id, centroid, extents, len, points ]]  
            
            geom = [] 
            for poly in self.gr_sort:
                geom.append( shp_ply(self.cvt_3d_to_2d(poly[4])) )
            
            return geom

        #def shply_buffer(self, geom_obj):
        #        #buff = buffer(Point(10, 10), 2, quad_segs=2)
        #        buff = buffer(poly, 2, quad_segs=2)
        #        
        #        return buff

        ##-------------------------------## 
        def cvt_shapely_grsort(self, shapely_geom, zval=0):
            xx, yy = shapely_geom.exterior.coords.xy

            x = xx.tolist()
            y = yy.tolist()

            ply = [] 

            for i in range(0,len(x)):
                #2d or 3D ??
                #ply.append( (x, y) )
                ply.append( (x[i], y[i], zval) )

            self.insert_gr_sort(ply)  


    ##-------------------------------## 
    ##-------------------------------## 
    def _set_cam_properties_mop(self, rh, ch, cdpi, cmax):
        self.rh = rh           
        self.ch = ch            
        self.cdpi = cdpi       
        self.cmax = cmax          
    
    ##-------------------------------## 
    def _set_extents_mop(self, bbox):
        """ set global extents for generating data 
            based on PIL coordinate which is [left, top, right, bottom] 
            [minx, miny, maxx, maxy]  
        """
        
        self.sort_minx = bbox[0]
        self.sort_miny = bbox[1]         
        self.sort_maxx = bbox[2]
        self.sort_maxy = bbox[3]

    ##-------------------------------##       
    def _omit_ids_mop(self, ids=None, span=None, unique=True, nth=None):
        
        #DEBUG - not fully working - see indexer docs - nths with crash if you try it 
        pop = pop3d()
        ids = pop.indexer(ids, span, unique, nth)
        
        print(" ### omit ids: ", ids)

        self.omit_ids = ids

    ##-------------------------------## 
    def _make_outline(self, iterations):
        """ 
            DEBUG - NEED TO ADD 1/2 DIA OF CUTTING TOOL 
            WE NEED DIALATE AND ERODE TO MAKE THIS WORK 
            first attempt at a Z operation - spiral down a polygon 

        """
        
        #spirals = [] 

        for idx, sort in enumerate(self.gr_sort):

            newply = []            
            ply = sort[4]

            for i, c in enumerate(range(iterations)):
                depth = self.ch-(i*self.cdpi)
                for pt in ply:
                    if depth < self.cmax:
                        newply.append( (pt[0] ,pt[1], depth) )
                    if depth>self.cmax:
                        print("cut too deep ")

            #spirals.append(newply)
            self.gr_sort[idx][4] = newply  

    ##-------------------------------## 
    def show_setup_mop(self):
        print('##--------------------------##')        
        print("retract height %s"%self.rh)
        print("cut height     %s"%self.ch)
        print("cut max        %s"%self.cmax)
        print("cut cdpi       %s"%self.cdpi)

        print('original extents  %s %s %s %s '%(self.orig_minx, self.orig_miny, self.orig_maxx, self.orig_maxy))
        print('sorted extents    %s %s %s %s '%(self.sort_minx, self.sort_miny, self.sort_maxx, self.sort_maxy)) 

    ##-------------------------------##
    def export_3dprint(self, filename):
        """   
            DEBUG - NOT DONE  
        """

        print("# exporting 3d printer NGC file ", filename)

        self._calculate_3dprint()

        print("gr_sort buffer has %s polys in it. "%(len(self.gr_sort)))
        fobj = open( filename,"w") #encoding='utf-8'
        for line in self.outfile: 
            fobj.write(line+'\n')
        fobj.close()

    ##-------------------------------##
    def _calculate_3dprint(self, do_retracts=True, doround=True):
        """ 
             DEBUGGY 
        """
        splat = 5 
        stax = 5
        zincr = .5

        pl = 2 #numeric rounding places 
        lastpt = (0,0,0)

        self.outfile.append('( exported with _calculate_3dprint )')
        self.outfile.append('(linear scale set to %s of internal coordinates)'%self.global_scale )
        self.outfile.append('  ')

        self.outfile.append('M73 P0 R19')
        self.outfile.append('M201 X2500 Y2500 Z800 E2500 ; sets maximum accelerations, mm/sec^2')
        self.outfile.append('M203 X300 Y250 Z8 E80 ; sets maximum feedrates, mm / sec')
        self.outfile.append('M204 S2500 T2500 ; sets acceleration (S) and retract acceleration (R), mm/sec^2')
        self.outfile.append('M205 X15.00 Y10.00 Z2.00 E10.00 ; sets the jerk limits, mm/sec')
        self.outfile.append('M205 S0 T0 ; sets the minimum extruding and travel feed rate, mm/sec')
        self.outfile.append(';TYPE:Custom')
        self.outfile.append('  ')

        self.outfile.append('G90 ; use absolute coordinates')
        self.outfile.append('M83 ; extruder relative mode')
        self.outfile.append('M104 S215 ; set extruder temp')
        self.outfile.append('M140 S40 ; set bed temp')
        self.outfile.append('M190 S40 ; wait for bed temp')
        self.outfile.append('M109 S215 ; wait for extruder temp')
        self.outfile.append('  ')

        self.outfile.append('G28                      ; move X/Y/Z to min endstops')
        self.outfile.append('G1 Z0.28                 ; lift nozzle a bit ')
        self.outfile.append('G92 E0 ')
        self.outfile.append('  ')

        self.outfile.append('G1 Y3 F1800              ; zero the extruded length ')
        self.outfile.append('G1 X60  E25 F500         ; Extrude 25mm of filament in a 5cm line. ')
        self.outfile.append('G92 E0                   ; zero the extruded length again ')
        self.outfile.append('G1 E-2 F500              ; Retract a little ')
        self.outfile.append('G1 X70 F4000             ; Quickly wipe away from the filament line')
        self.outfile.append('M117                  ')
        self.outfile.append('  ')

        self.outfile.append('G21     ; set units to millimeters')
        self.outfile.append('G90     ; use absolute coordinates')
        #self.outfile.append('M82     ; use absolute distances for extrusion')
        self.outfile.append('M83     ; use relative distances for extrusion')

        self.outfile.append('G92 E0  ; Filament gcode')
        self.outfile.append('M107    ;fan off / LAYER_CHANGE')

        #self.outfile.append(';Z:0.28')
        #self.outfile.append(';HEIGHT:0.28')
        #self.outfile.append('; BEFORE_LAYER_CHANGE 0 @ 0.28mm')
        self.outfile.append('G1 E-2 F4800')
        self.outfile.append('G92 E0')
        self.outfile.append('G1 Z.28 F7200')
        
        #self.outfile.append('; AFTER_LAYER_CHANGE 0 @ 0.28mm')
        self.outfile.append('G1 X70.004 Y91.359')
        self.outfile.append('G1 E2 F2000') #extrusion feed rate 
        self.outfile.append('M204 S2000')
        self.outfile.append('  ')

        ##-----------------------------------------##

        #move to origin  
        
        #self.outfile.append('g0 x%s y%s z%s f30'% ( self.hp[0], self.hp[1], self.rh) )   #rapid move to 0 
        #self.ngc_to_obj.append( ( self.hp[0], self.hp[1], self.rh) ) 

   
        ##-----------------------------------------##
        self.outfile.append('  ')
        self.outfile.append('(exporting filled polygons )')

        ##-----------------------------------------##
        # build the gcode up with simple linear movements 

        self.outfile.append('  ')
        self.outfile.append('(exporting graphic polygons )')
        
        ##------------------------------


        #self.outfile.append( 'G1 E0') #extruder off                
        #self.outfile.append( 'G1 E-2 F500 ') #retract filament little 
        #self.outfile.append( 'G1 Z%s'%(self.rh ) )  #retract head to set height 

        # graphical polygons - build the gcode up with simple linear movements
        for si in range(stax):
            for rc,row in enumerate(self.gr_sort):

                # preprocess 
                export_ply = True 

                if row[0] in self.omit_ids:
                    print("# omitting polygon ID %s"%row[0])
                    export_ply = False

                #modified to export sorted data - easy peasy  
                gr_poly = row[4]

                ##grab a polygon off the stack and print it 
                if len(gr_poly) and export_ply:
                    self.outfile.append( ' ' )                
                    self.outfile.append('(exporting new polygon %s)'%rc)

                    # no formatting (full precision)
                    if doround:                
                        pt1=(round(gr_poly[0][0],pl) ,round(gr_poly[0][1],pl)) 
                    else:
                        pt1 = gr_poly[0]
           

                    ## iterate points in each polygon 
                    for i,pt in enumerate(gr_poly):
                        if doround:
                            tmp = ( round(pt[0],pl), round(pt[1],pl), round(pt[2],pl) ) 
                            pt = tmp  
                        self.outfile.append( 'G1 Z%s'%(pt[2]+(si*zincr) ) ) 
                        self.outfile.append( 'G1 X%s Y%s E%s'%( pt[0], pt[1], splat ) )  

 
        self.outfile.append( ' ' )
        self.outfile.append( ' ' )

        self.outfile.append(';WIPE_END')
        self.outfile.append('G92 E0')
        self.outfile.append('M107')
        self.outfile.append(';TYPE:Custom')
        self.outfile.append('; Filament-specific end gcode ')
        self.outfile.append(';END gcode for filament')
        self.outfile.append('M104 S0                                    ; Extruder off ')
        self.outfile.append('M140 S0                                    ; Heatbed off ')
        self.outfile.append('M107                                       ; Fan off ')
        self.outfile.append('G91                                        ; relative positioning ')
        self.outfile.append('G1 E-5 F3000  ')
        self.outfile.append('G1 Z+0.3 F3000                             ; lift print head ')
        self.outfile.append('G28 X0  F3000')
        self.outfile.append('M84                                        ; disable stepper motors')
        self.outfile.append('M73 P100 R0')

        self.outfile.append('  ')

        ##-----------------------------------------##
        # rapid move at end 
        self.outfile.append('m2') #program end

    ##-------------------------------##
    """
    def _calculate_paths3d(self, do_retracts=True, doround=True):
        pl = 6 #numeric rounding places 
        lastpt = (0,0,0)

        self.outfile.append('( exported with _calculate_paths3d )')
        self.outfile.append('(linear scale set to %s of internal coordinates)'%self.global_scale )
        self.outfile.append('  ')

        self.outfile.append('g20')                  #inches for unit 
        
        ##-----------------------------------------##

        #move to origin  
        self.outfile.append('g0 x%s y%s z%s f30'% ( self.hp[0], self.hp[1], self.rh) )   #rapid move to 0 
        self.ngc_to_obj.append( ( self.hp[0], self.hp[1], self.rh) ) 

        if do_retracts:
            self.outfile.append('g0z%s'% ( self.rh ) )  #retract in between cuts
   
        ##-----------------------------------------##
        self.outfile.append('  ')
        self.outfile.append('(exporting filled polygons )')

        ##-----------------------------------------##
        # build the gcode up with simple linear movements 

        self.outfile.append('  ')
        self.outfile.append('(exporting graphic polygons )')
        
        ##------------------------------

        # graphical polygons - build the gcode up with simple linear movements
        for row in self.gr_sort:
            
            # preprocess 
            export_ply = True 

            if row[0] in self.omit_ids:
                print("# omitting polygon ID %s"%row[0])
                export_ply = False

            #modified to export sorted data - easy peasy  
            gr_poly = row[4]

            ##--

            if len(gr_poly) and export_ply:
                self.outfile.append('(exporting new polygon )')
                
                ##-- 
                #DEBUG - need to sort out clean points - want to run as close to final export 
                #it looses precision 

                # no formatting (full precision)
                if doround:                
                    pt1=(round(gr_poly[0][0],pl) ,round(gr_poly[0][1],pl)) 
                else:
                    pt1 = gr_poly[0]

                ## first point at retract height 
                if do_retracts:
                    #move to first point RH 
                    self.outfile.append('x%s y%s z%s'% (  pt1[0] , pt1[1], self.rh ) )               
                    self.ngc_to_obj.append( ( pt1[0], pt1[1], self.rh ) )             

                ## iterate points in polygon 
                self.outfile.append( 'G1' )
                for i,pt in enumerate(gr_poly):
                    if doround:
                        tmp = ( round(pt[0],pl), round(pt[1],pl), round(pt[2],pl) ) 
                        pt = tmp  

                    self.outfile.append( 'x%s y%s z%s'%( pt[0], pt[1],  pt[2] ) )
                    self.ngc_to_obj.append( (pt[0], pt[1],  pt[2]) )                   
                
                self.outfile.append( 'G0' )

                if do_retracts:
                    if doround:
                        gpt=(round(gr_poly[0][0],pl) ,round(gr_poly[0][1],pl)) 
                    else:
                        gpt=gr_poly[0]

                    self.ngc_to_obj.append( (gpt[0], gpt[1], self.rh)  )
                    self.outfile.append( 'x%s y%s z%s'%( gpt[0], gpt[1], self.rh) )

                    #### retract in between cuts
                    self.outfile.append('g0z%s'% ( self.rh ) )  

                self.outfile.append('  ')

        ##-----------------------------------------##
        # rapid move at end 
        self.outfile.append('m2') #program end
    """

    ##-------------------------------##
    def export_centroids_mop(self, folder, name):
        """ export a centroid for each polygon as a geojson file """

        features = []

        #[[id, centroid, extents, len, points ]]
        for i,s in enumerate(self.gr_sort):
            features.append(Feature(geometry=Point((s[1][0],s[1][1])), properties={"id":i, "len" : len(s[4]) } ))

        feature_collection = FeatureCollection(features)
        with open('%s/%s_ply_cntrs.json'%(folder,name), 'w') as f:
            dump(feature_collection, f)

    ##-------------------------------##
    def export_extents_mop(self, name, folder):
        """ export the data extents as a geojson polyline """

        features = []

        #[[id, centroid, extents, len, points ]]
        for i,s in enumerate(self.gr_sort):
            #    features.append(Feature(geometry=gjln(coordinates=self.extents_fr_bbox(s[2])), properties={"id":i, "len" : len(s[4]) } ))
            coords = self.extents_fr_bbox(s[2], periodic=True)

            features.append(Feature(geometry=gjln(coordinates=coords), 
                                 properties={"id" : 0 
                                            }
                                 ) 
                         )

        #sort extents  
        coords = self.extents_fr_bbox([self.sort_minx,self.sort_miny,self.sort_maxx,self.sort_maxy], periodic=True)
        features.append(Feature(geometry=gjln(coordinates=coords), 
                                properties={"id" : 0 
                                           }
                                ) 
                        )
                        

        feature_collection = FeatureCollection(features)
        with open('%s/%s_ply_xtntx.json'%(folder,name), 'w') as f:
            dump(feature_collection, f)

    ##-------------------------------##
    def save_line_obj_mop(self, name):
        """ dump a bunch of 3d points as a giant line to visualize in gnolmec 
            
            do_retracts iterates the gr_polys buffer (not working yet)
            otherwise it iterates the ngc_to_obj buffer 

        """   
        
        #print("exporting %s polygons "%len(self.gr_polys))

        # make a series of connected lines from points 
        self.points = self.clean_pts_str(self.ngc_to_obj)

        for i,pt in enumerate(self.ngc_to_obj):
            if i>0 and i<len(self.ngc_to_obj):
                poly = (i,i+1)
                self.polygons.append( poly ) 

        #self.scale_pts(self.scale)
        self.save(name)

    ##-------------------------------##
    def export_extents_ngc_mop(self, rh, ch, cdpi, cmax, filename, do3d=False):
        
        #self.gl_extents()

        tempbuffer = self.gr_sort
        self.gr_sort = []

        #self.prim_triangle('z', (0,0,0), (0,0,0) )
        # [[id, centroid, extents, len, points ]] 
        #self.gr_sort.append([0,0,0,0,self.points]) 
        
        pts = self.calc_square_diag((self.sort_minx,self.sort_miny ),
                                   (self.sort_maxx,self.sort_maxy), add_zaxis=True ) 
        pts.append( (self.sort_minx, self.sort_miny, 0) )
        self.gr_sort.append([0,0,0,0,pts])

        if do3d==True:
            self._calculate_paths3d()
        else:
            self._calculate_paths2d()
        
        print("gr_sort buffer has %s polys in it. "%(len(self.gr_sort)))
        fobj = open( filename,"w") #encoding='utf-8'
        for line in self.outfile: 
            fobj.write(line+'\n')
        fobj.close()

        self.gr_sort = tempbuffer

    ##-------------------------------##
    def export_ngc_mop(self, rh, ch, cdpi, cmax, filename, do3d=False, do_retracts=False):

        if ch==None:
            print("# exporting NGC file - cut height disabled (3d) ", filename)
        else:
            print("# exporting NGC file ", filename)
            self.ch = ch          # cut height (top, start of cut)
        
        self.rh = rh          # retract height 
        self.cdpi = cdpi      # cut depth per iteration on Z axis
        self.cmax = cmax      # maximum cut depth on Z axis 

        if do3d==True:
            self._calculate_paths3d(do_retracts=do_retracts, doround=True)
        else:
            self._calculate_paths2d(do_retracts=do_retracts, doround=True)
        
        print("gr_sort buffer has %s polys in it. "%(len(self.gr_sort)))
        fobj = open( filename,"w") #encoding='utf-8'
        for line in self.outfile: 
            fobj.write(line+'\n')
        fobj.close()

    ##-------------------------------##
    def make_grid(self, distance, folder, xcuts, ycuts, bbox=None):
        """ chop a square into smaller squares 

        """

        if bbox:
            self.tesl._set_extents(bbox) 
        else:
            self.tesl._set_extents([self.sort_minx, self.sort_miny, self.sort_maxx, self.sort_maxy]) 

        self.tesl.build_2d_cells(xcuts, ycuts)
        

        #temporary export of geom I HAVE TO SEE THIS BEFORE I GO TO BED!
        features = []

        ##---
        # add attrs to derived DAG graph nodes 
        for c in self.tesl.nodes:
            cen = (c.coord_x+(c.width/2), c.coord_y+(c.height/2))
            c.addattr('centroid', cen )
            c.addattr('width' , c.width )
            c.addattr('height', c.height )

            if len(self.gr_sort):
                for sort in self.gr_sort:
                    # [[id, centroid, extents, len, points ]]
                    dist = self.mu.calc_line_length(cen[0], cen[1], sort[1][0], sort[1][1])
                    if float(dist) < float(distance):
                        
                        #DEBUG TODO - there is a problem of the same polygon getting added more than once 
                        #make sure we dont export multiple times 
                        # (check the number of points, then check thr actual points) 

                        tmp_pts=[cen, sort[1]]
                        features.append(Feature(geometry=gjln(coordinates=tmp_pts)) 
                                       )

        feature_collection = FeatureCollection(features)
        with open('%s/_distances.json'%(folder), 'w') as f:
            dump(feature_collection, f)


##-----------------------------------------##

class cam_op(cnc_op):
    """ 
        This is analagous to a fusion360 (chained) CAM OP 

        True 3D milling operations  
        the highest level object for CAM and milling 


        TODO - MULTIHEAD DISTRIBUTED OPERATIONS !!!
           define a rgion of work (presumably just half) 
           jobs get run in parallel !!




    """

    def __init__(self):
        super().__init__()  
        self.tesl = tessellator() 
        self.gcparser = gcode_object()

    def obj_to_wkt(self):
        print(self.points)

    ##-------------------------------## 
    #def hemisphere_slice(self, origin, normal, numdivs, path, infile, axis='y'):
    #    pass 
        

    ##-------------------------------## 
    def tm_section_test(self, origin, normal, numdivs, path, infile, axis='y'):

        if type(origin)==vec3:
            origin = origin.aspt
                    
        if type(normal)==vec3:
            normal = normal.aspt

        mesh = trimesh.load_mesh('%s/%s'%(path,infile))

        slice2d = mesh.section(
            plane_normal=normal,
            plane_origin=origin,

        )
        
        if slice2d==None:
            print('### tm_section_test NO POLYGOINS FOUND \n\n')
            return None 

        geom, m44 = slice2d.to_planar()
        pts = [poly for poly in geom.polygons_full]



        """
        # if we wanted to take a bunch of parallel slices, like for a 3D printer
        # we can do that easily with the section_multiplane method
        # we're going to slice the mesh into evenly spaced chunks along z
        # this takes the (2,3) bounding box and slices it into [minz, maxz]
        z_extents = mesh.bounds[:,2]
        # slice every .125 model units (eg, inches)
        z_levels  = np.arange(*z_extents, step=.125)

        # find a bunch of parallel cross sections
        sections = mesh.section_multiplane(plane_origin=mesh.bounds[0], 
                                           plane_normal=[0,0,1], 
                                           heights=z_levels)
        sections

        # summing the array of Path2D objects will put all of the curves
        # into one Path2D object, which we can plot easily
        combined = np.sum(sections)
        combined.show()

        # if we want to intersect a line with this 2D polygon, we can use shapely methods
        polygon = slice_2D.polygons_full[0]
        # intersect line with one of the polygons
        hits = polygon.intersection(LineString([[-4,-1], [3,0]]))
        # check what class the intersection returned
        hits.__class__

        # we can plot the intersection (red) and our original geometry(black and green)
        ax = plt.gca()
        for h in hits.geoms:
            ax.plot(*h.xy, color='r')
        slice_2D.show()

        # the medial axis is available for closed Path2D objects
        (slice_2D + slice_2D.medial_axis()).show()
        """

        
        #poly_union = shapely.geometry.MultiPolygon([poly for poly in slice_2D.polygons_full])

        return pts,m44

    ##-------------------------------##     
    def kdag_merge_points(self, pts):
        """ experiment with DG graphs in an attempt to merge points 
        """
        
        import networkx as nx
        G = nx.Graph()

    ##-------------------------------##     
    def nx_merge_points(self, pts):
        """ experiment with DG graphs in an attempt to merge points 
        """
        
        import networkx as nx
        G = nx.Graph()

    ##-------------------------------## 
    def boundary(mesh, close_paths=True):
        #https://stackoverflow.com/questions/76435070/how-do-i-use-python-trimesh-to-get-boundary-vertex-indices
        # Set of all edges and of boundary edges (those that appear only once).
        edge_set = set()
        boundary_edges = set()

        # Iterate over all edges, as tuples in the form (i, j) (sorted with i < j to remove ambiguities).
        # For each edge, three cases are possible:
        # 1. The edge has never been visited before. In this case, we can add it to the edge set and as a boundary
        #    candidate as well.
        # 2. The edge had already been visited once. We want to keep it into the set of all edges but remove it from the
        #    boundary set.
        # 3. The edge had already been visited at least twice. This is generally an indication that there is an issue with
        #    the mesh. More precisely, it is not a manifold, and boundaries are not closed-loops.
        for e in map(tuple, mesh.edges_sorted):
            if e not in edge_set:
                edge_set.add(e)
                boundary_edges.add(e)
            elif e in boundary_edges:
                boundary_edges.remove(e)
            else:
                raise RuntimeError(f"The mesh is not a manifold: edge {e} appears more than twice.")

        # Given all boundary vertices, we create a simple dictionary that tells who are their neighbours.
        neighbours = defaultdict(lambda: [])
        for v1, v2 in boundary_edges:
            neighbours[v1].append(v2)
            neighbours[v2].append(v1)

        # We now look for all boundary paths by "extracting" one loop at a time. After obtaining a path, we remove its edges
        # from the "boundary_edges" set. The algorithm terminates when all edges have been used.
        boundary_paths = []

        while len(boundary_edges) > 0:
            # Given the set of remaining boundary edges, get one of them and use it to start the current boundary path.
            # In the sequel, v_previous and v_current represent the edge that we are currently processing.
            v_previous, v_current = next(iter(boundary_edges))
            boundary_vertices = [v_previous]

            # Keep iterating until we close the current boundary curve (the "next" vertex is the same as the first one).
            while v_current != boundary_vertices[0]:
                # We grow the path by adding the vertex "v_current".
                boundary_vertices.append(v_current)

                # We now check which is the next vertex to visit.
                v1, v2 = neighbours[v_current]
                if v1 != v_previous:
                    v_current, v_previous = v1, v_current
                elif v2 != v_previous:
                    v_current, v_previous = v2, v_current
                else:
                    # This line should be un-reachable. I am keeping it only to detect bugs in case I made a mistake when
                    # designing the algorithm.
                    raise RuntimeError(f"Next vertices to visit ({v1=}, {v2=}) are both equal to {v_previous=}.")

            # Close the path (by repeating the first vertex) if needed.
            if close_paths:
                boundary_vertices.append(boundary_vertices[0])

            # "Convert" the vertices from indices to actual Cartesian coordinates.
            boundary_paths.append(mesh.vertices[boundary_vertices])

            # Remove all boundary edges that were added to the last path.
            boundary_edges = set(e for e in boundary_edges if e[0] not in boundary_vertices)

        # Return the list of boundary paths.
        return boundary_paths

    ##-------------------------------## 
    def tm_meshplane_test(self, origin, normal, numdivs, path, infile, axis='y'):

        out = [] 

        if type(origin)==vec3:
            origin = origin.aspt
        if type(normal)==vec3:
            normal = normal.aspt

        mesh = trimesh.load_mesh('%s/%s'%(path,infile))
        #mesh = trimesh.creation.icosphere()
        
        bbox = mesh.bounds.tolist()

        ############
        out = trimesh.intersections.mesh_plane(mesh, normal, origin, return_faces=False, local_faces=None, cached_dots=None)

        print(type(out))

        """
        pts = []  
        tmp = out.tolist()
        for pt in tmp: 
            pts.append( pt )
        
        #self.contiguous_segs_to_poly(pts)
        """ 

        return out 

    ##-------------------------------## 
    """
    # FROM https://lukeparry.uk/tag/slicing/

    # Generate a list of slices for each layer across the entire mesh
    zSlices = np.linspace(-5, 5, 200)
    k = 1000 # number of slices
    zMin, zMax = myMesh.bounds[:,2]
    zBox = np.linspace(zMin, zMax, k)


    tris = myMesh.triangles

    # Obtain the min and maximum Z values across the entire mesh
    zVals = tris[:, :, 2]
    triMin = np.min(zVals, axis=1)
    triMax = np.max(zVals, axis=1)


    # Below is a manual approach to sorting and collecting the triangles across each layer 

    if False:

        triSortIdx = np.argsort(triMinMax[1,:])
        triMinMaxSort = triMinMax[:,triSortIdx]

        startTime = time.time()

        sortTris = []

        iSects2 = []
        for i in range(len(zBox)):
            minInside = zBox[i].reshape(-1,1) > triMinMax[0, :].reshape(1, -1)
            maxInside = zBox[i].reshape(-1,1) < triMinMax[1, :].reshape(1, -1)
            iSects2.append(np.argwhere((minInside & maxInside).ravel()))

        print('endTime array base', time.time() - startTime)

    # An alteratnvie more compact way is to use binary search operator available within
    # numpy.searchsorted. This locates the bottom and top layer position 

    minIdx = np.searchsorted(zBox, triMin, side='left')
    maxIdx = np.searchsorted(zBox, triMax, side='left')

    # Attach the corresponding presorted triangles into the 
    iSects = [[] for i in range(len(zBox))]

    # The iterative part for assigning potential triangles for intersection on a triangle are performed here
    # Note: Process is very inefficient in native Python O(n*k) 
    for i in range(len(minIdx)):

        startLayer = minIdx[i]
        endLayer = maxIdx[i]
        for layer in iSects[startLayer:endLayer+1]:
            layer.append(i)


    #########################
    The pre-processing has been completed, now the final slincg may be complete. 
    Performing some micro-optimisations, further refactoring may be done to adapt the code previously present in Trimesh.

    The process of slicing or intersection, is simple. 
    Suprisingly, there are a no obvious references for this process. 
    Slicing a triangular mesh, relies firstly computing the potential intersection of each edge of a faceted mesh, 
    in this case it is generalised for an abritrary plane. Firstly, the vector between the mesh vertices 
     and the plane origin is calculated – red lines in the diagram. 

     The dot product is taken with the slicing plane normal nn. 
     The sign of the dot product indicates if the point lies above or below the plane – zero uniquely is the intersection.

    plane_origin = np.array([0,0,0])
    plane_normal = np.array([0,0,1])

    vertex_dots = np.dot(myMesh.vertices - plane_origin, plane_normal)

    """

    ##-------------------------------## 
    # def cross_section(mesh, plane_origin=[0,0,0], plane_normal=[1,0,0]):
    #     slice_ = mesh.section(plane_origin=plane_origin, 
    #                           plane_normal=plane_normal)
    #     # transformation matrix for to_planar 
    #      to_2D = trimesh.geometry.align_vectors(plane_normal, [0,0,-1])
    #     
    #     slice_2D, to_3D = slice_.to_planar(to_2D = to_2D)
    #     return slice_2D, to_3D

    ##---------------------------------------------
    
 
    def ltp(self, linestrings):
        """
        ltp = linestrings_to_polygon

        Return valid polygon for unordered set of 3D linestrings in numpy array (n,2,3)
        
        """
        n_vertices, dim = linestrings.shape[0] * linestrings.shape[1], linestrings.shape[2]
        vertices = linestrings.reshape(n_vertices, dim)
        uniq, index, inv = np.unique(vertices.round(decimals=4), return_index=True, return_inverse=True, axis=0)
        assert uniq.shape[0] == vertices.shape[0] / 2
        linestrings_indices = inv.reshape(inv.shape[0]//2, 2)
        
        G = nx.from_edgelist([(a,b) for a,b in linestrings_indices])
        polygon = shapely.geometry.shp_ply(vertices[index][[a for a, b in nx.find_cycle(G)]])
        assert polygon.is_valid
        return polygon
    
  

    # if you do trimesh.load_path(sections) It will give you a Path3D object 
    # (or if you want to call the specific function, it's trimesh.path.io.misc.lines_to_path). 
    # If you transform the sections onto the plane either manually or through path.to_planar, 
    # It will then have path.polygons_full constructed for you.
    # As for the result of triangulate_polygon, that looks correct- the vertices are 2D because 
    # it's coming from a planar polygon, you just need to np.column_stack it with zeros to use it in a 3D mesh.
    
    ##---------------------------------------------
    def tm_multiplane_test(self, heights, origin, normal, numdivs, path, infile, axis='y'):

        mesh = trimesh.load_mesh('%s/%s'%(path,infile))
        
        if type(origin)==vec3:
            origin = origin.aspt
        if type(normal)==vec3:
            normal = normal.aspt
    
        # get a single cross section of the mesh
        slice = mesh.section(plane_origin=mesh.centroid, 
                             plane_normal=[0,0,1])

        # the section will be in the original mesh frame
        #slice.show()
         
        # we can move the 3D curve to a Path2D object easily
        slice_2D, to_3D = slice.to_planar()

        if origin == None: 
            print('#WARNING tm_multiplane_test: origin not specified - using object bounds ')
            origin = mesh.bounds[0]
         
        # find a bunch of parallel cross sections
        sections = mesh.section_multiplane(plane_origin=origin, 
                                           plane_normal=normal, 
                                           heights=heights)

        cleansec = []
        for sec in sections:
            if sec:
                cleansec.append(sec)

        print('##tm_multiplane_test: sliced %s polygons '%len(cleansec) )

        pts2d = cleansec[0].vertices

        # summing the array of Path2D objects will put all of the curves into one Path2D
        combined = np.sum(cleansec)
        #combined.show()


        # # if we want to intersect a line with this 2D polygon, we can use shapely methods
        # polygon = slice_2D.polygons_full[0]
        # # intersect line with one of the polygons
        # hits = polygon.intersection(shp_ln([[-4,-1], [3,0]]))
        # # check what class the intersection returned
        # hits.__class__
        ###
        # # we can plot the intersection (red) and our original geometry(black and green)
        # ax = plt.gca()
        # for h in hits.geoms:
        #     ax.plot(*h.xy, color='r')
        # slice_2D.show()
        ###
        # the medial axis is available for closed Path2D objects
        #(slice_2D + slice_2D.medial_axis()).show()
        ###
        # #returns: 
        # # lines ((m,) sequence of (n, 2, 2) float) – Lines in space for m planes
        # # to_3D ((m, 4, 4) float) – Transform to move each section back to 3D
        # # face_index ((m,) sequence of (n,) int) – Indexes of mesh.faces for each segment

        """
        bbox = mesh.bounds.tolist()

        pts3d = combined.to_3D()
        #pts3d = combined.polygons_full
        
        #convert the Tracked array to NP array 
        npts = pts3d.vertices.view(np.ndarray)

        ptlist = [pts3d.vertices.tolist()] 
        
        """
        #ptlist = self.cvt_2d_to_3d(pts2d.tolist())

        """
        vflo = vectorflow()
        vflo.gr_polys = [ptlist]
        vflo._sort()
    
        #vflo.export_geojson_polygon(path, 'trimesh.json') 
        vflo.export_geojson_lines (path, 'trimesh.json') 
        vflo.export_geojson_lines(path, 'heed.json')
        vflo.cvt_grpoly_obj3d(objtype='singlepoly')
        vflo.save('heed.obj')
        
        #return trimesh.path.segments.clean(out) 
        #ee = trimesh.path.path.Path(out)
        """
        
        return slice.vertices 

    ##-----------------------------------
    def batch_ray_intersect(self, step, stacks, spokes, outname, axis='y'):
        """ 
           DEBUG - WAY TOO SLOW to use for real work
           DEBUG - only works with triangular gemoetry (you can run the triangulate command first )

           it had to iterate each face to check the intersections 

           probably only useful as a curiosity or toy 

           shoot a bunch of rays down from above in a circle to get the contours of an object
           raycast only works with triangles

           ARGS:
               step   - diameter change per iteration 
               stacks - number of concentric rings 
               spokes - number of spokes per ring 

        """
        if self.numtris==0:
            print('## batch_ray_intersect - NO TRIANGLES ')
            return None 

        DEBUG = True 
        MODE = 'radial' #default

        #MODE = 'radial'

        
        axis = axis.lower() 

        numx = 8 
        numy = 8

        # height to "shoot" from 
        #debug calc from 3d extents 
        dist = 1  
        
        #debug calc from 3d extents
        raylength = 2 

        flipnormal = False

        n = object3d()
        r = object3d()

        extents = self.calc_3d_bbox()
        cen = self.centroid() 
       
        stack_pts = [] 
        
        #--------------        
        # RECTANGULAR GRID
        if MODE =='grid':
            #set extents for teselator - calc them from geom 
            # [min_x, min_y, min_z, max_x, max_y, max_z ] 
            # [minx, miny, maxx, maxy] 
            if axis =='x':
                self.tesl._set_extents( self.calc_2d_bbox( 'x' ) )
            if axis =='y':
                self.tesl._set_extents( self.calc_2d_bbox( 'y' ) )
            if axis =='z':
                self.tesl._set_extents( self.calc_2d_bbox( 'z' ) )

            # build the grid using tesselator  
            self.tesl.build_2d_cells(numx,numy)
            for c in self.tesl.nodes:
                cen = (c.coord_x+(c.width/2), c.coord_y+(c.height/2))
                stack_pts.append((c.coord_x, c.coord_y, c.coord_z)) 
            stack_pts = [stack_pts]

         
        #--------------
        #RADIAL 
        if MODE =='radial':
            for s in range(1,stacks):
                #stack_pts.append( self.calc_circle(pos=(cen[0],cen[1],cen[2]), rot=(0,0,0), dia=step*s, axis='y', periodic=True, spokes=spokes) )
                stack_pts.append( self.calc_circle(pos=(0,0,0), rot=(0,0,0), dia=step*s, axis=axis, periodic=True, spokes=spokes) )



        curve_geom = [] 

        for ring in stack_pts:
            newcurve = [] 
            for pt in ring:
                
                # X axis
                if axis =='x':
                    ray = (vec3(dist, pt[1], pt[2]), vec3(-raylength, 0, 0))
                    if DEBUG:
                        r.one_vec_to_obj(ray[1], pos=ray[0])
                # Y axis
                if axis =='y':
                    ray = (vec3(pt[0], dist, pt[2]), vec3(0, -raylength, 0))
                    if DEBUG:
                        r.one_vec_to_obj(ray[1], pos=ray[0])
                # Z axis
                if axis =='z':
                    ray = (vec3(pt[0], pt[1], dist), vec3(0, 0, -raylength))
                    if DEBUG:
                        r.one_vec_to_obj(ray[1], pos=ray[0])

                hits = self.ray_hit( ray[0], ray[1], flipnormal)
                
                #n.prim_locator(pos=pt, size=.1)
                #n.one_vec_to_obj(ray[1], ray[0])
                
                if hits:
                    print('####### HITS ')
                    print(hits)

                for h in hits:
                    newcurve.append(h[1])
                    #n.prim_locator(h[1], size=.1)
                    
                    # g = self.get_face_geom(h[0], reindex=True, geom=None)
                    # n.insert(g)
                curve_geom.append(newcurve)

        for curve in curve_geom:
            if len(curve):
                #n.linegeom_fr_points(curve, periodic=True )
                n.linegeom_fr_points(curve, periodic=False )
        
        if DEBUG:
            r.save("%s_isect_rays.obj"%outname)
        
        n.save("%s_isects.obj"%outname)

    ##-----------------------------------
    def waterline(self, radius, spokes, heights, outname):
        """
            DEBUG - WAY TOO SLOW to use for real work
            it had to iterate each face to check the intersections 

            DEBUG - not done , may not ever work right
                    it works on convex objects but will probably fail on dimpled topology 

            it should return a list of lists of points 
            for not it builds and OBJ of intersection points and faces that intersect 

            do a radial raycast from world oringin facing inwards to get poly FIDs in a ring 
             - get the points from those to form a new polygon/plane (waterline) 
            
            ARGS:
                radius   - how far out 
                spokes   - how many rays to fire 
                heights  - how many rings and where to fire them 


        """

        extents = self.calc_3d_bbox()
        cen = self.centroid()   #pt = xzy
        dim = self.dimensions   # xyz - 3 floats 

        output = [] 

        ray_positions = [] 

        curve_geom =[]

        n = object3d()

        numheights = len(heights)

        for i,h in enumerate(heights):
            print("calculating ring %s of %s"%(i,numheights))
            ray_positions = self.calc_circle(pos=(0,h,0), rot=(0,0,0), dia=radius*2, axis='y', periodic=True, spokes=spokes) 

            
            newcurve =[]          
            for ray in ray_positions:
                aim = vec3()
                aim =  vec3(0,0,0) - vec3(ray) 

                hits = self.ray_hit( ray, aim, fastmode=True)
                if hits:
                    newcurve.append(hits[0][1])

                for h in hits:
                    #n.prim_locator(h[1], size=.1)
                    #g = self.get_face_geom(h[0], reindex=True, geom=None)
                    #n.insert(g)
                    pass 

                curve_geom.append(newcurve)

        #print(curve_geom)

        #print(ray_positions)
        for curve in curve_geom:
            if len(curve):
                #n.linegeom_fr_points(curve, periodic=True )
                n.linegeom_fr_points(curve, periodic=False )

        n.save("%s_waterline.obj"%outname)  

    ##-----------------------------------
    #def project_image(self, img_curves):
    #    pass 

    ##-----------------------------------
    def zigzag_intersect(self, size_x, size_y, rows, cols, overstep=0):
        """
           UNFINISHED 
           shoot a bunch of rays down from above in a square to get the contours of an object
        """
        
        # height to "shoot" from 
        yval = 5 

        n = object3d()

        extents = self.calc_3d_bbox()
        cen = self.centroid()   #pt = xzy
        dim = self.dimensions   # xyz - 3 floats 

        stack_pts = [] 

        for s in range(1,stacks):
            stack_pts.append( self.calc_circle(pos=(0,0,0), rot=(0,0,0), dia=step*s, axis='y', periodic=True, spokes=spokes) )

        
        curve_geom = [] 
        
        """
        for ring in stack_pts:
            newcurve = [] 
            for pt in ring:
                ray = (vec3(pt[0], yval, pt[2]), vec3(0, -1, 0))
                hits = self.ray_hit( ray[0], ray[1])
                print(hits)
                #n.prim_locator(pos=pt, size=.1)
                #n.one_vec_to_obj(ray[1], ray[0])
                for h in hits:
                    newcurve.append(h[1])
                    #n.prim_locator(h[1], size=.1)
                    # g = self.get_face_geom(h[0], reindex=True, geom=None)
                    ###WHY??## self.exprt_ply_idx = 1
                    # n.insert(g)
                curve_geom.append(newcurve)
        """ 

        for curve in curve_geom:
            if len(curve):
                n.linegeom_fr_points(curve, periodic=True )

        n.save("lotsa_locators.obj")


         

    ##-----------------------------------
    def zigzag_on_quad(self, fid, num):
        """ extract points on a 4 sided face 
            does not use raycasting - gets the edges of a face and sets points along those
        """

        out = [] 
        tmp = self.get_face_geom(fid, reindex=True) #returns [fidx, pts] 

        tmp2 = self.get_face_pts(fid)
        
        cuts = [] 

        print(tmp2)

        # walk the 4 edged as 2 point pairs - calc the in betweens 
        for i,pt in enumerate(tmp2):
            ep1 = tmp2[i]  
            ep2 = tmp2[i-1]
            cuts.append(self.locate_pt_along3d(ep1, ep2, num) )         

        #self.tesl.from_square_outline(cuts)

        # self.tesl.from_square_outline([
        #                                ['a','b','c','d'], 
        #                                ['e','f','g','h'], 
        #                                ['i','j','k','l'], 
        #                                ['m','n','o','p']
        #                               ])         

        # self.tesl.from_square_outline([
        #                             ['a','b','c','d','e','f'], 
        #                             ['g','h','i','j','k','l'], 
        #                             ['m','n','o','p','q','r'], 
        #                             ['s','t','u','v','w','x']
        #                            ]) 
 
        
        o = object3d()
        #o.insert(tmp) 
        o.insert(self) 

        grid = self.tesl.from_square_outline(cuts)
        
        #o.vectorlist_to_obj(grid) 
        o.pts_to_linesegment(grid)

        print(grid)

        # width = len(cuts) 
        # height = len(cuts[0])
        # for x in range(width):
        #     for y in range(height):
        #         o.prim_locator(cuts[x][y], size=.1)

        o.save('zigzag.obj')

        return cuts

    ##-----------------------------------
    def poly_bisect(self, outfile):
        """
           test of polygon raycast code  
           if face is a quad you could get the two "vertical" egdes as vectors 
        """

        # project_pt
        # locate_pt_along3d
        pop = object3d()

        flip = False 
        flipray = False 


        pop.prim_triangle('z',(0,0,-2),(0,70,0))

        pts = pop.points
        if flip:
            vc1 = vec3(pts[2])
            vc2 = vec3(pts[1])
            vc3 = vec3(pts[0])
        else:
            vc1 = vec3(pts[0])
            vc2 = vec3(pts[1])
            vc3 = vec3(pts[2])

        triangle = (vc1, vc2, vc3)

        cen = vec3() 
        cen.insert( pop.centroid(triangle))


        if flipray:
            ray = (vec3(0,0,1), vec3(0, 0, -1))
        else: 
            ray = (vec3(0,.5,0), vec3(-.2, 0, -3))

        test = vec3() 
        
        #result = test.poly_intersect(ray, triangle)

        result = test.ray_tri_intersect(ray[0], ray[1], vc1, vc2, vc3)
        
        o = object3d()

        #ray origin
        o.prim_locator(ray[0], size=.1)
        #ray vector 
        o.one_vec_to_obj(ray[1], ray[0], arrowhead=True)

        #polygon geom 
        o.insert_polygons(plyids=[(1,2,3)], points=triangle)

        if result:
            #hit location
            o.prim_locator((result[1]), size=.1)
            #polygon normal 
            #o.pts_to_linesegment([cen, cen+(result[2])] )

        #o.save(outfile)
        print(result)

    ##-----------------------------------


    def delaunay(self, object, height):
        #if you could get outline at a Z value - (and spiral) - you have working cam 
        pass 


    ##-----------------------------------    
    def dialate(self):
        pass 

    ##-- 

    def erode(self):
        pass 

    ##--    

    def scanlines(self, plyid):
        """ 
        DEBUG WIP 

        SEE scan_line_tool() in obj_to_cnc for a similar thing 

        run a 3d scanline across a polygon at any angle 

        start with a triangle 
        get the normal vector 
        build a grid of normal vecotrs orthoganal to the face 
        fire a ray from each one 


        return hits in the following order:
           1_______2              1_______2
          /        \             /        \   
        3/__________\4         6/__________\3 
         \           /          \           /     
         5\_________/6          5\_________/4 

        """
        #self.triangulate()

        if plyid is not None:
            print(self.get_face_geom(plyid))

   
         

    ##-----------------------------------
    def face_sprial(self):
       """ recursive erode->scanline -> repeat == spiral  

       """

       pass 

##-----------------------------------------##



##############################################################

## CLIPPER LIBRARY TESTS 

"""
def clipper_test1():
    import pyclipper

    subj = (
        ((180, 200), (260, 200), (260, 150), (180, 150)),
        ((215, 160), (230, 190), (200, 190))
    )
    clip = ((190, 210), (240, 210), (240, 130), (190, 130))

     
    #convert 2d into 3d 
    obj = object3d()

    pts = obj.cvt_2d_to_3d(subj[0])
    obj.linegeom_fr_points(pts) 

    pts = obj.cvt_2d_to_3d(subj[1])
    obj.linegeom_fr_points(pts) 


    pts = obj.cvt_2d_to_3d(clip)
    obj.linegeom_fr_points(pts) 

    #run clipper
    pc = pyclipper.Pyclipper()
    pc.AddPath(clip, pyclipper.PT_CLIP, True)
    pc.AddPaths(subj, pyclipper.PT_SUBJECT, True)

    solution = pc.Execute(pyclipper.CT_INTERSECTION, pyclipper.PFT_EVENODD, pyclipper.PFT_EVENODD)

    for s in solution:
        pts = obj.cvt_2d_to_3d(s)
        obj.linegeom_fr_points(pts) 

    obj.save('fooker.obj')


def clipper_test2():
    import pyclipper

    subj = ((348, 257), (364, 148), (362, 148), (326, 241), (295, 219), (258, 88), (440, 129), (370, 196), (372, 275))
 
    pco = pyclipper.PyclipperOffset()
    pco.AddPath(subj, pyclipper.JT_ROUND, pyclipper.ET_CLOSEDPOLYGON)

    solution = pco.Execute(-5)

     
    #convert 2d into 3d 
    obj = object3d()

    pts = obj.cvt_2d_to_3d(subj)
    obj.linegeom_fr_points(pts) 


    for s in solution:
        pts = obj.cvt_2d_to_3d(s)
        obj.linegeom_fr_points(pts) 

    obj.scale_pts((.1,.1,.1))
    obj.save('fooker.obj')

"""    


