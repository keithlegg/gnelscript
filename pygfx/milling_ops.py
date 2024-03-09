# ***** BEGIN GPL LICENSE BLOCK *****
#
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# ***** END GPL LICENSE BLOCK *****



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


GUI
laser pin?



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

#from pygfx.point_ops import polygon_operator
#from pygfx.math_ops import vec3 

import os 

from gnelscript import SHAPELY_IS_LOADED

if SHAPELY_IS_LOADED:
    from shapely import Point, LineString, Polygon


from gnelscript.pygfx.obj3d import object3d
from gnelscript.pygfx.obj2d import object2d

from gnelscript.pygfx.point_ops import *
from gnelscript.pygfx.math_ops import *

from gnelscript.pygfx.gis_vector_ops import *

from gnelscript.pygfx.grid_ops import tessellator


#from pygfx.gcode.bridgeport import  parser_commands
from gnelscript.pygfx.gcode.linuxcnc import  parser_commands


import geojson
from geojson import Point, Feature, LineString, FeatureCollection, dump




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
    ## first stab of gcode translator - got put on hold 
    ## kicad ops has the first working simple gcode exporter
    
    def __init__(self):
        super().__init__()  

        self.DEBUG_MODE = True
        # self.linear_units = 'in'
        # self.z_axis  ?? 

        #comments get dumped into this array (full or partial line)
        self.commented = [] # [index, string] 
        
        self.param_names  = [] # [name, value ]
        self.param_values = [] # [name, value ]
        # #swappable dialects of gcode commands 
        # self.dialect = parser_commands

        self.coord_words = ['N','G', 'X','Y','Z','U','V','W','I','J','K','R','P','F'] # ,'A','B','C','D']
        
        # simulated position of the cutting head
        self.POSX = 0
        self.POSY = 0
        self.POSZ = 0

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
    def milling_test_linear(self):
        pop = point_operator_3d()
        #calc_circle(self, pos=(0,0,0), rot=(0,0,0), dia=1, axis='z', periodic=True, spokes=23):
        circ_pts = pop.calc_circle()


        # n099    (This is a test plot nc program to be run on backplot)
        # n100    (Author Ray Henry 10-Feb-2000)
        # n101    g20
        # n102    g0 x0 y0 z0 f30
        # n103    x1 y1(start xy circle)
        # n104    g17 g02 i.5 j.5
        # n106    g0 z.1 (add xy lettering)
        # n107    y1.75
        # n108    z0
        # n109    g1 y1.25 x1.4
        # n110    y1.5 x1.2
        # n111    y1.25 x1
        # n112    y1.75 x1.4
        self.outfile.append('g20')                  #inches for unit 
        self.outfile.append('g0 x0 y0 z0 f30')      #rapid move to 0 
        for i,pt in enumerate(circ_pts):
            if i==0:
                self.outfile.append( 'g1 x%s y%s'%(pt[0], pt[1]) )
            else:    
                self.outfile.append( 'g1 x%s y%s'%(pt[0], pt[1]) )

        #rapid move at end 
        #self.outfile.append('g0z1')
        self.outfile.append('g0 x0 y0 z0 f30')      #rapid move to 0 

        #program end 
        self.outfile.append('%')
        self.outfile.append('m2')


    def show_data(self):
        print("## ## ## ## ## ## ## ## ## ##")

        for s in self.segments:
            print(s[1])

    def contains_coord_words(self, string): 
        # check for any known coordinate words
        has_coords = False 
        for cw in self.coord_words:
            if cw in string:
                has_coords = True
        return  has_coords       

    def which_coord_words(self, string):
        words = []
        for cw2 in self.coord_words:
            if cw2 in string:
                words.append(cw2)
        return words

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


    def save_gcode(self, filename):
        f = open( filename,"w", encoding='utf-8')
        
        for seg in self.segments:
            #f.write(seg[])
            print( seg[0], seg[2])


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
                gcc = gcode_op()
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

##---------------------------##

class gcode_op(gcode_object, vectorflow ):
    """ copy of vectorflow()  
    """


    def __init__(self):
        super().__init__()  

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
        pop = point_operator_3d()
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
                depth = self.ch+(i*self.cdpi)
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
    def _calculate_paths3d(self, do_retracts=True, doround=True):
        """ 
             DEBUGGY 
        """
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

    ##-------------------------------##
    def _calculate_paths2d_mop(self, do_retracts=True):
        """ 
            this walks self.gr_poly and builds an OBJ and NGC file in self.outfile


            https://linuxcnc.org/docs/html/gcode/g-code.html
            Top 10 tasty GCODE commands:
                S - surface speed
                G20 (Use inch)
                G21 (Use mm)
                G90 (Set Absolute Coordinates)
                G0 -  Rapid Move
                G1 -  Linear Move
                M3, M4, M5  S ($)   Spindle Control
                M6 - 
                M9 - 
                M3 - 
                M2 - program end 

        """

        lastpt = (0,0,0)
        
        self.outfile.append('(exported with _calculate_paths2d )')
        self.outfile.append('(linear scale set to %s of internal coordinates)'%self.global_scale )
        self.outfile.append('  ')

        self.outfile.append('g20') # inches for unit 
        #self.outfile.append('G21')  # mm for unit 

        ##-----------------------------------------##

        #move to origin  
        self.outfile.append('g0 x%s y%s z%s f30'% ( self.hp[0], self.hp[1], self.rh) )   #rapid move to 0 
        self.ngc_to_obj.append( ( self.hp[0], self.hp[1], self.rh) ) 

        if do_retracts:
            self.outfile.append('g0z%s'% ( self.rh ) )  #retract in between cuts
   
        ##-----------------------------------------##
        self.outfile.append('  ')
        self.outfile.append('(exporting filled polygons )')

        # build the gcode up with simple linear movements 
        ##------------------------------

        #DEBUG - FILL_POLYS ARE A LAYOVER FROM KICAD STUFF - UNTESTED 
        for fill_poly in self.filled_polys:
            if do_retracts:
                pt1 = fill_poly[0] 
                self.outfile.append('G0 x%s y%s z%s'% (  pt1[0], pt1[1], self.rh ) )  #first point at retract height   
                self.ngc_to_obj.append( ( self.hp[0], self.hp[1], self.rh) ) 
            
            for i,pt in enumerate(fill_poly):
                self.outfile.append( 'G1 x%s y%s z%s'%( pt[0], pt[1], self.ch ) )
                self.ngc_to_obj.append( ( pt[0], pt[1], self.ch ) )             
                lastpt =( pt[0], pt[1], self.ch )

            self.outfile.append( 'G0 x%s y%s z%s'%(lastpt[0], lastpt[1], self.ch) )
            self.ngc_to_obj.append( (lastpt[0], lastpt[1], self.ch)   ) 
            if do_retracts:
                self.outfile.append('g0z%s'% ( self.rh ) )  #retract in between cuts
                self.ngc_to_obj.append( (lastpt[0], lastpt[1], self.rh)   ) 
                self.outfile.append('  ')

        ####
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
                #DEBUG - need to sort out clean points - want to run as close rto final xport 
                #it looses precision 

                # no formatting (full precision)
                pt1 = gr_poly[0]
                # do round to get rid of bad string formatting ( exponent in floats) 
                #pt1 = self.clean_pts_str(gr_poly[0])

                ##-- 

                #### first point at retract height 
                if do_retracts:
                    #move to first point RH 
                    self.outfile.append('x%s y%s z%s'% (  pt1[0] , pt1[1], self.rh ) )               
                    self.ngc_to_obj.append( ( pt1[0], pt1[1], self.rh ) )             

                ## iterate points in polygon 
                self.outfile.append( 'G1' )
                for i,pt in enumerate(gr_poly):
                    self.outfile.append( 'x%s y%s z%s'%( pt[0], pt[1], self.ch ) )
                    self.ngc_to_obj.append( (pt[0], pt[1], self.ch) )                   
                
                self.outfile.append( 'G0' )

                # move to last point at CH  
                #self.ngc_to_obj.append( (gr_poly[0][0], gr_poly[0][1], self.ch))   
                #self.outfile.append( 'x%s y%s z%s'%( (gr_poly[0][0], gr_poly[0][1], self.ch) ) )

                if do_retracts:
                    self.ngc_to_obj.append( (gr_poly[0][0], gr_poly[0][1], self.rh)  )
                    self.outfile.append( 'x%s y%s z%s'%( gr_poly[0][0], gr_poly[0][1], self.rh) )

                    #### retract in between cuts
                    self.outfile.append('g0z%s'% ( self.rh ) )  

                self.outfile.append('  ')

        ##-----------------------------------------##        
        # self.outfile.append('(exporting segments )')

        ##-----------------------------------------##
        # rapid move at end 
        self.outfile.append('m2') #program end

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
            #    features.append(Feature(geometry=LineString(coordinates=self.extents_fr_bbox(s[2])), properties={"id":i, "len" : len(s[4]) } ))
            coords = self.extents_fr_bbox(s[2], periodic=True)

            features.append(Feature(geometry=LineString(coordinates=coords), 
                                 properties={"id" : 0 
                                            }
                                 ) 
                         )

        #sort extents  
        coords = self.extents_fr_bbox([self.sort_minx,self.sort_miny,self.sort_maxx,self.sort_maxy], periodic=True)
        features.append(Feature(geometry=LineString(coordinates=coords), 
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
                        features.append(Feature(geometry=LineString(coordinates=tmp_pts)) 
                                       )

        feature_collection = FeatureCollection(features)
        with open('%s/_distances.json'%(folder), 'w') as f:
            dump(feature_collection, f)


##---------------------------##

class cam_op(gcode_op):
    
    def __init__(self):
        super().__init__()  
        self.tesl = tessellator() 

    def obj_to_wkt(self):
        print(self.points)

    if SHAPELY_IS_LOADED:
        def test(self):
            line = LineString([(2, 0), (2, 4), (3, 4)])
            print(line)

    ##-----------------------------------
    def batch_ray_intersect(self, step, stacks, spokes, outname, axis='y'):
        """ 
           DEBUG - WAY TOO SLOW to use for real work
           it had to iterate each face to check the intersections 

           probably only useful as a curiosity or toy 

           shoot a bunch of rays down from above in a circle to get the contours of an object
           raycast only works with triangles

           ARGS:
               step   - diameter change per iteration 
               stacks - number of concentric rings 
               spokes - number of spokes per ring 

        """
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

##------------------------------------------



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


