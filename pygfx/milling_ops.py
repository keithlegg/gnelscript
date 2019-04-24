#!/usr/local/bin/python3



import math 
import re 

#from pygfx.point_ops import polygon_operator
#from pygfx.math_ops import vec3 

import os 

from pygfx.obj3d import object3d




##------------------------------------------

## GCODE KNOWN DIALECTS 



#  X-Axis (table)  
#  Y-Axis (saddle)  
#  Z-Axis (quill)  

#-------------------------------------------

# Address of Coordinate Word
# Basic Axes       Meaning
# X                X-Axis Move Absolute or Incremental
# Y                Y-Axis Move Absolute or Incremental
# Z                Z-Axis Move Absolute or Incremental
# U                X-Axis Move Incremental
# V                Y-Axis Move Incremental
# W                Z-Axis Move Incremental
#
# Parameters for Circular Interpolation (G75)
# I,J,K            Absolute: Coordinates of center of circle.
#                  Incremental: Signed distance from arc start point to arc center
#
# Polar and Spherical Coordinates
#
# R                 Pole Radius  
# I,J,K             Absolute: Coordinates of center of pole.      
#                   Incremental: Signed distance from pole start point to center of pole.  
# A                 Angular Distance, (Longitude) Degrees Absolute or Incremental  
# B                 Angular Distance, (Longitude) Incremental  
# E                 Agular Distance, (Colatitude) Degrees Absolute

#-------------------------------------------


# Preparatory Codes 
# G-Code   Group     Function

# 0         1        Rapid Traverse
# 1         1        Linear Interpolation (Feed)
# 2         1        Circular Interpolation Clockwise
# 3         1        Circular Interpolation Counterclockwise
# 4         0        Dwell
# 8         11       Modal Deceleration Override Off
# 9         11       Modal Deceleration Override On
# 12        0        Helical Interpolation CW
# 13        0        Helical Interpolation CCW
# 17        2        XY Plane Selection
# 18        2        ZX Plane Selection
# 19        2        YZ Plane Selection
# 22        0        Circular Interpolation, Fillet Input CW
# 23        0        Circular Interpolation, Fillet Input CCW
# 30        3        Mirror Image Off
# 31        3        Mirror Image X On
# 32        3        Mirror Image Y On
# 40        4        Cutter Diameter Offset Off
# 41        4        Cutter Compensation Left
# 42        4        Cutter Compensation Right
# 44        5        Cutter Compensation, Normal Feedrate
# 45        5        Cutter Compensation, Modify Feedrate
# 48        12       Corner Rounding in Cutter Comp Off
# 49        12       Corner Rounding in Cutter Comp On
# 70        6        Input in Inch
# 71        6        Input in Millimeter
# 72        7        Transformation Off
# 73        7        Transformation/Rotation, Scaling
# 74        8        Multi-quadrant Circle Input Off
# 75        8        Multi-quadrant Circle Input On
# 77        1        Zig-Zag Mill Cycle
# 78        1        Pocket Mill Cycle
# 79        1        Bore Mill Cycle
# 80        1        Drill Cycle Off
# 81        1        Z Cycle, Drill (Feed In, Rapid Out)
# 82        1        Z Cycle, Spot Face (Feed In, Rapid Out)
# 83        1        Z Cycle, Deep Hole (Peck, Rapid Out)
# 84        1        Z Cycle, Tap (Feed In, Feed Out)
# 85        1        Z Cycle, Bore (Feed In, Feed Out)
# 86        1        Z Cycle, Bore (Feed In, Stop-Wait, Rapid Out
# 87        1        Z Cycle, Chip Break (Peck, Rapid Out)
# 89        1        Z Cycle, Bore (Feed In, Dwell, Feed Out)
# 90        9        Absolute Programming
# 91        9        Incremental Programming
# 92        0        Preset Part Programming Zero Point
# 94        13       Feedrate Per Minute Mode
# 95        13       Feed Per Spindle Revolution (pitch) mode.
# 96        10       Restore Base Part Program Coordinate System
# 97        10       Set Work Coordinate System
# 99        0        Deceleration Override
# 170       1        Outside Frame Mill
# 171       1        Inside Frame Mill
# 172       1        Pocket Frame Mill
# 173       1        Outside Face Mill
# 174       1        Inside Face Mill
# 175       1        Outside Circle Mill
# 176       1        Inside Circle Mill
# 177       1        Pocket Circle Mill
# 179       1        Slot Mill
# 181-189   1        Z Cycle (Same as G81-G89) Multi-Hole
# 191-199   1        Z Cycle (Same as G81-G89) Frame of Holes

##------------------------------------------
# Power On and Reset State of G-Codes

# G-Code   Group   Function
# 0        1       Rapid Traverse
# 8        11      Modal Deceleration Override Off
# 17       2       XY Plane Selection
# 30       3       Mirror Image Off
# 40       4       Cutter Compensation Off
# 45       5       Cutter Compensation, Modify Feedrate
# 49       12      Corner Rounding in Cutter Comp On
# 70/71    6       Input in Inch or Metric, Non-volatile
# 72       7       Transformation Off
# 75       8       Multi-quadrant Circle Input On
# 90       9       Absolute Programming
# 94       13      Feedrate Per Minute mode
# 96       10      Base Coordinate System



bridgeport_gcode =  {
    "M8": "coolant_on", #?
    "M9": "coolant_off", #?
    

    "M5": "spindle_off", #?

    # G1 Code sets the Linear Interpolation mode. The format of the G1 command is:
    # G1X__Y__Z__F__;
    # Where X__Y__Z__ defines the endpoint of the move to be made. Simultaneous XYZ motion along a linear (straight
    # line) path occurs at the feedrate defined by the F__ word.    

    # "G0" : "rapid_traverse", #  G0/G1 - R__I__J__A__; Move to Pt
    "G00": "rapid_traverse",

    # "G1" : "linear_interpolation", 
    "G01": "linear_interpolation",

    # "G2": "circular_interpolation",        
    "G02": "circular_interpolation",

    # "G4" : "dwell",
    "G40" : "dwell",

    "G12": "helical_1", # G12/G13 - A__Z__F__; Do Helix 
    "G13": "helical_2", 

    #Circular interpolation is effective in the XY, YZ or ZX planes by preparatory codes as follows:
    "G17": "circle_XY_plane",
    "G18": "circle_ZX_plane",
    "G19": "circle_YZ_plane",

    "G70": "inch_conversion",
    "G71": "metric_conversion",

    "G80": "drill_cycle",

    "G90": "abs_progr",  # data in a word represents the coordinate value of the point from part program zero.
    "G91": "incr_progr", # data in a word is distance along designated axis from the existing position to the desired position

    # canned milling cycles 
    # G170: Outside Frame Mill
    # G171: Inside Frame Mill
    # G172: Pocket Frame Mill
    # G173: Outside Face Mill
    # G174: Inside Face Mill
    # G175: Outside Circle Mill
    # G176: Inside Circle Mill
    # G177: Pocket Circle Mill
    # G179: Slot Mill

    "G175": "inside_mill_cycle",   # G175 X__Y__Z__R__Z__Z__P__F__P__F__F__;
    "G176": "outside_mill_circle", # G176 X__Y__Z__R__Z__Z__P__F__P__F__F__;

    # "X": "x_axis",
    # "Y": "y_axis",
    # "Z": "z_axis",  

    # X-Y Plane (G17)
    # G2X__Y__I__J__F__; CW
    # G3X__Y__I__J__F__; CCW
    
    # Z-X Plane (G18)
    # G2X__Z__I__K__F__; CW
    # G3X__Z__I__K__F__; CCW
    
    # Y-Z Plane (G 19)
    # G2Y__Z__J__K__F__; CW
    # G3Y__Z__J__K__F__; CCW

    #If G75 is active and G90 (absolute data input) is also active, I__, J__, K__ specifies the arc center coordinates with
    #respect to part program zero. If G91 (incremental data input) is active, I__, J__, K__ specifies the signed distance from
    #the arc startpoint to the center of the arc.
    
    #If G74 is active, I__, J__, K__ specifies the unsigned magnitude of the distance from the arc startpoint to the arc center
    #in G90 (Absolute) or G91 (Incremental) mode.

    # "I": "i_axis",
    # "J": "j_axis",
    # "K": "k_axis",  

    # "U": "u_axis",
    # "V": "v_axis",
    # "W": "w_axis", 

    #"M3": "end?", 
    "M30": "end",

    "E": "end",
}


##------------------------------------------



##------------------------------------------

class generate_gcode(object):
    def __init__(self):
        pass
        

##------------------------------------------


class gcode_to_polyline(object3d):

    def __init__(self):
        super().__init__()  

        self.DEBUG_MODE = True
                
        # self.linear_units = 'in'
        # self.z_axis  ?? 

        #comments get dumped into this array (full or partial line)
        self.commented = [] # [index, string] 

        #swappable dialects of gcode commands (Bridgeport is the only one for now)
        self.dialect = bridgeport_gcode

        self.coord_words = ['X','Y','Z','U','V','W','I','J','K','R','P','F'] # ,'A','B','C','D']
        
        # simulated position of the cutting head
        self.POSX = 0
        self.POSY = 0
        self.POSZ = 0

        self.segments  = [] # [index, command, xyz_pos] 




    def show_data(self):
        for s in self.segments:
            print(s[1])

    def save_3d_object(self, filename):
        """ dump the XYZ positions into a polyline object """

        obj   = object3d() # container for 3D object 
        data = []

        line_data = []
        for pt in self.segments:
            line_data.append( pt[1] )

        obj.pts_to_linesegment(line_data)
        obj.save( filename )

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

    def between_token_list(self, string, tokens):
        """ give a list of tokens, return a list of betweens 
            
            string = 'y987a123b541c307d999'
            tokens = ['a','b','c','d']
            out = self.between_token_list( string, tokens)

            out will be:
                 ['123', '541', '307']
        """

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

        tmp = re.split( match, string)
        out = []
        for t in tmp:
            tmp2 = t.replace('\n','')
            tmp2 = tmp2.replace(';','')
            out.append(tmp2)
        return out[1:]  


    def save_gcode(self, filename):
        f = open( filename,"w", encoding='utf-8')
        
        #for seg in self.segments:
        #    #f.write(seg[])
        #    print( seg[0], seg[2])


    def load_gcode(self, filename):

        #if os.path.lexists(filename) == 0:
        #    self.scribe("%s DOES NOT EXIST !! "%filename )
        #    #raise
            
        if os.path.lexists(filename):
            f = open( filename,"r", encoding='utf-8')
            contents = f.readlines()

            for lin in contents:

                #if self.DEBUG_MODE:
                #    print("#### LINE IS ", lin.replace("\n","") )

                #seperate the line index out and split from rest of command
                lindex = self.between_token_list(lin, ['N', 'G', 'M'])
                n_idx = lindex[0]   
                comm = lin[len(str(n_idx)):]   
                
                # print("########### ", lin , n_idx , comm)

                if comm:
                   
                    # ignore commented out lines 
                    if lin[0]==';' or lin[0]=='\'':  
                        self.commented.append(lin)

                    #scan if not commented     
                    else:
                        #check for comments on the line but not at start   
                        if ';' in lin:
                            tmp = lin.split(';')
                            lin = tmp[0]
                            if tmp[1] is not '\n':
                                self.commented.append(tmp[1])
                        

                          
                        # check for known commands 
                        for key in self.dialect :
                            if key in comm:
                                print(" COM FOUND at INDEX %s ! "%n_idx , key, '---', self.dialect [key] )

                           

                        # check for any known coordinate words
                        has_coords = self.contains_coord_words(comm) 

                        # if has coordinate words, determine which ones it has 
                        if has_coords is True:             

                            line_contains = ( self.which_coord_words(comm) )
                            output        = self.between_token_list(comm, line_contains)

                            #print( line_contains, output )
                            for i,token in enumerate(line_contains):
                                if token=='X':
                                    self.POSX = float(output[i])
                                if token=='Y':
                                    self.POSY = float(output[i])
                                if token=='Z':
                                    self.POSZ = float(output[i])

                    self.segments.append( [n_idx, [self.POSX, self.POSY, self.POSZ], comm ] )



