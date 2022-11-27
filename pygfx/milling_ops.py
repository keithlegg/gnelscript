#!/usr/local/bin/python3



import math 
import re 

#from pygfx.point_ops import polygon_operator
#from pygfx.math_ops import vec3 

import os 

from gnelscript.pygfx.obj3d import object3d

from gnelscript.pygfx.point_ops import *



#from pygfx.gcode.bridgeport import  parser_commands
from gnelscript.pygfx.gcode.linuxcnc import  parser_commands



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



##------------------------------------------

class generate_gcode(object):
    def __init__(self):
        self.coords = []
        self.outfile = []

    def savengc(self, filename):
        print(os.getcwd() )

        fobj = open( filename,"w") #encoding='utf-8'
        for line in self.outfile: 
            fobj.write(line+'\n')
        fobj.close()


    def lineartest(self):
        pop = point_operator()
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
        




##------------------------------------------


class gcode(object3d):

    def __init__(self):
        super().__init__()  

        self.DEBUG_MODE = True
        # self.linear_units = 'in'
        # self.z_axis  ?? 

        #comments get dumped into this array (full or partial line)
        self.commented = [] # [index, string] 
        
        self.param_names  = [] # [name, value ]
        self.param_values = [] # [name, value ]

        #swappable dialects of gcode commands 
        self.dialect = parser_commands

        self.coord_words = ['N','G', 'X','Y','Z','U','V','W','I','J','K','R','P','F'] # ,'A','B','C','D']
        
        # simulated position of the cutting head
        self.POSX = 0
        self.POSY = 0
        self.POSZ = 0

        self.segments  = [] # [index, command, xyz_pos] 




    def show_data(self):
        print("## ## ## ## ## ## ## ## ## ##")

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

    def parse_params(self, string):
        """ discover and retain any parameters in the file """
        if PARAM in string:
            tmp = self.between_token_list(string, [PARAM, '>']) 

            if len(tmp)==2:
                if tmp[0] not in self.param_names:
                    if '=' in tmp[1]:
                        tmp2 = tmp[1].split('=')
                        self.param_names.append(tmp[0])
                        self.param_values.append(tmp2[1].replace(' ',''))
                else:
                    #IT IS ALREADY PARSED!
                    #deal with this here or in another function??
                    pass

    def between_token_list(self, string, tokens, return_match=False):
        """ give a list of tokens, return a list of betweens 
            
            string = 'y987a123b541c307d999'
            tokens = ['a','b','c','d']
            out = self.between_token_list( string, tokens)

            out will be:
                 ['123', '541', '307']

            unless return_match is True,
                then out will be:
                    [ ['a',123'], ['b',541'], ['c',307'] ]

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
        
        #for seg in self.segments:
        #    #f.write(seg[])
        #    print( seg[0], seg[2])

    def load_gcode(self, filename):

        if os.path.lexists(filename) == 0:
            self.scribe("%s DOES NOT EXIST !! "%filename )
            #raise
            
        if os.path.lexists(filename):
            f = open( filename,"r", encoding='utf-8')
            contents = f.readlines()


            #scan entire file for parametrs first 
            for lin in contents:
                self.parse_params(lin)

            print("#### params ", self.param_names )

            cleaned_contents = []

            for lin in contents:

                #if self.DEBUG_MODE:
                #print("#### LINE IS ", lin.replace("\n","") )

                # ignore commented out lines 
                if lin[0]==COMMENT or lin[0]=='\'':  
                    self.commented.append(lin)

                else:


                    # seperate the line index out and split from rest of command
                    lindex = self.between_token_list(lin, self.coord_words, return_match=True)
                    #n_idx = lindex[0]   
                    #comm = lin[len(str(n_idx)):] 
                    if lindex:
                        #print('############### ', lindex)

                        # execute the parameters(variables) if they exist 
                        for i,tok in enumerate(lindex):
                            for j,p in enumerate(self.param_names):
                                if p in tok[1]:
                                    # substitute the value with the name 
                                    tok[1] = (tok[1].replace('#<%s>'%self.param_names[j], self.param_values[j]) )
                                    try:
                                        # DONT LOAD ANY GCODE TEXTFILES WITH EXECUTABLE PYTHON IN THEM!
                                        # THIS IS A SECURITY HOLE FOR WANKERS TO EXPLOIT 
                                        
                                        # eval the damn thing and capture the result 
                                        #print( tok[0], eval(tok[1])[0])
                                        pass
                                        
                                    except:
                                        pass    

                        print( lindex )

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





