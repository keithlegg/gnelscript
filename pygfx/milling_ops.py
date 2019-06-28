#!/usr/local/bin/python3



import math 
import re 

#from pygfx.point_ops import polygon_operator
#from pygfx.math_ops import vec3 

import os 

from pygfx.obj3d import object3d


#from pygfx.gcode.bridgeport import  parser_commands
from pygfx.gcode.linuxcnc import  parser_commands





##------------------------------------------

## GCODE KNOWN DIALECTS 



#  X-Axis (table)  
#  Y-Axis (saddle)  
#  Z-Axis (quill)  


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

        #swappable dialects of gcode commands 
        self.dialect = parser_commands

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



