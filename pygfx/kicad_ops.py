
import re 


#from pygfx.obj3d import object3d
from pygfx.milling_ops import gcode


## class kicad_line(object):
##     """container for a kicad pad"""    
##     def __init__(self):
##         self.width    = 0
##         self.start_xy = []
##         self.end_xy   = []
## 
## class kicad_pad(object):
##     """container for a kicad pad"""    
##     def __init__(self):
##         self.name   = ''
##         self.size   = 0
##         self.xcoord = 0
##         self.ycoord = 0
## 
## class kicad_module(object):
##     """container for a kicad module"""
##     def __init__(self):
##         self.name  = ''
##         self.lines = []
##         self.pads  = []



"""
    <- point_operator (contains math_util)
      <- polygon_operator 
        <- object3d 
          <- gcode (linuxcnc)
            <- pcbfile  (kicad)
"""


class pcbfile(gcode):

    def __init__(self):
        super().__init__()  

        self.parse_depth = 0
        self.module_depth = 0 

        self.file_contents = []
        
        self.oneup_entity = None 
        self.cur_entity = None 
        self.cur_module = None 

        self.modules = []


        self.known_entities = ['module','gr_line','segment']

        self.kicad_units = 'mils'
        self.object_units  = 'cm'

        #self.GLOBAL_SCALE = 0.3904 # unit conversion (cm to mills?)

    ##############
    def show_geom(self):
        #call the inherited polygon "show" method
        self.show()

    ##############
    def show_modules(self):
        for m in self.modules:
            print(m.name)

    def save_3d_obj(self, name):
       
       #make sure we have numeric data not strings
       #they are string because we parsed from a file 
       newpts = []

       print( self.points )

       for pt in self.points:
           print('## pt ',pt)  
           tmp = [] 
           for n in pt:
              tmp.append(float(n) )  
           newpts.append( tmp )

       self.points = newpts 

       self.scale_pts( ( -1.0,-1.0, 1.0) )  

       self.save(name)

    ##############
    def load_kicadpcb(self, filename):
        """ a parser that is not recursive, but clever enough to scan all the 
            file entities and know what module they are in 
        """ 
        
        var_module_name    = ''
        var_module_pos     = []
        var_module_lines   = []

        var_pad_name       = ''
        var_pad_size       = 0
        var_pad_xcoord     = 0
        var_pad_ycoord     = 0

        var_segment_start  = []
        var_segment_end    = []

        var_line_width     = 0
        var_line_start_xy  = []
        var_line_end_xy    = []


        f = open(filename, 'r')
        for line in f:
            if '(' in line or ')' in line:
                toked = line.split(' ') 
                for i,tok in enumerate(toked):
                    if tok !='':
                      
                        # if we encounter a "(" - go deeper into the parse
                        for ct in re.findall('\(',tok):
                            
                            # often, we only care about the entity immediately above current. store what it is 
                            if self.cur_entity is not None:
                                self.oneup_entity = self.cur_entity

                            self.cur_entity = self.scrub(tok[1:]) #name of entity immediately after "("

                            # if module depth is set we are parsing a module, if it matches depth we must be done traversing     
                            if self.parse_depth == self.module_depth and self.cur_module is not None:
                            

                                # make a new "module" container object to store what we found in file 
                                # make it when we exit because all the constituent pieces will be scanned at this point 
                                #new_mod = kicad_module()
                                #new_mod.name = self.cur_module
                                #self.modules.append( new_mod )
                                
                                #------ 
                                #step out of the module 
                                self.module_depth = 0
                                self.cur_module = None 
                                
                            #-------------------------------
                            # if module depth is non zero we are parsing a module  
                            # most of the parsing happens in this block 
                            if self.module_depth>0:
                               #  print("we are in a module %s"%self.cur_module) 


                                if tok[1:] == 'at' and self.oneup_entity == 'module':
                                    print("MODULE %s AT %s %s" %(self.cur_module, toked[i+1], toked[i+2]) ) 

                                if tok[1:] == 'at' and self.oneup_entity == 'pad':
                                    #new_pad = kicad_pad()
                                    print("PAD %s AT %s %s" %(self.cur_module, toked[i+1], toked[i+2]) ) 


                            #-------------------------------
                            # if we found a module - store the depth it was found at and the name     
                            if tok[1:] == 'module':
                                self.module_depth = self.parse_depth 
                                self.cur_module = self.scrub(toked[i+1])


                            #-------------------------------
                            # parsing objects outside of modules 

                            #-------------------
                            # arc parsing 
                            if tok[1:] == 'start' and self.oneup_entity == 'gr_arc':
                                print("gr_arc start", toked[i+1], toked[i+2]  ) 

                            #-------------------
                            # via parsing 
                            if tok[1:] == 'start' and self.oneup_entity == 'via':
                                print("via start", toked[i+1], toked[i+2]  ) 

                            #-------------------
                            # polygon parsing 
                            if tok[1:] == 'pts' and self.oneup_entity == 'gr_poly':
                                print("gr_poly start", toked[i+1], toked[i+2]  ) 

                            #-------------------
                            # circle parsing 
                            if tok[1:] == 'center' and self.oneup_entity == 'gr_circle':
                                print("gr_circle ", toked[i+1], toked[i+2]  ) 

                            #-------------------

                            # segment parsing 
                            if tok[1:] == 'start' and self.oneup_entity == 'segment':
                                print("segment start ", toked[i+1], toked[i+2]  ) 
                                var_segment_start  = [toked[i+1], toked[i+2]]
                            if tok[1:] == 'end' and self.oneup_entity == 'segment':
                                print("segment end ", toked[i+1], toked[i+2]  ) 
                                var_segment_end  = [toked[i+1], toked[i+2]]

                            if  var_segment_start and var_segment_end:
                                pts =[ (self.scrub(var_segment_start[0]), self.scrub(var_segment_start[1])  ,0) , 
                                       (self.scrub(var_segment_end[0])  , self.scrub(var_segment_end[1])    ,0) ]
                                poly = [(1,2)]
                                self.insert_polygons(poly, pts)   
                                
                                # reset for next line 
                                var_segment_start = None                              
                                var_segment_end = None  

                            #-------------------

                            # Line parsing  
                            if tok[1:] == 'start' and self.oneup_entity == 'gr_line':
                                #print("GR LINE start", toked[i+1], toked[i+2]  ) 
                                var_line_start_xy = [toked[i+1], toked[i+2]]
                            if tok[1:] == 'end' and self.oneup_entity == 'gr_line':
                                #print("GR LINE end", toked[i+1], toked[i+2]  ) 
                                var_line_end_xy   = [toked[i+1], toked[i+2]]
                            
                            if  var_line_start_xy and var_line_end_xy:
                                print("build a line from %s to %s"%(var_line_start_xy, var_line_end_xy) )
                                
                                pts =[ (self.scrub(var_line_start_xy[0]), self.scrub(var_line_start_xy[1]  ,0)) , 
                                       (self.scrub(var_line_end_xy[0])  , self.scrub(var_line_end_xy[1])   ,0) ]
                                poly = [(1,2)]
                                self.insert_polygons(poly, pts)   
                                
                                # reset for next line 
                                var_line_start_xy = None                              
                                var_line_end_xy = None  

                            #-------------------------------                                
                            self.parse_depth += 1

                        # if we encounter a ")" - go up one in the parse                            
                        for ct in re.findall('\)',tok):
                            

                            if self.oneup_entity is not None:
                                # if they are the same - we returned from an entity 
                                #if self.cur_entity != self.oneup_entity:
                                #    print("BACK OUT FROM %s INTO %s"% (self.cur_entity, self.oneup_entity ) )
                                
                                # go "oneup" in history      
                                self.cur_entity = self.oneup_entity 
                            
                            self.parse_depth -= 1
                  

        if self.parse_depth != 0:
            print("ERROR - UNEVEN NUMBER OF PARENTHESIS IN FILE! ")
            exit()            



            #self.file_contents.append(line)

    ##############
    #@property
    def scrub(self, inp):
        """ clean up parsed characters from kicad file """

        out = inp.lower()
        out = out.strip()
        out = out.replace(')','')
        out = out.replace(' ','')
        out = out.replace('-','')
        out = out.replace('\"','')
        out = out.replace('\'','')
        return out
    
    ##############    
    def do_shift(self, group):
        """ transform data operation after loading but before buidling 
            just for now flip Y negative, maybe more later  
        """
        #flip axis and normalize transforms
        #cmds.scale( 1, -1, 1, group, pivot=(0, 0, 0), absolute=True )
        #cmds.makeIdentity( group, apply=True, t=1, r=1, s=1, n=2 )

        #now perform a final scaling to ( mils from cm ?)
        #cmds.scale( self.GLOBAL_SCALE, self.GLOBAL_SCALE, self.GLOBAL_SCALE, group, pivot=(0, 0, 0), absolute=True )
        #cmds.makeIdentity( group, apply=True, t=1, r=1, s=1, n=2 )
        pass


    ## ##############
    ## def process(self):
    ##     """ 
    ##     depreciated - original parser code,  not really needed  
    ##     quick and dirty maya script to get kicad footprints into 3D 
    ##     see https://github.com/keithlegg/import_kicad_maya3d for more info 
    ##
    ##     take loaded text data and parse it looking for info to build pads and lines objects with """
    ##     #(module KL_ALPS_10KPOT
    ##     #part_grp = cmds.group( em=True, name='null1' )
    ##     for l in  self.file_contents:
    ##         linedata =  l.strip().split('(fp_line ')
    ##         sx = 0
    ##         sy = 0
    ##         ex = 0
    ##         ey = 0
    ##         padx = 0
    ##         pady = 0
    ##         dia_padx = 0 #use same value x and y 
    ##         ######################             
    ##         #build layers and assign names to them
    ##         #color like kicad 
    ##         ###################### 
    ##         #build pads and display as NURBS circles
    ##         if len(linedata)<=1:
    ##             parse_pads = linedata[0].split(' ')
    ##             if parse_pads[0]=='(pad':
    ##                 padx = self.scrub(parse_pads[5])
    ##                 pady = self.scrub(parse_pads[6])
    ##                 print( 'pad found at %s %s'%(padx,pady) )
    ##                 if parse_pads[7]=='(size':
    ##                     dia_padx = parse_pads[8]  
    ##                     #dia_pady = parse_pads[9]  
    ##         ###################### 
    ##         #import line segments as first degree NURBS curve segments
    ##         if len(linedata)>1:
    ##             vtxdata = linedata[1].split()
    ##             if vtxdata[0]=='(start':
    ##                 sx = self.scrub(vtxdata[1])
    ##                 sy = self.scrub(vtxdata[2])  
    ##             if vtxdata[3]=='(end':
    ##                 ex = self.scrub(vtxdata[4])
    ##                 ey = self.scrub(vtxdata[5]) 
    ##             #print( 'the line is %s %s %s %s '%(sx,sy,ex,ey))  
            


















