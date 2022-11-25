
import re 


#from pygfx.obj3d import object3d
from gnelscript.pygfx.milling_ops import gcode


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

        #self.points        = []  NOPE! - conflicts with 3D OBJECT  
        #self.polygons      = []  NOPE! - conflicts with 3D OBJECT

        self.parse_depth = 0
        self.file_contents = []
        
        self.oneup_entity = None 
        self.cur_entity = None 
        self.cur_module = None 
        self.cur_module_pos = None 

        self.modules = []

        self.known_entities = ['module','gr_line','segment']

        self.kicad_units   = 'mils' #inch to mm -> scale 0.0393701
        self.object_units  = 'cm'

        self.outfile = []
        
        # GEOM ARRAYS for export   
        self.ngc_to_obj    = []

        self.line_segments = []  #list of list of points 
        self.ki_polygons   = []  #list of list of points 
        self.filled_polys  = []  #list of list of points 
        self.gr_polys      = []  #list of list of points 
        self.modules       = []  #list of list of [name, [(),()] ] 

        self.gr_polygon_buffer   = [] # buffer to dump into array of arrays (gr)
        self.fill_polygon_buffer = [] # buffer to dump into array of arrays (filled_polys)
        self.polygon_buffer      = [] # buffer to dump into array of arrays (modules)

        self.segment_buffer = [] # buffer to dump into array of arrays (line_segments)
        self.module_buffer  = [] # buffer to dump into array of arrays (modules)
 
        self.parsing_fillpolygon = False # stored state of parser object type
        self.parsing_grpolygon   = False # stored state of parser object type 
        self.parsing_polygon     = False # stored state of parser object type 

        self.module_depth    = 0     # stored state of parser object type
        #self.parsing_segment = False # stored state of parser object type

        self.rh = 1.0     # retract height 
        self.ch = .5      # cut height 
        self.hp = (0,0,0) # home position 

        self.scale =  -0.0393701 #mm to inch

    ##############
    def show_geom(self):
        #call the inherited polygon "show" method
        self.show()

    ##############

    def bufferinfo(self):
        print('buffer sizes ')
        print('graphic polygons %s '%len(self.gr_polys     ) )
        print('filled  polygons %s '%len(self.filled_polys ) )
        print('        polygons %s '%len(self.ki_polygons     ) )                        

    ##############

    def showbuffers(self):
        print('graphic polygons %s '%self.gr_polys     )
        print('filled  polygons %s '%self.filled_polys )
        print('        polygons %s '%self.ki_polygons     ) 

    ##############

    #def insert_fill_poly(self, pts):

    ##############
    def export_ngc(self, filename):
        """ convert the raw points read from kicad into a usuable  path(s) with retract 
        """
        
        lastpt = (0,0,0)

        self.outfile = []
        
        self.outfile.append('(exported with gnelscript kicad_ops )')
        self.outfile.append('(linear scale set to %s of internal coordinates)'%self.scale )
        self.outfile.append('  ')

        self.outfile.append('g20')                  #inches for unit 
        
        #move to origin ? 
        self.outfile.append('g0 x%s y%s z%s f30'% ( self.hp[0], self.hp[1], self.rh) )   #rapid move to 0 
        self.ngc_to_obj.append( ( self.hp[0], self.hp[1], self.rh) ) 

        self.outfile.append('g0z%s'% ( self.rh ) )  #retract in between cuts
   
        ####
        self.outfile.append('  ')
        self.outfile.append('(exporting filled polygons )')

        # build the gcode up with simple linear movements 
        for fill_poly in self.filled_polys:

            pt1 = fill_poly[0] 
            self.outfile.append('x%s y%s z%s'% (  round(pt1[0]*self.scale,6) , round(pt1[1]*self.scale,6), self.rh ) )  #first point at retract height   
            self.ngc_to_obj.append( ( self.hp[0], self.hp[1], self.rh) ) 

            for i,pt in enumerate(fill_poly):
                self.outfile.append( 'x%s y%s z%s'%( round(pt[0]*self.scale,6) , round(pt[1]*self.scale,6), self.ch ) )
                self.ngc_to_obj.append( ( round(pt[0]*self.scale,6) , round(pt[1]*self.scale,6), self.ch ) )             
                lastpt =( round(pt[0]*self.scale,6) , round(pt[1]*self.scale,6), self.ch )


            self.outfile.append( 'x%s y%s z%s'%(lastpt[0], lastpt[1], self.ch) )
            self.ngc_to_obj.append( (lastpt[0], lastpt[1], self.ch)   ) 

            self.outfile.append('g0z%s'% ( self.rh ) )  #retract in between cuts
            self.ngc_to_obj.append( (lastpt[0], lastpt[1], self.rh)   ) 
            self.outfile.append('  ')

        ####
        self.outfile.append('  ')
        self.outfile.append('(exporting graphic polygons )')

        # build the gcode up with simple linear movements 
        for fill_poly in self.gr_polys:

            if len(fill_poly):
                pt1 = fill_poly[0] 

                self.outfile.append('x%s y%s z%s'% (  round(pt1[0]*self.scale,6) , round(pt1[1]*self.scale,6), self.rh ) )  #first point at retract height              
                self.ngc_to_obj.append( ( round(pt1[0]*self.scale,6) , round(pt1[1]*self.scale,6), self.rh ) )             
                
                for i,pt in enumerate(fill_poly):
                    self.outfile.append( 'x%s y%s z%s'%( round(pt[0]*self.scale,6) , round(pt[1]*self.scale,6), self.ch ) )
                    self.ngc_to_obj.append( ( round(pt[0]*self.scale,6) , round(pt[1]*self.scale,6), self.ch ) )                   
                    lastpt = ( round(pt1[0]*self.scale,6) , round(pt1[1]*self.scale,6), self.ch )

                self.ngc_to_obj.append( (round(fill_poly[0][0]*self.scale,6), round(fill_poly[0][1]*self.scale,6), self.ch)   ) 
                self.outfile.append( 'x%s y%s z%s'%( (round(fill_poly[0][0]*self.scale,6), round(fill_poly[0][1]*self.scale,6), self.ch) ) )

                self.ngc_to_obj.append( (round(fill_poly[0][0]*self.scale,6), round(fill_poly[0][1]*self.scale,6), self.rh)   )   
                self.outfile.append( 'x%s y%s z%s'%( (round(fill_poly[0][0]*self.scale,6), round(fill_poly[0][1]*self.scale,6), self.rh) ) )

                self.outfile.append('g0z%s'% ( self.rh ) )  #retract in between cuts

                self.outfile.append('  ')

        ####        
        # self.outfile.append('(exporting segments )')

        #### 
        # rapid move at end 
        self.outfile.append('m2') #program end

        #print('#####################')
        #print(self.ngc_to_obj)

        #########################


        fobj = open( filename,"w") #encoding='utf-8'
        for line in self.outfile: 
            fobj.write(line+'\n')
        fobj.close()

    ##############
    def show_modules(self):
        for m in self.modules:
            print(m.name)
    
    ##############
    def save_3d_obj(self, name):
        """ dump a bunch of 3d poonts as a giant line to visualize in gnolmec"""   
        self.points = self.ngc_to_obj

        # make a series of connected lines from points  
        for i,pt in enumerate(self.ngc_to_obj):
            if i>0 and i<len(self.ngc_to_obj):
                poly = (i,i+1)
                self.polygons.append( poly ) 

        
        #self.scale_pts(self.scale)
        self.save(name)


    def read_pcb(self,infile):
        #(gr_line (start 99.695 105.41) (end 119.38 149.86) (layer Dwgs.User) (width 0.15) (tstamp 6068B63B))
        #(gr_line (start 160.02 134.62) (end 99.695 105.41) (layer Dwgs.User) (width 0.15))
        #(gr_line (start 174.625 131.445) (end 160.02 134.62) (layer Dwgs.User) (width 0.15))
        #(gr_line (start 107.95 97.155) (end 174.625 131.445) (layer Dwgs.User) (width 0.15))
        #(gr_line (start 126.365 17.78) (end 107.95 97.155) (layer Dwgs.User) (width 0.15))
        #(gr_line (start 114.3 24.765) (end 126.365 17.78) (layer Dwgs.User) (width 0.15))
        #(gr_line (start 97.155 99.06) (end 114.3 24.765) (layer Dwgs.User) (width 0.15))
        #(gr_line (start 73.66 38.735) (end 97.155 99.06) (layer Dwgs.User) (width 0.15))
        #(gr_line (start 61.595 43.815) (end 73.66 38.735) (layer Dwgs.User) (width 0.15))
        #(gr_line (start 97.155 152.4) (end 61.595 43.815) (layer Dwgs.User) (width 0.15))
        #(gr_line (start 119.38 149.86) (end 97.155 152.4) (layer Dwgs.User) (width 0.15))
        #(via (at 139.7 76.835) (size 0.8) (drill 0.4) (layers F.Cu B.Cu) (net 0))
        #(via (at 153.035 88.9) (size 0.8) (drill 0.4) (layers F.Cu B.Cu) (net 0))
        #(via (at 149.225 99.695) (size 0.8) (drill 0.4) (layers F.Cu B.Cu) (net 0))
        pass


    ##############
    def load_kicadpcb(self, filename):
        """ a parser that is not recursive, but clever enough to scan all the 
            entities in the file that have coordinates and store then in a list with a type identifier 
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

        geometry = []

        for lc,line in enumerate(f):
            if '(' in line or ')' in line:
                newpoly =[]
                linetoked = line.split(' ') 
                cleanspaces = []
                for t in linetoked:
                    if t!='':
                        cleanspaces.append( self.scrub(t) )

                firstelem = self.scrub(cleanspaces[0]) 
                if firstelem!='':
                    newpoly.append(firstelem)

                for i,tok in enumerate(cleanspaces):
                    if tok !='':
                        if tok == 'xy':   
                            newpoly.append( (float(self.scrub(cleanspaces[i+1]))  , float(self.scrub(cleanspaces[i+2]) )) ) 

                        if tok == 'start':   
                            newpoly.append( (float(self.scrub(cleanspaces[i+1]))  , float(self.scrub(cleanspaces[i+2]) )) ) 
                        if tok == 'end':   
                            newpoly.append( (float(self.scrub(cleanspaces[i+1]))  , float(self.scrub(cleanspaces[i+2]) )) ) 

                        if tok == 'center':   
                            newpoly.append( (float(self.scrub(cleanspaces[i+1]))  , float(self.scrub(cleanspaces[i+2]) )) ) 

                if len(newpoly)>2:  
                    geometry.append(newpoly)
        ####################
        #not perfect - it will miss things that are longer than one line!
        
        print('######## num geometry elements read ', len(geometry ) )
        for p in geometry:
            tmp = []
            #print("### ", p )
            if len(p):
                if p[0] =='gr_poly':
                    tmp.extend(p[1:])
        
            self.gr_polys.append(tmp)
        

    ##############
    #@property
    def scrub(self, inp):
        """ clean up parsed characters from kicad file """

        out = inp.lower()
        out = out.strip()
        out = out.replace('(','')
        out = out.replace(')','')
        out = out.replace(' ','')
        out = out.replace('\"','')
        out = out.replace('\'','')
        return out

    ##############
    #@property
    def fl_scrub(self, inp):
        out = self.scrub(inp)
        return float(out)

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
            


















