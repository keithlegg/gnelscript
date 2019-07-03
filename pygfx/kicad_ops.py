
import re 


#from pygfx.obj3d import object3d
from pygfx.milling_ops import gcode_to_polyline




class import_footprint(gcode_to_polyline):

    def __init__(self):
        super().__init__()  

        self.parse_depth   = 0
        self.file_contents = []
        self.entities      = []

        self.known_entities = ['module','gr_line','segment']

        self.kicad_units = 'mils'
        self.object_units  = 'cm'

        #self.GLOBAL_SCALE = 0.3904 # unit conversion (cm to mills?)

    ##############
    def load_kicadpcb(self, filename):
        print('reading file %s'%filename)
        f = open(filename, 'r')
        for line in f:
            if '(' in line or ')' in line:
                #self.known_entities
                toked = line.split(' ') 
                for tok in toked:
                    if tok !='':
                      
                        for ct in re.findall('\(',tok):
                            if tok[1:] in self.known_entities:
                                print("### found one! ", tok )
                            self.parse_depth += 1
                        for ct in re.findall('\)',tok):
                            self.parse_depth -= 1
                  

        if self.parse_depth != 0:
            print("ERROR - UNEVEN NUMBER OF PARENTHESIS IN FILE! ")
            exit()            



            #self.file_contents.append(line)

    ##############
    #@property
    def scrub(self, inp):
        """ clean up parsed characters from kicad file """

        out = inp
        out = out.strip()
        out = out.replace(')','')
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


    ##############
    def process(self):
        """ take loaded text data and parse it looking for info to build pads and lines objects with """

        #(module KL_ALPS_10KPOT
        #part_grp = cmds.group( em=True, name='null1' )

        for l in  self.file_contents:
            linedata =  l.strip().split('(fp_line ')

            sx = 0
            sy = 0
            ex = 0
            ey = 0

            padx = 0
            pady = 0
            dia_padx = 0 #use same value x and y 

            ######################             
            #build layers and assign names to them
            #color like kicad 

            ###################### 
            #build pads and display as NURBS circles
            if len(linedata)<=1:
                parse_pads = linedata[0].split(' ')
                if parse_pads[0]=='(pad':
                    padx = self.scrub(parse_pads[5])
                    pady = self.scrub(parse_pads[6])

                    print( 'pad found at %s %s'%(padx,pady) )
                    
                    if parse_pads[7]=='(size':
                        dia_padx = parse_pads[8]  
                        #dia_pady = parse_pads[9]  

            ###################### 
            #import line segments as first degree NURBS curve segments
            if len(linedata)>1:
                vtxdata = linedata[1].split()
                
                #print linedata

                if vtxdata[0]=='(start':
                    sx = self.scrub(vtxdata[1])
                    sy = self.scrub(vtxdata[2])  
                if vtxdata[3]=='(end':
                    ex = self.scrub(vtxdata[4])
                    ey = self.scrub(vtxdata[5]) 
                #print( 'the line is %s %s %s %s '%(sx,sy,ex,ey))  
            


















