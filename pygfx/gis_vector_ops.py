import geojson



from gnelscript.pygfx.grid_ops import  *

from gnelscript.pygfx.obj3d import  *




class generic_ngc(object3d):
    """copy of kicad parser for experimenting  """

    def __init__(self):
        super().__init__()  
        self.gr_polys      = []  #list of list of points 
        self.scale =  -0.0393701 #mm to inch
        
        # GEOM ARRAYS for export   
        self.ngc_to_obj    = []
        self.filled_polys  = []  #list of list of points 

        self.rh = .1     # retract height 
        self.ch = .5      # cut height 
        self.hp = (0,0,0) # home position 

    ##-------------------------------##
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
    def save_3d_obj(self, name):
        """ dump a bunch of 3d pionts as a giant line to visualize in gnolmec"""   
        
        export_retracts = False
        
        print("exporting %s polygons "%len(self.gr_polys))


        if export_retracts == False:
            #this will make a pretty OBJ without retracts 
            
            idx = 1  
            for jpoly in self.gr_polys:
                self.points.extend(jpoly)

                for i,pt in enumerate(jpoly):
                    if i>0 and i<len(jpoly):
                        poly = (i,i+1)
                        self.polygons.append( poly ) 


        if export_retracts == True: 
            # make a series of connected lines from points 
            self.points = self.ngc_to_obj

            for i,pt in enumerate(self.ngc_to_obj):
                if i>0 and i<len(self.ngc_to_obj):
                    poly = (i,i+1)
                    self.polygons.append( poly ) 


        #self.scale_pts(self.scale)
        self.save(name)
        
    ##-------------------------------##
    def load_geojson(self, inputfile, zaxis, getone=None):
        """ parse a geojson file - store points in arrays in GR_POLY buffer """

        plyidx = 1
        ptidx = 1 

        geojson_txt = None

        with open(inputfile) as f:
            gj = geojson.load(f)
        features = gj['features'] 
 
        if getone:
            if getone>len(features):
                return "getone out of range "

        
        for i,f in enumerate(features):

            for coord in f.geometry.coordinates:
                ptidx = 1
                
                tmp_poly = [] 

                #no getone is all features 
                if getone ==None:
                    for pt in coord:
                        if type(pt[0])==float and type(pt[1])==float:
                            tmp_poly.append( (pt[0], pt[1], zaxis) )
                            ptidx += 1  
                #getone is single for now - specify a feature ID to pull 
                if getone:
                    if i == getone: 
                        for pt in coord:
                            if type(pt[0])==float and type(pt[1])==float:
                                tmp_poly.append( (pt[0], pt[1], zaxis) )
                                ptidx += 1                      

                #print("loaded %s points in polygon "%len(tmp_poly)) 
                self.gr_polys.append(tmp_poly) 
            plyidx += 1 
        
        print("loaded %s polygons from %s "%(plyidx,inputfile)) 

    ##-------------------------------##
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

    ##-------------------------------##
    def calulate_paths(self):
        """ 
            DEBUG - it seems that this just makes a single line with retract caclulated  
            convert the raw points read from kicad into a usuable path(s) with retract 
        """
        
        do_retract = False 

        lastpt = (0,0,0)

        self.outfile = []
        
        self.outfile.append('(exported with gnelscript kicad_ops )')
        self.outfile.append('(linear scale set to %s of internal coordinates)'%self.scale )
        self.outfile.append('  ')

        self.outfile.append('g20')                  #inches for unit 
        
        ##-----------------------------------------##

        #move to origin ? 
        self.outfile.append('g0 x%s y%s z%s f30'% ( self.hp[0], self.hp[1], self.rh) )   #rapid move to 0 
        self.ngc_to_obj.append( ( self.hp[0], self.hp[1], self.rh) ) 

        self.outfile.append('g0z%s'% ( self.rh ) )  #retract in between cuts
   
        ##-----------------------------------------##
        self.outfile.append('  ')
        self.outfile.append('(exporting filled polygons )')

        ##-----------------------------------------##
        # build the gcode up with simple linear movements 
        for fill_poly in self.filled_polys:

            if do_retract:
                pt1 = fill_poly[0] 
                self.outfile.append('x%s y%s z%s'% (  round(pt1[0]*self.scale,6) , round(pt1[1]*self.scale,6), self.rh ) )  #first point at retract height   
                self.ngc_to_obj.append( ( self.hp[0], self.hp[1], self.rh) ) 

            for i,pt in enumerate(fill_poly):
                self.outfile.append( 'x%s y%s z%s'%( round(pt[0]*self.scale,6) , round(pt[1]*self.scale,6), self.ch ) )
                self.ngc_to_obj.append( ( round(pt[0]*self.scale,6) , round(pt[1]*self.scale,6), self.ch ) )             
                lastpt =( round(pt[0]*self.scale,6) , round(pt[1]*self.scale,6), self.ch )


            self.outfile.append( 'x%s y%s z%s'%(lastpt[0], lastpt[1], self.ch) )
            self.ngc_to_obj.append( (lastpt[0], lastpt[1], self.ch)   ) 

            if do_retract:
                self.outfile.append('g0z%s'% ( self.rh ) )  #retract in between cuts
                self.ngc_to_obj.append( (lastpt[0], lastpt[1], self.rh)   ) 
                self.outfile.append('  ')

        ####
        self.outfile.append('  ')
        self.outfile.append('(exporting graphic polygons )')


        ##-----------------------------------------##
        # graphical polygons - build the gcode up with simple linear movements
        for gr_poly in self.gr_polys:

            if len(gr_poly):
                pt1 = gr_poly[0] 

                #### first point at retract height 
                if do_retract:
                    #move to first point RH 
                    self.outfile.append('x%s y%s z%s'% (  round(pt1[0]*self.scale,6) , round(pt1[1]*self.scale,6), self.rh ) )               
                    self.ngc_to_obj.append( ( round(pt1[0]*self.scale,6) , round(pt1[1]*self.scale,6), self.rh ) )             
                
                # iterate points in polygon 
                for i,pt in enumerate(gr_poly):
                    self.outfile.append( 'x%s y%s z%s'%( round(pt[0]*self.scale,6) , round(pt[1]*self.scale,6), self.ch ) )
                    self.ngc_to_obj.append( ( round(pt[0]*self.scale,6) , round(pt[1]*self.scale,6), self.ch ) )                   
                    #lastpt = ( round(pt1[0]*self.scale,6) , round(pt1[1]*self.scale,6), self.ch )

                # move to last point at CH  
                #self.ngc_to_obj.append( (round(gr_poly[0][0]*self.scale,6), round(gr_poly[0][1]*self.scale,6), self.ch)   ) 
                #self.outfile.append( 'x%s y%s z%s'%( (round(gr_poly[0][0]*self.scale,6), round(gr_poly[0][1]*self.scale,6), self.ch) ) )

                if do_retract:
                    self.ngc_to_obj.append( (round(gr_poly[0][0]*self.scale,6), round(gr_poly[0][1]*self.scale,6), self.rh)   )   
                    self.outfile.append( 'x%s y%s z%s'%( (round(gr_poly[0][0]*self.scale,6), round(gr_poly[0][1]*self.scale,6), self.rh) ) )

                    #### retract in between cuts
                    self.outfile.append('g0z%s'% ( self.rh ) )  

                self.outfile.append('  ')

        ##-----------------------------------------##        
        # self.outfile.append('(exporting segments )')

        ##-----------------------------------------##
        # rapid move at end 
        self.outfile.append('m2') #program end


    ##-------------------------------##
    def export_ngc(self, filename):
 
        #build the paths - follow all the paths and add retract in (as a simple test to get started) 
        self.calulate_paths() 


        fobj = open( filename,"w") #encoding='utf-8'
        for line in self.outfile: 
            fobj.write(line+'\n')
        fobj.close()


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

