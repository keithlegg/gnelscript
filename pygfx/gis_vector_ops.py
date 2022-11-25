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

        self.rh = 1.0     # retract height 
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
        """ dump a bunch of 3d poonts as a giant line to visualize in gnolmec"""   
        self.points = self.ngc_to_obj

        # make a series of connected lines from points  
        for i,pt in enumerate(self.ngc_to_obj):
            if i>0 and i<len(self.ngc_to_obj):
                poly = (i,i+1)
                self.polygons.append( poly ) 

        
        #self.scale_pts(self.scale)
        self.save(name)
        
    ##-------------------------------##
    def load_geojson(self, inputfile):

        plyidx = 1
        ptidx = 1 

        geojson_txt = None

        with open(inputfile) as f:
            gj = geojson.load(f)
        features = gj['features'] 

        zaxis = 0

        #obj = object3d()

        for f in features:
            geom  = [[],[]]

            for coord in f.geometry.coordinates:
                ptidx = 1
                
                tmp_poly = [] 

                for pt in coord:
                    if type(pt[0])==float and type(pt[1])==float:

                        print(pt[0], pt[1], ptidx, plyidx)
                        
                        tmp_poly.append( (pt[0], pt[1], zaxis) )


                        # #add new geom and auto increment the ids
                        #polys = [(1,2,3), (2,3,4) ]
                        #pts = [(1,1,1),(0,1,1),(-1,-1,1),(2,-2,1)]
                        #geom = obj.insert_polygons(polys, pts, geom=geom) 

                        ptidx += 1  
                
                self.gr_polys.append(tmp_poly) 
                # use insert to add geom to object 
                #obj.insert(geom) 


            plyidx += 1 


        # see what we have done, or not done 
        #obj.show() 
        #obj.save()

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
