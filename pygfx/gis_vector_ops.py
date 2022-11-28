import geojson



from gnelscript.pygfx.grid_ops import  *
from gnelscript.pygfx.obj3d import  *





class generic_ngc(object3d):
    """copy of kicad parser for experimenting  """

    def __init__(self):
        super().__init__()  
        self.loadbuffer    = []  #list of list of points 
        self.gr_polys      = []  #list of list of points 
        
        # GEOM ARRAYS for export   
        self.ngc_to_obj    = []
        self.filled_polys  = []  #list of list of points 


        self.global_scale =  0.0393701 #inch to mm 

        self.rh = 1.5     # retract height 
        self.ch = 1       # cut height 
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

    ## def insert_pts(self, pts):
    ##     """ insert points into gr_polys 
    ##         fix rounding errors 
    ##         fix anything else .... 
    ##     """
    ##     self.gr_polys.append(clean)

    ##-------------------------------##
    def save_line_obj(self, name):
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
    def cvt_grpoly_obj3d(self, index=None):
        """ 
            load gr_polys into standard 3d data so we can run ops on it 
        """   
        
        print("converting %s polygons "%len(self.gr_polys))

        idx = 0  
        lidx = 0 #last indexd 
        
        for ig, grpoly in enumerate(self.gr_polys):
            
            # if index 

            if index is None:
                polygons = []
                for ip, pt in enumerate(grpoly):
                    idx = ip+lidx
                    polygons.append( (idx, idx+1) ) #2 point polygon indeces

                
                self.polygons.append(polygons)
                
                self.points.extend(self.clean_points(grpoly))
                lidx = len(grpoly)

                #print(grpoly) 


        #print(self.rotation)
        self.show() 


    ##-------------------------------##
    def grply_inspect(self, index=None):
        """
            #return info about how many features, sub features, etc are in a file 
            
            args:
                index is positive int or iterable 

            todo:
                slice , index 
                get extents 
                get centroid 
                get as vectors 

        """
    
        print("##------------##")

        if index == None:
            print("gr_poly buffer has %s polygons "%len(self.gr_polys) )
            for i,p in enumerate(self.gr_polys):
                print("    feature %s has %s points"%(i, len(p)))
                print(p)

        else:
            pass 



    ##-------------------------------##
    def load_geojson(self, inputfile, zaxis, getfids=None, getids=None):
        """ parse a geojson file - store points in arrays in GR_POLY buffer """

        plyidx = 1
        ptidx = 1 

        geojson_txt = None

        with open(inputfile) as f:
            gj = geojson.load(f)
        features = gj['features'] 
 
        fixedids = [] 

        ##---------------

        #check that ids are in range  
        if getids:
            print("#load_geojson getone mode ", getids) 
            for id in getids:
                if id>=len(features):
                    print("id out of range %s"%id)
                else:
                    fixedids.append(id)
        #check that ids are in range  
        if getfids:
            print("#load_geojson getone mode ", getids) 
            for id in getids:
                if id>=len(features):
                    print("id out of range %s"%id)
                else:
                    fixedids.append(id)

        ##---------------
        # export all 
        if getfids is None:
            for i,f in enumerate(features):
                for coord in f.geometry.coordinates:
                    ptidx = 1
                    tmp_poly = [] 
                    for pt in coord:
                        if type(pt[0])==float and type(pt[1])==float:
                            tmp_poly.append( (pt[0], pt[1], zaxis) )
                            ptidx += 1  
                    #print("loaded %s points in polygon "%len(tmp_poly)) 
                    if tmp_poly:
                        self.loadbuffer.append(tmp_poly) 
                plyidx += 1 

        ##---------------        
        # cherry pick the features to export (below we also can pick the sub-features) 
        if getfids:
            for fid in getfids:
                for i,f in enumerate(features):
                    if fid==i:
                        for coord in f.geometry.coordinates:
                            ptidx = 1
                            tmp_poly = [] 
                            for pt in coord:
                                if type(pt[0])==float and type(pt[1])==float:
                                    tmp_poly.append( (pt[0], pt[1], zaxis) )
                                    ptidx += 1  
                            #print("loaded %s points in polygon "%len(tmp_poly)) 
                            if tmp_poly:
                                self.loadbuffer.append(tmp_poly) 
                        plyidx += 1 

        ##---------------
        #optional processing to get sub features by id
        
        # get all sub features
        if getids is None:
            self.gr_polys = self.loadbuffer 
        
        # pick out some sub features
        else:
            print("ids to load are is %s"%fixedids)
            for id in fixedids: 
                self.gr_polys.append(self.loadbuffer[id]) 

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
    def calculate_paths(self, do_retracts=True):
        """ 
            DEBUG - it seems that this just makes a single line with retract caclulated  
            convert the raw points read from kicad into a usuable path(s) with retract 

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

        self.outfile = []
        
        self.outfile.append('(exported with gnelscript kicad_ops )')
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

        ###################################################

        #DEBUG - FILL_POLYS ARE A LAYOVER FROM KICAD STUFF - UNTESTED 
        for fill_poly in self.filled_polys:
            if do_retracts:
                pt1 = fill_poly[0] 
                self.outfile.append('G0 x%s y%s z%s'% (  pt1[0], pt1[1], self.rh ) )  #first point at retract height   
                self.ngc_to_obj.append( ( self.hp[0], self.hp[1], self.rh) ) 
            
            for i,pt in enumerate(fill_poly):
                self.outfile.append( 'G1 x%s y%s z%s'%( pt[0], pt[1], self.ch ) )
                self.ngc_to_obj.append( ( pt[0], pt[1], self.ch ) )             
                lastpt =( round(pt[0]*self.global_scale,pl) , pt[1], self.ch )

            self.outfile.append( 'G0 x%s y%s z%s'%(lastpt[0], lastpt[1], self.ch) )
            self.ngc_to_obj.append( (lastpt[0], lastpt[1], self.ch)   ) 
            if do_retracts:
                self.outfile.append('g0z%s'% ( self.rh ) )  #retract in between cuts
                self.ngc_to_obj.append( (lastpt[0], lastpt[1], self.rh)   ) 
                self.outfile.append('  ')

        ####
        self.outfile.append('  ')
        self.outfile.append('(exporting graphic polygons )')
        
        ####################################################

        # graphical polygons - build the gcode up with simple linear movements
        for gr_poly in self.gr_polys:
            self.outfile.append('(exporting new polygon )')
            if len(gr_poly):
            
                
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
                    #lastpt = ( round(pt1[0]*self.global_scale,6) , round(pt1[1]*self.global_scale,6), self.ch )
                
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
    def export_ngc(self, filename):

        fobj = open( filename,"w") #encoding='utf-8'
        for line in self.outfile: 
            fobj.write(line+'\n')
        fobj.close()


    ##-------------------------------##
    def import_ngc(self, filename):
        """ DEBUG - NOT WORKING OR TESTED """
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

