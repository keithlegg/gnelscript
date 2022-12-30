import geojson
from geojson import Point, Feature, LineString, FeatureCollection, dump


from gnelscript.pygfx.math_ops import  *
from gnelscript.pygfx.grid_ops import  *

from gnelscript.pygfx.raster_ops import  *

from gnelscript.pygfx.render import simple_render

from gnelscript.pygfx.obj3d import  *
from gnelscript.pygfx.obj2d import  *





#########################################################################
#########################################################################



class generic_ngc(object3d):
    """copy of kicad parser for experimenting  """

    def __init__(self):
        super().__init__()  
        self.mu          = math_util()
        self.tesl        = teselator()
        self.pop2d       = object2d() 

        self.loadbuffer    = []  #list of list of points 
        self.gr_polys      = []  #list of list of points 

        self.gr_sort       = []  #list of [(pt), (pts)]  
        self.tmp_geom      = []    

        self.numexported = 0 

        # GEOM ARRAYS for export   
        self.ngc_to_obj    = []
        self.filled_polys  = []  #list of list of points 


        self.global_scale =  0.0393701 #inch to mm 

        self.rh = 1.5     # retract height 
        self.ch = 1       # cut height 
        self.hp = (0,0,0) # home position 

        self.total_minx = 0
        self.total_miny = 0         
        self.total_maxx = 0
        self.total_maxy = 0

    ##-------------------------------##
    # IDEA 

    #1 - DAG - b tree 
    #2 - parent geom to nodes - add matrix scenegraph  
    #3-  crazy f-ing animation potential - vector ala cyriak  

    ##-------------------------------##
    #TODO: 

    ## def import_polys(self, pts):

    ## def translate(self, pts):
    ##     global/ indexer ?
    ## def rotate(self, pts):
    ##     global/ indexer ?
    ## def scale(self, pts):
    ##     global/ indexer ?

    ## def remove_polys(self, pts):
    ##     indexer 
    ## def export_polys(self, pts):
    ##     indexer  


    ##-------------------------------##    
    def _set_extents(self, bbox):
        """ set global extents for generating data 
            based on PIL coordinate which is [left, top, right, bottom] 
            [minx, miny, maxx, maxy]  
        """
        
        self.total_minx = bbox[0]
        self.total_miny = bbox[1]         
        self.total_maxx = bbox[2]
        self.total_maxy = bbox[3]

    ##-------------------------------##
    def make_grid(self, folder, xcuts, ycuts, bbox=None):
        """ chop a square into smaller squares """

        if bbox:
            self.tesl._set_extents(bbox) 
        else:
            self.tesl._set_extents([self.total_minx, self.total_miny, self.total_maxx, self.total_maxy]) 

        self.tesl.build_2d_cells(xcuts, ycuts)

        
        #temporary export of geom I HAVE TO SEE THIS BEFORE I GO TO BED!
        features = []
        
        DIST_THRESH = 5.0 

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
                    if dist < DIST_THRESH:
                        # DEBUG CLEAN THIS UP 
                        #print(c.name, sort[0], dist)
                        #self.tmp_geom.append( [] )
                        
                        #DEBUG TODO - there is a problem of the same polygon getting added more than once 
                        #make sure we dont export multiple times 
                        # (check the number of points, then check thr actual points) 

                        tmp_pts=[cen, sort[1]]
                        features.append(Feature(geometry=LineString(coordinates=tmp_pts)) 
                                       )

        feature_collection = FeatureCollection(features)
        with open('%s/_distances.json'%(folder), 'w') as f:
            dump(feature_collection, f)

    ##-------------------------------##
    def mc_escher(self):
        """ DAG + point generator = tesselator 
            make a grid and draw points at each location 
        """
        for cell in self.tesl.nodes:
            cen = cell.getattrib('centroid')
            #print('# ', cell.name,' ', cen )
            
            pts = self.pop2d.calc_circle_2d(cen[0],cen[1], .75, periodic=True, spokes=3)
            cell.points.append( pts )

            pts2 = self.pop2d.batch_rotate_pts_2d( pts, cen, 180 )
            cell.points.append( pts2 )

            pts3 = self.pop2d.calc_circle_2d(cen[0],cen[1], 1.5, periodic=True, spokes=12)
            cell.points.append( pts3 )

    ##-------------------------------##
    def tess_vec_render(self, renderscale, path, objfile):
        """ make a grid and invoke render at each location 
            
            DEBUG - add mesh optimizer 
                - eliminate lines that overlap 
                - only render faces that face you (backface culling) 
                - eliminate lines too small/big 


        """

        # render as one line (with no breaks)
        single_line = False 
        
        ropr = simple_render() 
        obj = object3d()
        obj.load('%s/%s'%(path,objfile))

        #obj.show()

        total = len(self.tesl.nodes)

        for i,cell in enumerate(self.tesl.nodes):
            cen = cell.getattrib('centroid')
            print('#rendering  %s@%s %s of %s '%(cell.name, cen, i,total) )

            #clear cache each time 
            ropr.vec_fr_buffer.points = [] 
            ropr.vec_fr_buffer.polygons = [] 

            ## color, rx, ry, rz, thick, scale 
            ropr.render_obj((100,0,255), i*20, i*20, i*20,  1, renderscale, object3d=obj)

            #coords are in pixels - rather huge for a model 
            #ropr.vec_fr_buffer.scale_pts((.01,.01,.01))
        
            ropr.vec_fr_buffer.move_center()
          
            # used to test the vector render engine 
            #ropr.vec_fr_buffer.save('%s/test_%s.obj'%(path,i) )
    

            # draw as a single line without breaks       
            if single_line:
                pts = [] 
                for pt in ropr.vec_fr_buffer.points:
                    pts.append( (pt[0]+cen[0], pt[1]+cen[1] ) )  
                cell.points.append( pts )

            # draw as proper line segments         
            else:
                for ply in ropr.vec_fr_buffer.polygons:
                    pts = []
                    #in this setup it will only be two point polys - vector render only does that                    
                    if len(ply)==2:
                        # pts are zero indexed hence the -1
                        pt1 = ropr.vec_fr_buffer.points[ply[0]-1] 
                        pt2 = ropr.vec_fr_buffer.points[ply[1]-1] 
                        pts.append( (pt1[0]+cen[0], pt1[1]+cen[1]) )
                        pts.append( (pt2[0]+cen[0], pt2[1]+cen[1]) )
                    cell.points.append( pts )


    ##-------------------------------##
    def tess_objclone(self, objfile):
        """ make a grid and insert points frm an object at each location 
            make sure to center the object for best results
            use object3D.move_center()  
        """

        self.pop2d.load(objfile)
        self.pop2d.show()

        for cell in self.tesl.nodes:
            cen = cell.getattrib('centroid')
            #print('# ', cell.name,' ', cen )
        
            #DEBUG move goes into a broken loop if you try to use it - debug  
            #self.pop2d.move( cen[0],cen[1] )

            pts = []

            for pt in self.pop2d.points:
                pts.append( (pt[0]+cen[0], pt[1]+cen[1] ) )  

            cell.points.extend( pts )

    ##-------------------------------##
    def sample_data(self):
        """ DAG + point generator = tesselator """
        for cell in self.tesl.nodes:
            cen = cell.getattrib('centroid')
            #print('# ', cell.name,' ', cen )
            
            pts = self.pop2d.calc_circle_2d(cen[0],cen[1], .75, periodic=True, spokes=6)
            cell.points.extend( pts )

            pts2 = self.pop2d.batch_rotate_pts_2d( pts, cen, 180 )
            cell.points.extend( pts2 )




    ##-------------------------------##
    def show_buffers(self):

        #[[id, centroid, extents, len, points ]
        #for ply in self.gr_sort:
        #    print("%s %s %s "%(ply[0], ply[1], len(ply[2]) )) 
        print('#################')
        print('size gr_polys %s'%len(self.gr_polys))
        print('size gr_sort %s'%len(self.gr_sort))
    
    ##-------------------------------##
    def export_grid_gfx(self, name, folder ):
        #export cells as graphics  

        if self.total_minx==0 and self.total_miny==0 and self.total_maxx==0 and self.total_maxy==0:
            print("WARNING EXTENTS ARE ALL ZERO ")

        features = []
        for c in self.tesl.nodes:

            features.append(Feature(geometry=Point((c.coord_x+(c.width/2), c.coord_y+(c.height/2))), 
                                    properties={"id":c.name} 
                                   )
                        )
            #draw a square around the boundary of object 
            features.append(Feature(geometry=LineString(coordinates=c.boundary_pts), 
                          properties={"id" : 0 
                                     }
                          ) 
                  )
            #render multilines in points (nested arrays - broken up lines)            
            for ptgrp in c.points:
                features.append(Feature(geometry=LineString(coordinates=ptgrp), 
                              properties={"id" : 0 
                                         }
                              ) 
                      )

        feature_collection = FeatureCollection(features)
        with open('%s/%s_cells.json'%(folder, name), 'w') as f:
            dump(feature_collection, f)

    ##-------------------------------##
    def export_sorted_centroids(self, name, folder):
        """ export a centroid for each polygon as a geojson file """

        features = []

        #[[id, centroid, extents, len, points ]]
        for i,s in enumerate(self.gr_sort):
            features.append(Feature(geometry=Point((s[1][0],s[1][1])), properties={"id":i, "len" : len(s[4]) } ))

        feature_collection = FeatureCollection(features)
        with open('%s/%s_ply_cntrs.json'%(folder,name), 'w') as f:
            dump(feature_collection, f)

    ##-------------------------------##
    def export_sorted_extents(self, name, folder):
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

        #total extents  
        coords = self.extents_fr_bbox([self.total_minx,self.total_miny,self.total_maxx,self.total_maxy], periodic=True)
        features.append(Feature(geometry=LineString(coordinates=coords), 
                                properties={"id" : 0 
                                           }
                                ) 
                        )
                        

        feature_collection = FeatureCollection(features)
        with open('%s/%s_ply_xtntx.json'%(folder,name), 'w') as f:
            dump(feature_collection, f)

    ##-------------------------------##       
    def index_grsort(self):
        """ assemble data into [[id, centroid, extents, len, points ]] - put that in self.gr_sort   

            set total extents of all data while running  
        """

        #print("indexing sort buffer ")
        self.gr_sort   = []

        #print(" gr_polys buffer has %s polys in it "%len(self.gr_polys) )
        for ply in self.gr_polys:
            minx = 0
            miny = 0         
            maxx = 0
            maxy = 0

            #print('### len ', len(ply))
            for i,pt in enumerate(ply):
                if i == 0:
                    minx=pt[0]
                    maxx=pt[0]
                    miny=pt[1]
                    maxy=pt[1]

                if pt[0]<minx:
                    minx=pt[0]    
                if pt[0]>maxx:
                    maxx=pt[0]  
                if pt[1]<miny:
                    miny=pt[1]  
                if pt[1]>maxy:
                    maxy=pt[1] 
            
            if minx<self.total_minx:
                self.total_minx=minx
            if miny<self.total_miny:
                self.total_miny=miny
            if maxx>self.total_maxx:
                self.total_maxx=maxx
            if maxy>self.total_maxy:
                self.total_maxy=maxy

            self.gr_sort.append([i, self.centroid(ply) ,[minx,miny, maxx, maxy], len(ply), ply])

            #print("poly extents %s %s %s %s "%(minx, maxx, miny, maxy) )
            
        print("total extents %s %s %s %s "%(self.total_minx, self.total_maxx, self.total_miny, self.total_maxy) )

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
        
        print('## ## num geometry elements read ', len(geometry ) )
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

        print("gr poly buffer has %s polys in it. "%(len(self.gr_polys)))
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
        
        print('## ## num geometry elements read ', len(geometry ) )
        for p in geometry:
            tmp = []
            #print("### ", p )
            if len(p):
                if p[0] =='gr_poly':
                    tmp.extend(p[1:])
        
            self.gr_polys.append(tmp)

