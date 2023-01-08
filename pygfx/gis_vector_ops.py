import geojson
from geojson import Point, Feature, LineString, FeatureCollection, dump


from gnelscript.pygfx.math_ops import  *
from gnelscript.pygfx.grid_ops import  *
from gnelscript.pygfx.raster_ops import  *

from gnelscript.pygfx.milling_ops import  *
#from gnelscript.pygfx.kicad_ops import  *


from gnelscript.pygfx.obj3d import  *
from gnelscript.pygfx.obj2d import  *


from gnelscript.pygfx.render import simple_render






class vectorflow(object3d):
    """ copy of kicad parser for experimenting  
        not really GIS related but it DOES export geoJSON 

        It started out as s kicad ciles (graphic polygons and ???) to GCODE tool.
        It turned into a tool turn geojson into GCODE  
        then the dag_ops tesselator was added for GCODE optimization and I got derailed making MC escher style vector renders  


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

    """

    def __init__(self):
        super().__init__()  
        
        self.outfile = []

        self.mu          = math_util()
        self.tesl        = tessellator()
        self.pop2d       = object2d() 
        #self.kiparser    = pcbfile()


        # geometry buffers for JSON, NGC,sorting, processing, etc 
        self.gr_polys      = [] # list of list of points 
        self.gr_sort       = [] # [[id, centroid, extents, len, points ]]  
        self.ngc_buffer    = [] # list of list of points 
        self.jsonbuffer    = [] # list of list of points 

        # GEOM ARRAYS for export   
        self.ngc_to_obj    = [] # text buffer for obj 
        self.filled_polys  = [] # list of list of points 

        self.omit_ids      = [] #list of feature ids to NOT export (but leave in)

        self.rh = 2            # retract height 
        self.ch = .5           # cut height (top, start of cut)
        self.cdpi = .01        # cut depth per iteration on Z axis
        self.cmax = 1          # maximum cut depth on Z axis 

        self.hp = (0,0,0)              # home position 

        self.orig_minx = 0
        self.orig_miny = 0         
        self.orig_maxx = 0
        self.orig_maxy = 0

        self.sort_minx = 0
        self.sort_miny = 0         
        self.sort_maxx = 0
        self.sort_maxy = 0

        self.global_scale =  0.0393701 #NOT FULLY IMPLEMENTED - inch to mm 

    ##-------------------------------## 
    ##-------------------------------##   
    def _set_cam_properties(self, rh, ch, cdpi, cmax):
        self.rh = rh           
        self.ch = ch            
        self.cdpi = cdpi       
        self.cmax = cmax          


    def _set_extents(self, bbox):
        """ set global extents for generating data 
            based on PIL coordinate which is [left, top, right, bottom] 
            [minx, miny, maxx, maxy]  
        """
        
        self.sort_minx = bbox[0]
        self.sort_miny = bbox[1]         
        self.sort_maxx = bbox[2]
        self.sort_maxy = bbox[3]
    ##-------------------------------##       
    def _omit_ids(self, ids=None, span=None, unique=True, nth=None):
        
        #DEBUG - not fully working - see indexer docs - nths with crash if you try it 
        pop = point_operator()
        ids = pop.indexer(ids, span, unique, nth)
        
        print(" ### omit ids: ", ids)

        self.omit_ids = ids

    ##-------------------------------##       
    def _make_periodic(self):
        """ DEBUG - SHOULD OPERATE ON GR_SORT 
            if data has polygon that are not closed - add the first point to the end to close 
        """
        for ply in self.gr_polys:
            first = ply[0]
            ply.append(first)

    ##-------------------------------## 

    def _make_outline(self, iterations):
        """ 
            DEBUG - NEED TO ADD 1/2 DIA OF CUTTING TOOL 
            WE NEED DIALATE AND ERODE TO MAKE THIS WORK 
            first attempt at a Z operation - spiral down a polygon 

        """
        
        #spirals = [] 

        for idx, sort in enumerate(self.gr_sort):

            newply = []            
            ply = sort[4]

            for i, c in enumerate(range(iterations)):
                depth = self.ch+(i*self.cdpi)
                for pt in ply:
                    if depth < self.cmax:
                        newply.append( (pt[0] ,pt[1], depth) )
                    if depth>self.cmax:
                        print("cut too deep ")

            #spirals.append(newply)
            self.gr_sort[idx][4] = newply  

    ##-------------------------------##       
    def _sort(self):
        """ assemble data into [[id, centroid, extents, len, points ]] - put that in self.gr_sort   

            set extents of original data while running   
        """

        #print("indexing sort buffer ")
        self.gr_sort   = []

        #print(" gr_polys buffer has %s polys in it "%len(self.gr_polys) )
        for x,ply in enumerate(self.gr_polys):
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
            
            if minx<self.orig_minx:
                self.orig_minx=minx
            if miny<self.orig_miny:
                self.orig_miny=miny
            if maxx>self.orig_maxx:
                self.orig_maxx=maxx
            if maxy>self.orig_maxy:
                self.orig_maxy=maxy

            self.gr_sort.append([x, self.centroid(ply) ,[minx,miny, maxx, maxy], len(ply), ply])
        

        # run initial extents calc on sorted polys         
        self.gl_extents()

        print("orig data extents %s %s %s %s "%(self.orig_minx, self.orig_maxx, self.orig_miny, self.orig_maxy) )

    ##-------------------------------------------## 
    ##-------------------------------------------## 
    def show_setup(self):
        print("retract height %s"%self.rh)
        print("cut height     %s"%self.ch)
        print("cut max        %s"%self.cmax)
        print("cut cdpi       %s"%self.cdpi)

        print('original extents  %s %s %s %s '%(self.orig_minx, self.orig_miny, self.orig_maxx, self.orig_maxy))
        print('sorted extents    %s %s %s %s '%(self.sort_minx, self.sort_miny, self.sort_maxx, self.sort_maxy)) 



    ##-------------------------------##

    def show_buffers(self, sid=None):
        """
        
            show evrything meaningful we can about our data 
            id is optional 
            - if passed use it to look up gr_poly and gr_sort, etc  

        """
        
        print('#################')
        if sid is None:
            #[[id, centroid, extents, len, points ]
            #for ply in self.gr_sort:
            #    print("%s %s %s "%(ply[0], ply[1], len(ply[2]) )) 

            print('size gr_polys     %s'%len(self.gr_polys))
            print('size gr_sort      %s'%len(self.gr_sort))

        else:
            grp = self.gr_polys[sid]
            grs = self.gr_sort[sid]
            print('gr_polys id %s size %s  '%( sid, len(grp) ) )
            print('gr_sort  id %s json id:%s extents %s len %s size: %s '%( sid, grs[0], grs[1], grs[2], grs[3] ) )


    ##---------------------- 
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

    ##-------------------------------------------## 
    ##-------------------------------------------## 
    def gl_extents(self):
        """ global extents DEBUG - NOT TESTED

            we only work on gr_sort - gr_poly is a copy pf the orignial data

            if you run sort() - it automatically sets extents while it is sorting  
            if you want to re-run, use this 
        """

        # we only work on gr_sort - gr_poly is a copy pf the orignial data 
        # [[id, centroid, extents, len, points ]]
        for row in self.gr_sort:
            ply = row[4] 

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
            
            if minx<self.sort_minx:
                self.sort_minx=minx
            if miny<self.sort_miny:
                self.sort_miny=miny
            if maxx>self.sort_maxx:
                self.sort_maxx=maxx
            if maxy>self.sort_maxy:
                self.sort_maxy=maxy

    ##-------------------------------------------## 
    def gl_centroid(self):
        """ global centroid - DEBUG - NOT TESTED
            we only work on gr_sort - gr_poly is a copy pf the orignial data
 
        """

        self.gl_extents()
        
        width  = abs(self.sort_maxx - self.sort_minx) 
        height = abs(self.sort_maxy - self.sort_miny)         
        
        return ( self.sort_maxx-(width/2), self.sort_maxy-(height/2) )

    ##-------------------------------------------## 
    def gl_move(self, pos):
        """
            DEBUG - NOT TESTED  
            global move the entire dataset in 2d (ignore Z axis) 
            this allows 3D moves but you probably want 2D - zero the Z axis for pos or just send it 2 coords
        """
        print("global move ", pos) 

        pop = polygon_operator()

        for i,row in enumerate(self.gr_sort):
            if len(pos)==2:
                self.gr_sort[i][4] = pop.trs_points(self.gr_sort[i][4], (-pos[0], -pos[1]) )
            if len(pos)==3:
                self.gr_sort[i][4] = pop.trs_points(self.gr_sort[i][4], (-pos[0], -pos[1], -pos[0] ))

        # dont forget to recalculate extents 
        self.gl_extents()

    ##-------------------------------------------## 
    def gl_rotate(self, rot):
        """
            DEBUG - NOT TESTED  
            global rotate the entire dataset in 2d  
            this allows 3D moves but you probably want 2D - zero the Z axis for pos or just send it 2 coords
        """
        print("global rotate ", rot) 

        pop = polygon_operator()

        for i,row in enumerate(self.gr_sort):
                # TRS all in one 
                #self.gr_sort[i][4] = pop.trs_points(self.gr_sort[i][4], rotate=(-rot[0], -rot[1], -rot[2] ))

                # ROUNDING WAS CAUSING GLITCHES IN THE NCVIEWER APP - MAY OR MAY NOT BE AN ISSUE 
                self.gr_sort[i][4] = pop.rotate_pts(rot=(-rot[0], -rot[1], -rot[2] ), pts=self.gr_sort[i][4], doround=True)

        # dont forget to recalculate extents 
        self.gl_extents()

    ##-------------------------------------------## 
    def gl_scale(self, scale):
        """
            DEBUG - NOT TESTED  
            global rotate the entire dataset in 2d  
            this allows 3D moves but you probably want 2D - zero the Z axis for pos or just send it 2 coords
        """
        print("global scale ", scale) 

        pop = polygon_operator()

        for i,row in enumerate(self.gr_sort):
            if len(scale)==2:
                self.gr_sort[i][4] = pop.trs_points(self.gr_sort[i][4], scale=(-scale[0], -scale[1]) )
            if len(scale)==3:
                self.gr_sort[i][4] = pop.trs_points(self.gr_sort[i][4], scale=(-scale[0], -scale[1], -scale[0] ))

        # dont forget to recalculate extents 
        self.gl_extents()

    ##-------------------------------------------## 
    def gl_move_center(self):
        """
            we only work on gr_sort - gr_poly is a copy pf the orignial data
        """
        cen = self.gl_centroid()

        pop = polygon_operator()

        for i,row in enumerate(self.gr_sort):
            self.gr_sort[i][4] = pop.trs_points(self.gr_sort[i][4], (-cen[0], -cen[1], 0 ))

        # dont forget to recalculate extents 
        self.gl_extents()
 
    ##-------------------------------##
    def _scrub(self, inp):
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
    def _calculate_paths3d(self, do_retracts=True, doround=True):
        """ 
             DEBUGGY 
        """
        pl = 6 #numeric rounding places 
        lastpt = (0,0,0)

        self.outfile.append('(exported with _calculate_paths3d() )')
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
        ##------------------------------

        # #DEBUG - FILL_POLYS ARE A LAYOVER FROM KICAD STUFF - UNTESTED 
        # for fill_poly in self.filled_polys:
        #     if do_retracts:
        #         pt1 = fill_poly[0] 
        #         self.outfile.append('G0 x%s y%s z%s'% (  pt1[0], pt1[1], self.rh ) )  #first point at retract height   
        #         self.ngc_to_obj.append( ( self.hp[0], self.hp[1], self.rh) ) 
        #     for i,pt in enumerate(fill_poly):
        #         self.outfile.append( 'G1 x%s y%s z%s'%( pt[0], pt[1], pt[2] ) )
        #         self.ngc_to_obj.append( ( pt[0], pt[1], pt[2] ) )             
        #         lastpt =( round(pt[0]*self.global_scale,pl) , pt[1], pt[2] )
        #     self.outfile.append( 'G0 x%s y%s z%s'%(lastpt[0], lastpt[1], lastpt[2]) )
        #     self.ngc_to_obj.append( (lastpt[0], lastpt[1], lastpt[2])   ) 
        #     if do_retracts:
        #         self.outfile.append('g0z%s'% ( self.rh ) )  #retract in between cuts
        #         self.ngc_to_obj.append( (lastpt[0], lastpt[1], self.rh)   ) 
        #         self.outfile.append('  ')

        ####
        self.outfile.append('  ')
        self.outfile.append('(exporting graphic polygons )')
        
        ##------------------------------

        # graphical polygons - build the gcode up with simple linear movements
        for row in self.gr_sort:
            
            # preprocess 
            export_ply = True 

            if row[0] in self.omit_ids:
                print("# omitting polygon ID %s"%row[0])
                export_ply = False

            #modified to export sorted data - easy peasy  
            gr_poly = row[4]

            ##--

            if len(gr_poly) and export_ply:
                self.outfile.append('(exporting new polygon )')
                
                ##-- 
                #DEBUG - need to sort out clean points - want to run as close to final export 
                #it looses precision 

                # no formatting (full precision)
                if doround:                
                    pt1=(round(gr_poly[0][0],pl) ,round(gr_poly[0][1],pl)) 
                else:
                    pt1 = gr_poly[0]

                ## first point at retract height 
                if do_retracts:
                    #move to first point RH 
                    self.outfile.append('x%s y%s z%s'% (  pt1[0] , pt1[1], self.rh ) )               
                    self.ngc_to_obj.append( ( pt1[0], pt1[1], self.rh ) )             

                ## iterate points in polygon 
                self.outfile.append( 'G1' )
                for i,pt in enumerate(gr_poly):
                    if doround:
                        tmp = ( round(pt[0],pl), round(pt[1],pl), round(pt[2],pl) ) 
                        pt = tmp  

                    self.outfile.append( 'x%s y%s z%s'%( pt[0], pt[1],  pt[2] ) )
                    self.ngc_to_obj.append( (pt[0], pt[1],  pt[2]) )                   
                
                self.outfile.append( 'G0' )

                if do_retracts:
                    if doround:
                        gpt=(round(gr_poly[0][0],pl) ,round(gr_poly[0][1],pl)) 
                    else:
                        gpt=gr_poly[0]

                    self.ngc_to_obj.append( (gpt[0], gpt[1], self.rh)  )
                    self.outfile.append( 'x%s y%s z%s'%( gpt[0], gpt[1], self.rh) )

                    #### retract in between cuts
                    self.outfile.append('g0z%s'% ( self.rh ) )  

                self.outfile.append('  ')

        ##-----------------------------------------##        
        # self.outfile.append('(exporting segments )')

        ##-----------------------------------------##
        # rapid move at end 
        self.outfile.append('m2') #program end


    ##-------------------------------##
    def _calculate_paths2d(self, do_retracts=True):
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
        
        self.outfile.append('(exported with _calculate_paths2d() )')
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
        ##------------------------------

        #DEBUG - FILL_POLYS ARE A LAYOVER FROM KICAD STUFF - UNTESTED 
        for fill_poly in self.filled_polys:
            if do_retracts:
                pt1 = fill_poly[0] 
                self.outfile.append('G0 x%s y%s z%s'% (  pt1[0], pt1[1], self.rh ) )  #first point at retract height   
                self.ngc_to_obj.append( ( self.hp[0], self.hp[1], self.rh) ) 
            
            for i,pt in enumerate(fill_poly):
                self.outfile.append( 'G1 x%s y%s z%s'%( pt[0], pt[1], self.ch ) )
                self.ngc_to_obj.append( ( pt[0], pt[1], self.ch ) )             
                lastpt =( pt[0], pt[1], self.ch )

            self.outfile.append( 'G0 x%s y%s z%s'%(lastpt[0], lastpt[1], self.ch) )
            self.ngc_to_obj.append( (lastpt[0], lastpt[1], self.ch)   ) 
            if do_retracts:
                self.outfile.append('g0z%s'% ( self.rh ) )  #retract in between cuts
                self.ngc_to_obj.append( (lastpt[0], lastpt[1], self.rh)   ) 
                self.outfile.append('  ')

        ####
        self.outfile.append('  ')
        self.outfile.append('(exporting graphic polygons )')
        
        ##------------------------------

        # graphical polygons - build the gcode up with simple linear movements
        for row in self.gr_sort:
            
            # preprocess 
            export_ply = True 

            if row[0] in self.omit_ids:
                print("# omitting polygon ID %s"%row[0])
                export_ply = False

            #modified to export sorted data - easy peasy  
            gr_poly = row[4]
            
            ##--

            if len(gr_poly) and export_ply:
                self.outfile.append('(exporting new polygon )')

                
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
    ##-------------------------------##
    def export_grid_gfx(self, name, folder ):
        #export cells as graphics  
        
        self.gl_extents()

        if self.sort_minx==0 and self.sort_miny==0 and self.sort_maxx==0 and self.sort_maxy==0:
            raise ValueError('export_grid_gfx - extents are not set') 

        features = []
        for c in self.tesl.nodes:
            #cell centroids  
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

        #sort extents  
        coords = self.extents_fr_bbox([self.sort_minx,self.sort_miny,self.sort_maxx,self.sort_maxy], periodic=True)
        features.append(Feature(geometry=LineString(coordinates=coords), 
                                properties={"id" : 0 
                                           }
                                ) 
                        )
                        

        feature_collection = FeatureCollection(features)
        with open('%s/%s_ply_xtntx.json'%(folder,name), 'w') as f:
            dump(feature_collection, f)

    ##-------------------------------##
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

    ##-----------
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
    ##-------------------------------##
    def load_geojson(self, inputfile, zaxis, getfids=None, getids=None):
        """ parse a geojson file - store points in arrays in GR_POLY buffer 

            Args:
            
                inputfile - the file to load, duh. 
                
                zaxis - value to set false zaxis to (2d to 3d)  

                getfids - feature ids to load 
                
                SUB FEATURE NEEDS WORK - DEBUG ONLY WORKS WITH ONE FEATURE
                getids - sub-feature ids to load -  


        """

        print('loading geojson file %s'%inputfile)

        plyidx = 1
        ptidx = 1 

        geojson_txt = None

        with open(inputfile) as f:
            gj = geojson.load(f)
        features = gj['features'] 
 
        fixedids = [] 

        ##---------------
        ##---------------

        #check that ids are in range  
        if getfids:
            print("#load_geojson getone mode ", getfids) 
            for id in getfids:
                if id>=len(features):
                    print("id out of range %s"%id)
                else:
                    fixedids.append(id)

            print("DEBUG IDS TO IMPORT ARE ", fixedids )

        #check that ids are in range  
        if getids:
            print("#load_geojson getone mode ", getids) 
            for id in getids:
                if id>=len(features):
                    print("id out of range %s"%id)
                else:
                    fixedids.append(id)

            print("DEBUG IDS TO IMPORT ARE ", fixedids )

        ##---------------
        ##---------------
        # import all data ... 
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
                        self.jsonbuffer.append(tmp_poly) 
                plyidx += 1 

        ##---------------        
        # or just cherry pick some features to import (below we also can pick the sub-features) 
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
                                self.jsonbuffer.append(tmp_poly) 
                        plyidx += 1 

        ##---------------
        # optional processing to get sub features by id 
        # DEBUG seems sketchy - test this more 
        
        # get all sub features
        if getids is None:
            self.gr_polys = self.jsonbuffer 
        
        # pick out some sub features
        else:
            print("ids to load are is %s"%fixedids)
            for id in fixedids: 
                self.gr_polys.append(self.jsonbuffer[id]) 

        ##---------------

        print('cloning data, sorting and calculating extents ')
        #make a copy of the data to work with - leaving the original in case we want to do something with it  
        self._sort()
        
        print("loaded %s polygons from %s "%(plyidx,inputfile)) 

    ##-------------------------------##
    def export_extents_ngc(self, rh, ch, cdpi, cmax, filename, do3d=False):
        
        #self.gl_extents()

        tempbuffer = self.gr_sort
        self.gr_sort = []

        #self.prim_triangle('z', (0,0,0), (0,0,0) )
        # [[id, centroid, extents, len, points ]] 
        #self.gr_sort.append([0,0,0,0,self.points]) 
        
        pts = self.calc_square_diag((self.sort_minx,self.sort_miny ),
                                   (self.sort_maxx,self.sort_maxy), add_zaxis=True ) 
        pts.append( (self.sort_minx, self.sort_miny, 0) )
        self.gr_sort.append([0,0,0,0,pts])

        if do3d==True:
            self._calculate_paths3d()
        else:
            self._calculate_paths2d()
        
        print("gr_sort buffer has %s polys in it. "%(len(self.gr_sort)))
        fobj = open( filename,"w") #encoding='utf-8'
        for line in self.outfile: 
            fobj.write(line+'\n')
        fobj.close()

        self.gr_sort = tempbuffer

    ##-------------------------------##
    def export_ngc(self, rh, ch, cdpi, cmax, filename, do3d=False):
        print("# exporting NGC file ", filename)

        self.rh = rh          # retract height 
        self.ch = ch          # cut height (top, start of cut)
        self.cdpi = cdpi      # cut depth per iteration on Z axis
        self.cmax = cmax      # maximum cut depth on Z axis 

        if do3d==True:
            self._calculate_paths3d()
        else:
            self._calculate_paths2d()
        
        print("gr_sort buffer has %s polys in it. "%(len(self.gr_sort)))
        fobj = open( filename,"w") #encoding='utf-8'
        for line in self.outfile: 
            fobj.write(line+'\n')
        fobj.close()


    ##-------------------------------##
    ##-------------------------------##
    def make_grid(self, distance, folder, xcuts, ycuts, bbox=None):
        """ chop a square into smaller squares 

        """

        if bbox:
            self.tesl._set_extents(bbox) 
        else:
            self.tesl._set_extents([self.sort_minx, self.sort_miny, self.sort_maxx, self.sort_maxy]) 

        self.tesl.build_2d_cells(xcuts, ycuts)
        

        #temporary export of geom I HAVE TO SEE THIS BEFORE I GO TO BED!
        features = []

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
                    if float(dist) < float(distance):
                        
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

        # how deep in the parentheticals are we?
        depth = 0 
        newpoly = []

        for lc,line in enumerate(f):
            if '(' in line or ')' in line:
    
                #okay we are in a parenthetical - count the depth to break up the gr_polys
                for x in range(line.count('(')):
                    depth+=1
                for x in range(line.count(')')):
                    depth-=1
                
                ## ---

                linetoked = line.split(' ') 
                cleanspaces = []
                for t in linetoked:
                    if t!='':
                        cleanspaces.append( self._scrub(t) )
                
                ## ---
                if len(cleanspaces[0]):
                    #print(cleanspaces)
                    coord = [] 

                    # walk the cleaned line tokens, we know depth, we know if there was a gr_poly  
                    for i,tok in enumerate(cleanspaces):
                        if tok !='':
                            if tok=='gr_poly' and depth==2:
                                if len(newpoly):
                                    self.gr_polys.append(newpoly)
                                    newpoly = []

                            if depth==3:
                                if tok == 'xy':   
                                    xcoord = float( self._scrub(cleanspaces[i+1]) )
                                    ycoord = float( self._scrub(cleanspaces[i+2]) )

                                    coord.append( (xcoord, ycoord, self.ch ) )

                                # if tok == 'start':   
                                #     coord.append( (float(self._scrub(cleanspaces[i+1]))  , float(self._scrub(cleanspaces[i+2]) )) ) 
                                # if tok == 'end':   
                                #     coord.append( (float(self._scrub(cleanspaces[i+1]))  , float(self._scrub(cleanspaces[i+2]) )) ) 
                                # if tok == 'center':   
                                #     coord.append( (float(self._scrub(cleanspaces[i+1]))  , float(self._scrub(cleanspaces[i+2]) )) ) 
                        
                    #print(coord)    
                    if coord:
                        newpoly.append(coord[0] )
            
        #if data is in buffer - keep it 
        if len(newpoly):
            self.gr_polys.append(newpoly)

        ##----
        # copy all gr_polys to gr_sort to process
        self._sort()
















##------------------------------##
##------------------------------##            
##------------------------------##




def arc_to_degree(NS, degrees, minutes, seconds, EW):
    """ arc minutes to decimal degree ( example n50d0'02"e ) """
    
    outdegrees = 0.0

    if NS =='n':
        outdegrees = degrees
        outdegrees = outdegrees + (minutes*.0166667) #1/60
        outdegrees = outdegrees + (seconds*.0166667*.0166667) #1/60
    if NS =='s':
        outdegrees = 180.0
        outdegrees = outdegrees + degrees
        outdegrees = outdegrees + (minutes*.0166667) #1/60
        outdegrees = outdegrees + (seconds*.0166667*.0166667) #1/60
    if EW =='w' and NS =='s':
        outdegrees = outdegrees * -1
    if EW =='e' and NS =='n':
        outdegrees = outdegrees * -1
  
    return outdegrees

