
import os
import math 

from gnelscript import NUMPY_IS_LOADED




from gnelscript.pygfx.math_ops import math_util as mu
from gnelscript.pygfx.math_ops import matrix33, matrix44, vec2, vec3  

from gnelscript.pygfx.point_ops import pop3d





class polygon_operator(pop3d):
    """ polygon operator - should be called GEOM operator 

        3D model and tools to work on polygons 
       
        -------------------------------------------

        all code is experimental! 

        GEOM and GROUPS should really be another object type

        I was torn between simplicity of code and keeping a flat structure
        My reasoning was that if you have a format made of arrays of numbers 
        it offers more flexiblity for importing from other languages, formats, text files, etc  

        This tool was built as a toy to play with. It has potential to be a really awesome tool 
        but that would probably require a big think and a third rewrite. 

        -------------------------------------------

        Known problems:

            None of the tools "stiches" faces together.
            They obey rules of good topology, but all faces are unconnected

            - Extrude does not "clean up" the old geom it extrudes, it simply adds more on top 
            - Multi face triangulate will not work on a concave polygon, it needs to be symetrical 
            - GEOM format is a mess, it was designed to be an "OBJ file in memory"
              the indexing works, but is sloppy and the more you try to sub-select at a time the worse it gets 
            - EDGE GEOM does not obey the rules, it groups points in pairs, and is all wacky 
            


        --------------------------------------------
            Concepts 
        --------------------------------------------

        SUBSELECT              = slice (start, end ), ids = [LIST OF PTS]
  
        --------------------------------------------

        POLY                   = [ FID, FID, ... ]
        GEOM (just points)     = [ [(pt),(pt),(pt),(pt),... ] ]
        GEOM                   = [ [(poly),(poly),..], [(pt),(pt),(pt),(pt),... ] ]
        GEOM extended          = [ [(poly),..], [(pt),... ] , [normal...], [uv,...], [color,...]   ]

        --------------------------------------------
        
        PTGRP                  = [ [ID,(PT)], [ID,(PT)] ]

        FACEGRP (not implemented)          = [ [ID,(POLY)], [ID,(POLY)], .. ]

        --------------------------------------------

        NAMING 

        ids  = generic 
        fac  = face, polygon 
        fids = face ids 
        pids = process ids 

        -------------------------------------------
    """

    def __init__(self):
        super().__init__()  

        # geometry properties 
        self.uv_points       = []   # UV single point coordinates
        self.uv_polys        = []   # UV face indices  
        self.normals         = []   # face_normals

        self.points          = []    # list of tuples of XYZ points per vertex         -  [(x,y,z), (x,y,z)]  
        self.polygons        = []    # list of tuples for 2 or more vertex connections -  [(1,2,5,8) , (1,2)] 

        #self.point_normals  = []
        #self.point_colors   = []
        #self.line_colors    = []

        # render properties embedded in geometry 
        self.linecolors = None       # iterable of colors for lines 
        self.linecolor  = (0,240,00)
        self.vtxcolor   = (0,240,00)

        # "unavoidable side effect" variables - when reindexing between ops - we need to keep this
        self.exprt_ply_idx   = 1     # obj is NOT zero indexed
        #self.exprt_pnt_idx   = 0    # pts ARE zero indexed (everything BUT .obj face idx's are)

 
    @property
    def lastpt(self):
        """ get the highest indexed point geom in this object """
        return self.points[len(self.points)-1]

    @property
    def lastfid(self):
        """ get the highest indexed face index in this object """
        return self.polygons[len(self.polygons)-1]

    @property
    def numfids(self):
        """ get the total number of face indices in this object (-1 because of 0 index)"""
        return len(self.polygons)-1

    @property
    def numpts(self):
        """ get the total number of points in this object (-1 because of 0 index)"""        
        return len(self.points)-1


    ##-------------------------------------------## 
    def dump_pts_csv(self, pts, outfile):

        buffer = ''

        for pt in pts:
            print(pt)

            buffer = buffer+'%s,%s,%s\n'%(pt[0], pt[1], pt[2])

        #save it to disk now
        fobj = open( outfile,"w") #encoding='utf-8'
        fobj.write(buffer)
        fobj.close()
 


    ##-------------------------------------------## 
    def contiguous_segs_to_poly (self, linesegs ):
        """ DEBUG NOT DONE 
            take an array of [[start,end],..] and connect it to a single polygon
            this only works if the points are mathematically a match 

            (we could perhaps add some distance matching later )

            ARGS: 
                linesegs [[XYZ,XYZ], [XYZ,XYZ],.. ]

            RETURNS:
                polygon  [XYZ,XYZ,XYZ]                
        """
    
        #DEBUG
        pass 
        
    
    ##-------------------------------------------## 
    def clean_pts_str(self, pts=None):
        """ 
            returns clamped strings from float, should work with 2D or 3D points 

            destructive - may loose precision - only run for exporting 
            get rid of the pesky exponect/scientific notation on large/small floats 

            usage:
                pts = x.clean_pts_str([(-2.692345105448116e-09, 10.692345105448116e-09)])
                print(pts)

        """
        
        cleaned = []
        if pts==None:
            #print(self.points)
            for pt in self.points:
                cleaned.append( self._csp_str(pt) )
            self.points = cleaned
        else:
            #print(pts)
            for pt in pts:
                cleaned.append( self._csp_str(pt) )            
            return cleaned
    ##-------------------------------------------##  
    def _scribe(self, str):
        print(str)

    ##-------------------------------------------##           
    def _flush(self):
        """ set all geometry to a clean state """

        self.points          = [] 
        self.polygons        = []      
        self.normals         = []
        self.face_uvs        = []


    ##-------------------------------------------##
    def move(self, x,y,z ):
        # ( pts, translate=(0,0,0), rotate=(0,0,0), scale=(1,1,1) ):
        self.points = self.trs_points(self.points, (x,y,z) )

    ##-------------------------------------------##
    def rotate(self, x,y,z ):
        # ( pts, translate=(0,0,0), rotate=(0,0,0), scale=(1,1,1) ):
        self.points = self.trs_points(self.points,translate=(0,0,0), rotate=(x,y,z), scale=(1,1,1) )

    ##-------------------------------------------##  
    def z_sort(self, reverse=False):
        """ return a sorted list of polygons by average Z coordinate 
            lowest z value is (closest?) to camera 
            this is pretty much a crappy hack for the renderer 
            
            DEBUG!
                it gives the appearance of faces not being on top of each other
                Will be replaced by a visible faces algorithm someday 
                Be aware, the rendered renders EVERY polygon even when they are not seen
        """
        out = [];

        tmp = [] 

        #build a list of [mean z, [points, polygon] ]
        #sort by mean z value, and return the sorted list 
        for p in self.polygons:
            tri = []
            for idx in p:
                tri.append(self.points[idx-1])
            
            
            mean = self.triangle_mean_z(tri)
            #mean = self.centroid(tri)[2]

            if not isinstance(mean,float):
                print('error in  z_sort - mean is NOT float ', mean )
                mean = 0.0 
            if isinstance(mean,float):
                # print('## mean ',  mean )
                tmp.append( (mean, p) )
        
        def custom_sort(t):
            return t[0]
        
        #L = [("Alice", 25), ("Bob", 20), ("Alex", 5)]
        #L.sort(key=custom_sort)
        #print(L)

        #THIS IS BUGGY - ATTEMPTING MY OWN SORT 
        if reverse:
            tmp.sort(key=custom_sort, reverse=True)
        else:
            tmp.sort(key=custom_sort)
        for t in tmp:
            out.append(t[1])    
        
        self.polygons = out
        #return out

    ##-------------------------------------------##  
    def inspect_geom(self, geom):
        """ analyze a GEOM object and see how it is constructed """

        if geom == None:
            print('## inspect: error no geometry ')
            return None 

        print('## geom has %s top level items '%len(geom) )
        print('## geom has %s polygons        '%len(geom[0]) )
        print('## geom has %s points          '%len(geom[1]) )

        for i, poly in enumerate(geom[0]):
            print (' poly %s is type %s and has %s items '% (i,type(poly), len(poly)) )

    ##-------------------------------------------##  
    def verify_geom(self, geom, verbose=True):
        """ check that indices are within range for data 
            
            WELL FORMED GEOM IS: 

            [ 
                [(FIDS), (FIDS)], 
                [(PT),(PT),(PT)] 
            ]
 
        """
        
        if geom is None:
            if verbose:
                print('## debug verify geom is none ')
            return False 

        if not isinstance(geom[0],list) and not isinstance(geom[0],tuple):
            if verbose:
                print('## debug verify - bad data type first top elements ')
                print( type(geom[0]) ) 
            return False 

        if not isinstance(geom[1],list) and not isinstance(geom[1],tuple):
            if verbose:
                print('## debug verify - bad data type first top elements ')
                print( type(geom[1]) ) 
            return False 

        if geom[0] is None or geom[1] is None:
            if verbose:
                print('## debug verify geom element is none ')
            return False 

        if len(geom)!=2:
            if verbose:
                print('## debug verify geom does not have 2 items ')            
            return False 
 
        for ply in geom[0]:
            for pid in ply:
                if int(pid)>len(geom[1]):
                    if verbose:
                        print('## debug verify geom contains invalid index ')    
                    return False       
        return True 

    ##-------------------------------------------##  
    def _reindex_ply(self, f_idx, offset):
        """ only used with obj.append() 

            take a tuple of face indexes and re-index to offset+value
            not used yet, but may be a good thing to have 
        """ 

        out_face = [] 
        for i in f_idx:
           out_face.append(i+offset)
        return tuple(out_face)


    ##-------------------------------------------## 
    
    def triangle_mean_z(self, triangle):
        """ center of a triangle in 3d - used for render z sort 
            replace with centroid() ??  
        """
        z1 = triangle[0][2]
        z2 = triangle[1][2]
        z3 = triangle[2][2]
        return (z1+z2+z3)/3

    ##-------------------------------------------## 
    def move_center(self):
        """
           relocate an object so its centroid is at world origin 

            usage:
                pop3 = object3d()
                pop3.load('foo.obj')
                pop3.move_center()
                pop3.save('bar.obj')
        """                
        c = self.centroid()
        self.move(-c[0], -c[1], -c[2])

    ##-------------------------------------------##  
    @property
    def dimensions(self):
        """ X,Y,Z distance across object 

            usage:

                pop3 = object3d()
                pop3.load('file,obj')
                c = pop3.dimensions
                print(c)

        """

        bb = self.calc_3d_bbox()
        #[min_x, min_y, min_z, max_x, max_y, max_z ]
        x = abs(bb[3]-bb[0])
        y = abs(bb[4]-bb[2])
        z = abs(bb[5]-bb[3])
        
        return[x,y,z]

    ##-------------------------------------------## 
    def calc_2d_bbox(self, axis='z', pts=None, aspts=False):
        """ UNTESTED 
            derive a 2D BBOX from the built in 3D BBOX 
            you need to tell it which axis to generate on 
        """

        # 3d format is [0 min_x, 1 min_y, 2 min_z, 3 max_x, 4 max_y, 5 max_z ]
        if pts==None:
            bbox = self.calc_3d_bbox()
        else:
            bbox = self.calc_3d_bbox(pts=pts)


        if axis=='x':
            #miny, minz , maxy, maxz
            if aspts:
                return [(bbox[1],bbox[2]),
                        (bbox[1],bbox[5]),
                        (bbox[4],bbox[5]),
                        (bbox[4],bbox[2])]
            else:    
                return [bbox[1], bbox[2], bbox[4], bbox[5]]
        if axis=='y':
            #minx, minz , maxx, maxz
            if aspts:
                return [(bbox[0],bbox[2]),
                        (bbox[0],bbox[5]),
                        (bbox[3],bbox[5]),
                        (bbox[3],bbox[2])]
            else:             
                return [bbox[0], bbox[2], bbox[3], bbox[5]]       
        if axis=='z':
            #minx, miny , maxx, maxy
            if aspts:
                return [(bbox[0],bbox[1]),
                        (bbox[0],bbox[4]),
                        (bbox[3],bbox[4]),
                        (bbox[3],bbox[1])]
            else:             
                return [bbox[0], bbox[1], bbox[3], bbox[4]]

    ##-------------------------------------------##
    def calc_3d_bbox(self, pts=None, ptgrp=None, facgrp=None):
        """ 
            DEBUG - crashes if self.points is empty! 
                    made attempt to work with 2D data - treat Z as 0 , maybe a bad idea 

            get 3D bounds of points (or 3d object) 

            returns [min_x, min_y, min_z, max_x, max_y, max_z ] 

            usage:

                pop3 = object3d()
                pop3.load('file,obj')
                c = pop3.calc_3d_bbox()
                print(c)

        """

        if pts is None and ptgrp is None and facgrp is None:
            pts = self.points

        #grab a point from data to initialize
        min_x=pts[0][0]
        min_y=pts[0][1]
        if len(pts)==3:
            min_z=pts[0][2]
        else:
            min_z = 0 

        max_x=pts[0][0]
        max_y=pts[0][1]
        if len(pts)==3:        
            max_z=pts[0][2]
        else:
            max_z = 0 

        for pt in pts:
            if pt[0]<min_x:
                min_x=pt[0]
            if pt[0]>max_x:
                max_x=pt[0]
            
            if pt[1]<min_y:
                min_y=pt[1]
            if pt[1]>max_y:
                max_y=pt[1]
            
            if len(pt)==3:
                if pt[2]<min_z:
                    min_z=pt[2]
                if pt[2]>max_z:
                    max_z=pt[2]

        return [min_x, min_y, min_z, max_x, max_y, max_z ]

    ##-------------------------------------------##   
    def centroid(self, pts=None):
        """ get 3D center of a list of points 
            (average of a list of XYZ points) 

            usage:

                pop3 = object3d()
                pts = [(-1,0,0),(1,0,0)]
                c = pop3.centroid(pts)
                print(c)

        """

        if pts is None:
            pts = self.points

        if pts is None:
            print('warning centroid .. no points loaded')
            return None  
            
        ptsx = []
        ptsy = []
        ptsz = []
        
        if len(pts[0])==2:
            print('# centroid: error data appears to not be 3d')
            return None 
            
        if pts is not None:
            for pt in pts:
                ptsx.append(pt[0])
                ptsy.append(pt[1])
                ptsz.append(pt[2])  

            #average them 
            x = sum(ptsx)/len(ptsx)
            y = sum(ptsy)/len(ptsy)
            z = sum(ptsz)/len(ptsz)
            #return vec3(x,y,z)
            return [x,y,z]
        
        #hmm should we return None or 0,0,0 if no data... 
        return None  

    ##-------------------------------------------## 
    def three_vec3_to_normal(self, v1, v2, v3, unitlen=False):
        """ take 3 vec3 objects and return a face normal 

            usage:
                pop3 = object3d()

                v1 = vec3(-1,0,0)
                v2 = vec3(1,1,0)
                v3 = vec3(1,0,0)

                n = pop3.three_vec3_to_normal(v1, v2, v3)

        """
        # make sure type(v1,v2,v3) is vec3! 

        # calculate the face normal  
        a = v1 - v2
        b = v1 - v3

        if unitlen:
            f_nrml = a.cross(b).normal
        else:    
            f_nrml = a.cross(b)          
        
        return f_nrml 

    ##-------------------------------------------##  
    def calc_tripoly_normal(self, three_pts, unitlen):
        """  create a normal vector (vec3) from 3 points that represent a polygon  
             uses the internal function that requires vectors instead of points
             it allows "raw" point data to interface to that
             
        """

        v1=vec3();v2=vec3();v3=vec3(); 

        v1.insert( three_pts[0] )
        v2.insert( three_pts[1] )
        v3.insert( three_pts[2] )  

        return self.three_vec3_to_normal(v1, v2, v3, unitlen=unitlen)

    ##-------------------------------------------## 
    def any_pt_is_near(self, pt_list, pt2, dist ):
        """ run pt_is_near() for a list of points, 
            exiting if any are within distance
        """
        
        for pt in pt_list:
            if self.pt_is_near(pt, pt2, dist): 
                return True   
        return False

    ##-------------------------------------------## 
    def pt_is_near(self, pt1, pt2, dist ):
        """ compare two 3D points and return True if they are within 
            the specified distance to each other 

            pt1 = (0,0,0)
            pt2 = (2,2,2)
            o = object3d()
            out = o.pt_is_near(pt1,pt2,5)
            print(out)

        """
        # convert them to vec3 objects for the built in tools 
        if type(pt1)is tuple or type(pt1) is list:
            pt1vec  = vec3( pt1[0], pt1[1], pt1[2])
        else:
            pt1vec = pt1 

        if type(pt1)is tuple or type(pt1) is list:
             pt2vec  = vec3( pt2[0], pt2[1], pt2[2])        
        else: 
             pt2vec = pt2

        # get the vector between two points              
        b = pt1vec.between(pt2vec)

        if b.length <= dist:
            return True   
        return False

    ##-------------------------------------------## 
    ##-------------------------------------------##  
    # selection and inspection tools
    def select_by_location(self, select_type, pt_two, dist):
        """ UNFINISHED - seems to be selecting the wrong faces 
            subselct is wonky in a bunch of ways DEBUG 

            select by angle 
            direction to other things 
            select by distance to other objects, points , etc 

        """
        
        #pt2vec  = vec3( pt_two[0], pt_two[1], pt_two[2] )
        
        near_fids = [] 
  
        # not done yet 
        # if select_type == 'points':
        #     # for pnt in self.points:
        #     #     # convert the point to a vec3 object 
        #     pass


        if select_type == 'polygons':

            for fid in range(len(self.polygons)):
                pts = self.get_face_pts(fid)
                if self.any_pt_is_near(pts, pt_two, dist):
                    near_fids.append(fid)                   

            return near_fids

        return None
    ##-------------------------------------------##
    def ray_hit(self, ray_orgin, ray_vector, fastmode=False, flipnormal=False):
        """ ONLY WORKS FOR TRIANGLES  
         
            Iterates all polygons (tris) return [[fid, hitlocation]]

            fastmode will stop after it finds a single hit. If off it will find all hits. 

            return a polygon id for a ray intersect 
         
            usage:
                pop = object3d()
                pop.load('3d_obj/cube.obj')
                pop.triangulate()
                ray = (vec3(0,0,-1), vec3(0, 0, 2))
                hits = pop.ray_hit( ray[0],  ray[1])
                print(hits)  
            
                ##----------

                pop = object3d()
                pop.load('3d_obj/cube.obj')
                pop.triangulate()
                pop.rotate(45,45,45)
                ray = (vec3(0,0,-3), vec3(0, .1, 1))
                hits = pop.ray_hit( ray[0],  ray[1])
                x = object3d()
                x.one_vec_to_obj(ray[1], ray[0])
                for h in hits:
                    x.prim_locator(h[1], size=.1)
                    g = pop.get_face_geom(h[0], reindex=True, geom=None)
                    pop.exprt_ply_idx = 1
                    x.insert(g)
                x.save('ray_hits.obj')
                pop.save('ray_obj.obj')

        """

        hits = []
        test = vec3()

        for x in range(self.numfids):
            geom = self.get_face_geom(x)

            if len(geom[1])>2:
                
                if flipnormal:
                    # right or left handed coords - IDK 
                    vc1=vec3(geom[1][0]) 
                    vc2=vec3(geom[1][1]) 
                    vc3=vec3(geom[1][2])
                else: 
                    # right or left handed coords - IDK 
                    vc1=vec3(geom[1][2]) 
                    vc2=vec3(geom[1][1]) 
                    vc3=vec3(geom[1][0])

                result = test.ray_tri_intersect(ray_orgin, ray_vector, vc1, vc2, vc3)
                if result:
                    if fastmode:
                        return [[x,result]]
                    else:    
                        hits.append([x,result])

        # annoying - we must reset this so other get_face_geom will work 
        self.exprt_ply_idx = 1

        return hits
    
    ##-------------------------------------------## 
    # def raster_lines_poly(self, idx):
    #     """generate a series of vector lines that represent a 
    #        series of rasters across a polygon
    #     """
             
    ##-------------------------------------------## 
    def geom_to_ptgrp(self, geom):
        """ convert one weird data type into another 
        """

        out = []
        for i,g in enumerate(geom):
            out.append( [1, g[1][i]] )
        return out     

    ##-------------------------------------------##  
    def sub_select_geom(self, span=None, ids=None, reindex=False):
        """ 
            make work with xform_pts, rotate_pts, scale_pts 

            # DEBUG - This works for now, but... 
            # the only way to do this right is make get_face_geom() smarter
            # it needs to understand slicing, and not select the same points twice while iterating faces 
            # for example, do a subselect on a cube IDS 1-6 to see what I mean 
            # it will give you 24 points instead of 8 


            span - tuple of (start,end)  
            ids   - list of single ids 

            quick select chunks of geometry to feed into other tools: 
   
            get one or more faces as a new object 
            specify a list of ids, or a range
        """
       
        out_poly = []
        out_pts  = []

        # reset this when exporting with reindex 
        # it needs to be stored outside the function
        # this allows multiple polygons to be exported with an incrementing "relative" index  
        # relative, to each sub select 
        # other function can auto increment, thus allowing polygons reordering in chunks 
        
        if span:
            if span[1]=='n' or span[1]=='N' or span[1]>self.numfids:
                span = (span[0], self.numfids )

        pids = self.indexer(span=span, ids=ids)

        # print('## debug, indexer_geom pids : ', pids) 

        for i in pids:
         
            # DEBUG - This works for now, but... 
            # the only way to do this right is make get_face_geom() smarter
            # it needs to have slicing and not select the same thing twice 
            geom = self.get_face_geom(i, reindex=reindex)
            out_poly.append(geom[0])
            for pt in geom[1]:
                out_pts.append(pt)               

        if len(out_poly[0])==0:
            return None

        return ( out_poly, out_pts )

    ##-------------------------------------------##  
    def get_pt_ids(self, fids=None):
        """ lookup point indices (poly) from a list of face IDs 


        """
        out_poly = [] 
        for i in fids:
             tmp = self.get_face_geom(i, reindex=False)
             if tmp:
                 out_poly.append(tmp[0])
             if tmp == None:
                print('## error looking up face ids for ID %s '%i)
                return None
        
        return out_poly 

    ##-------------------------------------------##   
    def get_face_pts(self, fid, asvec3=False):
        """ lookup and return the point geometry of a face from a face ID
            for fancier features look at get_face_geom() 
            
        """

        tmp = []

        if fid<0 or fid > len(self.polygons)-1:
            print('# get_face_pts- bad face index : %s'%fid)
            return None

        for v_id in self.polygons[fid]:
            if asvec3:
                tmp.append(vec3(self.points[v_id-1]))
            else:    
                tmp.append(self.points[v_id-1])

        return tmp

    ##-------------------------------------------##   
    def pts_to_ptgrp(self, pts):
        """ 
            pts is [pt, pt, pt  ]
            ptgrp is [ [id, pt], [id, pt] ]

        """

        out = []
        # print('############# ptgrp ', len(ptgrp), '  ', len(self.points) )  

        for i,p in enumerate(pts):
            out.append( [i,p] )  
        return out   

    ##-------------------------------------------##   
    def insert_pt_grp(self, ptgrp):
        """ ptgrp is [ [id, (pt)], [id, (pt)] ]
            this will overwrite internal geometry with point group geometry 

        """

        out = []
        
        # print('############# ptgrp ', len(ptgrp), '  ', len(self.points) )  
         

        for p in ptgrp:
            #print('### ptgrp pt ', p , len(self.points))
            # debug this runs out of range on last point 
            # KEEP ZERO INDEX for all but OBJ file? debug 
            self.points[p[0]] = p[1]
        return out    

    ##-------------------------------------------##  
    def append_pt_grp(self, ptgrp):
        """ same as a list, but compatible with indexed point groups
            ptgrp is [ [id, pt], [id, pt] ]
        """
        out = []
        for p in ptgrp:
            self.points.append(p[1]) 
        return out 

    ##-------------------------------------------##  
    def get_pt_grp(self, span=None, ids=None):
        """ gets a point group, 
            a point group is a list of  

            [ [ID, (X,Y,Z)], ... ] 

            we can process them and put them back where 
            we found them in the model, allowing for modification of partial or whole objects 

            if nothing is specified, get all the points from self

            usage:


        """

        out = []

        if span is None and ids is None:
            for i,p in enumerate( self.points ):
                out.append([i,p])  
            return out

        else:
            pids = self.indexer( span=span, ids=ids)

            for p in pids:
                out.append( [p, self.points[p]] )
            return out      

    ##-------------------------------------------##  
    def get_geom_edges(self, geom ):
        """ 
            GEOM is [[IDS], [PTS] ]


            takes a geom object and returns another geom of the edges 
            it does this by iterating poly indices in groups of 2 

        """
        DEBUG = False  

        if type(geom)is not list and type(geom)is not tuple:
            raise ValueError('input must be geom type (list/tuple) not %s'%type(geom)) 

        out_edge_ids = []
        out_edge_pts = []

        polys = geom[0][0] 

        for ply in polys:
            
            if DEBUG:
                print('### POLY     ', ply   )
                print('### POINTS   ', geom[1] , type(geom[1]) )
            
            for idx in range(len(ply)):
                # iterate by two and store segments
                cidx = ply[idx-1]
                out_edge_ids.append((  cidx            ,ply[idx]          )) # poly index
                out_edge_pts.append((  geom[1][cidx-2] , geom[1][cidx-1]  )) # point index

        return [out_edge_ids, out_edge_pts]

    ##-------------------------------------------##   
    def get_face_edges(self, fid, reindex=False, geom=None):
        """ UNTESTED 
            return [[VTX_IDS], [VTX_PTS]]

        """
        if geom is None:
            geom = self.get_face_geom(fid, reindex=True)  # [[poly idx], [pt data]] 
        return self.get_geom_edges(geom)

    ##-------------------------------------------##  
    def get_face_geom(self, fids,  reindex=False, geom=None):
        """ lookup and return the polygon indices and points for a single polygon 

            
            ANNOYING - YOU MUST RESET THIS GLOBAL BEFORE USING 
            self.exprt_ply_idx


            reindex - if True  - renumber the new polygon indices startring at 1, 
                      if False - retain the original numbering 

            geom - act on a geom obj passed in, or on self

            usage: 

                o = object3d()
                o.load('3d_obj/cube.obj')
                g = o.get_face_geom(3, reindex=True, geom=None)

                x = object3d()
                x.insert(g)
                x.save('newobj.obj')



        """

        # validate inputs 
        #if fids<1 or fids > self.numply:
        #    print('# get_face_geom- bad face index : %s'%fids)
        #    return None

        # decide what the input is, fallback on self.poly/self.points  
        if geom is None:
            polygr  = self.polygons
            pointgr = self.points
        else:
            polygr  = geom[0]  #faces and polys common indexed    
            pointgr = geom[1]        


        if self.verify_geom( [polygr, pointgr] ) is False:
            return None 
        
        #if type(fids) is list:
        if type(fids) is int:
            fids = [fids]

        tmp_pts = []
        out_geom = [[],[]]

    
        for f_id in fids: 
            reindex_id = []             
            if f_id<len(polygr):  
                for v_id in polygr[f_id]:
                    # keep a count of points stored to use as new index
                    reindex_id.append(int(self.exprt_ply_idx ))
                    # store points that are indexed in geom 
                    tmp_pts.append(pointgr[v_id-1]) #data is NOT zero index but all else IS 
                    self.exprt_ply_idx +=1

                # geom is always [ [(poly),..], [(point),(point),...]  ]
                if reindex is False:
                    out_geom[0].append(polygr[f_id])
                if reindex is True:
                    out_geom[0].append(tuple(reindex_id))
            
            out_geom[1] = tmp_pts
        return out_geom

    ##-------------------------------------------##   
    def get_face_normal(self, fid=None, unitlen=False ):
        """ lookup a face(s) and calculate a face normal(s) for it  
            
            fid IS NOT zero indexed ... we need to fix this inconsistency 
            DEBUG i think the idea is that OBJ files are not ... all functions should do either


            only tested for 3 or 4 sided polygon 
            also returns the center position of a face

            returns vec3 type 

            DEBUG - the need for a standardized interface for slice, fid lookup is very apparent 
        """
        
        if fid == None:
            print("## error get_face_normal - need face id(s) to get normal")
            return None 
        if fid>len(self.polygons):
            print("## error get_face_normal - face id is too high")
            return None 

        if isinstance(fid, int):
            fids = [fid]
        if isinstance(fid,list):
            fids = fid
 
        out = []

        for f in fids:    
            tmp = self.get_face_geom(f) #returns [fidx, pts] 

            poly = tmp[0][0] #poly = face indices  

            nrmlvec = self.calc_tripoly_normal( (self.points[poly[0]-1],
                                                 self.points[poly[1]-1],
                                                 self.points[poly[2]-1]),
                                                 unitlen )
              
            out.append( nrmlvec )

        if isinstance(fid, int):
            return out[0]
        else:
            return out   

    ##-------------------------------------------##         
    def get_face_centroid(self, fid):
        """ lookup a face by id, return a 3d point representing the center
            wont work on nasty concave bad topology
            if the face is roundish, symetical, etc, it will work okay 
        """
        pts = self.get_face_pts(fid)
        out = self.centroid(pts) 
        return out
    ##-------------------------------------------## 
    ##-------------------------------------------##  
    # operators that modify geometry data and/or build new geom 

    ##-------------------------------------------##  
    def apply_matrix_ptgrp(self, ptgrp, m33=None, m44=None):
        """ same as apply_matrix_pts() but for point groups 

            batch mutliply a point group by a matrix 
            used for translate, rotate, and scaling. 
        """

        pt_ids = []
        pt_data = []

        tmp_buffer = [] 
        for pt in ptgrp:
            #print('## debug apply matrix ', pt )

            pt_ids.append( pt[0])
            pt_data.append(pt[1])

        # apply the transform here
        for pt in pt_data:  
            if m33 is not None:
                tmp_buffer.append( m33 * pt )
            if m44 is not None:
                tmp_buffer.append( m44 * pt )

        out = []
        for i,p in enumerate(tmp_buffer):
            out.append( [pt_ids[i],p ] ) 
        return out

    ##-------------------------------------------##   
    def scale_pts(self, amt, pts=None, ptgrp=None ):

        # build a scale matrix 
        
        amtx = 1.0
        amty = 1.0
        amtz = 1.0
                        
        if isinstance(amt,tuple) or isinstance(amt,list):
            amtx = amt[0]
            amty = amt[1]
            amtz = amt[2]
        else:
            amtx = amt
            amty = amt
            amtz = amt            

        sc_m33 = self.m33.identity
        sc_m33[0]  = amtx
        sc_m33[4]  = amty
        sc_m33[8]  = amtz    

        # sc_m44 = self.m44.identity
        # sc_m44[0]  = amt
        # sc_m44[5]  = amt
        # sc_m44[10] = amt        
    
        ################################################
        if pts and ptgrp is None: 
            return self.apply_matrix_pts(pts, m33=sc_m33)  # m44=sc_m44 

        if pts is None and ptgrp is None:
            # no args gets all the points of this object 
            ptgrp = self.get_pt_grp()    
  
        scaled = self.apply_matrix_ptgrp(ptgrp, m33=sc_m33)  # m44=sc_m44
        self.insert_pt_grp(scaled)

    ##-------------------------------------------##   
    def rotate_pts(self, rot, pts=None, ptgrp=None, doround=False):
        """  
            rotate some points with a matrix
            works on a pointgroup, or a list of points 
            if none specified, apply to entire object (self.points)

        """
        # construct a rotation matrix from euler angles 
        # angles are entered in degrees XYZ

        rx=rot[0]
        ry=rot[1]
        rz=rot[2]

        # degree to radian function 
        dtr = self.mu.dtr

        # build rotationY (see diagram above) 
        y_matrix =  self.m44.identity
        y_matrix[0]  =  math.cos(dtr( ry ))
        y_matrix[2]  = -math.sin(dtr( ry ))
        y_matrix[8]  =  math.sin(dtr( ry ))
        y_matrix[10] =  math.cos(dtr( ry ))
              
        # build rotationZ (see diagram above) 
        z_matrix    =  self.m44.identity
        z_matrix[0] =  math.cos(dtr( rz ))
        z_matrix[1] =  math.sin(dtr( rz ))
        z_matrix[4] = -math.sin(dtr( rz ))
        z_matrix[5] =  math.cos(dtr( rz ))
        tmp_matr = y_matrix * z_matrix

        # build rotationX (see diagram above) 
        x_matrix =  self.m44.identity
        x_matrix[5]  =   math.cos(dtr( rx )) 
        x_matrix[6]  =   math.sin(dtr( rx )) 
        x_matrix[9]  =  -math.sin(dtr( rx ))
        x_matrix[10] =   math.cos(dtr( rx ))

        rot_matrix = self.m44.identity
        rot_matrix = x_matrix * tmp_matr   
       
        ################################################
        # if points passed in but no point group operate on pts
        if pts is not None and ptgrp is None: 
            if doround:
                return self.apply_matrix_pts_round(4, pts,  m44=rot_matrix)  
            else:
                return self.apply_matrix_pts(pts,  m44=rot_matrix) 
            

        # if neither is specified, apply to whole object 
        if pts is None and ptgrp is None:
            ptgrp = self.get_pt_grp()    

        # if point group is passed put output into that 
        rotated = self.apply_matrix_ptgrp(ptgrp, m44=rot_matrix) 
        
        self.insert_pt_grp(rotated)
        #self.append_pt_grp(rotated)

    ##-------------------------------------------## 
    def xform_pts(self, pos, pts=None, ptgrp=None):
        """ 
            BROKEN - FIX 
            shift points without using a matrix 
            if no points are specified - apply to whole object 
        """

        if pts is not None and ptgrp is None: 
            out = []
            for pt in pts:
                x = pt[0] + pos[0]
                y = pt[1] + pos[1]
                z = pt[2] + pos[2]
                out.append( (x,y,z) )
            return out     


        if pts is None and ptgrp is None: 
            # no args gets all the points of this object 
            ptgrp = self.get_pt_grp()    

        if ptgrp:
            pt_ids  = []
            pt_data = []
            tmp_buffer = [] 

            for ptg in ptgrp:
                pt_ids.append(  ptg[0] )
                pt_data.append( ptg[1] )
            
            for i,pt in enumerate(pt_data):  
                x = pt[0] + pos[0]
                y = pt[1] + pos[1]
                z = pt[2] + pos[2]
                tmp_buffer.append( [i+1,(x,y,z)] )

            #print(tmp_buffer)    
            self.insert_pt_grp(tmp_buffer)

        return None 

    ##-------------------------------------------##          
    def insert_line(self, pts):
        """ insert tuple into THIS object as a 3d line 
            auto add it as another 2 point polygon  
        """

        lastidx = len(self.points)

        # DEBUG make sure num points is divisible by two, or ignore last one? 

        #look at first point and assume all data is similar
        #works with 2D data only 
        if len(pts[0])==2:
            #print("insert_line: data appears to be 2D")
            for i,pt in enumerate(pts):
                self.points.append( (pts[i][0], pts[i][1], 0) )
                if i%2==0:
                    self.polygons.append( [lastidx+1, lastidx+2] ) 


        if len(pts[0])==3:
            self.points.extend(pts)

    ##-------------------------------------------##  
    def linegeom_fr_points(self, pts, color=(100,0,100), periodic=False, singlepoly=False ):
        """ create renderable lines from array of 3D pts 
        """
        lptidx = self.numpts
        
        buffer = []
        for i in range(len(pts)):
            if i>0:
                pt1 = pts[i-1]
                pt2 = pts[i]
               
                self.points.append( (pt1[0], pt1[1], pt1[2], color[0], color[1], color[2]) ); lptidx+=1
                self.points.append( (pt2[0], pt2[1], pt2[2], color[0], color[1], color[2]) ); lptidx+=1 
                # each segment becomes a two point polygon 
                if singlepoly == False: 
                    self.polygons.append([lptidx-1, lptidx])
                # make all segmants as one big polygon    
                if singlepoly:
                    buffer.append(lptidx-1)
                    buffer.append(lptidx)
        
        # make all segmants as one big polygon
        if singlepoly:
            self.polygons.append(buffer)

        if periodic:
            pt1 = pts[0]
            pt2 = pts[len(pts)-1]            
            self.points.append( (pt1[0], pt1[1], pt1[2], color[0], color[1], color[2]) ); lptidx+=1
            self.points.append( (pt2[0], pt2[1], pt2[2], color[0], color[1], color[2]) ); lptidx+=1 
            self.polygons.append([lptidx-1, lptidx])

    ##-------------------------------------------## 
    ##-------------------------------------------## 

    def push_insert_poly(self, ids):
        """
            push a stack of ids by reindexing
            insert directly into self 
        """
        n = self.numpts 
        for p in ids:
            t=[]
            for idx in p:
                t.append(idx+n)
            #print('new polygon ',n, t)
            self.polygons.append(t)
   

    def push_plys(self, ids):
        """
            added as convenience for insert_polygons()
            push a stack of ids by reindexing
            array in <-> array out 

        """

        n = self.numpts 
        out = []
        for p in ids:
            t = []
            for idx in p:
                t.append(idx+n)
            out.append(t)
        return out     

    #def push_pts(self, ids):


    def insert_polygons(self, plyids=None, points=None, geom_obj=None, 
                              incrone=False, pos=None, ans=True ):
        """    
            Insert new geometry into an object
                you can pass a geom object to insert into 
                if you dont specify it goes into self 
                you can merge it with existing points or add new points 
            data is NOT ZERO INDEXED 
            you can pass in zero indexed data if you iuse incrone option 
            - if you pass no plyids it will exit 
            - if you pass plyids with no points it will merge polygons into existing points
            -if you pass plyids with new points it will ALWAYS be a new shell regardless of asn 
            it will attempt to merge into existing geometry first unless ans/as new shell is True 
            plyids, points    - use these or geom, but not both at same time 
            ans/ asnew_shell  - reindex the points and append, else keep the same indices
            geom_obj          - geom to insert into, instead of object.polygons, object.points
                                if true, will return the geom object when done 
            incrone         - offset IDS by +1 (zero index nightmare)
            pos             - offset point posititons  
            if points are 2D - automatically insert in 3D on the 0 Z axis
        
            if you pass plyids with no points and "ans True" it makes a copy of existing points 
            as new shell True will reindex new polygons 
            any points automatically are treated as new shell 
            as new shell False is only for adding polygons to existing geom 
        """

        if points is None:
            points=[]

        if not plyids and not geom_obj:
            raise ValueError("insert_polygons: no new polygons to insert")    
        
        #if plyids and points:
        #    ans=True 

        #if geom_obj and plyids:
        #    raise ValueError("insert_polygons: use geom obj or polyids but not both at same time") 
        #if geom_obj and points:
        #    raise ValueError("insert_polygons: use geom obj or polyids but not both at same time") 
        

             
                      
        #############################              
        #debug checks for data -  make sure it is good first 
        
        # no empty point data 
        if points is None and geom_obj is None and self.numpts==0 :
            raise ValueError("insert_polygons: no existing or new point data to work with")

        n = self.numpts 
        for ply in plyids:
            for fid in ply:
                #non numeric type
                if not isinstance(fid, int):
                    raise ValueError('## insert_polygons, bad data for index: %s'%fid)
                # no indexes less than 0 with incrone                
                if incrone:
                    if fid<0:
                        raise ValueError("insert_polygons: poly id too low: %s"%fid)  
                # no indexes less than 1 without incrone 
                else:    
                    if fid<1:
                        raise ValueError("insert_polygons: poly id too low: %s"%fid)  
                # no indexes larger than dataset 
                #if fid>n and fid>len(points):
                #    print("##### ",len(points)-1, n-1)
                #    raise ValueError("insert_polygons: poly id too high: %s"%fid)                    

        # no duplicate indexes in same polygon
        # no duplicate indexes in different order

        if points:
            #if vec3 data was passed in convert to tuple
            if type(points[0]) is vec3:
                tmp = []
                for pv in points: 
                    tmp.append(pv.aspt)
                points = tmp
        #############################
        tmp_pts = [] 
        tmp_plys = [] 

        #############################

        # option to offset the points in space 
        if pos:
            tmp=[]
            for pt in points:
                tmp.append((pt[0]+pos[0], pt[1]+pos[1], pt[2]+pos[2])) 
            points=tmp

        #############################
        # option to increment polyids by one automatically to solve zero index problem 
        newplyids = []
        if incrone:
            for pid in plyids:
                tmp = []
                for pt in pid:
                    tmp.append(int(pt)+1)
                newplyids.append(tmp)
            plyids = newplyids

        ############################# 
        # append polygons (using existing point buffer)
        if ans is False:   
            for poly in plyids:
                plytmp = []      
                for idx in poly:
                    plytmp.append(idx)  

                # do the insert operation                    
                if geom_obj is None: 
                    tmp_plys.append( tuple(plytmp) ) 
                else:
                    geom_obj[0].append( tuple(plytmp) )

        ################## 
        # do the insert operation
        if ans:
            if geom_obj is None: 
                tmp_pts.extend(points)
                tmp_plys.extend(self.push_plys(plyids))
            else:
                geom_obj[0].extend(plyids)  
                geom_obj[1].extend(points)    
        if not ans:
            if geom_obj is None: 
                tmp_pts= points
                tmp_plys = plyids
            else:
                geom_obj[0]= plyids  
                geom_obj[1]= points                  
        ##-----
        # merge the new geom into self 
        if geom_obj:
            return geom_obj
        else:
            self.points.extend(tmp_pts)
            self.polygons.extend(tmp_plys)


    ##-------------------------------------------## 
    def add_poly_frpts(self, pts):
        """ build a new triangle and auto-step the fids """
        npts = self.numpts 
        self.points.extend(pts)
        self.polygons.append( [npts+1,npts+2,npts+3] )
    
    ##-------------------------------------------## 
    def add_quad_frpts(self, pts):
        """ build a new 4 sided polygon and auto-step the fids """
        npts = self.numpts 
        self.points.extend(pts)
        self.polygons.append( [npts+1,npts+2,npts+3, npts+4])

    ##-------------------------------------------##   
    ##-------------------------------------------## 
    # MODELING TOOLS 

    ##-------------------------------------------##  
    def revolve_points(self, numdivs, axis, pts):
        """ simple lathe function 

        """

        
        step = int(360/numdivs)
 
        out = []
        if axis == 'y':
            for r in range(1, 360, step):
                #self.rotate_pts( (0,r,0), ptgrp=self.pts_to_ptgrp(pts) )
                out.append(self.rotate_pts( (0,r,0), pts=pts) )

        return out
        #return  self.modulo(numdivs, self.points) 

    ##-------------------------------------------## 
    def lathe(self, pts, num, axis='y'):
        """ spin a set of 3d points 360 degrees and make a renderable surface 
            needs to have the same num U and V to work

            usage: 
            
            # simple example 
            
                obj = object3d()
                pts = [(.1,.1,0),(1,1,0),(2,2,0),(3,3,0)]
                obj.lathe(pts, 4)


            # using bezier curve function 

                obj = object3d()
                num = 23
                start = (1 ,  0, 0)
                ctrl1 = (.5,  0, 0)
                ctrl2 = ( 0, .5, 0)
                end   = (0 ,  1, 0)
                curve = obj.cubic_bezier(num, start, ctrl1, ctrl2, end)
                obj.lathe(curve, num)
                obj.save('lathe.obj')

        """

        # use readable indices for testing iterator
        # pt_grid = [ ['a','b','c','d'],
        #             ['e','f','g','h'],
        #             ['i','j','k','l'],
        #             ['m','n','o','p'] ]
        #  print('##################\n\n')
        # self.print_grid(pt_grid)

        pt_grid = self.revolve_points( num, axis, pts )

        #print( pt_grid )
        viewhulls = False 

        #view hulls for debugging
        if viewhulls:
            for n in range(num):
                rows = self.get_grid_column( pt_grid , n)
                cols = self.get_grid_row( pt_grid    , n)
                self.linegeom_fr_points( rows )
                self.linegeom_fr_points( cols )                
       
        for u in range(num):
            for v in range(num):
                if u>0 and v>0:
                    tri1 = []
                    tri1.append( pt_grid[u][v]    )
                    tri1.append( pt_grid[u-1][v-1])
                    tri1.append( pt_grid[u][v-1]  )
                    self.add_poly_frpts(tri1)
                    tri2 = []
                    tri2.append( pt_grid[u-1][v]  )
                    tri2.append( pt_grid[u-1][v-1])
                    tri2.append( pt_grid[u][v]    )
                    self.add_poly_frpts(tri2)
                
                #if last row connect back to the first     
                if u==num -1:
                    if v>0:
                        tri2 = []
                        tri2.append( pt_grid[u][v-1] )
                        tri2.append( pt_grid[0][v] )
                        tri2.append( pt_grid[u][v]   )
                        self.add_poly_frpts(tri2)
                    if v<num-1:
                        tri1 = []
                        tri1.append( pt_grid[0][v+1] )
                        tri1.append( pt_grid[u][v]   )
                        tri1.append( pt_grid[0][v]   )
                        self.add_poly_frpts(tri1)

    ##-------------------------------------------## 
    def lathe2(self, pts, num):
        """ attempt to build quads instead of tris"""
        pt_grid = self.revolve_points( num, 'y', pts )

        for u in range(num):
            for v in range(num):
                if u>0 and v>0:
                    tri1 = []
                    tri1.append( pt_grid[u][v]    )
                    tri1.append( pt_grid[u][v-1]  )
                    tri1.append( pt_grid[u-1][v-1])
                    tri1.append( pt_grid[u-1][v]  )
                    self.add_quad_frpts(tri1)
                    #self.linegeom_fr_points(tri1, color=(100,0,100) )
                
                #if last row connect back to the first     
                if u==num -1:
                    tri1.append( pt_grid[u][v]    )
                    tri1.append( pt_grid[u][0]  )
                    tri1.append( pt_grid[u-1][0])
                    tri1.append( pt_grid[u-1][v]  )
                    self.add_quad_frpts(tri1)

    ##-------------------------------------------##  
    def extrude_face(self, f_id, distance):
        """ 
            DEBUG ITS BROKEN 
            proof of concept - but it makes bad topology - the edges and cap polys dont share verts

            obj = object3d()
            obj.load('3d_obj/cube.obj')
            obj.extrude_face(2, 2)
            obj.extrude_face(3, 2)  <--- CRASHES IF YOU RUN TWICE 
            obj.save('ext.obj')

        """
        DEBUG = False 

        geom  = self.sub_select_geom(ids=[f_id] , reindex=True)
        if geom == None:
            raise ValueError('no geometry found for FID %s'%f_id)

        nrml = self.get_face_normal(fid=f_id, unitlen=True) 
        
        nrml = nrml * distance 
        # edge selection iterates a polygons points 2 at a time, 
        # and forms each pair of points into a new line segment 
        s_edges = self.get_geom_edges(geom)  

        # move the face up and build the walls connecting to it 
        moved = self.xform_pts( nrml, pts=geom[1])
        e_edges = self.get_geom_edges([geom[0],moved]) 

        if DEBUG:
            print(geom)
            print(s_edges)
  
        # "wall" polygons, geometry connecting the new poly to the old  
        # iterate one set of edges assuming they both have the same number 
        for w in e_edges[0]:
            wall_poly = []
            #this is the old "broken" extrude that works, but uses a paired point data struct
            wall_poly.extend(s_edges[1][w[0]-1]) # bottom half of quad polygon 
            wall_poly.extend(e_edges[1][w[0]-1]) # top half of quad polygon  
            
            # stitch the 4 points into a quad polygon                 
            self.insert_polygons( plyids=[(1,2,4,3)], points=wall_poly, ans=True) 

        # transformed face along normal (cap polygon) 
        self.insert_polygons(plyids=geom[0][0], points=moved, ans=True) 
 


    ##-------------------------------------------##  
    def copy_sop(self, slice=None, ids=None, reindex=False, offset=(0,1,0), rot=(0,0,0), num=2, distance=2):
        """ UNFINISHED ,  mimmic the copy SOP in Houdini 
             
            offset normal per face would be slick           
        """


        pids = self.indexer( span=slice, ids=ids) 

        geom     = self.sub_select_geom( ids=pids, reindex=True )
        tmpnrmls = self.get_face_normal(fid=pids, unitlen=True) 

        #print("#### DEBUG ", tmpnrmls , ids )

        for i in range(num):
            for j in range(len(tmpnrmls)):
                
                # experimental transform on surface normal  
                f_nrml = tmpnrmls[j]*distance #normal vector * magnitude 
                amtx = f_nrml[0]
                amty = f_nrml[1]
                amtz = f_nrml[2]
 
                # absolute transform in world coordinates 
                # amtx = offset[0]
                # amty = offset[1]
                # amtz = offset[2]

                ox = amtx * i  
                oy = amty * i
                oz = amtz * i

                newpts = self.xform_pts((ox,oy,oz), geom[1] )
 
                ############# 
                # DEBUG - this seems not right, grinds to a halt on 20+ polygons
                self.insert_polygons(geom[0], newpts  ) 

    ##-------------------------------------------##  
    def radial_triangulate_face(self, fid, offset=None, as_new_obj=False ):
        """ put a vertex at the center of polygon 
            then form triangles in a circle 
            for N sided polygons 

 
            as_new_obj - replace object OR append to it 
            offset     - optional spatial offset for new center point 
                         added so I could turn a circle into an arrow :)   
        """
        out_polys = []
        out_pts   = []

        #if fid >= len(self.polygons):
        #    print('## error radial_triangulate_face - bad face index ')
        #    return None 

        tmp = self.get_face_geom(fid)
        poly = tmp[0][0]

        fac_pts = []
        for ptidx in poly:

            # build a list of points that make up polygon
            fac_pts.append(self.points[ptidx-1]) #not zero indexed?

        # calculate the center of each polygon
        # this will be added as a new point, center of radial mesh      
        fac_ctr = self.centroid(fac_pts)
        
        # offset is a spatial transform of the radial center point
        #DEBUG TODO - option to move along face normal!!  
        if offset is not None:
            fac_ctr[0] = fac_ctr[0]+offset[0]
            fac_ctr[1] = fac_ctr[1]+offset[1]
            fac_ctr[2] = fac_ctr[2]+offset[2]

        # start a new polygon dataset to append later 
        out_pts.append(tuple(fac_ctr))
        # add all the old points to our new dataset, plus our new center point   
        out_pts.extend(self.points)

        # iterate by two and connect to new radial center     
        for i in range(int(len(poly))):
            out_polys.append( (1, poly[i-1]+1, poly[i]+1 ) ) 
        
        if as_new_obj:
            self.points   = out_pts
            self.polygons = out_polys
        else:    
            self.insert_polygons(out_polys, out_pts)

    ##-------------------------------------------## 
    def radial_triangulate_obj(self, as_new_obj=False, offset=None, norm_mult=None ):
        """ 
            put a vertex at the center of polygon 
            then form triangles in a circle 
            for N sided polygons 

 
            as_new_obj - replace object OR append to it 
            offset     - optional spatial offset for new center point 
                         added so I could turn a circle into an arrow :)   
        """
        
        centers = [] 
        normals = [] 
       
        out_polys = []
        out_pts   = []

        for i,poly in enumerate(self.polygons):
            n = self.get_face_normal(i)
            if norm_mult is None:
                n = n.normal
            else:    
                n = n*norm_mult
              
            # calculate the center of each polygon
            # this will be added as a new point, center of radial mesh 
            fac_ctr = self.get_face_centroid(i)
            # offset is a spatial transform of the radial center point
            if offset is not None:
                if offset == 'normal':
                    fac_ctr[0] = fac_ctr[0]+n[0]
                    fac_ctr[1] = fac_ctr[1]+n[1]
                    fac_ctr[2] = fac_ctr[2]+n[2]
                elif offset == 'flipnormal':
                    fac_ctr[0] = fac_ctr[0]-n[0]
                    fac_ctr[1] = fac_ctr[1]-n[1]
                    fac_ctr[2] = fac_ctr[2]-n[2]
                else:    
                    fac_ctr[0] = fac_ctr[0]+offset[0]
                    fac_ctr[1] = fac_ctr[1]+offset[1]
                    fac_ctr[2] = fac_ctr[2]+offset[2]

            #centers.append(fac_ctr)
            #normals.append(n)
            # add new point to attach geom to 
            out_pts.append(fac_ctr)
            
            # iterate by two and connect to new radial center     
            for j in range(len(poly)):
                newply = (len(self.points)+i+1, poly[j-1] , poly[j] )
                out_polys.append( newply ) 

        if as_new_obj:
            #DEBUG - THIS OVERWRITES EACH TIME 
            #self.points = out_pts
            #self.polygons = out_polys
            pass 
        else:    
            #print(out_pts)
            #print(len(out_polys)) 
            self.insert_polygons(ans=False, plyids=out_polys, points=out_pts )


 
        #for c in centers:
        #    self.prim_locator(c) 
       



    ##-------------------------------------------## 
    def triangulate(self, force=False, offset=(0,0,0)):
        """ 
            Only works on 3 or 4 sided polygons. 3 are passed unchanged, 4 are triangulated  
            turn a quad into two triangles 
            return a new 3d object (possibly with more new polygons, all triangles) 
        """
        out_polys = []

        if force is True:
            self.radial_triangulate_obj(offset=offset)

        else: 
            for poly in self.polygons:
                num_vtx = len(poly)

                if num_vtx==3:
                    out_polys.append(poly)

                elif num_vtx==4:
                    v1=poly[0];v2=poly[1];
                    v3=poly[2];v4=poly[3]; 
                    out_polys.append( (v1,v3,v4) ) 
                    out_polys.append( (v1,v2,v3) )              
     
                else:
                    print('# experimental N sided triangulation for %s sided poly '%num_vtx)
                    self.radial_triangulate_obj(offset=(0,0,0))

                # overwrite old data 
                self.polygons = out_polys

    ##-------------------------------------------##  
    def poly_loft(self, obj2, as_new_obj=True):
        """ UNFINISHED 
            assume two profiles have been passed in with identical polygon ordering 
            connect them togther with new side wall polygons 
        """
        print("### debug begin loft op ")

        if len(self.polygons) != len(obj2.polygons):
            print("## debug ## objects must have the same num polys! ")
            return None 


        new_pts = []
        new_plys = [] 

        print("### poly counts for both objects %s %s"%( len(self.points), len(obj2.points) )  )
        for i,poly in enumerate(self.polygons):

            # iterate the first object 
            for pidx in poly:
                new_pts.append(self.points[pidx-1]) # first half of new poly quad

            # iterate the second object (should match first)             
            for pidx in obj2.polygons[i]:
                new_pts.append(obj2.points[pidx-1]) # second half of new poly quad
            
            num = len(new_pts)
            print( obj2.polygons[i] )
                #new_plys.append( (num-3, num-2, num ) )                   # assemble the new polygon indices 
                #new_plys.append( (1,2,3,4) )            # assemble the new polygon indices 
                #new_plys.append( (num-3, num-2, num ) ) # assemble the new polygon indices 

        
        print('#####################')
        #print(new_plys)
        #print(len(new_pts))

        if as_new_obj:
            self.points   = new_pts
            self.polygons = new_plys
        else:    
            self.insert_polygons(new_plys, new_pts)

    ##-------------------------------------------## 
    def slice_axis(self, low, high, axis='y'):
        """ extract a section of polygons from a model based on position in space 
            it should catch any polygon with a point within the range... I think 
        """

        sliced=[]
        exists = []

    
        def add_once(data):
            if data[0] not in exists:
                sliced.append(data)                
                exists.append(data[0])

        for fid,ply in enumerate(self.polygons):
            tmply=[]
            for pid in ply:
                pt = self.points[pid-1]
                tmply.append(pt)

            for a in tmply:
                if axis=='x':
                    if a[0]>low and a[0]< high: 
                        add_once([fid, [ply ,tmply]])
 
                if axis=='y':
                    if a[1]>low and a[1]< high:  
                        add_once([fid, [ply ,tmply]])
                          
                if axis=='z':
                    if a[2]>low and a[2]< high:  
                        add_once([fid,[ply ,tmply]])
       
        return sliced

    ##-------------------------------------------##   
    ##-------------------------------------------## 
    #file IO / mesh analysis, etc 

    def _repair(self):
        """ UNFINISHED 
            walk internal data and fix any bad data found  (empty point tuples, etc) 
        """

        fix = []

        for i,pt in enumerate(self.points):
            
            if pt is None:
                self._scribe("pt idx %s is None"%i )

            elif type(pt)==vec3:
                fix.append( (pt.x, pt.y, pt.z) )

            elif len(pt)==0:
                self._scribe('found bad data - empty vertex')
             
            #CASE     
            #elif len(pt)==0:
            #    self._scribe('found bad data ')
            
            #CASE
            #

            else:
                fix.append( pt )
        self.points = fix

    ##-------------------------------------------## 
    #load/dump numbered point caches and reload - very powerful idea!
    # def import_pts(self)

    # def import_faces(self)

    ##-------------------------------------------## 
    def load(self, filename, doflush=True):
        """ 
            DEBUG - DOES NOT CLEAR BUFFERS FIRST!!
            so if you load two models, the points - polygons will be merged and have bad topology

            ALSO NEEDS TO DO A FILE CHECK FOR VALID DATA 

            DOC: 
                load a wavefront OBJ file into model's point/poly memory 
                you can save shape or cache data 

                doflush clears out all memory 
                if you dont flush it will attempt to fuse existing geometry with loaded 
        """

        #if doflush is True:
        #    self._flush() 

        if os.path.lexists(filename) == 0:
            self._scribe("%s DOES NOT EXIST !! "%filename )
            #raise
            
        if os.path.lexists(filename):
            f = open( filename,"r", encoding='utf-8')
            contents = f.readlines()
            
            #if data exists it, FIDS need to be offset 
            fs = self.numpts 
            
            print('######num ', fs)
             
            for x in contents :
                #lines = x
                nonewline = x.split('\n')
                tok =  nonewline[0].split(" ") 
                if tok[0]!='#':
                    ###
                    #THIS NONSENSE IS TO CLEAN UP ERRANT SPACES IN FILE 
                    clndat = []
                    for f in tok:
                        if(f!='' ):
                            clndat.append(f) 

                    #weak attempt to clean up the textfile a little                    
                    tok=clndat    
                    numtok = len(tok)               
                    if tok: 
                        # VERTICIES 
                        if tok[0]=='v':
                            self.points.append( (float(tok[1]), float(tok[2]), float(tok[3]) ) ) 

                        # LINES
                        if tok[0]=='l':
                            ## LINE IMPORT IS UNTESTED !
                            fids = tok[1:] #remove the first item (letter f )
                            polyline = []
                            for fid in fids:
                                polyline.append(int(fid))   
                            self.polygons.append( polyline )                         

                        # FACES
                        if tok[0]=='f':
                            
                            fids = tok[1:] #remove the first item (letter f )
                            
                            poly    = []
                            uv_poly = []



                            for fid in fids:
                                
                                ## DEAL WITH THIS STUFF - '47//1'
                                if '/' in fid:
                                    
                                    # Vertex texture coordinate indices
                                    # f v1/vt1 v2/vt2 v3/vt3

                                    # Vertex normal indices
                                    # f v1/vt1/vn1 v2/vt2/vn2 v3/vt3/vn3

                                    # Vertex normal indices without texture coordinate indices
                                    # f v1//vn1 v2//vn2 v3//vn3

                                    tmp = fid.split('/')
                                    if len(tmp):
                                        # first slash delineated integer is face ID
                                        poly.append(int(tmp[0])+fs)    
                                        # second slash delineated integer is face UV
                                        if tmp[1]: 
                                            uv_poly.append(int(tmp[1])+fs)
                                        
                                        # third item is vertex normal 
                                        #if tmp[2]: 


                                else:    
                                    poly.append(int(fid)+fs)   

                            self.polygons.append( poly )
                            self.uv_polys.append(uv_poly) 
                             
                        # NORMALS
                        if tok[0]=='vn':
                            # debug - count the data being loaded and do sanity checks,
                            # for example, does the number of normals match the faces 
                            
                            #print('normal found     ', (tok[1],tok[2],tok[3]) )
                            
                            self.normals.append( ( float(tok[1]), float(tok[2]), float(tok[3]) ) )    

                        # UV's
                        if tok[0]=='vt':
                            #print('texture UV found ', tok)
                            
                            #print(tok)

                            if len(tok)==3: 
                                self.uv_points.append( (float(tok[1]), float(tok[2])) )  

                            """
                            v  0.000000 2.000000 0.000000
                            v  0.000000 0.000000 0.000000
                            v  2.000000 0.000000 0.000000
                            v  2.000000 2.000000 0.000000
                            vt 0.000000 1.000000 0.000000
                            vt 0.000000 0.000000 0.000000
                            vt 1.000000 0.000000 0.000000
                            vt 1.000000 1.000000 0.000000
                            # 4 vertices
                            usemtl wood
                            # The first number is the point,
                            # then the slash,
                            # and the second is the texture point
                            f 1/1 2/2 3/3 4/4
                            # 1 element
                            """

                        # What the hell are "Parameter space vertices"?                        
                        #if tok[0]=='vp':
                        #    print('Parameter space vertices found ', tok)
                               
        # DEBUG put a file integrity check here 
        # I had one with a bad point index 

    ##-------------------------------------------## 

    def save(self, filename, as_lines=False):
        """ format the data so blender (or anything else) can read it  
            this will save points and polygons as a wavefront OBJ file 
            you can save shape or cache data. 
        """    
        
        ##################
        self._repair()#optional but this will help find/fix problems later


        buf = [] #array of strings to be written out as the OBJ file

        buf.append("# Created by Magic Mirror render toy.")
        buf.append("# Keith Legg - December 2015.")        
        buf.append("# version2   - November 2018.\n")

        buf.append('\n# Define the vertices')

        #DEBUG - PUT MORE ERROR CHECKING ON VERTS, I HAD SOME BAD DATA GET THROUGH 
        #EX: - v 1 2 3) (4,5,6)
        for p in self.points:
            if len(p) == 2:
                buf.append('v %s %s %s'%( p[0], p[1], 0) ) #attempt to work with 2D data. BAD IDEA? 
            elif len(p) == 3:
                buf.append('v %s %s %s'%( p[0], p[1], p[2]) ) #x y z components 
            elif len(p) == 6:
                buf.append('v %s %s %s %s %s %s'%( p[0], p[1], p[2], p[3], p[4], p[5]) ) #x y z components                 
            else:
                print('## object save - bad vertex coordinate ', p )
                return None 

        buf.append('\n# Define the polygon geometry')
        buf.append('# No UV or normals at this time')
        
        for ply in self.polygons:
            plybuf = ''
            for f in ply:
                #DEBUG - hack to deal with tuple - probably a bad idea 
                if type(f)==tuple or type(f)==list:
                    print('#DEBUG tuple point index hack! ')
                    for idx in f:
                        plybuf = plybuf +(' %s'%str(int(idx)) ) 
                    print(plybuf)   
                else:    
                    #plybuf = plybuf +(' %s'%str(int(f)+1) ) #add one because OBJ is NOT zero indexed
                    plybuf = plybuf +(' %s'%str(int(f)) ) 

            if as_lines:
                # save as lines
                buf.append('l %s'%plybuf)
            else:
                # save as polygons 
                # format for f command is : f position_id/texture_coordinates_id/normal_id 
                # for now we leave off slashes and leave other two blank   
                buf.append('f %s'%plybuf)
 
        buf.append('\n')

        ################################### 
        #Our filebuffer is an array, we need a string so flatten it 
        output = ''
        for s in buf:
            output=output+s+'\n'

        self._scribe('### file "%s" saved' % filename)

        #save it to disk now
        fobj = open( filename,"w") #encoding='utf-8'
        fobj.write(output)
        fobj.close()



