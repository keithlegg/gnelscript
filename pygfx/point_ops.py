#!/usr/local/bin/python3

import os
import math 





from pygfx.math_ops import math_util as mu
from pygfx.math_ops import NUMPY_IS_LOADED, matrix33, matrix44, vec2, vec3  




if NUMPY_IS_LOADED:
    # print(' ## debug - loading numpy module in point ops. ')
    import numpy as np  
else:
    print(' ## debug - numpy module disabled in point ops. ')




class point_operator(object):
    """ what became of the original point generator 
        deal with raw points, not geom, ptgrps, facgrps. JUST points 
    """

    def __init__(self):
        self.mu   = mu()
        self.m33  = matrix33()      
        self.m44  = matrix44()  # used for point rotations, and more? 
        self.vec2     = vec2()     
        self.vec3     = vec3()      

    def tuple_pop(self, listTuples, tup):
        """ take a list of tuples, remove one by filtering it out and return the rest back """
        out = []
        for t in listTuples:
            if(t != tup):
                out.append(t)
        return out
            
    def add_margin_bbox(self, bbox, size):
        """ return center (x,y) from two diagonal coordinates """
        
        out = []
        out.append( bbox[0]-size  ) 
        out.append( bbox[1]-size  )
        out.append( bbox[2]+size  )
        out.append( bbox[3]+size  )
        return out

    def extents_fr_bbox(self, bbox, offset=None):
        """ return center (x,y) from a tuple of 4 numbers (PIL bbox [left, top, right, bottom]) 
            offset (int) adds a margin to the size of the page edges 
        """
        
        out = []
        if not offset:
            out.append( (bbox[0], bbox[1]) ) #top left
            out.append( (bbox[2], bbox[1]) ) #top right
            out.append( (bbox[2], bbox[3]) ) #bottom right
            out.append( (bbox[0], bbox[3]) ) #bottom left
        
        #increase the margins, (pixels are indexed with top left at 0,0)
        if offset:
            out.append( (bbox[0]-offset , bbox[1]-offset ) ) #top left
            out.append( (bbox[2]-offset , bbox[1]+offset ) ) #top right
            out.append( (bbox[2]+offset , bbox[3]+offset ) ) #bottom right
            out.append( (bbox[0]-offset , bbox[3]-offset ) ) #bottom left

        return out

    def closest_to_axis(self, points, val, axis):
        """ return the closest point to X OR Y """
        
        tmpsort = []

        for pt in points:
            if axis=='x':
                sp = ( val ,pt[1]  )
            if axis=='y':
                sp = ( pt[0]  ,val )

            tmpsort.append( (self.calc_line_length(pt[0], pt[1], sp[0], sp[1] ), pt) )
        
        ###
        tmpsort.sort()
        return tmpsort[0][1]

    ############################################### 
    def calc_square_diag(self, tl, br):
        """ creates 4 coordinates representing an extent box from a diagonal coordinate pair """ 
        out =[]  
        out.append( (tl[0], tl[1])  ) #tl
        out.append( (br[0], tl[1])  ) #tr
        out.append( (br[0], br[1])  ) #br
        out.append( (tl[0], br[1])  ) #bl

        return out

    ############################################### 
    def calc_square_pt(self, size, origin=None ):
        """ UNTESTED 
            DEBUG - this is unclamped  
            calc the XY coodrinates for a square  
            (top left, top right , bottom right, bottom left ) 
        """
        
        out =[]  

        if origin:
            out.append( ( origin[0]-size, origin[1]+size)  ) #tl
            out.append( ( origin[0]+size, origin[1]+size)  ) #tr
            out.append( ( origin[0]+size, origin[1]-size)  ) #br
            out.append( ( origin[0]-size, origin[1]-size)  ) #bl
        else:
            out.append( ( -size,  size)  )
            out.append( (  size,  size)  )
            out.append( (  size, -size)  )
            out.append( ( -size, -size)  )
        return out

    ###############################################     
    def calc_bbox_pt(self, size, origin=None ):
        """ DEBUG - this is unclamped!
            calc extents for a square (left, upper, right, and lower)
        """
                
        out =[]  
        
        if origin:
            out.append(  origin[0]-size  ) #west
            out.append(  origin[1]-size  ) #north
            out.append(  origin[0]+size  ) #east
            out.append(  origin[1]+size  ) #south
        else:
            out.append(  -size  ) #west
            out.append(  -size  ) #north
            out.append(  size  ) #east
            out.append(  size  ) #south
        return out

    ############################################### 
    def calc_circle(self, pos=(0,0,0), dia=1, axis='z', periodic=True, spokes=23):
        """ spokes = num spokes """

        px=0
        py=0
        pz=0 

        ox = pos[0]
        oy = pos[1]
        oz = pos[2]

        out = []
        
        degs = 360/spokes
        dit = 0
        
        dtr = self.mu.dtr

        for i in range(spokes):
            if axis=='x':
                py = oy + (math.sin(dtr(dit))*dia) 
                pz = oz + (math.cos(dtr(dit))*dia) 
                out.append( (px, py, pz))
            if axis=='y':
                px = ox + (math.sin(dtr(dit))*dia) 
                pz = oz + (math.cos(dtr(dit))*dia) 
                out.append( (px, py, pz))
            if axis=='z':
                px = ox + (math.sin(dtr(dit))*dia) 
                py = oy + (math.cos(dtr(dit))*dia) 
                out.append( (px, py, pz))                
            dit+=degs
        
        if periodic:
            out.append( out[0] )

        return out
 
    ############################################### 
    def sort_3_distances(self, mode, coords):
        """ take 3 XY points (triangle) and get the distances between all 3 
            return the sorted distances represented as two coordinate pairs
        """
        out = []
        tmp = []

        tmp.append((self.calc_line_length(coords[0][0], coords[0][1], coords[1][0], coords[1][1]), coords[0], coords[1]))
        tmp.append((self.calc_line_length(coords[1][0], coords[1][1], coords[2][0], coords[2][1]), coords[2], coords[1]))
        tmp.append((self.calc_line_length(coords[2][0], coords[2][1], coords[0][0], coords[0][1]), coords[2], coords[0]))

        tmp.sort()

        if mode=='shortest':
            out.append(tmp[0][1]);out.append(tmp[0][2])

        if mode=='middle':
            out.append(tmp[1][1]);out.append(tmp[1][2])

        if mode=='longest':
            out.append(tmp[2][1]);out.append( tmp[2][2])

        #we have the data but we need to keep the order intact 
        newout = [0,0,0]#create a work buffer with 3 elements
        
        for x in out:
            count = 0
            for f in coords:
                if x == f:
                    newout[count]=x
                count += 1
        
        #remove empty elements
        tmp = []
        for y in newout:
            if y:
                tmp.append(y)
        return tmp


    ############################################### 
    def get_mean_z(self, triangle):
        z1 = triangle[0][2]
        z2 = triangle[1][2]
        z3 = triangle[2][2]
        return (z1+z2+z3)/3

    ###############################################  
    def cvt_2d_to_3d(self, points):
        """ convert 2d points into 3d by adding an empty z axis """

        newpts = []
        for pt in points:
            newpts.append( (pt[0], pt[1], 0)   )
        return newpts


    ############################################### 
    def locate_pt_along3d(self, x1, y1, z1, x2, y2, z2, num):
        """
            given two 3D points, return a series of N number connecting points in 3D 
        """

        pts_created = []
        fpos=(x1, y1, z1)
        spos=(x2, y2, z2)
         
        for n in range(num):
            npos = [3]

            npos[0]     = spos[0]+(((fpos[0]-spos[0])/(num+1))*(n+1))
            npos.append(  spos[1]+(((fpos[1]-spos[1])/(num+1))*(n+1))  )
            npos.append(  spos[2]+(((fpos[2]-spos[2])/(num+1))*(n+1))  )

            pts_created.append( (npos[0], npos[1], npos[2]) )
        return pts_created

###############################################
class polygon_operator(point_operator):
    """ polygon operator - should be called GEOM operator 

        3D model and tools to work on polygons 
       

        -------------------------------------------
            Concepts 
 
        SUBSELECT        = slice (start, end ), ids = [LIST OF PTS]

        GEOM             = [ [], [] ]
        GEOM             = [ [(poly),(poly),..], [(pt),(pt),(pt),(pt),... ] ]
        GEOM (extended)  = [ [(poly),(poly),..], [(pt),(pt),(pt),(pt),... ] ]

        PTGRP            = [ [ID,(PT)], [ID,(PT)] ]
        PTGRP (extended) = [ [ID,(PT)], [ID,(PT)] ,   [COLOR], [NORMAL], [META] ]        

        FACEGRP = [ [ID,(POLY)], [ID,(POLY)], .. ]

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
        self.points          = []    # list of tuples of XYZ points per vertex         -  [(x,y,z), (x,y,z)]  
        self.polygons        = []    # list of tuples for 2 or more vertex connections -  [(1,2,5,8) , (1,2)] 
        self.face_normals    = []
        #self.point_normals  = []
        #self.point_colors   = []
        #self.line_colors    = []

        # render properties embedded in geometry 
        self.linecolors = None       # iterable of colors for lines 
        self.linecolor  = (0,240,00)
        self.vtxcolor   = (0,240,00)

        # "unavoidable side effect" variables 
        self.exprt_ply_idx   = 1     # obj is NOT zero indexed
        #self.exprt_pnt_idx   = 0    # pts ARE zero indexed (everything BUT .obj face idx's are)


    def flush(self):
        """ set all geometry to a clean state """

        self.points       = [] 
        self.polygons     = []      
        self.face_normals = []  

        self.exprt_ply_idx   = 1 #obj is NOT zero indexed
        #self.exprt_pnt_idx   = 0 #pts ARE zero indexed (everything BUT .obj face idx's are)


    ############################################### 
    #def grow_selection(self, ptgrp, facgrp, u_num, v_num):

    ############################################### 
    #def weld_edges(self, obj):

    ############################################### 
    #def scan_shells(self, obj):
    #    """ look for all the chunks of geometry that are not connected """

    ###############################################  
    #def poly_seperate(self, obj):
    #    """ check for geometry that is not connected, if any found, break it off """

    ###############################################  
    #def get_edge_centroid(self, f_id , e_id):
    #    pass
    
    ############################################### 
    #def get_obj_centroid(self):
    #    pass

    ###############################################  
    ###############################################  
    ###############################################  
    # stupid little helpers, knick knacks, tin toys 

    def scribe(self, str):
        print(str)

    ###############################################     
    def calc_bbox(self, prgrp=None, facgrp=None ):
        """ UNFINISHED  
            get the boudning area of an object or face(s)
        """
        maxx = 0
        maxy = 0
        maxz = 0
        
        for p in self.points:
            print(p)

    ############################################### 
    def _reindex_ply(self, f_idx, offset):
        """ UNTESTED 
            take a tuple of face indexes and re-index to offset+value
            not used yet, but may be a good thing to have 
        """ 

        out_face = [] 
        for i in f_idx:
           out_face.append(i+offset)
        return tuple(out_face)

    ###############################################
    def three_vec3_to_normal(self, v1, v2, v3, unitlen=False):
        """ take 3 vec3 objects and return a face normal """

        # calculate the face normal  
        a = v1 - v2;b = v1 - v3;
        if unitlen:
            f_nrml = a.cross(b).normal
        else:    
            f_nrml = a.cross(b)          
        
        return f_nrml 

    ############################################### 
    def sub_select(self, slice=None, ids=None):
        """ interface to get a list of IDS for whatever
        """
        pids = []
  
        # block the same face ID from being selected twice
        # there may be cases where you WANT it twice 
        optimize = True  

        ####
        # insert list of individual ids first 
        if ids:
            if optimize:
                for i in ids:
                    if i not in pids:
                        pids.append(i)
            else: 
                pids.extend(ids)            
        #### 
        # then do the slice
        if slice:
            tids = []
            # insert slice IDs (start-end range) 
            for i in range(slice[0], slice[1]+1):
                if optimize:
                    if i not in pids:
                        tids.append(i)
                else: 
                    tids.append(i)
      
            pids.extend(tids)

        return pids 

    ############################################### 
    def geom_to_ptgrp(self, geom):
        out = []
        for i,g in enumerate(geom):
            out.append( [1, g[1][i]] )
        return out     

    ###############################################  
    def sub_select_geom(self, slice=None, ids=None, reindex=False):
        """ 
            make work with xform_pts, rotate_pts, scale_pts 

            # DEBUG - This works for now, but... 
            # the only way to do this right is make get_face_geom() smarter
            # it needs to understand slicing, and not select the same points twice while iterating faces 
            # for example, do a subselect on a cube IDS 1-6 to see what I mean 
            # it will give you 24 points instead of 8 


            slice - tuple of (start,end)  
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
        self.exprt_ply_idx = 1

        pids = self.sub_select(slice=slice, ids=ids)

        #print('## debug, test of optimize  ', pids) 

        for i in pids:
         
            # DEBUG - This works for now, but... 
            # the only way to do this right is make get_face_geom() smarter
            # it needs to have slicing and not select the same thing twice 
            geom = self.get_face_geom(i, reindex=reindex)
            out_poly.append(geom[0][0])
            for pt in geom[1]:
                out_pts.append(pt)               

        return ( out_poly, out_pts )

    ###############################################  
    def get_pt_ids(self, fids=None):
        """ lookup point indices from a listf of face IDs 
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

    ###############################################  
    def get_face_pts(self, fid):
        """ lookup and return the point geometry of a face from a face ID
            for fancier features look at get_face_geom() 
            
        """

        tmp = []

        if fid<0 or fid > len(self.polygons)-1:
            print('# show_poly- bad face index : %s'%fid)
            return None

        for v_id in self.polygons[fid]:
            tmp.append(self.points[v_id-1])

        return tmp

    ###############################################  
    def insert_pt_grp(self, ptgrp):
        """ ptgrp is [ [id, pt], [id, pt] ]
            this will overwrite internal geometry with point group geometry 

        """

        out = []
        for p in ptgrp:
            self.points[p[0]] = p[1]
        return out    

    ###############################################  
    def get_pt_grp(self, slice=None, ids=None):
        """ gets a point group, 
            a point group is a list of  

            [ [ID, (X,Y,Z)], ... ] 
            
            we can process them and put them back where 
            we found them in the model, allowing for modification of partial or whole objects 


        """

        out = []

        if slice is None and ids is None:
            for i,p in enumerate( self.points ):
                out.append([i,p])  
            return out

        else:
            pids = self.sub_select( slice=slice, ids=ids)

            for p in pids:
                out.append( [p, self.points[p]] )
            return out      

    ###############################################  
    def get_face_group(self, slice=None, ids=None):
        """ get a face group, a list of faces and IDS so 
            we can process them and put them back 

            data format [ [ID, face] ]

        """
        fids = self.sub_select( slice=slice, ids=ids)


        pass

    ###############################################  
    def get_geom_edges(self, geom ):
        """ only works on a single polygon geom object 
            DEBUG - make work with multiple 
        """

        out_edge_ids = []
        out_edge_pts = []
         
        ply = geom[0][0]

        for idx in ply:
            #print('### POLY DATA   ', ply , type(ply) )
            #print('### POINTS LEN  ', ply , type(geom[1]) )

            # iterate by two and store segments
            out_edge_ids.append((  ply[idx-2]           ,ply[idx-1]              )) # poly index
            out_edge_pts.append((  geom[1][ply[idx-1]-2], geom[1][ply[idx-1]-1]  )) # point index

        return [out_edge_ids, out_edge_pts]

    ###############################################  
    def get_face_edges(self, fid, reindex=False, geom=None):
        """ UNTESTED 
            return [[VTX_IDS], [VTX_PTS]]
        """
        if geom is None:
            self.exprt_ply_idx = 1
            geom = self.get_face_geom(fid, reindex=True)  # [[poly idx], [pt data]] 
        return self.get_geom_edges(geom)

    ###############################################  
    def get_face_edges2(self, fid, reindex=False):
        """ UNTESTED 
            return [[VTX_IDS], [VTX_PTS]]
        """
        tmp = self.get_face_geom(fid)  # [poly idx, pt data] 

        out_edge_pts = []
        out_edge_ids = []

        poly = tmp[0]

        # iterate by two and connect to new radial center   
        # thanks to pythons negative index, this works a treat  
        for i in range(len(poly)):
            out_edge_ids.append( (poly[i-1], poly[i] ) ) 
            out_edge_pts.append( (self.points[poly[i-1]-1], self.points[poly[i]-1] ) ) 

        return [out_edge_ids, out_edge_pts]

    ############################################### 
    def insert_polygons(self, plyids, points, asnew_shell=True, geom=None):
        """  
             Insert NEW geometry into this object

             asnew_shell  - reindex the points and append, else keep the same indices
             geom         - geom to insert into, instead of object.polygons, object.points
                            if true, will return the geom object when done 

        """

        #if isinstance(points, vec3):
        #    self.points.extend(points)
        
        ######### 
        # append polygons 
        for poly in plyids:
            plytmp = []      
            for idx in poly:
                if not isinstance(idx, int):
                    print('## insert_polygons, bad data for index ')
                    return None 

                if asnew_shell is True:
                    plytmp.append(idx+self.numpts) # add the poly index to current count    
                else:
                    plytmp.append(idx)  

            # do the insert operation                    
            if geom is None: 
                self.polygons.append( tuple(plytmp) ) 
            else:
                geom[0].append( tuple(plytmp) )

        #########        
        # append points,  only if new geom - just use python extend 
        if asnew_shell is True:
            if isinstance(points, tuple) or isinstance(points, list):
                # do the insert operation
                if geom is None:
                    self.points.extend(points)
                else:
                    geom[1].extend(points)
        
        if geom is not None:
            return geom

    ###############################################  
    def extrude_face(self, f_id, distance):
 
        geom  = self.sub_select_geom(ids=[f_id] , reindex=True)
        nrmls = self.get_face_normal(fid=f_id, unitlen=True) 

        nrmls = nrmls * distance 
        #s_edges = self.get_face_edges(f_id) 
        s_edges = self.get_geom_edges(geom)  
        moved = self.xform_pts( nrmls, pts=geom[1])
        e_edges = self.get_geom_edges([geom[0],moved]) 

        # wall polygons 
        # iterate one set of edges assuming they both have the same number 
        for w in e_edges[0]:
            wall_poly = []
            wall_poly.extend(s_edges[1][w[0]-1]) # bottom half of quad polygon 
            wall_poly.extend(e_edges[1][w[0]-1]) # top half of quad polygon                  
            self.insert_polygons( [(1,2,4,3)], wall_poly, asnew_shell=True) 

        # transformed face along normal (cap polygon) 
        self.insert_polygons(geom[0], moved, asnew_shell=True) 

    ###############################################  
    def select_by_location(self, reindex=False):
        """ UNFINISHED
            select by angle 
            direction to other things 
            select by distance to other objects, points , etc 

        """
        pass

    ###############################################  
    def copy_sop(self, slice=None, ids=None, reindex=False, offset=(0,1,0), rot=(0,0,0), num=2, distance=2):
        """ UNFINISHED ,  mimmic the copy SOP in Houdini 
             
            offset normal per face would be slick           
        """


        pids = self.sub_select( slice=slice, ids=ids) 

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

    ###############################################  
    def inspect(self, geom):
        """ analyze a GEOM object and see how it is constructed """

        if geom == None:
            print('## inspect: error no geometry ')
            return None 

        print('## geom has %s top level items '%len(geom) )
        print('## geom has %s polygons        '%len(geom[0]) )
        print('## geom has %s points          '%len(geom[1]) )

        for i, poly in enumerate(geom[0]):
            print (' poly %s is type %s and has %s items '% (i,type(poly), len(poly)) )

    ###############################################  
    def verify(self, geom):
        """ check that indices are within range for data 
            
            WELL FORMED GEOM IS: 

            [ 
                [(FIDS), (FIDS)], 
                [(PT),(PT),(PT)] 
            ]
 
        """
        
        verbose = True

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

    ###############################################  
    def get_face_geom(self, fid,  reindex=False, geom=None):
        """ lookup and return the polygon indices and points for a single polygon 

            reindex - if True  - renumber the new polygon indices startring at 1, 
                      if False - retain the oringial numbering 

            geom - act on a geom obj passed in, or on self
        """

        tmp_pts = []
       
        
        # validate inputs 
        #if fid<1 or fid > self.numply:
        #    print('# get_face_geom- bad face index : %s'%fid)
        #    return None

        # decide what the input is, fallback on self.poly/self.points  
        if geom is None:
            polygr  = self.polygons
            pointgr = self.points
        else:
            polygr  = geom[0]  #faces and polys common indexed    
            pointgr = geom[1]        

        # get the data we want         
        reindex_id = [] 

        # print('## $$$$ fid %s data %s '% ( fid,   polygr[fid-1]) ) 
        if self.verify([polygr,pointgr]) is False:
            return None 

        for v_id in polygr[fid-1]:
            # keep a count of points stored to use as new index
            reindex_id.append(int(self.exprt_ply_idx ))
            # store points that are indexed in geom 
            tmp_pts.append(pointgr[v_id-1]) #data is NOT zero index but all else IS 
            self.exprt_ply_idx +=1

        # geom is always [ [(poly),..], [(point),(point),...]  ]
        if reindex is False:
            return [[polygr[fid-1]]     , tmp_pts]
        if reindex is True:
            return [[tuple(reindex_id)] , tmp_pts]

    ###############################################  
    def get_face_normal(self, fid=None, unitlen=False ):
        """ lookup a face(s) and calculate a face normal(s) for it  
            only tested for 3 or 4 sided polygon 
            also returns the center position of a face

            returns vec3 type 

            DEBUG - the need for a standardized interface for slice, fid lookup is very apparent 
        """

        if fid == None:
            print("## error - need face id(s) to get normal")
            return None 

        if isinstance(fid, int):
            fids = [fid]
        if isinstance(fid,list):
            fids = fid
 
        out = []

        for f in fids:    
            # create a vec3 for each vertex (3 or 4 sided polys)
            v1=vec3();v2=vec3()
            v3=vec3();v4=vec3()

            tmp = self.get_face_geom(f) #returns [fidx, pts] 
            poly = tmp[0][0] #poly = face indices  

            v1.insert( self.points[poly[0]-1] )
            v2.insert( self.points[poly[1]-1] )
            v3.insert( self.points[poly[2]-1] )                
            out.append( self.three_vec3_to_normal(v1, v2, v3, unitlen=unitlen) )

        if isinstance(fid, int):
            return out[0]
        else:
            return out   

    ###############################################        
    def get_face_centroid(self, fid):
        pts = self.get_face_pts(fid)
        return self.poly_centroid(pts) 

    ###############################################   
    def poly_centroid(self, pts):
        """ get 3D center of object (average XYZ point) """

        ptsx = []
        ptsy = []
        ptsz = []
        for pt in pts:
            ptsx.append(pt[0])
            ptsy.append(pt[1])
            ptsz.append(pt[2])            
        #average them 
        x = sum(ptsx)/len(ptsx)
        y = sum(ptsy)/len(ptsy)
        z = sum(ptsz)/len(ptsz)
        return [x,y,z]

    ###############################################  
    def z_sort(self, reverse=False):
        """ return a sorted list of polygons by average Z coordinate 
            lowest z value is (closest?) to camera 
        """
        out = [];tmp = []

        #build a list of [mean z, [points, polygon] ]
        #sort by mean z value, and return the sorted list 
        for p in self.polygons:
            tri = []
            for idx in p:
                tri.append(self.points[idx-1])
            mean = self.get_mean_z(tri)
            tmp.append( (mean, p) )
        
        if reverse:
            tmp.sort(reverse=True)
        else:
            tmp.sort()

        for t in tmp:
            out.append(t[1])    
        
        self.polygons = out
        #return out

    ############################################### 
    def apply_matrix_pts(self, pts, m33=None, m44=None):
        """ batch mutliply points by a matrix 
            used for translate, rotate, and scaling. 
            
            Can be used for many other experiments 

        """
      
        tmp_buffer = [] 

        # apply the transform here
        for pt in pts:  

            if m33 is not None:
                tmp_buffer.append( m33 * pt )
            if m44 is not None:
                tmp_buffer.append( m44 * pt )

        return tmp_buffer
    
    ############################################### 
    def apply_matrix_ptgrp(self, ptgrp, m33=None, m44=None):
        """ batch mutliply a point group by a matrix 
            used for translate, rotate, and scaling. 
            
            Can be used for many other experiments 

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

    ###############################################  
    def scale_pts(self, amt, pts=None, ptgrp=None ):

        # build a scale matrix 

        sc_m33 = self.m33.identity
        sc_m33[0]  = amt
        sc_m33[4]  = amt
        sc_m33[8]  = amt    

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

    ###############################################  
    def rotate_pts(self, rot, pts=None, ptgrp=None):
        #self.repair() # may fix or find problems 

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
        if pts is not None and ptgrp is None: 
            return self.apply_matrix_pts(pts,  m44=rot_matrix)  

        if pts is None and ptgrp is None:
            ptgrp = self.get_pt_grp()    

        rotated = self.apply_matrix_ptgrp(ptgrp, m44=rot_matrix) 
        self.insert_pt_grp(rotated)

    ############################################### 
    def xform_pts(self, pos, pts=None, ptgrp=None):
        """ shift points without using a matrix 
            if no points are specified - apply to whole object 
        """

        ################################################
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
            

        self.insert_pt_grp(tmp_buffer)

    ############################################### 
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
        fac_ctr = self.poly_centroid(fac_pts)
        
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

    ############################################### 
    def radial_triangulate_obj(self, as_new_obj=False, offset=None ):
        """ put a vertex at the center of polygon 
            then form triangles in a circle 
            for N sided polygons 

 
            as_new_obj - replace object OR append to it 
            offset     - optional spatial offset for new center point 
                         added so I could turn a circle into an arrow :)   
        """
        

        out_polys = []
        out_pts   = []

        for poly in self.polygons:
            fac_pts = []
            for ptidx in poly:

                # build a list of points that make up polygon
                fac_pts.append(self.points[ptidx-1]) #not zero indexed?

            # calculate the center of each polygon
            # this will be added as a new point, center of radial mesh      
            fac_ctr = self.poly_centroid(fac_pts)
            
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
            self.points = out_pts
            self.polygons = out_polys
        else:    
            self.insert_polygons(out_polys, out_pts)

    ############################################### 
    def triangulate(self, force=False, offset=(0,0,0)):
        """ 
            Only works on 3 or 4 sided polygons. 3 are passed unchanged, 4 are triangulated  
            turn a quad into two triangles 
            return a new 3d object (possibly with more new polygons, all triangles) 
        """
        out_polys = []

        #print("### DEBUG TRIANGULATE CALLED ", len(self.polygons) )

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

    ############################################### 
    def poly_loft(self, obj2, as_new_obj=True):
        """ assume two profiles have been passed in with identical polygon ordering 
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


    ###############################################  
    ############################################### 
    ############################################### 
    #file IO / mesh analysis, etc 

    def repair(self):
        """ UNFINISHED 
            walk internal data and fix any bad data found  (empty point tuples, etc) 
        """

        fix = []

        for pt in self.points:

            #iftype ==  <class 'tuple'>
            if len(pt)==0:
                self.scribe('found bad data - empty vertex')
             
            #CASE     
            #elif len(pt)==0:
            #    self.scribe('found bad data ')
            
            #CASE
            #

            else:
                fix.append( pt )
        self.points = fix

    ############################################### 

    #load/dump numbered point caches and reload - very powerful idea!
    def load(self, filename, doflush=True):
        """ 
            DEBUG - DOES NOT CLEAR BUFFERS FIRST!!
            so if you load two models, the points - polygons will be merged and have bad topology

            load a wavefront OBJ file into model's point/poly memory 
            you can save shape or cache data 


            doflush clears out all memory 
            if you dont flush it will attempt to fuse existing geometry with loaded 
        """

        #if doflush is True:
        #    self.flush() 

        if os.path.lexists(filename) == 0:
            self.scribe("%s DOES NOT EXIST !! "%filename )
            #raise
            
        if os.path.lexists(filename):
            f = open( filename,"r", encoding='utf-8')
            contents = f.readlines()
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
                    if len(clndat)==5 or len(clndat)==4:
                        tok=clndat    
                        numtok = len(tok)               
                        #THIS NONSENSE IS TO CLEAN UP ERRANT SPACES IN FILE 
  
                        #VERTICIES 
                        if tok[0]=='v':
                            self.points.append( (float(tok[1]), float(tok[2]), float(tok[3]) ) ) 

                        #FACES
                        if tok[0]=='f':
                           # print('face found ', tok[1],tok[2],tok[3],tok[4])
                           b1 = tok[1].split('/')
                           b2 = tok[2].split('/')
                           b3 = tok[3].split('/')

                           #four sided polygons   
                           if numtok==4:
                               poly = ( int(b1[0]),  int(b2[0]), int(b3[0])  )
                               
                           #three sided polygons      
                           if numtok==5:
                               b4 = tok[4].split('/')
                               poly = ( int(b1[0]),  int(b2[0]), int(b3[0]), int(b4[0]) )


                           self.polygons.append( poly )

                #NORMALS 
                #if tok[0]=='vn':
                #    print('normal found ', tok)
                
                #UV's
                #if tok[0]=='vt':
                #    print('texture UV found ', tok)

    ###############################################  

    def save(self, filename, as_lines=False):
        """ format the data so blender (or anything else) can read it  
            this will save points and polygons as a wavefront OBJ file 
            you can save shape or cache data. 
        """    
        
        ##################
        self.repair()#optional but this will help find/fix problems later


        buf = [] #array of strings to be written out as the OBJ file

        buf.append("# Created by Magic Mirror render toy.")
        buf.append("# Keith Legg - December 2015.")        
        buf.append("# version2   - November 2018.\n")

        buf.append('\n# Define the vertices')

        #DEBUG - PUT MORE ERROR CHECKING ON VERTS, I HAD SOME BAD DATA GET THROUGH 
        #EX: - v 1 2 3) (4,5,6)
        for p in self.points:
            if len(p) != 3:
                print('## object save - bad vertex coordinate ', p )
                return None 

            buf.append('v %s %s %s'%( p[0], p[1], p[2]) ) #x y z components 
        
        buf.append('\n# Define the polygon geometry')
        buf.append('# No UV or normals at this time')
        for ply in self.polygons:
            plybuf = ''
            for f in ply:
                #plybuf = plybuf +(' %s'%str(int(f)+1) ) #add one because OBJ is NOT zero indexed
                plybuf = plybuf +(' %s'%str(int(f)) ) #add one because OBJ is NOT zero indexed

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

        self.scribe('###file "%s" saved' % filename)

        #save it to disk now
        fobj = open( filename,"w") #encoding='utf-8'
        fobj.write(output)
        fobj.close()


