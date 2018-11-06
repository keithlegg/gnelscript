#!/usr/local/bin/python3

import os
import math 


import numpy as np  


from pygfx.math_ops import math_util as mu
from pygfx.math_ops import vec2  
from pygfx.math_ops import vec3  
from pygfx.math_ops import matrix33 
from pygfx.math_ops import matrix44


class point_operator(object):
    """ what became of the original point generator 

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

    def calc_square_diag(self, tl, br):
        """ creates 4 coordinates representing an extent box from a diagonal coordinate pair """ 
        out =[]  
        out.append( (tl[0], tl[1])  ) #tl
        out.append( (br[0], tl[1])  ) #tr
        out.append( (br[0], br[1])  ) #br
        out.append( (tl[0], br[1])  ) #bl

        return out

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

    def calc_circle(self, origin=(0,0,0), dia=1, axis='z', periodic=True, spokes=23):
        """ spokes = num spokes """

        px=0
        py=0
        pz=0 

        ox = origin[0]
        oy = origin[1]
        oz = origin[2]

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
class polygon_operator(point_operator):
    """ polygon operator - 3D model and tools to work on polygons """

    def __init__(self):
        super().__init__()  

        self.points          = []    # list of tuples of XYZ points per vertex         -  [(x,y,z), (x,y,z)]  
        self.polygons        = []    # list of tuples for 2 or more vertex connections -  [(1,2,5,8) , (1,2)] 
        self.face_normals    = []

        self.exprt_ply_idx   = 1 #obj is NOT zero indexed
        self.exprt_pnt_idx   = 0 #pts ARE zero indexed (everything BUT .obj face idx's are)
        
        #self.point_normals  = []
        #self.point_colors   = []
        #self.line_colors    = []

        #render properties
        self.linecolors = None #iterable of colors for lines 
        self.linecolor  = (0,240,00)
        self.vtxcolor   = (0,240,00)


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
    def scribe(self, str):
        print(str)

    ############################################### 
    def _reindex_ply(self, f_idx, offset):
        """ take a tuple of face indexes and re-index to offset+value""" 

        out_face = [] 
        for i in f_idx:
           out_face.append(i+offset)
        return tuple(out_face)

    ############################################### 
    def insert_polygons(self, plyids, points, autoincrement=True):
        """  
             DEBUG IDEAS TO WORK ON : 
                 - single polygons, multiple polygons 
                 - tuple, list, vector data for points 
                 - add to existing, or replace it all  
                 - operate on an external array, or internal object polygons 

        """

        #if isinstance(points, vec3):
        #    self.points.extend(points)

        #########
        # add poly idx  - indexing and reindexing becomes a thing. 

        for poly in plyids:
            plytmp = []      
            for idx in poly:
                if autoincrement is True:
                    plytmp.append(idx+self.numpts) # add the poly index to current count    
                else:
                    plytmp.append(idx)  

            self.polygons.append( tuple(plytmp) ) 

        #########        
        # add points - easy peasy, just use extend 

        #if pass_array is None:
        if isinstance(points, tuple) or isinstance(points, list):
            self.points.extend(points)

    ############################################### 
    def get_mean_z(self, triangle):
        z1 = triangle[0][2]
        z2 = triangle[1][2]
        z3 = triangle[2][2]
        return (z1+z2+z3)/3

    ###############################################     
    def calc_bbox(self, object=None, fids=None ):
        """ UNFINISHED  
            get the boudning area of an object or face(s)
        """
        maxx = 0
        maxy = 0
        maxz = 0
        
        for p in self.points:
            print(p)

    ###############################################  
    def get_poly_geom(self, slice=None, ids=None, reindex=False):
        """ 
            slice - tuple of (start,end)  
            ids   - list of single ids 
            
            get one or more faces as a new object 
            specify a list of ids, or a range
        """
        out_poly = []
        out_pts  = []

        self.exprt_ply_idx = 1 #reset this when exporting with reindex 
        #self.exprt_pnt_idx = 1  

        if slice:
            # start-end id range 
            for i in range(slice[0], slice[1]):
                #print(i)
                tmp = self.get_face_data(i, reindex=reindex)
                out_poly.append(tmp[0])
                for pt in tmp[1]:
                    out_pts.append(pt)
    
        if ids:
            # list of specific ids 
            for i in ids:
                 tmp = self.get_face_data(i, reindex=reindex)
                 out_poly.append(tmp[0])
                 for pt in tmp[1]:
                    out_pts.append(pt)                

        return ( out_poly, out_pts )


    ###############################################  
    def get_face_pts(self, fid):
        """ lookup and return the points of a polygon in this object 
            for fancier features look at get_face_data() 
            
        """

        tmp = []

        if fid<0 or fid > len(self.polygons)-1:
            print('# show_poly- bad face index : %s'%fid)
            return None

        for v_id in self.polygons[fid]:
            tmp.append(self.points[v_id-1])

        return tmp

    ###############################################  
    def get_face_data(self, fid,  reindex=False):
        """ lookup and return the polygon indices and points or a single polygon 
            same as get_face_pts() , but this will get the indices and points

            reindex = reorder the indices, effectively making a new object 
                    if you reindex , set self.exprt_ply_idx to 1 first 
                    jankey, but it works 
        """

        tmp_pts = []

        if fid<0 or fid > len(self.polygons)-1:
            print('# show_poly- bad face index : %s'%fid)
            return None
        
        reindex_id = [] 

        for v_id in self.polygons[fid]:
            reindex_id.append(int(self.exprt_ply_idx ))
            tmp_pts.append(self.points[v_id-1])
            self.exprt_ply_idx +=1

        if reindex is False:
            return (self.polygons[fid], tmp_pts)
        if reindex is True:
            return (tuple(reindex_id), tmp_pts)

    ###############################################

    def three_vec3_to_normal(self, v1, v2, v3):
        """ take 3 vec3 objects and return a face normal """

        #secondary tweaks to the normal data 
        scale     = 1
        normalize = False  #make each face normal unit length 

        # calculate the face normal  
        a = v1 - v2;b = v1 - v3;
        if normalize:
            f_nrml = a.cross(b).normal*scale
        else:    
            f_nrml = a.cross(b)*scale         
        
        return f_nrml 

    ###############################################  
    def get_face_normal(self, fid=None):
        """ lookup a face and calulate a face normal for it  
            only tested for 3 or 4 sided polygon 
            also returns the center position of a face

            returns ???
        """

        if fid == None:
            print("## error - need a face id to get normal")
            return None 

        # create a vec3 for each vertex (3 or 4 sided polys)
        v1=vec3();v2=vec3()
        v3=vec3();v4=vec3()

        tmp = self.get_face_data(fid) #returns [fidx, pts] 
        f = tmp[0] #poly = face indices  

        v1.insert( self.points[f[0]-1] )
        v2.insert( self.points[f[1]-1] )
        v3.insert( self.points[f[2]-1] )                
        f_nrml = self.three_vec3_to_normal(v1, v2, v3)

        return f_nrml   


    ###############################################  
    def get_face_edges(self, fid, reindex=False):
        """ UNTESTED 
            return [[VTX_IDS], [VTX_PTS]]
        """
        tmp = self.get_face_data(fid)  # [poly idx, pt data] 

        out_edge_pts = []
        out_edge_ids = []

        poly = tmp[0]

        # iterate by two and connect to new radial center   
        # thanks to pythons negative index, this works a treat  
        for i in range(len(poly)):
            if reindex is False:
                out_edge_ids.append( (poly[i-1], poly[i] ) ) 
            if reindex is True:
                #out_edge_ids.append( (poly[i-1], poly[i] ) )            
                pass

            out_edge_pts.append( (self.points[poly[i-1]-1], self.points[poly[i]-1] ) ) 

        return [out_edge_ids, out_edge_pts]

    ###############################################        
    def get_face_centroid(self, fid):
        pts = self.get_face_pts(fid)
        return self.poly_centroid(pts) 

    ###############################################  
    def extrude_face(self, f_id):
        """ UNFINISHED """
        #edges = self.get_face_edges
        #for e in edge:
        #    build_poly(e)
        #ETC
        pass 

    ###############################################  
    """ 
    # replaced by poly_centroid
    def triangle_centroid(self, triangle):
        # all 3 x coordinates 
        x1 = triangle[0][0]
        x2 = triangle[1][0]
        x3 = triangle[2][0]  

        # all 3 y coordinates         
        y1 = triangle[0][1]
        y2 = triangle[1][1]
        y3 = triangle[2][1]
        
        # all 3 z coordinates 
        z1 = triangle[0][2]
        z2 = triangle[1][2]
        z3 = triangle[2][2]
        
        #average them 
        x= (x1+x2+x3)/3
        y= (y1+y2+y3)/3
        z= (z1+z2+z3)/3

        return [x,y,z]
    """  

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
    def apply_transforms(self):
        """ update the point data to reflect internal matricies """

        tmp_buffer = []
        for pvec in self.points:  
            tmp_buffer.append( self.m44 * pvec )
        self.points= tmp_buffer
        
    ###############################################  
    def move_pts(self, offset=(0,0,0)):
        """ transform POINTS not object - actually changes geometry"""
        shifted = []
        for pt in self.points:
            shifted.append( (pt[0]+offset[0], pt[1]+offset[1], pt[2]+offset[2] )  ) 
        self.points = shifted

    ###############################################  
    def scale_pts(self, offset=(0,0,0)):
        """ transform POINTS not object - actually changes geometry"""        
        shifted = []
        for pt in self.points:
            shifted.append( (pt[0]*offset[0], pt[1]*offset[1], pt[2]*offset[2] )  ) 
        self.points = shifted

    ###############################################  
    def rotate_pts(self, rot=(0,0,0) ):
        """ 
            transform POINTS not object - actually changes geometry        
            derived from pointgen2d method rotate_mat4 
            this simply rotates a 4X4 matrix to be attached to a 3D object
        """
        
        #self.repair()#may fix or find problems 

        rx=rot[0];ry=rot[1];rz=rot[2];

        dtr = self.mu.dtr

        #build rotationY (see diagram above) 
        y_matrix =  self.m44.identity
        y_matrix[0]  =  math.cos(dtr( ry ))
        y_matrix[2]  = -math.sin(dtr( ry ))
        y_matrix[8]  =  math.sin(dtr( ry ))
        y_matrix[10] =  math.cos(dtr( ry ))
              
        #build rotationZ (see diagram above) 
        z_matrix    =  self.m44.identity
        z_matrix[0] =  math.cos(dtr( rz ))
        z_matrix[1] =  math.sin(dtr( rz ))
        z_matrix[4] = -math.sin(dtr( rz ))
        z_matrix[5] =  math.cos(dtr( rz ))
        tmp_matr = y_matrix * z_matrix

        #build rotationX (see diagram above) 
        x_matrix =  self.m44.identity
        x_matrix[5]  =   math.cos(dtr( rx )) 
        x_matrix[6]  =   math.sin(dtr( rx )) 
        x_matrix[9]  =  -math.sin(dtr( rx ))
        x_matrix[10] =   math.cos(dtr( rx ))
        self.m44 = x_matrix * tmp_matr   
        
        self.apply_transforms()

    ############################################### 
    def xform_pts(self, amt=(0,0,0)):
        """ shift points without using a matrix """
        tmp = []
        for pt in self.points:  
            x = pt[0]+ amt[0]
            y = pt[1]+ amt[1]
            z = pt[2]+ amt[2]
            tmp.append( (x,y,z) )
        self.points = tmp

    ############################################### 
    def radial_triangulate_face(self, offset=None ):
        """ put a vertex at the center of polygon 
            then form triangles in a circle 
            for N sided polygons 

 
            as_new_obj - replace object OR append to it 
            offset     - optional spatial offset for new center point 
                         added so I could turn a circle into an arrow :)   
        """
        out_polys = []
        out_pts   = []

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
    def cvt_2d_to_3d(self, points):
        """ convert 2d points into 3d by adding an empty z axis """

        newpts = []
        for pt in points:
            newpts.append( (pt[0], pt[1], 0)   )
        return newpts

    ############################################### 
    def locate_pt_along3d(self, x1, y1, z1, x2, y2, z2, num):
        """
            overloaded method for 3D 
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

    #load/dump numbered point caches and reload - very powerful idea!
    def load_obj(self, filename):
        """ 
            DEBUG - DOES NOT CLEAR BUFFERS FIRST!!
            so if you load two models, the points - polygons will be merged and have bad topology

            load a wavefront OBJ file into model's point/poly memory 
            you can save shape or cache data 
        """

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

    def save_obj(self, filename, as_lines=False):
        """ format the data so blender (or anything else) can read it  
            this will save points and polygons as a wavefront OBJ file 
            you can save shape or cache data. 
        """    
        
        ##################
        self.repair()#optional but this will help find/fix problems later


        buf = [] #array of strings to be written out as the OBJ file

        buf.append("# Created by Magic Mirror render toy.")
        buf.append("# Keith Legg - December 2015.\n\n")        
        buf.append("# version2   - November 2018.\n\n")

        buf.append('\n#these define the vertecies')

        #DEBUG - PUT MORE ERROR CHECKING ON VERTS, I HAD SOME BAD DATA GET THROUGH 
        #EX: - v 1 2 3) (4,5,6)
        for p in self.points:
             buf.append('v %s %s %s'%( p[0], p[1], p[2]) ) #x y z components 
        
        buf.append('\n#these define the polygon geometry')
        for ply in self.polygons:
            plybuf = ''
            for f in ply:
                #plybuf = plybuf +(' %s'%str(int(f)+1) ) #add one because OBJ is NOT zero indexed
                plybuf = plybuf +(' %s'%str(int(f)) ) #add one because OBJ is NOT zero indexed

            if as_lines:#save as lines
                buf.append('l %s'%plybuf)
            else:#save as polygons 
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


