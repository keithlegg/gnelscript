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

        self.points    = []             #list of tuples of XYZ points per vertex         -  [(x,y,z), (x,y,z)]  
        self.polygons  = []             #list of tuples for 2 or more vertex connections -  [(1,2,5,8) , (1,2)] 

        #render properties
        self.linecolors = None #iterable of colors for lines 
        self.linecolor  = (0,240,00)
        self.vtxcolor   = (0,240,00)

    ############################################### 
    def get_mean_z(self, triangle):
        z1 = triangle[0][2]
        z2 = triangle[1][2]
        z3 = triangle[2][2]
        return (z1+z2+z3)/3

    ###############################################     
    def calc_bbox(self, object):
        """ UNFINISHED """
        maxx = 0
        maxy = 0
        maxz = 0
        
        for p in self.points:
            print(p)

    ###############################################  
    def get_face_verts(self, fid):
        """ lookup and return the constituent points of a polygon """
        tmp = []

        if fid<0 or fid > len(self.polygons)-1:
            print('#show_poly- bad face index : %s'%fid)
            return None

        for v_id in self.polygons[fid]:
            tmp.append(self.points[v_id-1])

        return tmp

    ###############################################  
    def get_face_normal(self, fid):
        """ UNFINISHED """
        tmp = self.get_face_verts(fid) 

        data = []
        data.append('\n############################')
        #data.append(str(tmp) )
        data.append('############################\n')

        #for d in data:
        #    print(d)  

    ###############################################  
    def get_face_edges(self, fid):
        """ UNFINISHED 
            return boundary edge segments as arrays of lines in 3D
            for extruding  
        """
        tmp = self.get_face_verts(fid) 

        data = []
        data.append('\n############################')
        #data.append(str(tmp) )
        data.append('############################\n')

        #for d in data:
        #    print(d) 

    ###############################################  
    def get_face_centroid(self, fid):
        """ UNFINISHED 
            return center of a polygon by ID 

            for extruding  
        """
        tmp = self.poly_centroid(self.points) 

        data = []
        data.append('\n############################')
        #data.append(str(tmp) )
        data.append('############################\n')

        #for d in data:
        #    print(d)  

    ###############################################  
    def show_poly(self, id):
        """ lookup and return the constituent points of a polygon """
        tmp = self.get_face_verts(id)   
        data = []
        data.append('\n############################')
        data.append(str(tmp) )
        data.append('############################\n')

        for d in data:
            print(d)

    ###############################################  
    def extrude_face(self, f_id):
        """ UNFINISHED """
        #edges = self.get_face_edges
        #for e in edge:
        #    build_poly(e)
        #ETC
        pass 

    ###############################################   
    def triangle_centroid(self, triangle):
        """ get 3D center of object (average XYZ point) """

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
        #return(0,0,0)

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
    def scribe(self, str):
        print(str)

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
    def radial_triangulate(self, as_new_obj=False, offset=None ):
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
            out_pts.append(fac_ctr)
            # add all the old points to our new dataset, plus our new center point   
            out_pts.extend(self.points)

            # iterate by two and connect to new radial center     
            for i in range(int(len(poly))):
                out_polys.append( (1, poly[i-1]+1, poly[i]+1 ) ) 
         
        if as_new_obj:
            self.points = out_pts
            self.polygons = out_polys
        else:    
            # STUPID BUG ALERT, you have to do the indecies first
            self._insert_poly_idxs(out_polys)
            self._insert_points(out_pts)


    ############################################### 
    def poly_loft(self, obj2, as_new_obj=False):
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

            #print("obj 1 ## ", poly )
            for pidx in poly:
                print('## ', self.points[pidx-1] )
                new_pts.append(self.points[pidx-1]) # first half of new poly quad

            print("obj 2 ## ", obj2.polygons[i] )
            for pidx in obj2.polygons[i]:
                print('## ', obj2.points[pidx-1] )
                new_pts.append(obj2.points[pidx-1]) # second half of new poly quad
                num = len(new_pts)
                new_plys.append( (num-3, num-2, num ) )                   # assemble the new polygon indices 
                #new_plys.append( (1,2,3,4) )                   # assemble the new polygon indices 

        
        print('#####################')
        print(new_plys)
        print(len(new_pts))

        self._insert_poly_idxs(obj2.polygons)
        self._insert_points(obj2.points)

        if True:#as_new_obj:
            self.points = new_pts
            self.polygons = new_plys
        else:    
            # STUPID BUG ALERT, you have to do the indecies first
            self._insert_poly_idxs(new_plys)
            self._insert_points(new_pts)

    ############################################### 
    def triangulate(self):
        """ 
            Only works on 3 or 4 sided polygons. 3 are passed unchanged, 4 are triangulated  
            turn a quad into two triangles 
            return a new 3d object (possibly with more new polygons, all triangles) 
        """
        out_polys = []

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
                self.radial_triangulate(offset=(0,0,0))

        # overwrite old data 
        self.polygons = out_polys



    ###############################################  
    
    def repair(self):
        """ walk internal data and fix any bad data found  (empty point tuples, etc) """

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

        buf.append('\n#these define the vertecies')
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

###############################################

class object3d(polygon_operator):

    def __init__(self):
        super().__init__()  

        self.face_normals    = []
        #self.point_normals  = []

        #self.point_colors   = []
        #self.line_colors    = []


        #self.geom_history = []
 
        self.rot       = [0,0,0]
        self.pos       = [0,0,0]
        self.scale     = [1,1,1]


    ############################################### 

    def flush(self):
        self.points       = [] # inherited
        self.polygons     = [] # inherited     
        self.face_normals = []  

        self.rot          = [0,0,0]
        self.pos          = [0,0,0]
        self.scale        = [1,1,1]

    ############################################### 

    def copy(self):
        new = type(self)()
        new.points   = self.points
        new.polygons = self.polygons  
        return new

    ############################################### 
    def _reindex_ply(self, f_idx, offset):
        """ take a tuple of face indexes and re-index to offset+value""" 

        out_face = [] 
        for i in f_idx:
           out_face.append(i+offset)
        return tuple(out_face)
        


    ############################################### 

    def append(self, otherobj):
        """ add another object to this one 
            points can be added, paying attention to the index 
            values. 

            polys have to be re-indexed to match the higher index values
        """
        n = len(self.points)

        for pt in otherobj.points:
            self.points.append(pt)

        for ply in otherobj.polygons:
            self.polygons.append(self._reindex_ply(ply,n))  

    ###############################################  
    def show(self):
        data = []
 
        tris = 0
        quads = 0
        other = 0        
        for p in self.polygons:
            if len(p)==3:
                tris+=1
            elif len(p)==4:
                quads+=1  
            else:
                other+=1

        data.append('\n############################')
        data.append('  position      : %s %s %s'%( self.pos[0], self.pos[1], self.pos[2]) )
        data.append('  rotation      : %s %s %s'%( self.rot[0], self.rot[1], self.rot[2]) )
        data.append('  scale         : %s %s %s'%( self.scale[0], self.scale[1], self.scale[2]) )
        data.append(' --------------------------- ' )        
        data.append('  num face normals   : %s' %(  len(self.face_normals) ) )#unimplemented
        data.append('  num verts          : %s' %(  len(self.points) ) )
        data.append('  num polygons       : %s' %(  len(self.polygons) ) )
        data.append(' --------------------------- ' )         
        data.append('  num triangles      : %s' %(   tris) )
        data.append('  num quads          : %s' %(   quads) )  
        data.append('  num other          : %s' %(   other) ) 
        data.append('############################\n')

        for d in data:
            print(d)


    def insert(self, obj):
        """ insert an objects geometry into this object 
        """

        if isinstance(obj, object3d):
            #print("########## DEBUG ", obj.points )
            #print("########## DEBUG ", obj.polygons )

            # STUPID BUG ALERT, you have to do the indecies first
            self._insert_poly_idxs(obj.polygons)
            self._insert_points(obj.points)


    ############################################### 
    def _insert_points(self, pt_array):
        """ append points to internal """

        #print("insert points count got called ", len(pt_array) )

        #for p in pt_array:
        self.points.extend(pt_array)

    ############################################### 
    def _insert_poly_idxs(self, idx_array):
        """ 
            if idx_array is a tuple, assume its a single polygon

            if idx_array is a list, assume it is a list of tuples
            representing multiple polygons
        """
        
        #print("insert polygon count got called ", idx_array, type(idx_array) )

        # add this to the indexes, so they get added "on top" of existing polygons 
        n_ply = len(self.polygons)

        #if any polygons are loaded, auto increment the face indecies  
        if n_ply>0:
            n = len(self.points) 
        else:
            n = 0

        # single polygon 
        if isinstance(idx_array, tuple):
            plytmp = []      
            for x in idx_array:
                plytmp.append(x+n) #add the poly index to current count            
            self.polygons.append( tuple(plytmp) )   

        # list of multiple polygons 
        if isinstance(idx_array, list):
            for tup in idx_array:
                
                plytmp = []      
                for x in tup:
                    plytmp.append(x+n) #add the poly index to current count            
                self.polygons.append( tuple(plytmp) ) 
    

    ############################################### 
    def convert_pts_vec3(self):
        """ return a list of pygfx vec3 objects for each vertex """

        vectrx = []
        #load each point into a pygfx.vec3 object 
        for pt in self.points:
            v = vec3().insert(pt)
            ##
            # tests of vec3 stuff 
            # v = v*1.1 
            #v = v.normal 

            ##
            vectrx.append(v.copy(vtype='tuple'))
        
        #tmp = object3d()
        #tmp.vectorlist_to_obj(vectrx)
        #tmp.save_obj('vtx_vectrz.obj')
        return vectrx

    ############################################### 
    def calc_face_normals(self, f_index=None):
        """ iterate each vertex by face assignment and convert 
            to vec3 
        """
        
        vectrx = []

        #secondary tweaks to the normal data 
        scale     = 1
        normalize = False  #make each face normal unit length 

        out_face_normals = [] 

        # iterate each face and convert eart vertex into a vec3 
        for f in self.polygons:
            
            f_nrml = None 

            # create a vec3 for each vertex (3 or 4 sided polys)
            v1=vec3();v2=vec3()
            v3=vec3();v4=vec3()

            # load each point into a pygfx.vec3 object 
            # 3 sided polys
            if len(f) == 3:
                v1.insert( self.points[f[0]-1] )
                v2.insert( self.points[f[1]-1] )
                v3.insert( self.points[f[2]-1] )
                 
                #calculate the face normal  
                a = v1 - v2;b = v1 - v3;
                if normalize:
                    f_nrml = a.cross(b).normal*scale
                else:    
                    f_nrml = a.cross(b)*scale 
                #get the position of the face 
                fpos = self.poly_centroid([v1,v2,v3])
                out_face_normals.append( (f_nrml,fpos) ) #store vec3, position

            # 4 sided polys     
            if len(f) == 4:
                v1.insert( self.points[f[0]-1] )
                v2.insert( self.points[f[1]-1] )
                v3.insert( self.points[f[2]-1] )                
                #v4.insert( self.points[f[3]-1] )

                #calculate the face normal  
                a = v1 - v2;b = v1 - v3;
                if normalize:
                    f_nrml = a.cross(b).normal*scale
                else:    
                    f_nrml = a.cross(b)*scale  
                #get the position of the face 
                fpos = self.poly_centroid([v1,v2,v3])

                out_face_normals.append( (f_nrml,fpos) ) #store vec3, position

            #store it in object for later use 
            self.face_normals.append(f_nrml) 

        return out_face_normals

    ############################################### 
    def test_vector_thingy(self, f_index=None):
        """ iterate each vertex by face assignment and convert 
            to vec3 
        """
        
        vectrx = []

        out_face_normals = [] 

        # iterate each face and convert eart vertex into a vec3 
        for f in self.polygons:
            
            # create a vec3 for each vertex (3 or 4 sided polys)
            v1=vec3();v2=vec3()
            v3=vec3();v4=vec3()

            # load each point into a pygfx.vec3 object 
            # 3 sided polys
            if len(f) == 3:
                v1.insert( self.points[f[0]-1] )
                v2.insert( self.points[f[1]-1] )
                v3.insert( self.points[f[2]-1] )

 
                a = v1 - v2;
                b = v1 - v3;
                out_face_normals.append(a.cross(b) )


                #print( self.mu.rtd( v1.angle_between(v2) ) )
                #out_face_normals.append( v3.vector_mean([v1,v2,v3] ).normal * 2 )
                #to get angle??
                #av = v3.average_angle( (1,0,0), [v1,v2,v3])
                #print( ' #average of multiple vecotrs  is ?? ', av )
                               

            # 4 sided polys     
            if len(f) == 4:
                v1.insert( self.points[f[0]-1] )
                v2.insert( self.points[f[1]-1] )
                v3.insert( self.points[f[2]-1] )                
                v4.insert( self.points[f[3]-1] )


                a = v1 - v2;
                b = v1 - v3;
                out_face_normals.append(a.cross(b) )

            ##
            #vectrx.append(v.copy(vtype='tuple'))
        
        #tmp = object3d()
        #tmp.vectorlist_to_obj(vectrx)
        #tmp.save_obj('vtx_vectrz.obj')
        #return vectrx
        return out_face_normals

    ###############################################  
    def prim_cone(self, axis='y', pos=(0,0,0), rot=(0,0,0), size=1):
        ## """ Not done yet - this makes a cone """

        spokes = 8
        #circpts = self.calc_circle(0, 0, size, spokes, False)   #2D data
        ## self.points = self.cvt_2d_to_3d(circpts)            #converted to 3D
        ## tmp = []
        ## for x in range(spokes):
        ##     tmp.append(x)
        ## tmp.append(0) #make periodic - insert start at end to close the cicle
        ## self.polygons = [tuple(tmp)]
        
        self.rotate_pts( rot )
        self.xform_pts( pos )


    ############################################### 
    def prim_sphere(self, pos=(0,0,0), rot=(0,0,0), radius=1):
        
        #UNFINISHED 

        #icosahedron  - from http://www.songho.ca/opengl/gl_sphere.html 

        #// constants
        PI = 3.1415926;
        H_ANGLE = PI/ 180*72;       # 72 degree = 360 / 5
        V_ANGLE = math.atan(1.0/2); # elevation = 26.565 degree

        i1 = 0
        i2 = 0                            # indices
        z  = 0
        xy = 0                            # coords
        hAngle1 = -PI / 2 - H_ANGLE / 2   # start from -126 deg at 1st row
        hAngle2 = -PI / 2                 # start from -90 deg at 2nd row

        # the first top vertex at (0, 0, r)
        self.points.append( (0,0,radius) )
        fid = len(self.points)-1

        faces = [] 

        # compute 10 vertices at 1st and 2nd rows
        for i in range(1,8):
            n = len(self.points) # add to this index each time

            z  = radius * math.sin(V_ANGLE)  # elevaton
            xy = radius * math.cos(V_ANGLE)  # length on XY plane

            vtmp1 = [];vtmp2 = []

            vtmp1.append( xy * math.cos(hAngle1)  )# x
            vtmp2.append( xy * math.cos(hAngle2)  )

            vtmp1.append( xy * math.sin(hAngle1)  )# y
            vtmp2.append( xy * math.sin(hAngle2)  )

            vtmp1.append(  z                      )# z
            vtmp2.append( -z                      )

            #// next horizontal angles
            hAngle1 += H_ANGLE
            hAngle2 += H_ANGLE
            
            self.points.append(tuple(vtmp1))
            self.points.append(tuple(vtmp2))


            if i>1:
                self.polygons.append( (n, n+1, n+2) )
            if i<7:
                self.polygons.append( (n+1, n+2, n+3) )
        
        # the last bottom vertex at (0, 0, -r)
        self.points.append( (0,0,-radius) )
               
        #end caps
        #self.points.append( (fid,fid+1,fid+2) )
        

        self.rotate_pts( rot )
        self.xform_pts( pos )

    ###############################################  
    def prim_locator(self, pos=(0,0,0), rot=(0,0,0), size=1):
   
        pts = [
               # x axis indicator  
               (pos[0],       pos[1],        pos[2]),
               (pos[0]+size,  pos[1],        pos[2]), 
               # y axis indicator
               (pos[0],       pos[1],        pos[2]),
               (pos[0],       pos[1]+size,   pos[2]), 
               # z axis indicator
               (pos[0],       pos[1],        pos[2]),
               (pos[0],       pos[1],        pos[2]+size)                
              ]

        plyidx = [(1,2),
                  (3,4),
                  (5,6)
                 ]
        
        self.linecolors = [
                 (255,0,0),         
                 (0,255,0), 
                 (0,0,255)                  
        ]

        self._insert_poly_idxs(plyidx) 
        self._insert_points(pts)  

        self.rotate_pts( rot )
        self.xform_pts( pos )

    ###############################################  
    def prim_locator_xyz(self, pos=(0,0,0), rot=(0,0,0), size=1):
       
        fs = size * .1  # font_size variable 
        fd = size * 1.1 # font distance (from origin)

        pts = [
               # x axis indicator  
               (pos[0]     ,     pos[1],        pos[2]),
               (pos[0]+size,     pos[1],        pos[2]), 

               # letter "X" on the X/Z plane  
               (pos[0]+fd+fs ,  pos[1],   pos[2]-fs  ), 
               (pos[0]+fd-fs ,  pos[1],   pos[2]+fs  ),
               
               (pos[0]+fd-fs ,  pos[1],   pos[2]-fs  ), 
               (pos[0]+fd+fs ,  pos[1],   pos[2]+fs  ),

               # y axis indicator
               (pos[0],       pos[1]     ,   pos[2]),
               (pos[0],       pos[1]+size,   pos[2]), 

               # letter "Y" on the Y/Z plane  
               (pos[0]-fs,    pos[1]+fd+(fs*3)  ,  pos[2]),
               (pos[0],       pos[1]+fd+(fs*2)  ,  pos[2]), 
               (pos[0]+fs,    pos[1]+fd+(fs*3)  ,  pos[2]),

               (pos[0],       pos[1]+fd         ,  pos[2]), 
               (pos[0],       pos[1]+fd+(fs*2)  ,  pos[2]),

               # z axis indicator
               (pos[0],       pos[1],        pos[2]),
               (pos[0],       pos[1],        pos[2]+size),  

               # letter "Z" on the Z/X plane  
               (pos[0]+(fs*1),  pos[1],   pos[2]+fd+(fs*1.5) ), 
               (pos[0]-(fs*1),  pos[1],   pos[2]+fd+(fs*1.5) ),
               (pos[0]+(fs*1),  pos[1],   pos[2]+fd        ), 
               (pos[0]-(fs*1),  pos[1],   pos[2]+fd        )
              ]
      
        
        plyidx = [ (1,2),
                   (3,4),

                   (5,6),

                   (7,8),
                   (9,10,11),                   
                   (12,13),
                   (14,15),
                   (16,17,18,19)                   
                 ]


        self._insert_poly_idxs(plyidx) 
        self._insert_points(pts)
        self.rotate_pts( rot )
        self.xform_pts( pos )

    ###############################################  
    def prim_line(self, axis='y', pos=(0,0,0), rot=(0,0,0), size=1):
        """ 3d lines, 2 point polygons """

        if axis=='x':
            pts =[ (-size,0,0), (size,0,0) ]
        if axis=='y':
            pts =[ (0,-size,0), (0,size,0) ]
        if axis=='z':
            pts =[ (0,0,-size) , (0,0,size) ]

        plyidx = (1,2)
       
        self._insert_poly_idxs(plyidx) 
        self._insert_points(pts)
        self.rotate_pts( rot )
        self.xform_pts( pos )

    ###############################################  
    def prim_triangle(self, axis='y', pos=(0,0,0), rot=(0,0,0), size=1):
        """ single polygon operations (that can be stacked together ?) """

        if axis=='x':
            pts =  [(0,-size,0), (0,0,size), (0,size,0) ]
        if axis=='y':
            pts =  [(0,0,-size), (size,0,0), (0,0,size) ]
        if axis=='z':
            pts =  [(-size,0,0), (0,size,0), (size,0,0) ]

        plyidx = (1,2,3)
       
        self._insert_poly_idxs(plyidx) 
        self._insert_points(pts)
        self.rotate_pts( rot )
        self.xform_pts( pos )

    ###############################################  
    def prim_quad(self, axis='y', pos=(0,0,0), rot=(0,0,0), size=1):
        """ single polygon operations (that can be stacked together ?) """
        
        if axis == 'x':
            pts = [(0,-size,-size), (0,-size,size), (0,size,size), (0,size,-size) ] #X AXIS

        if axis == 'y':
            pts = [(-size,0,-size), (-size,0,size), (size,0,size), (size,0,-size) ] #Y AXIS
            
        if axis == 'z':
            pts = [(-size,-size,0), (-size,size,0), (size,size,0), (size,-size,0) ] #Z AXIS

        plyidx    = (1,2,3,4)
       
        self._insert_poly_idxs(plyidx) 
        self._insert_points(pts)
        self.rotate_pts( rot )
        self.xform_pts( pos )

    ###############################################  
    def prim_circle(self, axis='z', pos=(0,0,0), rot=(0,0,0), size=1, spokes = 5):
        """ UNFINSIHED single polygon operations  """    

        # ARGS: calc_circle_2d( x_orig, y_orig, dia, periodic, spokes=23):
        #pts2d = self.calc_circle_2d(1, 0, size, False, spokes)   #2D data
        #for pt in self.cvt_2d_to_3d(pts2d):  #converted to 3D
        #    pts.append(pt)

        pts    = []
        plyidx = []
    
        # calc_circle( origin=(0,0,0), dia=1, axis='z', periodic=True, spokes=23):
        pts = self.calc_circle( (0,0,0), 1, axis, True, spokes )
        
        #we add one because calc_circle_2d returns zero indexed data but OBJ is NOT zero indexed        
        for x in range(spokes):
            plyidx.append(x+1) 

        #print(len(pts))
        #print(plyidx)

        self._insert_poly_idxs(tuple(plyidx))
        self._insert_points(pts)

    ############################################### 
    def prim_cube(self, linecolor=None, pos=(0,0,0), rot=(0,0,0), size=1):
        """ single polygon operations (that can be stacked togteher ?) """
        pts = [];plybfr = []
                
        # by adding position argument we can create it in different places
        pts.append( (-size+pos[0], -size+pos[1], size+pos[2])  ) #vertex 0
        pts.append( (-size+pos[0],  size+pos[1], size+pos[2])  ) #vertex 1
        pts.append( ( size+pos[0],  size+pos[1], size+pos[2])  ) #vertex 2  
        pts.append( ( size+pos[0], -size+pos[1], size+pos[2])  ) #vertex 3
        # notice these next 4 are the same coordinates with a negative Z instead!
        pts.append( (-size+pos[0], -size+pos[1], -size+pos[2])  ) #vertex 4
        pts.append( (-size+pos[0],  size+pos[1], -size+pos[2])  ) #vertex 5
        pts.append( ( size+pos[0],  size+pos[1], -size+pos[2])  ) #vertex 6  
        pts.append( ( size+pos[0], -size+pos[1], -size+pos[2])  ) #vertex 7

        # plot the connections between the points that will form polygons
        plybfr.append( (1,2,3,4) ) #polygon 0  #front
        plybfr.append( (5,8,7,6) ) #polygon 1  #top
        plybfr.append( (1,5,6,2) ) #polygon 2  #back
        plybfr.append( (2,6,7,3) ) #polygon 3  #right
        plybfr.append( (3,7,8,4) ) #polygon 4  #bottom
        plybfr.append( (5,1,4,8) ) #polygon 4  #bottom

        if linecolor:
            if not self.linecolors:
                self.linecolors = []
            for i in range(18):
                self.linecolors.append(linecolor)        

        self._insert_poly_idxs(plybfr) 
        self._insert_points(pts)
        self.rotate_pts( rot )
        self.xform_pts( pos )

    ############################################### 
    def one_vec_to_obj(self, r3, pos=None):
        """ single vector into a renderable 3D line """
        
        if pos:
            pts = [
                   (pos[0]       , pos[1]      , pos[2]       ),
                   (pos[0]+r3[0] , pos[1]+r3[1], pos[2]+r3[2] ),                   
                  ]

        if not pos:    
            pts = [
                   (0    , 0    , 0    ),
                   (r3[0], r3[1], r3[2]), 
                  ]

        n = len(self.points) # add this number to the indexes in case of existing geom 

        plyidx = [(n+1,n+2)]

        #append points to internal 
        for p in pts:
            self.points.append(p)
        for vec in plyidx:    
            self.polygons.append( vec )  

    ############################################### 
    def two_vecs_to_obj(self, r3_1, r3_2):
        """ a vector between two other vectors 
            probably not useful, but interesting 
        """
        
        pts = [
               (r3_1[0], r3_1[1], r3_1[2]),
               (r3_2[0], r3_2[1], r3_2[2]), 
              ]

        n = len(self.points) # add this number to the indexes in case of existing geom 

        plyidx = [(n+1,n+2)]

        #append points to internal 
        for p in pts:
            self.points.append(p)
        for vec in plyidx:    
            self.polygons.append( vec )  
 
    ############################################### 
    def vectorlist_to_obj(self, vecs):
        """ take a list of vectors and convert it to renderable geometry 
        
            vecs:
            - can be a list of single values (vec3)
              or a list 2 two values         (vec,position)


        """
        for v in vecs:
            if len(v) == 1:
                self.one_vec_to_obj(v) 
            if len(v) == 2:
                self.one_vec_to_obj(v[0], v[1])                 
                
        # experiment to make a line bewteen each 2 points 
        #for i,v in enumerate(vecs):
        #    if i>0:    
        #        self.two_vecs_to_obj(vecs[i-1], v) 

###############################################

class point_operator_2d(object):
    def __init__(self):
        self.mu   = mu()

    """
    def rotate_points_shifted(self, points, oldpivot, newpivot, angle, doOffset=False, doRound=False):
        #old 2d rotate function  
        rotated_fids =  self.batch_rotate_pts( points, oldpivot, angle , doRound)

        deltax = newpivot[0]-oldpivot[0]
        deltay = newpivot[1]-oldpivot[1]

        newfids = []
        for pt in rotated_fids:
            if not doOffset:
                newfids.append( (pt[0]+deltax, pt[1]+deltay) )
            if doOffset:
                newfids.append( (pt[0]+(deltax+doOffset[0]), pt[1]+(deltay+doOffset[1])  ) )

        return newfids 
    """

    def calc_line(self, x1, y1, x2, y2):        
        """ bresenham's algorithm 
            from http://www.roguebasin.com/index.php?title=Bresenham%27s_Line_Algorithm
        """
        x1 = int(x1);y1 = int(y1)
        x2 = int(x2);y2 = int(y2)

        points = []
        issteep = abs(y2-y1) > abs(x2-x1)
        if issteep:
            x1, y1 = y1, x1
            x2, y2 = y2, x2
        rev = False
        if x1 > x2:
            x1, x2 = x2, x1
            y1, y2 = y2, y1
            rev = True
        deltax = x2 - x1
        deltay = abs(y2-y1)
        error = int(deltax / 2)
        y = y1
        ystep = None
        if y1 < y2:
            ystep = 1
        else:
            ystep = -1
        for x in range(x1, x2 + 1):
            if issteep:
                points.append((y, x))
            else:
                points.append((x, y))
            error -= deltay
            if error < 0:
                y += ystep
                error += deltax
        # Reverse the list if the coordinates were reversed
        if rev:
            points.reverse()
        return points

    def aggregate_extents(self, points, offset=False):
        """ needs at least three points to work 
            take any number of points - (3 or more tuples of (x,y) ),  and find the min/max of all of them as a whole
        """
       
        minx = points[0][0];maxx = 0
        miny = points[0][1];maxy = 0
        if not offset:
            for p in points:
                if p[0]<=minx:
                    minx = p[0]
                if p[0]>=maxx:
                    maxx = p[0]
                if p[1]<=miny:
                    miny = p[1]
                if p[1]>=maxy:
                    maxy = p[1]
        if offset:
            for p in points:
                if p[0]<=minx:
                    minx = p[0]-offset
                if p[0]>=maxx:
                    maxx = p[0]+offset
                if p[1]<=miny:
                    miny = p[1]-offset
                if p[1]>=maxy:
                    maxy = p[1]+offset

        return [minx, miny, maxx, maxy ]
        
    def rotate_point_2d(self, point, pivot, angle, doRound=False):
        """ helper to rotate a single point in 2 dimensions """
        
        x1 = point[0] - pivot[0]  #x
        y1 = point[1] - pivot[1]  #y  

        dtr = self.dtr

        a = x1 * math.cos(dtr(angle)) - y1 * math.sin(dtr(angle))
        b = x1 * math.sin(dtr(angle)) + y1 * math.cos(dtr(angle))
                
        if not doRound:
            return (a + pivot[0], b + pivot[1]);
        if doRound:
            return ( round(a + pivot[0]), round(b + pivot[1]) );

    def batch_transform_pts_2d(self, pts, new_coord, doround=False):
        """ UNTESTED ? - takes a list of tuples and returns a shifted list of tuples , based on a tuple of (x,y) shift"""

        out = []
        for pt in pts:
            if not doround:
                out.append( (pt[0]+ new_coord[0], pt[1]+ new_coord[1] ) ) 
            if doround:
                out.append( ( int(pt[0]+ new_coord[0]), int(pt[1]+ new_coord[1]) ) )         
        return out

    def batch_rotate_pts_2d(self, pts, pivot, angle, doround=False):
        out = []
        for pt in pts:
            out.append( self.rotate_point_2d(pt, pivot, angle, doround) )
        return out

    def calc_circle_2d(self, x_orig, y_orig, dia, periodic=True, spokes=23):
        """ spokes = num spokes """

        plot_x = 0;plot_y = 0;
        out = []
        
        degs = 360/spokes
        dit = 0
        
        dtr = self.mu.dtr

        for i in range(spokes):
            plot_x = x_orig + (math.sin(dtr(dit))*dia) 
            plot_y = y_orig + (math.cos(dtr(dit))*dia) 
            out.append( (plot_x, plot_y))
             
            dit+=degs
        
        if periodic:
             out.append( out[0] )

        return out


    def pt_offset_to_line( self, vector, pt): # x3,y3 is the point
        """
           given a vector and a point, return a vector between them 
           (you can get the length of that to calculate distance)
        """
        x1 = vector[0][0] ; y1 = vector[0][1]
        x2 = vector[1][0] ; y2 = vector[1][1]
        x3 = pt[0]        ; y3 = pt[1]

        px = x2-x1
        py = y2-y1
        
        notzero = px*px + py*py
        u =  ((x3 - x1) * px + (y3 - y1) * py) / float(notzero)

        if u > 1:
            u = 1
        elif u < 0:
            u = 0

        x = x1 + u * px
        y = y1 + u * py
        dx = x - x3
        dy = y - y3

        return [ pt , (pt[0]+dx, pt[1]+dy) ]
        
        ###################
        #if you want to get distance
        #return math.sqrt(dx*dx + dy*dy)

    def extract_pt_vector(self, vector, offset):
        """
          given a vector [(x,y), (x,y)] and an offset - return a point (x,y)

          same as vector class project pt along a vector - even beyond its extent
        """

        a=(vector[0][0], vector[0][1])
        b=(vector[1][0], vector[1][1])
        
        nX = b[0] - a[0];nY = b[1] - a[1]
        distX = pow( (a[0] - b[0] ) , 2.0 ) 
        distY = pow( (a[1] - b[1] ) , 2.0 ) 
        vecLength = math.sqrt(distX + distY )
        # normalized vector  
        calcX = nX / vecLength
        calcY = nY / vecLength
        # project point along vector with offset (can use negative too)
        ptX = b[0] + (calcX * offset)
        ptY = b[1] + (calcY * offset)
        return (ptX, ptY)

    def locate_pt_along(self, x1, y1, x2, y2, num):
        """ taken from my old 3d character rigging code (with the Z axis removed)
            this was used to place vertebra along a spine
        """

        pts_created = []
        fpos=(x1, y1)
        spos=(x2, y2)
         
        for n in range(num):
            npos = [2]
            npos[0]    =spos[0]+(((fpos[0]-spos[0])/(num+1))*(n+1))
            npos.append(spos[1]+(((fpos[1]-spos[1])/(num+1))*(n+1)) )
            pts_created.append( (npos[0], npos[1]) )
        return pts_created

    def calc_mid(self, pt1, pt2, num=1, doround=False):
        """ get the midpoint of a 2d vector """

        out = self.locate_pt_along( pt1[0], pt1[1], pt2[0], pt2[1], num )
        if doround:
            return ( int(out[0][0]), int(out[0][1]) )
        if not doround:
            return ( out[0][0], out[0][1] )
                
    def old_calc_theta_vert(self, bot_xy, top_xy):
        try:
            a = bot_xy[1]-top_xy[1]  
            o = bot_xy[0]-top_xy[0]
            #if a==0 or o ==0:
            #    r = self.rtd(math.atan(o)) 
            #else:
            r = self.rtd(math.atan(o/a)) 
            return r
        except:
            print('error calc angle - div by 0?') 
            return None 

    def calc_theta_vert(self, start_xy, end_xy, corner, mode='relative'):
        """
           convert a 2 point line segment into a meaningful rotation value in degrees
           outputs a 0-360(?) degree rotation with 0/360 facing left from center , 180 right from center
           converts the line into a right triangle, calculates theta,
               absolute mode - projects the angle into a quadrant of a 0-360 degree circle
               relative mode - number of degrees of current rotation to get back to a normal page orientation

            instead of attempting more complex math, the rotation is filtered though a series of steps 
            that use the corner fiducial as a variable to determine which end is proper +Y axis.
            There are 4 cases, positive and negative (i.e.  rotated to the left , or the right for each of the 4 orintations) 
            the four states represent four page orientations (i.e.  up, down left,right)   
        """


        #get corner to build a right triangle
        a_x = end_xy[0]-start_xy[0]  
        a_y = end_xy[1]-start_xy[1]

        r = 0
        #relative offset (depending on order of start-end)
        if a_x != 0 and a_y !=0:
            r = ( self.rtd(math.atan(a_x/a_y))  )

        #make it positive
        if r <0:
            r = -r

        ################################
        ### if you want absolute rotation - "projected" into a circle to get 0-360
        #this is with 0/360 pointing to the left of center
        #DEBUG - ABSOLUTE MODE IS UNFINISHED/UNTESTED 
        if mode=='absolute':
            #this works with positive rotation ??debug            
            if end_xy[1]<start_xy[1]:
                if end_xy[0]<start_xy[0]:
                    r+=270
                if end_xy[0]>start_xy[0]:
                    r+=180       
            if end_xy[1]>start_xy[1]:
                  if end_xy[0]>start_xy[0]:
                    r+=90  

        ################################
        ###  if you want relative offset 
        if mode=='relative':        
            isnegative = False
            #first we need to determine if offset is positive or negative 
            if corner=='bl':
                if end_xy[0]<start_xy[0]:
                    isnegative = True
            if corner=='tr':
                if end_xy[0]>start_xy[0]:
                    isnegative = True            
            if corner=='br' :
                if end_xy[1]>start_xy[1]:
                    isnegative = True
            if corner== 'tl':
                if end_xy[1]<start_xy[1]:
                    isnegative = True

            #deal with positive rotation 
            if not isnegative:
                if end_xy[1]<start_xy[1]:
                    if end_xy[0]<start_xy[0]:
                        r=180-r
                    if end_xy[0]>start_xy[0]:
                        r=r+180       
                if end_xy[1]>start_xy[1]:
                      if end_xy[0]<start_xy[0]:
                        r=r  
                      if end_xy[0]>start_xy[0]:
                        r=360-r  

            #deal with negative rotation
            if isnegative:
                if end_xy[1]<start_xy[1]:
                    if end_xy[0]<start_xy[0]:
                        r=180-r
                    if end_xy[0]>start_xy[0]:
                        r=r+180       
                if end_xy[1]>start_xy[1]:
                    if end_xy[0]<start_xy[0]:
                        r=r   
                    if end_xy[0]>start_xy[0]:
                        r=360-r  

        #####
        if r == 0 or r == -0:
            """
              at this point we are done - unless the rotation came back as 0 or -0.
              why does it come back zero you may ask? because it uses a right triangle to calculate the theta - if the page is 
              too perfect it gets a line instead of a triangle and we cant form an anlge from it. In that case we fall back on 
              the self.corner to calculate the rotation.
            """
            #if it is zero - we can safely assume page is axis aligned - therefore we simply look at self.corner and Bob's yer uncle.
            if corner=='tl':
                r  = -90
            if corner=='tr':
                r  = 0                
            if corner=='br':
                r  = 90                
            if corner=='bl':       
                r = 180  
                                                           
        return r

