
import os
import math 

from gnelscript import NUMPY_IS_LOADED, GEOJSON_IS_LOADED, SHAPELY_IS_LOADED, \
                       NETWORKX_IS_LOADED, TRIMESH_IS_LOADED


from gnelscript.pygfx.math_ops import math_util as mu
from gnelscript.pygfx.math_ops import matrix33, matrix44, vec2, vec3  

##-----------------------------##

if GEOJSON_IS_LOADED:
    from geojson import dump  
    from geojson import Point as gjpt
    from geojson import Polygon as gjply
    from geojson import Feature as gjftr
    from geojson import LineString as gjln
    from geojson import FeatureCollection as gjfc

if SHAPELY_IS_LOADED:
    from shapely import buffer, BufferCapStyle, BufferJoinStyle
    from shapely import Point as shp_pt
    from shapely import Polygon as shp_ply
    from shapely import LineString as shp_ln
    #from shapely import Feature as shp_ftr
    #from shapely import FeatureCollection as shp_fc

if NETWORKX_IS_LOADED:
    import networkx as nx

if TRIMESH_IS_LOADED:
    import trimesh

if NUMPY_IS_LOADED:
    import numpy as np  


##-----------------------------##
##-----------------------------## 


"""

    defined objects:

    BBOX - 2D left, top, right, bottom on Z axis 
    
    BBOX3 - 3D left, top , right, bottom, fromt, back on Z axis 

    GEOM  -  [[[face ids]], [vertices] ]
        simplest way to define a 3d model. It is two arrays, face ids and vertices
        can be used to write translators to other formats  
        geom type is a 3d model basically obj format in memory 

        BEWARE THE ZERO INDEX HEADACHE 
        OBJS ARE NOT SERO INDEXED 
        EVERTYTHING ELSE SHOULD BE 

        # HOW TO ITERATE GEOM OBJECTS (POLY IDS)  
        for ply in o[0]:
            for i in range(len(ply)):
                print(ply[i])

        # HOW TO ITERATE GEOM OBJECTS WITH POLYIDS WITH ZERO INDEX   
        for ply in o[0]:
            for i in range(len(ply)):
                print(o[1][ply[i]-1])


    POINT GROUPS -  [ [ID, (X,Y,Z)], ... ] 
        
        POINT GROUPS SHOULD BE ZERO INDEXED!
        I think the choice was made NOT to zero index because OBJ is not zero indexec 
        give this a big think 

        add check in all functions to ensure this - fix all fucntions!
        a point group is another type of container (DEBUG partially implemented)
        it is only points with a an ID for each point
        it is a way to work with partial objects and not loose the ID of each point 

 
    TODO:
        upgrade to these to full objects instead of lists 
        make operators work universally on GEOM, OBJECT, PTGRPS, ETC 
        SOLVE THE STUPID ZERO INDEX ISSUE 

        class ptgrp(object):
            def __init__(self):
                ## store in zero index add function to get NON ZERO INDEXED
                self.pts = [ ] #[idx, (pt), ..] 

            def non_zero(self)
                pass 
            def centroid(self)
                pass 
            def bbox(self)
                pass  

 

"""






def printgeom(geom):
    for i,fid in enumerate(geom[0]):
        print( " f_id %s - %s"%(i,fid))
    for pt in geom[1]:
        print( " pt %s"%str(pt) )



##-------------------------------------------##
##-------------------------------------------## 
## PLAYING WITH TRIMESH - LOOKS VERY INTERESTING 

#print(dir(trimesh))

# trimesh.util.append_faces(vertices_seq, faces_seq)

# trimesh.load_path(file_obj, file_type=None, **kwargs)

#edges_unique


class tm_pop3(object):
    def __init__(self):
        self.mesh = None

    def load(self, infile):
        print('#tm_pop3 loading %s'%infile)
        self.mesh = trimesh.load_mesh(infile) 

    def move(self):
        #trimesh.transform_points(points, matrix, translate=True)
        print()

    def convex_hull(self, pts, outfile):
        pc = trimesh.PointCloud(vertices=pts)
        cvh = pc.convex_hull
        cvh.export(outfile)





##-------------------------------------------##
##-------------------------------------------##  

class pop3d(object):
    """ point_operator_3d() 
        what became of the original point generator 

        deal with raw points, not geom, ptgrps, facgrps. JUST points 
    """

    def __init__(self):
        self.mu     = mu()
        self.m33    = matrix33()      
        self.m44    = matrix44()  # used for point rotations, and more? 
        self.vec2   = vec2()     
        self.vec3   = vec3()      


    ##-------------------------------------------##
    #def mirror(self, origin, axis):
        """ mirror a set of points """

    ##-------------------------------------------## 
    #def grow_selection(self, ptgrp, facgrp, u_num, v_num):


    ##-------------------------------------------## 
    #def merge_all_pts(self):

    #def merge_pts_line(self):
    #attempt to form a series of pts into a single line 
    
    #def merge_segs_line(self):
    #attempt to form a series of line segs into a single line 

    #def weld_edges(self, obj):

    ##-------------------------------------------##  
    #def scan_shells(self, obj):
    #    """ look for all the chunks of geometry that are not connected """

    ##-------------------------------------------##   
    #def poly_seperate(self, obj):
    #    """ check for geometry that is not connected, if any found, break it off """

    ##-------------------------------------------##   
    #def get_edge_centroid(self, f_id , e_id):
    #    pass


    ##-------------------------------------------##
    def copy_rotate(self, points, pos=None, num=4, axis='y'):
        """DEBUG UNTESTED 
           for tesselating 
           take a group of points - duplicate them, apply a rotation and repeat 
       
        """

        out = [] 

        tmp =[] 
        
        num=num+1 
        
        angle = 360/num 

        for i in range(num):
            tmp = points 
            
            m33= matrix33()

            if axis=='x':
                m33.from_euler([angle*i,0       ,0])
            if axis=='y':
                m33.from_euler([0      ,angle*i ,0])
            if axis=='z':
                m33.from_euler([0      ,0       ,angle*i])

            tmp = self.apply_matrix_pts(pts=points, m33=m33)
            
            #optional move  
            if pos is not None:
                #tmp = self.trs_points(out, translate=pos) 

                tmp2 = []
                for pt in tmp: 
                    x = pt[0] + pos[0]
                    y = pt[1] + pos[1]
                    z = pt[2] + pos[2]
                    tmp2.append( (x,y,z) )
                out.extend(tmp2)

            

        return out 


    ##-------------------------------------------## 
    def pt_in(self, pt, pts):
        for sample in pts:
            if str(pt)==str(sample):
                return True 
        return False            


    #def pt_is_close(self, pt, pts):
    #   check to see if points are near each other 


    ##-------------------------------------------## 
    def get_pt_extreme(self, pts, mode='idx'):
        """ 
            get the farthest extents of a point group 
            (the most +X,-X,+Y, etc )

            return the indeces to the points relatove to the array
            
            DEBUG - not sure how to handle planar objects where an axis gets the same for both 

        """

        minx=0;miny=0;minz=0 
        maxx=0;maxy=0;maxz=0 
    
        minx_id=0;maxx_id=0
        miny_id=0;maxy_id=0 
        minz_id=0;maxz_id=0 

        #print(" gr_polys buffer has %s polys in it "%len(self.gr_polys) )
        for ii,pt in enumerate(pts):
            if ii == 0:
                minx=pt[0]
                maxx=pt[0]
                miny=pt[1]
                maxy=pt[1]
                minz=pt[2]
                maxz=pt[2]

            if pt[0]<minx:
                minx=pt[0]
                minx_id=ii    
            if pt[0]>maxx:
                maxx=pt[0]  
                maxx_id=ii                  
            if pt[1]<miny:
                miny=pt[1]  
                miny_id=ii                  
            if pt[1]>maxy:
                maxy=pt[1]
                maxy_id=ii                   
            if pt[2]<minz:
                minz=pt[2]  
                minz_id=ii
            if pt[2]>maxz:
                maxz=pt[2] 
                maxz_id=ii

        if mode=='idx':
            return [minx_id, maxx_id, miny_id, maxy_id, minz_id, maxz_id] 
        
        if mode=='scalar':
            return [minx   , maxx   , miny   , maxy   , minz   , maxz ]

        if mode=='points':
            return [pts[minx_id], pts[maxx_id], pts[miny_id], pts[maxy_id], pts[minz_id], pts[maxz_id] ]

    ##-------------------------------------------## 
    def _csp_str(self, pt):
        """ clean single point - tuple of 2 or 3 floats
            destructive command - cleans floats for export but will loose precision 
            only run when object is exported 

             clean single point to 8 places of precision  
               -remove scientific notation (exponents)

        """
        if len(pt)==2:
            return (f'{pt[0]:.8f}',f'{pt[1]:.8f}')
        else:    
            return (f'{pt[0]:.8f}',f'{pt[1]:.8f}',f'{pt[2]:.8f}')
    
    ##-------------------------------------------## 
    def cf(self, pts):
        """ clean float data 
            destructive command - cleans floats for export but will loose precision 
            only run when object is exported 
               -remove scientific notation (exponents)
        """
        outpts = []
        
        for pt in pts:
            outpts.append( f'{pt:.8f}' )
        return outputs

    ##-------------------------------------------## 
    def apply_matrix_pts_round(self, pl, pts, m33=None, m44=None):
        """ UNTESTED 
            DO ROUNDING AND batch mutliply points by a matrix 
        """
        tmp_buffer = [] 

        # apply the transform here
        for pt in pts:  
            if m33 is not None:
                tmp = m44*pt                
                rpt = (round(tmp[0],pl),
                       round(tmp[1],pl),
                       round(tmp[2],pl)
                       )  
                tmp_buffer.append( rpt )

            if m44 is not None:
                tmp = m44*pt
                rpt = (round(tmp[0],pl),
                       round(tmp[1],pl),
                       round(tmp[2],pl)
                       )  
                tmp_buffer.append( rpt )
        return tmp_buffer

    ##-------------------------------------------## 
    def apply_matrix_pts(self, pts, m33=None, m44=None):
        """ 
            batch mutliply points by a matrix 
        """
  
        tmp_buffer = [] 

        # apply the transform here
        for pt in pts: 
            if m33 is not None:
                tmp_buffer.append( m33 * pt )
            if m44 is not None:
                tmp_buffer.append( m44 * pt )

        return tmp_buffer

    ##-------------------------------------------## 
    def tuple_pop(self, listups, tup):
        """ take a list of tuples, remove one by filtering it out and return the rest back 

            tups = [(1,1), (2,2), (2,3), (42,23)]
            ids = pop3.tuple_pop(tups, (1,1) )
            print(ids)

            output: [(2, 2), (2, 3), (42, 23)]

        """

        out = []
        for t in listups:
            if(t != tup):
                out.append(t)
        return out

    ##-------------------------------------------## 
    def test_data_grid(self, width, height, divs):
        """ 
            usage
                pop3 = pop3d()
                ids = pop3.test_data_grid( 5,5,1)
                print(ids)

                output= [ [0, 1, 2, 3, 4], 
                          [0, 1, 2, 3, 4], 
                          [0, 1, 2, 3, 4], 
                          [0, 1, 2, 3, 4], 
                          [0, 1, 2, 3, 4] ]

        """

        out = []

        va = []

        for u in range(0, width, divs):
            va = []
            for v in range(0, height, divs):
                va.append(v)
                 
            out.append(va)     

        return out

    ##-------------------------------------------##         
    def indexer(self, ids=None, span=None, unique=True, nth=None):
        """ generates a fancy list of positive ints (indeces)

            the idea is to additively build up a list of ids 
            you can start with an optional list of ids,
            then add in a range.  
 
            ids    - list of ids (optional data to start with)

            span   - batch add numbers in a [start, end] [start,end] 
                     choose to count by Nths, single or list of them 
                     negative Nths remove 

            unique is True by default - it guarantees each id is unique
                    - if off, an index will be repeated 

            nth - removes every N index
                  skips over N indices while iterating. Outputs two lists ,
                  the "goods" and the "bads" in order of  [ goods, rejects ] 

            ------------------------------------------

            usage:

            idz = range(1,20)
            ids = pop3.indexer( ids=idz, span=None, unique=True, nth=None)
            ouputs: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]


            ids = [41,52,75]
            span = [2,15]
            ids = pop3.indexer( ids=ids, span=span, unique=True, nth=None)
            ouputs: [41, 52, 75, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]


            ids = [41,52,75]
            span = [2,15]
            ids = pop3.indexer( ids=ids, span=span, unique=True, nth=5)
            outputs: [[41, 4, 9, 14], [52, 75, 2, 3, 5, 6, 7, 8, 10, 11, 12, 13, 15]]


        """

        pids = []
 

        # insert list of individual ids first 
        if ids:
            if unique:
                for i in ids:
                    if type(i) is int:    
                        if i not in pids:
                            pids.append(i)
            else: 
                for i in ids:
                    if type(i) is int:                 
                        pids.extend(ids)            
        #### 
        # then do the span
        if span and len(span)>1:
            tids = []
            # insert span IDs (start-end range) 
            for i in range(span[0], span[1]+1):
                if unique:
                    if i not in pids:
                        tids.append(i)
                else: 
                    tids.append(i)
            pids.extend(tids)

        ####        
        # iterate that on Nth 
        #not just for span - operate on IDS too ??
        if nth:
            if nth<len(pids):
                new_pids = []
                rejects = []
                if type(nth)==int:
                    for i,pd in enumerate(pids):
                        if i%nth==0:
                            new_pids.append(pd)
                        else:
                            rejects.append(pd)

                pids = [new_pids,rejects] 
        

        return pids 

    ##-------------------------------------------##  
    def chunker(self, pts, mod):
        """
            iterate a "mass" of points by N and separate into N sized groups 
            omits leftovers that dont fit 

            usage:
                pop3 = pop3d()
                pts = [1,2,3,4,5,6,7,8,9,0]
                out = pop3.chunker(pts, 7) 
                output: [[2, 3, 4, 5, 6, 7, 8]]

        """

        outarrays = []
        newgrp = []
        for i,p in enumerate(pts):
           if i%mod==0:
               newgrp.append(p)
               if len(newgrp)>1:
                   outarrays.append(newgrp)
               newgrp=[] 
           else:
               newgrp.append(p)
        
        return outarrays

    ##-------------------------------------------##
    def print_gridinfo(self, grid_array):
        """ 

            usage:
                pop3 = pop3d()
                ids = pop3.test_data_grid( 5,5,1)
                pop3.print_gridinfo(ids)

        """
        print("number of columns " , len(grid_array) )
        for column in grid_array:
            print("number of rows " , len(column) )
    ##-------------------------------------------## 
    def print_grid(self, grid_array):
        """
            usage:        
                pop3 = pop3d()
                ids = pop3.test_data_grid( 5,5,1)
                pop3.print_grid(ids)   

        """
        for column in grid_array:
            #print("number of rows " , len(column) )
            print( column)   
    ##-------------------------------------------##        
    def get_grid_column(self, grid_array , colidx):
        """
            usage:        
                pop3 = pop3d()
                ids = pop3.test_data_grid( 5,5,1)
                out = pop3.get_grid_column(ids , 2)
                print(out)

        """
    
        out = []
        for u,column in enumerate(grid_array):
            for v,row in enumerate(column):
                if v == colidx:
                    out.append(row)  

        return out 

    ##-------------------------------------------##
    def get_grid_row(self, grid_array , rowidx):
        """
            usage:
                pop3 = pop3d()
                ids = pop3.test_data_grid( 5,5,1)
                out = pop3.get_grid_row(ids , 2)
                print(out)       
        """
        #print("number of columns " , len(grid_array) )
        
        out = []
        for u,column in enumerate(grid_array):
            for v,row in enumerate(column):
                #print("number of rows " , len(column) )
                if u == rowidx:
                    out.append(row)  

        return out 

    ##-------------------------------------------##
    def cubic_bezier2d(self, draw_steps, start, end, ctrl1, ctrl2):
        """  2D spline 
            pretty much what it says it does 
        """
        out = []
            
        for i in range( draw_steps):
            t = i / draw_steps
            tt = t * t
            ttt = tt * t
            u = 1 - t
            uu = u * u
            uuu = uu * u

            x = uuu * start[0];
            x += 3 * uu * t * ctrl1[0]
            x += 3 * u * tt * ctrl2[0]
            x += ttt * end[0]

            y = uuu * start[1]
            y += 3 * uu * t * ctrl1[1]
            y += 3 * u * tt * ctrl2[1]
            y += ttt * end[1]


            out.append( [x, y] )

        return out
    
    ##-------------------------------------------##
    def cubic_bezier(self, draw_steps, start, ctrl1, ctrl2, end):
        """  3D spline 

            usage:
                obj = object3d()
                num = 23
                start = (1 ,  0, 0)
                ctrl1 = (.5,  0, 0)
                ctrl2 = ( 0, .5, 0)
                end   = (0 ,  1, 0)
                curve = obj.cubic_bezier(num, start, ctrl1, ctrl2, end)
                obj.lathe(curve, num)
                obj.save('rotateowl.obj')  

        """
        
        out = []
       
        
        for i in range( draw_steps):
            t = i / draw_steps
            tt = t * t
            ttt = tt * t
            u = 1 - t
            uu = u * u
            uuu = uu * u

            x = uuu * start[0];
            x += 3 * uu * t * ctrl1[0]
            x += 3 * u * tt * ctrl2[0]
            x += ttt * end[0]

            y = uuu * start[1]
            y += 3 * uu * t * ctrl1[1]
            y += 3 * u * tt * ctrl2[1]
            y += ttt * end[1]

            z = uuu * start[2]
            z += 3 * uu * t * ctrl1[2]
            z += 3 * u * tt * ctrl2[2]
            z += ttt * end[2]

            out.append( [x, y, z] )

        return out

    ##-------------------------------------------##
    def draw_splines(self, num, curves, drawctrls=False, drawhulls=False, singlepoly=False):
        """ render a spline made of multiple cubic bezier curves
            
            ARGS:
                curves - array of [start, ctrl1 , ctrl2, end ]

            usage:

                obj = object3d()

                start = (1 ,  0, 0)
                ctrl1 = (.5,  0, 0)
                ctrl2 = ( 0, .5, 0)
                end   = (0 ,  1, 0)
                kurve = [start, ctrl1, ctrl2, end]

                start  = (0   ,  1   , 0)
                ctrl1  = (1.5 ,  1.0   , 0)
                ctrl2  = (1   ,  1.5 , 0)
                end    = (2   ,  2   , 0)
                kurve2 =[start, ctrl1, ctrl2, end]
                
                curves = [kurve, kurve2] 
                
                obj.draw_splines( 10, curves, drawctrls=True, drawhulls=True)

                obj.save('splines.obj')                

        """
        size = .1

        for c in curves:    
            if drawctrls:
                self.prim_locator_color(pos=c[0] , rot=(0,0,0), size=size) # start 
                self.prim_locator_color(pos=c[1] , rot=(0,0,0), size=size) # ctrl1
                self.prim_locator_color(pos=c[2] , rot=(0,0,0), size=size) # ctrl2
                self.prim_locator_color(pos=c[3] , rot=(0,0,0), size=size) # end
            if drawhulls:
                self.linegeom_fr_points( [c[0], c[1], c[2], c[3]], color=(0,255,0) ) 
            curvepts = self.cubic_bezier(num, c[0], c[1], c[2], c[3])

            if singlepoly:
                self.linegeom_fr_points( curvepts, singlepoly=True ) 
            else:
                self.linegeom_fr_points( curvepts ) 

 

    ##-------------------------------------------##            
    def bbox_buffer_2d(self, bbox, size):
        """ was called add_margin_bbox 

            BBOX is 2D [-x,-y, +x, +y ]

            return center (x,y) from two diagonal coordinates 
            assuming a bbox is [left, top, right, bottom ] it "grows" the size of it 

        """
        
        out = []
        out.append( bbox[0]-size  ) 
        out.append( bbox[1]-size  )
        out.append( bbox[2]+size  )
        out.append( bbox[3]+size  )
        return out

    ##-------------------------------------------##
    def extents_fr_bbox(self, bbox, offset=None, periodic=False):
        """ return 2D pts geom from a 2D bbox  
            
            args:
               bbox   - iterable of 4 numbers (PIL bbox [left, top, right, bottom]) 
               offset - (int) adds a margin to the size of the page edges 

        """
        
        out = []
        if not offset:
            out.append( (bbox[0], bbox[1]) ) #top left
            out.append( (bbox[2], bbox[1]) ) #top right
            out.append( (bbox[2], bbox[3]) ) #bottom right
            out.append( (bbox[0], bbox[3]) ) #bottom left
            if periodic:
                out.append( (bbox[0], bbox[1]) ) #top left                

        #increase the margins, (pixels are indexed with top left at 0,0)
        if offset:
            out.append( (bbox[0]-offset , bbox[1]-offset ) ) #top left
            out.append( (bbox[2]-offset , bbox[1]+offset ) ) #top right
            out.append( (bbox[2]+offset , bbox[3]+offset ) ) #bottom right
            out.append( (bbox[0]-offset , bbox[3]-offset ) ) #bottom left
            if periodic:
                out.append( (bbox[0]-offset , bbox[1]-offset ) ) #top left

        return out


    ##-------------------------------------------## 
    def calc_square_diag(self, tl, br, add_zaxis=False):
        """ creates 4 coordinates representing a 2D extent box from a diagonal coordinate pair 
            
            usage:

            bbox = self.calc_3d_bbox()
            pts = self.calc_square_diag((bbox[0],bbox[1]),(bbox[3],bbox[4]), add_zaxis=True)
         
        """

        out =[]  
        if add_zaxis:
            out.append( (tl[0], tl[1], 0)  ) #tl
            out.append( (br[0], tl[1], 0)  ) #tr
            out.append( (br[0], br[1], 0)  ) #br
            out.append( (tl[0], br[1], 0)  ) #bl
        else:
            out.append( (tl[0], tl[1])  ) #tl
            out.append( (br[0], tl[1])  ) #tr
            out.append( (br[0], br[1])  ) #br
            out.append( (tl[0], br[1])  ) #bl            
        return out

    ##-------------------------------------------##
    def calc_2d_square_pt(self, size, origin=None ):
        """ UNTESTED 
            DEBUG - this is unclamped  
            calc the XY coodrinates for a square  
            (top left, top right , bottom right, bottom left ) 

            usage: 
                pop3 = pop3d()
                sq = pop3.calc_2d_square_pt(1,[0,0])
                print(sq)
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

    ##-------------------------------------------##    
    def calc_2d_bbox_pts(self, size, origin=None ):
        """ 
            BBOX is 2D on Z axis 

            DEBUG - this is unclamped!
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

    ##-------------------------------------------##
    def calc_icosahedron(self, pos=(0,0,0), rot=(0,0,0), size=1):
        """ DEBUG POS and ROT not working yet 

        
        """

        radius = size/2

        #tmp = object3d()
        tmp_pts =[]

        #// constants
        PI = 3.1415926;
        H_ANGLE = PI/ 180*72;       # 72 degree = 360 / 5
        V_ANGLE = math.atan(1.0/2); # elevation = 26.565 degree

        z  = 0
        xy = 0                            # coords
        hAngle1 = -PI / 2 - H_ANGLE / 2   # start from -126 deg at 1st row
        hAngle2 = -PI / 2                 # start from -90 deg at 2nd row

        # the first top vertex at (0, 0, r)
        tmp_pts.append( (0,0,radius) )

        # compute 10 vertices at 1st and 2nd rows
        for i in range(1,7):
            n = len(tmp_pts) # add to this index each time

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
            
            tmp_pts.append(tuple(vtmp1))
            tmp_pts.append(tuple(vtmp2))
        
        # the last bottom vertex at (0, 0, -r)
        tmp_pts.append( (0,0,-radius) )

        #tmp_pts = tmp.xform_pts(pos=pos, pts=tmp_pts)

        return tmp_pts

    ##-------------------------------------------##
    def calc_circle(self, pos=(0,0,0), rot=(0,0,0), dia=1, axis='z', periodic=True, spokes=23):
        """ spokes = num spokes 
            
            point generators build raw points, not renderable geometry    
        """

        px=0
        py=0
        pz=0 

        ox = 0
        oy = 0
        oz = 0

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
        
        #out = self.xform_pts( pos, out) 
        #out = self.rotate_pts(rot, out) 
        if len(out)==0:
            print('calc_circle: error - no points (bad axis?)')
            raise ValueError('no data') 

        out = self.trs_points(out, translate=pos, rotate=rot)
        if periodic:
            out.append( out[0] )

        return out

    ##-------------------------------------------##  
    def cvt_3d_to_2d(self, points, asvec3=False):
        """ inverse of 2d to 3d - throw out z axis 

        """

        if type(points) is tuple:
            newpts = []
            for pt in points:
                if asvec3:
                    newpts.append(vec3(pt[0], pt[1]))  
                else:
                    newpts.append( (pt[0], pt[1])) 
            return newpts

        else:
            newpts = []
            for pt in points:
                if asvec3:
                    newpts.append( vec3(pt[0], pt[1])   )                
                else:
                    newpts.append( (pt[0], pt[1])   )
            return newpts
        return None
    ##-------------------------------------------##  
    def cvt_2d_to_3d(self, points, asvec3=False, zval=0):
        """ convert a single (tuple) or multiple (list) 2d points 
            into 3d by adding an empty z axis 
            
            this assumes tuple is a single point 

        """

        if type(points) is tuple:
            newpts = []
            for pt in points:
                if asvec3:
                    newpts.append(vec3(pt[0], pt[1], zval))  
                else:
                    newpts.append( (pt[0], pt[1], zval)) 
            return newpts

        else:
            newpts = []
            for pt in points:
                if asvec3:
                    newpts.append( vec3(pt[0], pt[1], zval)   )                
                else:
                    newpts.append( (pt[0], pt[1], zval)   )
            return newpts
        return None
 
     ##-------------------------------------------##  
    def cvt_vec3(self, points):
        """ convert a single (tuple) or multiple (list) into vec3 objects
            this assumes tuple is a single point 

        """

        if type(points) is tuple:
            return vec3(pt[0], pt[1], pt[2])                             
        else:
            newpts = []
            for pt in points:
                newpts.append( vec3(pt[0], pt[1], pt[2])   )                
            return newpts
        return None

    ##-------------------------------------------##
    def locate_pt_along3d(self, fpos, spos, num):
        """
            given two 3D points, return a series of N number connecting points in 3D 

            usage:

        """

        pts_created = []
        #fpos=(x1, y1, z1)
        #spos=(x2, y2, z2)
         
        for n in range(num):
            npos = [3]

            npos[0]     = spos[0]+(((fpos[0]-spos[0])/(num+1))*(n+1))
            npos.append(  spos[1]+(((fpos[1]-spos[1])/(num+1))*(n+1))  )
            npos.append(  spos[2]+(((fpos[2]-spos[2])/(num+1))*(n+1))  )

            pts_created.append( (npos[0], npos[1], npos[2]) )
        return pts_created

    ##-------------------------------------------##
    def trs_points(self, pts, translate=(0,0,0), rotate=(0,0,0), scale=(1,1,1) ):
        #UNTESTED 

        ##--------------
        #rotate

        rx=rotate[0]
        ry=rotate[1]
        rz=rotate[2]

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
 
        pts = self.apply_matrix_pts (pts, m44=rot_matrix) 

        ##--------------
        #translate 
        tmp = []
        for pt in pts: 
            x = pt[0] + translate[0]
            y = pt[1] + translate[1]
            z = pt[2] + translate[2]
            tmp.append( (x,y,z) )

        pts = tmp 

        ##--------------
        #scale 
    
        # build a scale matrix 
        sc_m33 = self.m33.identity
        sc_m33[0]  = scale[0]
        sc_m33[4]  = scale[1]
        sc_m33[8]  = scale[2]    
         
        ################################################
        pts = self.apply_matrix_pts(pts, m33=sc_m33)  # m44=sc_m44 
        
        return pts 






