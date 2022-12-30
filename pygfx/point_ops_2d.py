
import math 
import os 

from gnelscript import NUMPY_IS_LOADED

from gnelscript.pygfx.math_ops import math_util as mu
from gnelscript.pygfx.math_ops import NUMPY_IS_LOADED, matrix22, matrix33, vec2, vec3  


class point_operator_2d(object):
    def __init__(self):
        self.mu   = mu()
        self.dtr  = self.mu.dtr
        
        self.m22  = matrix22()
        self.m33  = matrix33()      
        
        self.vec2     = vec2()     
        self.vec3     = vec3()  

        #self.m44  = matrix44()   

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


    ##-------------------------------------------##
    def extents_fr_bbox(self, bbox, offset=None, periodic=False):
        """ return pt geom from a bbox  
            
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
    def calc_bbox(self, size, origin=None ):
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

    ##-------------------------------------------##  
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

    ##-------------------------------------------##          
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

    ##-------------------------------------------##  
    def batch_transform_pts_2d(self, pts, new_coord, doround=False):
        """ UNTESTED ? - takes a list of tuples and returns a shifted list of tuples , based on a tuple of (x,y) shift"""

        out = []
        for pt in pts:
            if not doround:
                out.append( (pt[0]+ new_coord[0], pt[1]+ new_coord[1] ) ) 
            if doround:
                out.append( ( int(pt[0]+ new_coord[0]), int(pt[1]+ new_coord[1]) ) )         
        return out

    ##-------------------------------------------##  
    def batch_rotate_pts_2d(self, pts, pivot, angle, doround=False):
        out = []
        for pt in pts:
            out.append( self.rotate_point_2d(pt, pivot, angle, doround) )
        return out

    ##-------------------------------------------##  
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

    ##-------------------------------------------##  
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

    ##-------------------------------------------##  
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

    ##-------------------------------------------##  
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

    ##-------------------------------------------##  
    def calc_mid(self, pt1, pt2, num=1, doround=False):
        """ get the midpoint of a 2d vector """

        out = self.locate_pt_along( pt1[0], pt1[1], pt2[0], pt2[1], num )
        if doround:
            return ( int(out[0][0]), int(out[0][1]) )
        if not doround:
            return ( out[0][0], out[0][1] )
                
    ##-------------------------------------------##                  
    def old_calc_theta_vert(self, bot_xy, top_xy):
        try:
            a = bot_xy[1]-top_xy[1]  
            o = bot_xy[0]-top_xy[0]
            if a==0 or o ==0:
                r = self.rtd(math.atan(o)) 
            else:
                r = self.rtd(math.atan(o/a)) 
            return r
        except:
            print('error calc angle - div by 0?') 
            return None 

    ##-------------------------------------------##  
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



class polygon_operator_2d(point_operator_2d):

    def __init__(self):
        super().__init__()  

        self.points          = []    # list of tuples of XYZ points per vertex         -  [(x,y,z), (x,y,z)]  
        self.polygons        = []    # list of tuples for 2 or more vertex connections -  [(1,2,5,8) , (1,2)] 

        self.vec_buffer      = []  # vector work buffer - scratch area to store vectors for operations 

        # render properties embedded in geometry 
        self.linecolors = None       # iterable of colors for lines 
        self.linecolor  = (0,240,00)
        self.vtxcolor   = (0,240,00)

        # "unavoidable side effect" variables 
        self.exprt_ply_idx   = 1     # obj is NOT zero indexed
        #self.exprt_pnt_idx   = 0    # pts ARE zero indexed (everything BUT .obj face idx's are)      

    ##-------------------------------------------## 
    def reset(self):
        self.rot          = [0,0]
        self.pos          = [0,0]
        self.scale        = [1,1]

    ##-------------------------------------------## 
    def apply_matrix_pts(self, pts, m22=None, m33=None):
        """ 
            DEBUG BAD INTERFACE! 
             
            batch mutliply points by a matrix 
            used for translate, rotate, and scaling. 
            
            Can be used for many other things as well.  

        """
      
        tmp_buffer = [] 

        # apply the transform here
        for pt in pts: 

            #print("###### PTS ", type(pt) , pt) 
            #https://web.cse.ohio-state.edu/~shen.94/681/Site/Slides_files/transformation_review.pdf

            if m33 is not None:
                tmp_buffer.append( m33 * (pt[0], pt[1], 1) )
            if m33 is not None:
                tmp_buffer.append( m33 * (pt[0], pt[1], 1) )
            #if m44 is not None:
            #    tmp_buffer.append( m44 * pt )

        return tmp_buffer

    ##-------------------------------------------##
    def trs_points(self, pts, translate=(0,0), rotate=(0,0), scale=(1,1) ):
        """ 
            DEBUG - NOT WORKING YET 

            To combine rotation and translation in one operation one extra dimension is needed than the model requires. 
            For planar things this is 3 components and for spatial things this is 4 components. 
            The operators take 3 components and return 3 components requiring 3x3 matrices.


            GO READ: 
            https://web.cse.ohio-state.edu/~shen.94/681/Site/Slides_files/transformation_review.pdf

        """

        ##--------------
        #rotate

        rx=rotate[0]
        ry=rotate[1]

        # degree to radian function 
        dtr = self.mu.dtr

        # build rotationY (see diagram above) 
        y_matrix =  self.m33.identity
        #y_matrix[0]  =  math.cos(dtr( ry ))
        #y_matrix[2]  = -math.sin(dtr( ry ))
        #y_matrix[8]  =  math.sin(dtr( ry ))
        #y_matrix[10] =  math.cos(dtr( ry ))
              

        # build rotationX (see diagram above) 
        x_matrix =  self.m33.identity
        # x_matrix[5]  =   math.cos(dtr( rx )) 
        # x_matrix[6]  =   math.sin(dtr( rx )) 
        # x_matrix[9]  =  -math.sin(dtr( rx ))
        # x_matrix[10] =   math.cos(dtr( rx ))

        rot_matrix = self.m33.identity
        #rot_matrix = x_matrix * tmp_matr   
 
        pts = self.apply_matrix_pts (pts, m33=rot_matrix) 

        ##--------------
        #translate 
        tmp = []
        for pt in pts: 
            x = pt[0] + translate[0]
            y = pt[1] + translate[1]
            tmp.append( (x,y) )

        pts = tmp 

        ##--------------
        #scale 
    
        # build a scale matrix 
        sc_m33 = self.m33.identity
        # sc_m33[0]  = scale[0]
        # sc_m33[4]  = scale[1]
        # sc_m33[8]  = scale[2]    
         
        ################################################
        pts = self.apply_matrix_pts(pts, m33=sc_m33)  # m44=sc_m44 
        
        return pts 

    ##-------------------------------------------##
    def move(self, x,y ):
        
        # ( pts, translate=(0,0), rotate=(0,0), scale=(1,1) )

        self.points = self.trs_points(self.points, (x,y) )

    ##-------------------------------------------##  
    def scribe(self, str):
        print(str)

    ##-------------------------------------------## 
    ##-------------------------------------------##  

    def load(self, filename):

        """
            copied from pointgen 3d 
            load a 3d object - ignoring the third coordinate 
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
                                    tmp = fid.split('/')
                                    if len(tmp):
                                        # first slash delineated integer is face ID
                                        poly.append(int(tmp[0]))    
                                        # second slash delineated integer is face UV
                                        if tmp[1]: 
                                            uv_poly.append(int(tmp[1]))


                                else:    
                                    poly.append(int(fid))   

                            self.polygons.append( poly )
                            self.uv_polys.append(uv_poly) 

    ##-------------------------------------------##  

    def save(self, filename, as_lines=False):
        ## copied from pointgen 3d 

        buf = [] #array of strings to be written out as the OBJ file

        buf.append("# Created by Magic Mirror render toy.")
        buf.append("# Keith Legg - December 2015.")        
        buf.append("# version2   - November 2018.\n")

        buf.append('\n# Define the vertices')

        for p in self.points:
            if len(p) == 2:
                #add empty Z if 2 otherwise it becomes and error 
                p = (p[0],p[1],0)

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
                buf.append('f %s'%plybuf)
 
        buf.append('\n')

        ################################### 
        #Our filebuffer is an array, we need a string so flatten it 
        output = ''
        for s in buf:
            output=output+s+'\n'

        #save it to disk now
        fobj = open( filename,"w") #encoding='utf-8'
        fobj.write(output)
        fobj.close()


    ##-------------------------------------------## 
    


    ##-------------------------------------------## 




    ##-------------------------------------------## 
    


    ##-------------------------------------------## 




