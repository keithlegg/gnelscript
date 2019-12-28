

import math 
import os 

from gnelscript.pygfx.math_ops import math_util as mu


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
            if a==0 or o ==0:
                r = self.rtd(math.atan(o)) 
            else:
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


class object2d(point_operator_2d):
    """
        2.5D object - modeled after obj3d but much simpler
        
        I kept the 3D data but z is unused
        that allows loading and saving of OBJ files and other goodies 

        polygons and objects are combined in this single object  
    """

    def __init__(self):
        super().__init__()  

        self.uv_points       = []   # UV single point coordinates
        self.uv_polys        = []   # UV face indices  
        self.normals         = []   # face_normals

        self.points          = []    # list of tuples of XYZ points per vertex         -  [(x,y,z), (x,y,z)]  
        self.polygons        = []    # list of tuples for 2 or more vertex connections -  [(1,2,5,8) , (1,2)] 

        self.rot             = [0,0,0]
        self.pos             = [0,0,0]
        self.scale           = [1,1,0]

    def reset(self):
        self.rot          = [0,0,0]
        self.pos          = [0,0,0]
        self.scale        = [1,1,1]


    def prim_square(self,  pos=(0,0,0), rot=(0,0,0), size=1):
        pts = [] 

        #keep Z but leave it zero
        pts.append((-size, -size, 0))  
        pts.append((-size,  size, 0))  
        pts.append(( size,  size, 0)) 
        pts.append(( size, -size, 0)) 

        #not a good solution - this doesnt take into account existing geom 
        self.points.extend(pts)
        self.polygons.append( (1,2,3,4) )


    def prim_triangle(self, pos=(0,0,0), rot=(0,0,0), size=1):
        pts =  [(-size,0,0), (0,size,0), (size,0,0) ]
        poly = [(1,2,3)]

        #not a good solution - this doesnt take into account existing geom 
        self.points.extend(pts)
        self.polygons.append( (1,2,3) )


    # def prim_circle(self,  center=(0,0,0), dia=1):
    #     pass


    ###############################################  

    def load(self, filename):
        ## copied from pointgen 3d 

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

    ###############################################  

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



###################################################

class tnode(object):
    def __init__(self):
        self.x=0
        self.y=0
    
    def get(self):
        return (self.x, self.y)
    
    def set(self,x,y):
        self.x = x
        self.y = y


def fractal(tree, depth, maxdepth):

    growth = depth*3

    frwrd = 40-growth
    brfrwrd = 20+growth
    brlen = 30-growth

    #print('depth is %s'%depth)
    if depth==maxdepth:
        return None 

    else:    
        # #trunk
        # if depth==0:

        # else:

        last = tree[len(tree)-1]
        # print('## last tree xy is ', last.x, last.y )

        # advance forward 
        tn = tnode()
        tn.set(last.x, last.y+frwrd) 
        tree.append(tn)

        newlast = tree[len(tree)-1]

        # right side of tree
        tn = tnode()
        tn.set(newlast.x+brlen*2, newlast.y+brfrwrd) 
        tree.append(tn)

        # return to stalk 
        tn = tnode()
        tn.set(newlast.x, newlast.y) 
        tree.append(tn)

        # left side of tree
        tn = tnode()
        tn.set(newlast.x-brlen*2, newlast.y+brfrwrd) 
        tree.append(tn)

        # return to stalk 
        tn = tnode()
        tn.set(newlast.x, newlast.y) 
        tree.append(tn)

        depth += 1
        fractal(tree, depth, maxdepth)        


def tree_to_lines(tree):
    """ iterate a tree and get the points XY coords to draw"""

    out_pts = []
    ct = 0 
    for i,t in enumerate(tree):
        #print(' # node %s x:%s y:%s '%(ct, t.x ,t.y ) )
        if i>0:
            last = tree[i-1]    
            out_pts.append( ( (last.x, last.y), (t.x, t.y) ) )
        ct += 1

    return out_pts 


