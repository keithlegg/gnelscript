

import os 
import math 

#from gnelscript import NUMPY_IS_LOADED

#from gnelscript.pygfx.math_ops import math_util as mu
#from gnelscript.pygfx.math_ops import NUMPY_IS_LOADED, matrix22, matrix33, vec2, vec3  

from gnelscript.pygfx.point_ops_2d import pop2d




class polygon_operator_2d(pop2d):

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


    def ray_hit_2d(self):

        """
        # https://rosettacode.org/wiki/Ray-casting_algorithm#Python

        from collections import namedtuple
        from pprint import pprint as pp
        import sys

        Pt = namedtuple('Pt', 'x, y')               # Point
        Edge = namedtuple('Edge', 'a, b')           # Polygon edge from a to b
        Poly = namedtuple('Poly', 'name, edges')    # Polygon

        _eps = 0.00001
        _huge = sys.float_info.max
        _tiny = sys.float_info.min

        def rayintersectseg(p, edge):
            ''' takes a point p=Pt() and an edge of two endpoints a,b=Pt() of a line segment returns boolean
            '''
            a,b = edge
            if a.y > b.y:
                a,b = b,a
            if p.y == a.y or p.y == b.y:
                p = Pt(p.x, p.y + _eps)

            intersect = False

            if (p.y > b.y or p.y < a.y) or (
                p.x > max(a.x, b.x)):
                return False

            if p.x < min(a.x, b.x):
                intersect = True
            else:
                if abs(a.x - b.x) > _tiny:
                    m_red = (b.y - a.y) / float(b.x - a.x)
                else:
                    m_red = _huge
                if abs(a.x - p.x) > _tiny:
                    m_blue = (p.y - a.y) / float(p.x - a.x)
                else:
                    m_blue = _huge
                intersect = m_blue >= m_red
            return intersect

        def _odd(x): return x%2 == 1

        def ispointinside(p, poly):
            ln = len(poly)
            return _odd(sum(rayintersectseg(p, edge)
                            for edge in poly.edges ))

        def polypp(poly):
            print ("\n  Polygon(name='%s', edges=(" % poly.name)
            print ('   ', ',\n    '.join(str(e) for e in poly.edges) + '\n    ))')

        if __name__ == '__main__':
            polys = [
              Poly(name='square', edges=(
                Edge(a=Pt(x=0, y=0), b=Pt(x=10, y=0)),
                Edge(a=Pt(x=10, y=0), b=Pt(x=10, y=10)),
                Edge(a=Pt(x=10, y=10), b=Pt(x=0, y=10)),
                Edge(a=Pt(x=0, y=10), b=Pt(x=0, y=0))
                )),
              Poly(name='square_hole', edges=(
                Edge(a=Pt(x=0, y=0), b=Pt(x=10, y=0)),
                Edge(a=Pt(x=10, y=0), b=Pt(x=10, y=10)),
                Edge(a=Pt(x=10, y=10), b=Pt(x=0, y=10)),
                Edge(a=Pt(x=0, y=10), b=Pt(x=0, y=0)),
                Edge(a=Pt(x=2.5, y=2.5), b=Pt(x=7.5, y=2.5)),
                Edge(a=Pt(x=7.5, y=2.5), b=Pt(x=7.5, y=7.5)),
                Edge(a=Pt(x=7.5, y=7.5), b=Pt(x=2.5, y=7.5)),
                Edge(a=Pt(x=2.5, y=7.5), b=Pt(x=2.5, y=2.5))
                )),
              Poly(name='strange', edges=(
                Edge(a=Pt(x=0, y=0), b=Pt(x=2.5, y=2.5)),
                Edge(a=Pt(x=2.5, y=2.5), b=Pt(x=0, y=10)),
                Edge(a=Pt(x=0, y=10), b=Pt(x=2.5, y=7.5)),
                Edge(a=Pt(x=2.5, y=7.5), b=Pt(x=7.5, y=7.5)),
                Edge(a=Pt(x=7.5, y=7.5), b=Pt(x=10, y=10)),
                Edge(a=Pt(x=10, y=10), b=Pt(x=10, y=0)),
                Edge(a=Pt(x=10, y=0), b=Pt(x=2.5, y=2.5))
                )),
              Poly(name='exagon', edges=(
                Edge(a=Pt(x=3, y=0), b=Pt(x=7, y=0)),
                Edge(a=Pt(x=7, y=0), b=Pt(x=10, y=5)),
                Edge(a=Pt(x=10, y=5), b=Pt(x=7, y=10)),
                Edge(a=Pt(x=7, y=10), b=Pt(x=3, y=10)),
                Edge(a=Pt(x=3, y=10), b=Pt(x=0, y=5)),
                Edge(a=Pt(x=0, y=5), b=Pt(x=3, y=0))
                )),
              ]
            testpoints = (Pt(x=5, y=5), Pt(x=5, y=8),
                          Pt(x=-10, y=5), Pt(x=0, y=5),
                          Pt(x=10, y=5), Pt(x=8, y=5),
                          Pt(x=10, y=10))
            
            print ("\n TESTING WHETHER POINTS ARE WITHIN POLYGONS")
            for poly in polys:
                polypp(poly)
                print ('   ', '\t'.join("%s: %s" % (p, ispointinside(p, poly))
                                       for p in testpoints[:3]))
                print ('   ', '\t'.join("%s: %s" % (p, ispointinside(p, poly))
                                       for p in testpoints[3:6]))
                print ('   ', '\t'.join("%s: %s" % (p, ispointinside(p, poly))
                                       for p in testpoints[6:]))
        """    
        pass 
        
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