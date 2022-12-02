
import math 


from gnelscript.pygfx.point_ops import polygon_operator
from gnelscript.pygfx.math_ops import vec3 


###############################################

class object3d(polygon_operator):

    def __init__(self):
        super().__init__()  

        #self.geom_history = []
 
        self.rot       = [0,0,0]
        self.pos       = [0,0,0]
        self.scale     = [1,1,1]

        #self.uv_points       = []
        #self.uv_polys        = []
        #self.normals         = []

    def reset(self):
        self.rot          = [0,0,0]
        self.pos          = [0,0,0]
        self.scale        = [1,1,1]


    @property
    def rotation(self):
        return self.rot

    @property
    def position(self):
        return self.pos

    ############################################### 
    def copy(self):
        new = type(self)()
        new.points   = self.points
        new.polygons = self.polygons  
        return new

    ############################################### 
    def insert_2pp(self, obj):
        """ insert tuple as two point polygon (line) 
            
        """

        # if tuple or list assume its [polyidx, points]
        if isinstance(obj, tuple) or isinstance(obj, list):
            if replace is True:
                pass    
                #self.points=obj[0]
                #self.polygons=obj[1] 

            else:
                self.insert_polygons(obj[0], obj[1])

    ############################################### 
    def insert(self, obj, replace=False):
        """ insert an objects geometry into this object 
            
        """

        # if tuple or list assume its [polyidx, points]
        if isinstance(obj, tuple) or isinstance(obj, list):
            if replace is True:
                self.points=obj[0]
                self.polygons=obj[1]                 
            else:
                self.insert_polygons(obj[0], obj[1])

        if isinstance(obj, object3d):
            if replace is True:
                self.points=obj.points
                self.polygons=obj.polygons                 
            else:
                self.insert_polygons(obj.polygons, obj.points)

    ############################################### 
    def append(self, otherobj):
        """ add another object to this one 
            points can be added, paying attention to the index 
            values. 

            polys have to be re-indexed to match the higher index values
        """

        for pt in otherobj.points:
            self.points.append(pt)

        for ply in otherobj.polygons:
            self.polygons.append(self._reindex_ply(ply, self.numpts))  

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
        data.append('  num face normals   : %s' %  self.numfacnrml )   
        data.append('  num verts          : %s' %  self.numpts     )
        data.append('  num polygons       : %s' %  self.numply     )
        
        data.append('  num UV points      : %s' %  len(self.uv_points)    )
        data.append('  num UV faces       : %s' %  len(self.uv_polys)     )

        data.append(' --------------------------- ' )         
        data.append('  num triangles      : %s' %  tris  )
        data.append('  num quads          : %s' %  quads )  
        data.append('  num other          : %s' %  other ) 
        data.append('############################\n')

        for d in data:
            print(d)

    ############################################### 
    @property
    def numply(self):
        return len(self.polygons)

    @property
    def numpts(self):
        return len(self.points)   

    @property
    def numfacnrml(self):
        return len(self.normals)  

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
        #tmp.save('vtx_vectrz.obj')
        return vectrx

    ############################################### 
    def calc_face_normals(self):
        """ calculate the normals for each face of object
            only tested for 3 and 4 sided polygons   
        """
        
        vectrx = []

        #secondary tweaks to the normal data 
        scale     = 1
        normalize = False  #make each face normal unit length 

        cache_face_normals = [] 

        # iterate each face and convert eart vertex into a vec3 
        for idx in range(self.numply):
            f_nrml = self.get_face_normal(idx)
            #print('############ ', f_nrml)
            # store it in object for later use 
            self.normals.append(f_nrml) 

        return cache_face_normals

    ############################################### 
    def one_vec_to_arrow(self, r3, pos=None):
        """ single vector into a renderable 3D line 
            
            can be from world origin or from a point 

        """
        
        # get the mag of the vector 

        # build the geom for the model 

        # create a matrix from the vec??

        # multiply the obj by the matrix?


        pts = [
               (0    , 0    , 0    ),
               (r3[0], r3[1], r3[2]), 
              ]


        ############                              
        n = self.numpts # add this number to the indices in case of existing geom 
        plyidx = [(n+1,n+2)]
        #append points to internal 
        for p in pts:
            self.points.append(p)
        for vec in plyidx:    
            self.polygons.append( vec )  

    ############################################### 
    def one_vec_to_obj(self, r3, pos=None, arrowhead=False):
        """ single vector into a renderable 3D line 
            
            can be from world origin or from a point 

        """

        #if isinstance(r3,vec4):
        #    r3 = [r3[0]
        

        if pos is not None:
            pts = [
                   (pos[0]       , pos[1]      , pos[2]       ),
                   (pos[0]+r3[0] , pos[1]+r3[1], pos[2]+r3[2] ),                   
                  ]

        if pos is None:    
            if arrowhead is False:
                pts = [
                       (0    , 0    , 0    ),
                       (r3[0], r3[1], r3[2]), 
                      ]
            if arrowhead is True:
                pts = [
                       (0    , 0    , 0    ),
                       (r3[0], r3[1], r3[2]), 
                      ]

        ############                              
        n = self.numpts # add this number to the indices in case of existing geom 
        plyidx = [(n+1,n+2)]
        #append points to internal 
        for p in pts:
            self.points.append(p)
        for vec in plyidx:    
            self.polygons.append( vec )  

    ############################################### 
    """
    def vector_between_two_obj(self, r3_1, r3_2):
        # a vector between two other vectors 
        # probably not useful, but interesting 
        pts = [
               (r3_1[0], r3_1[1], r3_1[2]),
               (r3_2[0], r3_2[1], r3_2[2]), 
              ]

        n = self.numpts # add this number to the indexes in case of existing geom 
        plyidx = [(n+1,n+2)]
        #append points to internal 
        for p in pts:
            self.points.append(p)
        for vec in plyidx:    
            self.polygons.append( vec )  
    """
    ############################################### 

    #def edgegeom_to_vectorlist(self, geom):

    ############################################### 
    def pts_to_linesegment(self, pt_list):

        for i,pt in enumerate(pt_list):
            
            if pt == None:
                return None

            if i>0:
                n = self.numpts # add this number to the indices in case of existing geom 
                plyidx = [(n+1,n+2)]
                
                pts = [ pt_list[i-1], pt_list[i] ]

                #append points to internal 
                for p in pts:
                    self.points.append(p)
                for vec in plyidx:    
                    self.polygons.append( vec )                          

    ############################################### 
    def vectorlist_to_obj(self, vecs, pos=None):
        """ take a list of vectors and convert it to renderable geometry 
        
            vecs:
            - can be a list of single values (vec3)
              or a list 2 two values         (vec,position)
                  -- you can also pass a POS seperately , like in case of a vec3 type 

        """
        
        if not isinstance(vecs, list):
            print ("## error vectorlist_to_obj - only accepts list of vectors ")
            return None         

        for v in vecs:

            # somehow this gets recursive if None - exit just in case 
            if v == None:
                return None 

            if isinstance(v, vec3):
                self.one_vec_to_obj( (v[0],v[1],v[2]), pos  ) 

            if isinstance(v,tuple) or isinstance(v, list):    
                if len(v) == 1:
                    self.one_vec_to_obj(v, pos) 
                if len(v) == 2:
                    self.one_vec_to_obj(v[0], v[1])                 
               

    ###############################################  
    ###############################################  
    #        BUILT IN PRIMITIVE OBJECTS
    ###############################################  
    ###############################################  
 
    #def prim_cylinder(self, axis='z', pos=(0,0,0), rot=(0,0,0), dia=1, spokes = 9):
    #    pass  
    
    ############################################### 
    def prim_line(self, axis, pos, rot, size=1):
        """ 3d lines, 2 point polygons """

        if axis=='x':
            pts =[ (-size,0,0), (size,0,0) ]
        if axis=='y':
            pts =[ (0,-size,0), (0,size,0) ]
        if axis=='z':
            pts =[ (0,0,-size) , (0,0,size) ]

        poly = [(1,2)]
       
        self.insert_polygons(poly, pts)
        #self.rotate_pts( rot )
        #self.xform_pts( pos )

    ###############################################  
    def prim_triangle(self, axis, pos, rot, size=1):
        """ single polygon operations (that can be stacked together ?) """

        if axis=='x':
            pts =  [(0,-size,0), (0,0,size), (0,size,0) ]
        if axis=='y':
            pts =  [(0,0,-size), (size,0,0), (0,0,size) ]
        if axis=='z':
            pts =  [(-size,0,0), (0,size,0), (size,0,0) ]

        poly = [(1,2,3)]
       
        self.insert_polygons(poly, pts)
        #self.rotate_pts( rot )
        #self.xform_pts( pos )
        
    ###############################################  
    def prim_quad(self, axis='y', pos=(0,0,0), rot=(0,0,0), size=1):
        """ single polygon operations (that can be stacked together ?) """
        
        if axis == 'x':
            pts = [(0,-size,-size), (0,-size,size), (0,size,size), (0,size,-size) ] #X AXIS

        if axis == 'y':
            pts = [(-size,0,-size), (-size,0,size), (size,0,size), (size,0,-size) ] #Y AXIS
            
        if axis == 'z':
            pts = [(-size,-size,0), (-size,size,0), (size,size,0), (size,-size,0) ] #Z AXIS

        poly    = [(1,2,3,4)]
       
        self.insert_polygons(poly, pts)
        #self.rotate_pts( rot )
        #self.xform_pts( pos )

    ###############################################  
    def prim_circle(self, axis, pos=(0,0,0), rot=(0,0,0), dia=1, spokes=9):
        """ UNFINSIHED single polygon operations  """    

        # print('## debug prim circle ', axis , pos   )

        pts  = []
        poly = []

        # calc_circle( pos=(0,0,0), dia=1, axis='z', periodic=True, spokes=23):
        pts = self.calc_circle( dia=dia, axis=axis, periodic=True, spokes=spokes)
        
        # we add one because calc_circle_2d returns zero indexed data but OBJ is NOT zero indexed        
        for x in range(spokes):
            poly.append(x+1) 

        pts = self.xform_pts(  pos, pts)
        pts = self.rotate_pts( rot, pts)

        self.insert_polygons([tuple(poly)], pts)

    ###############################################
    def prim_cone(self, axis, pos=(0,0,0), rot=(0,0,0), dia=1, spokes=8):
        """ first prim tool to use other tools and prims 
            made so we can make an arrow prim 
        """

        # print("## debug pos cone ", pos )

 
        self.prim_circle(axis=axis, pos=pos, dia=dia, spokes=spokes)
        
        tiplen = dia*2

        #debug = use a normal instead of world coords
        #this currently will not work if the circle is rotated !! 
        if axis=='x':
            oset = (tiplen,0,0)
        if axis=='y':
            oset = (0,tiplen,0)            
        if axis=='z':
            oset = (0,0,tiplen) 

        self.radial_triangulate_face(0, offset=oset )

        #pts = self.xform_pts(  pos, pts)
        #pts = self.rotate_pts( rot, pts)

    ############################################### 
    #def prim_sphere2(self, pos=(0,0,0), rot=(0,0,0), size=1 ):

    ############################################### 
    def prim_sphere(self, pos, rot, size=1 ):
        
        #UNFINISHED 

        #radius is called size for uniformity in ARGS 
        radius = size

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
        fid = self.numpts-1

        faces = [] 

        # compute 10 vertices at 1st and 2nd rows
        for i in range(1,8):
            n = self.numpts # add to this index each time

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
        
        #pts = self.xform_pts(  pos, pts)
        #pts = self.rotate_pts( rot, pts)



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

        polys = [(1,2),
                  (3,4),
                  (5,6)
                 ]
        
        self.linecolors = [
                 (255,0,0),         
                 (0,255,0), 
                 (0,0,255)                  
        ]

        self.insert_polygons(polys, pts) 

        #self.rotate_pts( rot )
        #self.xform_pts( pos )
    ###############################################  
    def prim_locator_color(self, color=(255,0,0), pos=(0,0,0), rot=(0,0,0), size=1):
   
        pts = [
               # x axis indicator  
               (pos[0],       pos[1],        pos[2], color[0], color[1], color[2]),
               (pos[0]+size,  pos[1],        pos[2], color[0], color[1], color[2]), 
               # y axis indicator
               (pos[0],       pos[1],        pos[2], color[0], color[1], color[2]),
               (pos[0],       pos[1]+size,   pos[2], color[0], color[1], color[2]), 
               # z axis indicator
               (pos[0],       pos[1],        pos[2], color[0], color[1], color[2] ),
               (pos[0],       pos[1],        pos[2]+size, color[0], color[1], color[2])                
              ]

        polys = [(1,2),
                  (3,4),
                  (5,6)
                 ]
        
        self.linecolors = [
                 (255,0,0),         
                 (0,255,0), 
                 (0,0,255)                  
        ]

        self.insert_polygons(polys, pts) 

        #self.rotate_pts( rot )
        #self.xform_pts( pos )
    ###############################################  
    def prim_locator_xyz(self, pos, rot, size=1):
       
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
      
        
        polys = [ (1,2),
                   (3,4),

                   (5,6),

                   (7,8),
                   (9,10,11),                   
                   (12,13),
                   (14,15),
                   (16,17,18,19)                   
                 ]

        self.insert_polygons(polys, pts) 
        #self.rotate_pts( rot )
        #self.xform_pts( pos )

    ###############################################  
    def prim_arrow(self, axis='y', pos=(0,0,0), rot=(0,0,0), vec3=None, size=1, pivot='obj'): 
        """ fully 3D model of an arrow 
            will be used for visualizing vectors 
        """

        spokes  = 4  # num turns around axis 

        length  = 1

        dia      = .1        # extrude length is double this, or .2 
        shaftlen = length-.2 # cone is .2, that plus this = 1 
        
        tmpobj1 = object3d()

        if axis=='x':
            tmpobj1.prim_cone( axis=axis, pos=(shaftlen,0,0), dia=dia, spokes=spokes )
        if axis=='y':
            tmpobj1.prim_cone( axis=axis, pos=(0,shaftlen,0), dia=dia, spokes=spokes )        
        if axis=='z':
            tmpobj1.prim_cone( axis=axis, pos=(0,0,shaftlen), dia=dia, spokes=spokes )

        tmpobj = object3d()
        tmpobj.prim_circle( axis=axis, pos=(0,0,0), spokes=spokes , dia=dia/5)

        # normal is flipped wrong only on some axis  
        # look into why this happens - extrude used face normal its probably that 
        if axis=='y':
            tmpobj.extrude_face(0, distance=shaftlen)            
        else:
            tmpobj.extrude_face(0, distance=-shaftlen)

        tmpobj1.insert(tmpobj)
        
        #tmpobj1.points = self.rotate_pts( rot , tmpobj1.points )
        #tmpobj1.points = self.xform_pts( pos , tmpobj1.points)

        self.insert(tmpobj1) 

    ###############################################  
    #def prim_line_arrow(self, axis, pos, rot, size=1):
    #    """ UNFINSIHED make an arrow out of line geom without polygons """   
    #    pass

    ############################################### 
    def prim_cube(self, linecolor=(0,128,0), pos=(0,0,0), rot=(0,0,0), size=1, pivot='obj'):
        """ single polygon operations (that can be stacked togteher ?) """
        pts = [];plybfr = []
                
        # pivot = 'obj' #make object at origin, then move VS make in place       
        
        if pivot == 'world':  
            """ build it in place - making the pivot point at world zero """
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

        if pivot == 'obj':
            """ build it center of world THEN move it into place """

            # by adding position argument we can create it in different places
            pts.append( (-size, -size, size)  ) #vertex 0
            pts.append( (-size,  size, size)  ) #vertex 1
            pts.append( ( size,  size, size)  ) #vertex 2  
            pts.append( ( size, -size, size)  ) #vertex 3
            # notice these next 4 are the same coordinates with a negative Z instead!
            pts.append( (-size, -size, -size)  ) #vertex 4
            pts.append( (-size,  size, -size)  ) #vertex 5
            pts.append( ( size,  size, -size)  ) #vertex 6  
            pts.append( ( size, -size, -size)  ) #vertex 7


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

        self.insert_polygons(plybfr, pts)

        #self.rotate_pts( rot )
        #if pivot == 'obj':
        #    self.xform_pts( pos )




