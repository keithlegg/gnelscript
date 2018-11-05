
from pygfx.point_ops import polygon_operator

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
        data.append('  num face normals   : %s' %  self.numfacnrml )  # unimplemented
        data.append('  num verts          : %s' %  self.numpts     )
        data.append('  num polygons       : %s' %  self.numply     )
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
        return len(self.face_normals)  

    ############################################### 


    def insert(self, obj):
        """ insert an objects geometry into this object 
        """
        
        if isinstance(obj, object3d):
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
        n_ply = self.numply

        #if any polygons are loaded, auto increment the face indecies  
        if n_ply>0:
            n = self.numpts 
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
    def one_vec_to_obj(self, r3, pos=None):
        """ single vector into a renderable 3D line 
            
            can be from world origin or from a point 

        """
        
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

        n = self.numpts # add this number to the indexes in case of existing geom 
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

        n = self.numpts # add this number to the indexes in case of existing geom 

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
    ###############################################  
    #        BUILTIN PRIMITIVE OBJECTS
    ###############################################  
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
    def prim_circle(self, axis='z', pos=(0,0,0), rot=(0,0,0), size=1, spokes = 9):
        """ UNFINSIHED single polygon operations  """    

        pts    = []
        plyidx = []

        pts = self.calc_circle( pos, size, axis, True, spokes )
        
        #we add one because calc_circle_2d returns zero indexed data but OBJ is NOT zero indexed        
        for x in range(spokes):
            plyidx.append(x+1) 

        #print(len(pts))
        #print(plyidx)

        self._insert_poly_idxs(tuple(plyidx))
        self._insert_points(pts)

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
    def prim_arrow(self, axis='z', pos=(0,0,0), rot=(0,0,0), size=1, spokes = 5):
        """ UNFINSIHED single polygon operations  """   
        pass

    ############################################### 
    def prim_cube(self, linecolor=None, pos=(0,0,0), rot=(0,0,0), size=1):
        """ single polygon operations (that can be stacked togteher ?) """
        pts = [];plybfr = []
                
        pivot = 'obj' #make object at origin, then move VS make in place       
        
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

        self._insert_poly_idxs(plybfr) 
        self._insert_points(pts)

        self.rotate_pts( rot )
        if pivot == 'obj':
            self.xform_pts( pos )




