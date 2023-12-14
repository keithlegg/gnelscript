# ***** BEGIN GPL LICENSE BLOCK *****
#
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# ***** END GPL LICENSE BLOCK *****


from gnelscript import NUMPY_IS_LOADED

if NUMPY_IS_LOADED:
    import numpy as np             
    #from numpy.linalg import inv

import math 



from gnelscript.pygfx.point_ops import polygon_operator
from gnelscript.pygfx.math_ops import vec3, matrix33, matrix44  


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

   
    def __mul__(self, matrix):
        """ you can multiply an object by a matrix (or vector?) to transform- cool huh? """

        #if type(matrix)is vec3:
        #    self.points = self.apply_matrix_pts(pts=self.points, m33=matrix) 

        print(type(matrix)) 
        
        if type(matrix)is matrix33:
            self.points = self.apply_matrix_pts(pts=self.points, m33=matrix) 
        if type(matrix)is matrix44:
            self.points = self.apply_matrix_pts(pts=self.points, m44=matrix) 
        if NUMPY_IS_LOADED:
            if type(matrix)==np.ndarray:
                if matrix.size==9:
                    m33=matrix33()
                    self.points = self.apply_matrix_pts(pts=self.points, m33=m33.from_np(matrix)) 
                if matrix.size==16:
                    m444=matrix44()
                    self.points = self.apply_matrix_pts(pts=self.points, m44=m44.from_np(matrix)) 

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

    ##------------------------------------------------##
    def copy(self):
        new = type(self)()
        new.points   = self.points
        new.polygons = self.polygons  
        return new

    ##------------------------------------------------##
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

    ##------------------------------------------------##
    def insert(self, obj, replace=False, ans=False, pos=None):
        """ insert an object's geometry into this object 
            
            DEBUG - deafults to sharing the point indeces 
            ans/asnew_shell was added but seems to be broken 

        """
        #DEBUG - check for data in object first 

        # VEC3 TYPE - UNTESTED
        if isinstance(obj, vec3):
            if pos:            
                self.one_vec_to_obj(obj, pos=pos, arrowhead=False)
                #self.insert_line(obj)
                #self.linegeom_fr_points( [(0,0,0), tuple(obj.aspt)] )
                pass 
            else:            
                self.one_vec_to_obj(obj)

        # GEOMTYPE tuple or list assume its [polyidx, points]
        if isinstance(obj, tuple) or isinstance(obj, list):
            if replace is True:
                if pos:
                    tmp=[]
                    for pt in obj[0]:
                        tmp.append((pt[0]+pos[0], pt[1]+pos[1], pt[2]+pos[2])) 
                    self.points=tmp
                else:
                    self.points=obj[0]                    
                self.polygons=obj[1]                 
            else:
                self.insert_polygons(plyids=obj[0], points=obj[1])

        # object3d type  
        if isinstance(obj, object3d):
            if replace is True:
                if pos:
                    tmp=[]
                    for pt in obj[0]:
                        tmp.append((pt[0]+pos[0], pt[1]+pos[1], pt[2]+pos[2])) 
                    self.points=tmp
                else:
                    self.points=obj[0]  
            
                if len(obj.polygons):
                    self.polygons=obj.polygons                 
            else:
                if pos:
                    self.insert_polygons(plyids=obj.polygons, points=obj.points, pos=pos)
                else:   
                    self.insert_polygons(plyids=obj.polygons, points=obj.points)

    ##------------------------------------------------##
    def append(self, otherobj):
        """ add another object to this one 
            points can be added, paying attention to the index 
            values. 

            polys have to be re-indexed to match the higher index values
        """
        if type(otherobj)==self:
            for pt in otherobj.points:
                self.points.append(pt)

            for ply in otherobj.polygons:
                self.polygons.append(self._reindex_ply(ply, self.numpts))  

    ##------------------------------------------------## 
    def show(self):
        data = []
 
        tris = 0
        quads = 0
        lines = 0
        other = 0   

        for p in self.polygons:
            if len(p)==3:
                tris+=1
            elif len(p)==4:
                quads+=1  
            elif len(p)==2:
                lines+=1  
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
        data.append('  num lines          : %s' %  lines ) 
        data.append('  num other          : %s' %  other ) 
        data.append('############################\n')

        for d in data:
            print(d)

    ##-------------------------------------------
    @property
    def numply(self):
        return len(self.polygons)

    @property
    def numpts(self):
        return len(self.points)   

    @property
    def numfacnrml(self):
        return len(self.normals)  

    ##-------------------------------------------
    def cvt_pts_2d(self, axis='xy'):
        pts2d = [] 
        for pt in self.points: 
            if axis=='xy':
                pts2d.append((pt[0],pt[1]))
            if axis=='yz':
                pts2d.append((pt[1],pt[2]))
            if axis=='xz':
                pts2d.append((pt[0],pt[2]))
        return pts2d
            

    ##-------------------------------------------
    def cvt_pts_vec3(self):
        """ return a list of pygfx vec3 objects for each vertex """

        vectrx = []
        for pt in self.points:
            v = vec3().insert(pt)
            vectrx.append(v.copy(vtype='tuple'))
        return vectrx

    ##-------------------------------------------
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

    ##------------------------------------------------##
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

    ##------------------------------------------------##
    def one_2d_vec_to_obj(self, r2, pos=None, arrowhead=False):
        if pos==None:
            self.one_vec_to_obj((r2[0],r2[1],0) )
        else:
            self.one_vec_to_obj((r2[0],r2[1],0), pos=(pos[0],pos[1],0) )

    def one_vec_to_obj(self, r3, pos=None, arrowhead=False):
        """ single vector into a renderable 3D line 
            
            can be from world origin or from a point 

        """

        if type(r3)is tuple or type(r3)is list:
            r3 = vec3(r3)

        #calc_circle(self, pos=(0,0,0), rot=(0,0,0), dia=1, axis='z', periodic=True, spokes=4):
        #circle = self.calc_circle(periodic=True, spokes=4) 
        #print(circle)   
        
        if pos is not None:
            pts = [
                   (pos[0]       , pos[1]      , pos[2]       ),
                   (pos[0]+r3[0] , pos[1]+r3[1], pos[2]+r3[2] ),                   
                  ]

        if pos is None:    
            if arrowhead is False:
                pts = [
                       (0    , 0    , 0    ),
                       (r3[0], r3[1], r3[2]) 
                      ]
            if arrowhead is True:
                pts = [
                       (0    , 0    , 0    ),
                       (r3[0], r3[1], r3[2])

                      ]

        ############                              
        n = self.numpts # add this number to the indices in case of existing geom 
        
        if arrowhead is False:
            plyidx = [(n+1,n+2)]
        if arrowhead:
            plyidx = [(n+1,n+2)] 

        #append points to internal 
        for p in pts:
            self.points.append(p)
        for vec in plyidx:    
            self.polygons.append( vec )  

    ##------------------------------------------------##
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
    ##------------------------------------------------##

    #def edgegeom_to_vectorlist(self, geom):

    ##------------------------------------------------##
    def pts_to_linesegment(self, pt_list, periodic=False):
        """
           iterate a group of points and return insert redenrable geom geom into self 
        """

        for i,pt in enumerate(pt_list):
            if type(pt) is vec3:
                pt = pt.aspt

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

        if periodic:
            n = self.numpts # add this number to the indices in case of existing geom 
            self.polygons.append( (n, 1) )
        

    ##------------------------------------------------##
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

            # DEBUG somehow this gets recursive if None - exit just in case 
            if v == None:
                return None 

            if isinstance(v, vec3):
                self.one_vec_to_obj( (v[0],v[1],v[2]), pos  ) 

            if isinstance(v,tuple) or isinstance(v, list):    
                if len(v) == 1:
                    self.one_vec_to_obj(v, pos) 
                if len(v) == 2:
                    self.one_vec_to_obj(v[0], v[1])                 
               

    ##------------------------------------------------## 
    ##------------------------------------------------## 
    #        BUILT IN PRIMITIVE OBJECTS
    ##------------------------------------------------## 
    ##------------------------------------------------## 
 
    #def prim_cylinder(self, axis='z', pos=(0,0,0), rot=(0,0,0), dia=1, spokes = 9):
    #    pass  
    
    ##------------------------------------------------##
    def prim_line(self, axis='y', pos=(0,0,0), rot=(0,0,0), size=1):
        """ 3d lines, 2 point polygons """

        if axis=='x':
            pts =[ (-size,0,0), (size,0,0) ]
        if axis=='y':
            pts =[ (0,-size,0), (0,size,0) ]
        if axis=='z':
            pts =[ (0,0,-size) , (0,0,size) ]

        poly = [(1,2)]
       
        pts = self.rotate_pts( rot, pts )
        pts = self.xform_pts( pos, pts )
        self.insert_polygons(poly, pts)

    ##------------------------------------------------## 
    def prim_triangle(self, axis='z', pos=(0,0,0), rot=(0,0,0), size=1, flip=False):
        """ single polygon operations (that can be stacked together ?) """

        if axis=='x':
            pts =  [(0,-size,0), (0,0,size), (0,size,0) ]
        if axis=='y':
            pts =  [(0,0,-size), (size,0,0), (0,0,size) ]
        if axis=='z':
            pts =  [(-size,0,0), (0,size,0), (size,0,0) ]

        if flip:
            poly = [(3,2,1)]        
        else:
            poly = [(1,2,3)]
       

        pts = self.rotate_pts( rot, pts )
        pts = self.xform_pts( pos, pts )
        self.insert_polygons( poly, pts)

    ##------------------------------------------------##  
    def prim_quad(self, axis='y', pos=(0,0,0), rot=(0,0,0), size=1):
        """ single polygon operations (that can be stacked together ?) """
        
        if axis == 'x':
            pts = [(0,-size,-size), (0,-size,size), (0,size,size), (0,size,-size) ] #X AXIS

        if axis == 'y':
            pts = [(-size,0,-size), (-size,0,size), (size,0,size), (size,0,-size) ] #Y AXIS
            
        if axis == 'z':
            pts = [(-size,-size,0), (-size,size,0), (size,size,0), (size,-size,0) ] #Z AXIS

        poly    = [(1,2,3,4)]
       
        pts = self.rotate_pts( rot, pts )
        pts = self.xform_pts( pos, pts )
        self.insert_polygons(poly, pts)

    ##------------------------------------------------##  
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

    ##------------------------------------------------## 
    def prim_cone(self, axis='y', pos=(0,0,0), rot=(0,0,0), dia=1, spokes=8):
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

    ##------------------------------------------------## 
    #def prim_sphere2(self, pos=(0,0,0), rot=(0,0,0), size=1 ):
    """
    http://www.songho.ca/opengl/gl_sphere.html

    // clear memory of prev arrays
    std::vector<float>().swap(vertices);
    std::vector<float>().swap(normals);
    std::vector<float>().swap(texCoords);

    float x, y, z, xy;                              // vertex position
    float nx, ny, nz, lengthInv = 1.0f / radius;    // vertex normal
    float s, t;                                     // vertex texCoord

    float sectorStep = 2 * PI / sectorCount;
    float stackStep = PI / stackCount;
    float sectorAngle, stackAngle;

    for(int i = 0; i <= stackCount; ++i)
    {
        stackAngle = PI / 2 - i * stackStep;        // starting from pi/2 to -pi/2
        xy = radius * cosf(stackAngle);             // r * cos(u)
        z = radius * sinf(stackAngle);              // r * sin(u)

        // add (sectorCount+1) vertices per stack
        // the first and last vertices have same position and normal, but different tex coords
        for(int j = 0; j <= sectorCount; ++j)
        {
            sectorAngle = j * sectorStep;           // starting from 0 to 2pi

            // vertex position (x, y, z)
            x = xy * cosf(sectorAngle);             // r * cos(u) * cos(v)
            y = xy * sinf(sectorAngle);             // r * cos(u) * sin(v)
            vertices.push_back(x);
            vertices.push_back(y);
            vertices.push_back(z);

            // normalized vertex normal (nx, ny, nz)
            nx = x * lengthInv;
            ny = y * lengthInv;
            nz = z * lengthInv;
            normals.push_back(nx);
            normals.push_back(ny);
            normals.push_back(nz);

            // vertex tex coord (s, t) range between [0, 1]
            s = (float)j / sectorCount;
            t = (float)i / stackCount;
            texCoords.push_back(s);
            texCoords.push_back(t);
        }
    }
    """
    
    ##------------------------------------------------## 
    def _icosahedron(self, pos, radius, build_wire=False, build_geom=False):
        """ should be in point gen?
            builds the points but not polygons for an icosahedron
            there are some experimental geom functions here  
        """
    
        tmp = object3d()

        #// constants
        PI = 3.1415926;
        H_ANGLE = PI/ 180*72;       # 72 degree = 360 / 5
        V_ANGLE = math.atan(1.0/2); # elevation = 26.565 degree

        z  = 0
        xy = 0                            # coords
        hAngle1 = -PI / 2 - H_ANGLE / 2   # start from -126 deg at 1st row
        hAngle2 = -PI / 2                 # start from -90 deg at 2nd row

        # the first top vertex at (0, 0, r)
        tmp.points.append( (0,0,radius) )
        fid = tmp.numpts-1

        faces = [] 
        front_loop=[]
        back_loop=[]
        last_pt = None 

        # compute 10 vertices at 1st and 2nd rows
        for i in range(1,7):
            n = tmp.numpts # add to this index each time

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
            
            tmp.points.append(tuple(vtmp1))
            tmp.points.append(tuple(vtmp2))
            
            if n>1:
                if (n%2)-1==0:
                    front_loop.append(n)
                    back_loop.append(n-1)

            if build_wire:
                tmp.prim_locator(pos=vtmp1,size=.1)
                tmp.prim_locator(pos=vtmp2,size=.1)
                tmp.insert(vec3(vtmp2)-vec3(vtmp1), pos=vtmp1)  

            if build_geom:
                if i>1:
                    tmp.polygons.append( (n, n+1, n+2) )
                if i<6:
                    tmp.polygons.append( (n+1, n+2, n+3) )
        
        # the last bottom vertex at (0, 0, -r)
        tmp.points.append( (0,0,-radius) )
        last_pt = tmp.numpts 
 
        if build_geom:
            ## build the front endcaps
            for i,idx in enumerate(front_loop):
                #cool star pattern  
                #self.polygons.append( (idx, last_pt, front_loop[i-2]) )            
                tmp.polygons.append( (idx, last_pt, front_loop[i-1]) )    
            ## build the back endcaps
            for i,idx in enumerate(back_loop):
                tmp.polygons.append( (idx, 1, back_loop[i-1]) )  

        tmp.points = tmp.xform_pts(pos=pos, pts=tmp.points)

        return tmp

    ##------------------------------------------------## 
    def _subdiv(self, size):
        """
            find middle point of 2 vertices
            NOTE: new vertex must be resized, so the length is equal to the radius
            from http://www.songho.ca/opengl/gl_sphere.html 
        """

        def computeHalfVertex( radius, v1, v2):
            newV = []
            
            print("####### computeHalfVertex ", v1 , v2 )

            newV.append(v1[0] + v2[0])    #x
            newV.append(v1[1] + v2[1])    #y
            newV.append(v1[2] + v2[2])    #z
            scale = radius / math.sqrt(newV[0]*newV[0] + newV[1]*newV[1] + newV[2]*newV[2]);
            newV[0] *= scale
            newV[1] *= scale
            newV[2] *= scale
            
            print("####### computeHalfVertex ", newV)

            return newV

        # The subdivision algorithm is splitting the 3 edge lines of each triangle in half, 
        # then extruding the new middle point outward, 
        # so its length (the distance from the center) is the same as sphere's radius. 
        
        subdivision = 1

        tmpVertices = []
        tmpIndices = []

        v1 = .0
        v2 = .0
        v3 = .0
 
        # new vertex positions
        newV1 = vec3()
        newV2 = vec3()
        newV3 = vec3() 

        index = 0

        # iterate all subdivision levels
        for i in range(subdivision):
            # copy prev vertex/index arrays and clear
            tmpVertices = self.points
            tmpIndices = self.polygons

            self.points = []
            self.polygons = []
            
            index = 0;

            # perform subdivision for each triangle
            for j in range(0, len(tmpIndices), 3):

                ##--------------------------------##
                ## DEBUG

                # get 3 vertices of a triangle
                #v1 = tmpVertices[tmpIndices[j] * 3]
                #v2 = tmpVertices[tmpIndices[j + 1] * 3]
                #v3 = tmpVertices[tmpIndices[j + 2] * 3]

                v1 = tmpVertices[ tmpIndices[j][0] ]
                v2 = tmpVertices[ tmpIndices[j+1][1] ]
                v3 = tmpVertices[ tmpIndices[j+2][2] ]

                print('###############', j)
                print(v1)
                print(v2)
                print(v3)
                
                ## DEBUG
                ##--------------------------------##

                # compute 3 new vertices by spliting half on each edge
                #         v1       
                #        / \       
                # newV1 *---* newV3
                #      / \ / \     
                #    v2---*---v3   
                #       newV2      

                #newV1 = computeHalfVertex(size, v1, v2)
                #newV2 = computeHalfVertex(size, v2, v3)
                #newV3 = computeHalfVertex(size, v1, v3)

                # add 4 new triangles to vertex array
                self.points.append((v1,    newV1, newV3) )
                self.points.append((newV1, v2,    newV2) )
                self.points.append((newV1, newV2, newV3) )
                self.points.append((newV3, newV2, v3)    )

                # add indices of 4 new triangles
                self.polygons.append((index,   index+1, index+2) )
                self.polygons.append((index+3, index+4, index+5) )
                self.polygons.append((index+6, index+7, index+8) )
                self.polygons.append((index+9, index+10,index+11))
                index = index + 12  # next index
            
                self.points = tmpVertices
                self.polygons = tmpIndices

    ##------------------------------------------------## 

    def prim_sphere(self, pos=(0,0,0), rot=(0,0,0), size=1 ):
        """
            build a icosahedron : UNFINISHED - need to add smoothing/triangulate
            from http://www.songho.ca/opengl/gl_sphere.html 
        """
        obj = self._icosahedron(pos=pos, radius=size, build_geom=True)
        
        #DEBUG NOT WORKING YET # 
        #obj._subdiv(size)
        
        self.insert(obj)



     
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
    def prim_ray(self, pos=(0,0,0), vec3=(0,1,0) ): 
        """ derived from prim_arrow (unfinished)
            fully 3D model of an arrow on a vector 
             
        """

        spokes  = 4  # num turns around axis 

        length  = 1

        dia      = .1        # extrude length is double this, or .2 
        shaftlen = length-.2 # cone is .2, that plus this = 1 
        
        tmpobj1 = object3d()
        tmpobj1.prim_cone( axis='y', pos=(0,shaftlen,0), dia=dia, spokes=spokes )        
        
        tmpobj = object3d()
        tmpobj.prim_circle( axis='y', pos=(0,0,0), spokes=spokes , dia=dia/5)
        tmpobj.extrude_face(0, distance=shaftlen)

        tmpobj1.insert(tmpobj)


        self.insert(tmpobj1) 

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




