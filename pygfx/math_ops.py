#!/usr/local/bin/python3


import itertools 
import math
import os

import numpy as np             #for testing - remove later 
from scipy.linalg import expm  #for testing - remove later
from numpy.linalg import inv


DEG_TO_RAD = 0.0174532925 # degree = radian * (180 / PI) # PI = 3.14159265
RAD_TO_DEG = 57.29577951  # radian = degree * (PI/180)


###############################################
class math_util(object):    
    """ general math library - may contain some of the same functions 
        as the vecotr and matrix objects, but those will be implemented 
        to operate relative to themselves, these will work more generally

        for example self.dot_vec3() will be vec3.dot , ect 
    """

    def dot_vec3 (self, v1, v2):
         """ scalar - mag but not direction """
         return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2] 

    def cross_vec3 (self, v1, v2):
        return (v1[1]*v2[2]-v1[2]*v2[1], v1[2]*v2[0]-v1[0]*v2[2], v1[0]*v2[1]-v1[1]*v2[0])

    def dtr (self, deg ):
       return deg * DEG_TO_RAD

    def rtd (self, rad ):
       return rad * RAD_TO_DEG

    def dtr_vec3(self, invec):
        return ( self.dtr( invec[0] ),
                 self.dtr( invec[1] ),
                 self.dtr( invec[2] ),
               )

    def rtd_vec3(self, invec):
        return ( self.rtd( invec[0] ),
                 self.rtd( invec[1] ),
                 self.rtd( invec[2] ),
               )

    def vect_add(self, v1, v2):
        return [v1[0]+v2[0], v1[0]+v2[1], v1[2]+v2[2]]

    def mult_scalar (self, scalar, v):
        return [v[0]*scalar, v[1]*scalar, v[2]*scalar ]

    def mult_m33_vec3(self, m, v):
        """ multiplies a 3X3 matrix by a 3D vector - returns a vector tuple 
            NUMPY DOT is the same as multiplying 
        """
        
        outx = m[0] * v[0] + m[3] * v[1] + m[6] * v[2] 
        outy = m[1] * v[0] + m[4] * v[1] + m[7] * v[2] 
        outz = m[2] * v[0] + m[5] * v[1] + m[8] * v[2] 
          
        return  (outx, outy, outz)

    def mult_m44_vec3(self, m, v):
        """ multiplies a 4X4 matrix by a 3D vector - returns a vector tuple """
        
        outx = m[0]*v[0] + m[4]*v[1] + m[8]  * v[2]+m[12]
        outy = m[1]*v[0] + m[5]*v[1] + m[9]  * v[2]+m[13]
        outz = m[2]*v[0] + m[6]*v[1] + m[10] * v[2]+m[14]
          
        return  (outx, outy, outz)

    def calc_line_length(self, x1, y1, x2, y2):
        """ distance between two points 
            DEBUG -- make work in 2d and 3d!
            DEBUG -- merge with -  self.calc_line_length(pt1[0], pt1[1], pt2[0], pt2[1] )

        """
        return math.sqrt( ((x1-x2)**2)+ ((y1-y2)**2) )

    def mult_mat33(self, m, n):
        """ multiply two 3X3 matricies together """
        return [
                m[0]*n[0] + m[1]*n[3] + m[2]*n[6],
                m[0]*n[1] + m[1]*n[4] + m[2]*n[7],
                m[0]*n[2] + m[1]*n[5] + m[2]*n[8],
                m[3]*n[0] + m[4]*n[3] + m[5]*n[6],
                m[3]*n[1] + m[4]*n[4] + m[5]*n[7],
                m[3]*n[2] + m[4]*n[5] + m[5]*n[8],
                m[6]*n[0] + m[7]*n[3] + m[8]*n[6],
                m[6]*n[1] + m[7]*n[4] + m[8]*n[7],
                m[6]*n[2] + m[7]*n[5] + m[8]*n[8]   
               ]        

    def mult_mat44(self, m, n):
        """multiply two 4X4 matricies together """
        return [
                m[0]*n[0]  + m[1]*n[4]  + m[2]*n[8]   + m[3]*n[12],
                m[0]*n[1]  + m[1]*n[5]  + m[2]*n[9]   + m[3]*n[13],
                m[0]*n[2]  + m[1]*n[6]  + m[2]*n[10]  + m[3]*n[14],
                m[0]*n[3]  + m[1]*n[7]  + m[2]*n[11]  + m[3]*n[15],
                m[4]*n[0]  + m[5]*n[4]  + m[6]*n[8]   + m[7]*n[12],
                m[4]*n[1]  + m[5]*n[5]  + m[6]*n[9]   + m[7]*n[13],
                m[4]*n[2]  + m[5]*n[6]  + m[6]*n[10]  + m[7]*n[14],
                m[4]*n[3]  + m[5]*n[7]  + m[6]*n[11]  + m[7]*n[15],
                m[8]*n[0]  + m[9]*n[4]  + m[10]*n[8]  + m[11]*n[12],
                m[8]*n[1]  + m[9]*n[5]  + m[10]*n[9]  + m[11]*n[13],
                m[8]*n[2]  + m[9]*n[6]  + m[10]*n[10] + m[11]*n[14],
                m[8]*n[3]  + m[9]*n[7]  + m[10]*n[11] + m[11]*n[15],
                m[12]*n[0] + m[13]*n[4] + m[14]*n[8]  + m[15]*n[12],
                m[12]*n[1] + m[13]*n[5] + m[14]*n[9]  + m[15]*n[13],
                m[12]*n[2] + m[13]*n[6] + m[14]*n[10] + m[15]*n[14],
                m[12]*n[3] + m[13]*n[7] + m[14]*n[11] + m[15]*n[15]
               ]


    def normalize_vec3(self, in_vec):
        try:
           invLength = 1.0/math.sqrt(in_vec[0]*in_vec[0] + in_vec[1]*in_vec[1]+ in_vec[2]*in_vec[2])
           return (in_vec[0] *invLength, in_vec[1] * invLength,  in_vec[2] * invLength)
        except:
           print('normalize_vec3: divide by zero error.') 
           return [0,0,1]




###############################################
class vec2(object):    

    def __init__(self,x=0,y=0):
        self.x = x;self.y = y    

    def __repr__(self):
        return '(%s, %s)' % (self.x, self.y)

    def __abs__(self):
        return type(self)(abs(self.x), abs(self.y))

    def __add__(self, other):
        return type(self)(self.x + other.x, self.y + other.y)

    def __sub__(self, other):
        return type(self)(self.x - other.x, self.y - other.y)

    def __mul__(self, other):
        return type(self)(self.x * other.x, self.y * other.y)

    ## ## ##
    def div_scalar(self, scalar):
        self.x = self.x/scalar
        self.y = self.y/scalar
        return type(self)(self.x , self.y )

    def __truediv__(self, other):

        #check type and switch 

        #combine with self.div_scalar() 

        #untested - normalized? 
        #return type(self)(self.x / other.x, self.y / other.y)
        self.x = other.x/other.length()
        self.y = other.y/other.length()
        return type(self)(self.x, self.y)

    def __getitem__(self,index):
        if index==0:
            return self.x
        if index==1:
            return self.y

    @property
    def length(self):
        return math.sqrt(self.x*self.x + self.y*self.y)

    def dtr(self):
        """ degree to radian """ 
        return ( self.mu.dtr( self.x ),
                 self.mu.dtr( self.y ),
               )
    def rtd(self, invec):
        """ radian to degree """
        return ( self.mu.rtd( self.x ),
                 self.mu.rtd( self.y ),
               ) 

    def dot(self, other):
        """ dot product """
        return self.x * other.x + self.y * other.y

    def distance_to(self, other, doRound=False):
        val = math.hypot((self.x - other.x), (self.y - other.y))
        if not doRound:        
            return val
        if doRound: 
           return int(val)

    def project_pt(self, A, B, offset , doRound=False):
        nX = B.x - A.x;nY = B.y - A.y
        distX = pow( (A.x - B.x ) , 2.0 ) 
        distY = pow( (A.y - B.y ) , 2.0 ) 
        vecLength = math.sqrt(distX + distY )
        # normalized vector  
        calcX = nX / vecLength
        calcY = nY / vecLength
        # project point along vector with offset (can use negative too)
        ptX = B.x + (calcX * offset)
        ptY = B.y + (calcY * offset)
        if not doRound:
            return type(self)(ptX, ptY)
        if doRound:
            return type(self)(int(ptX), int(ptY) )

    def intersect(self, v1s, v1e, v2s, v2e ):
        """ intersect 2 lines in 2D 
            v1s - line1 start 
            v1e - line1 end 
            v2s - line2 start
            v2e - line2 end 
        """

        #start and end coords for two 2D lines
        p0_x = float(v1s[0])
        p0_y = float(v1s[1])
        p1_x = float(v1e[0])
        p1_y = float(v1e[1])
        p2_x = float(v2s[0])
        p2_y = float(v2s[1])
        p3_x = float(v2e[0])
        p3_y = float(v2e[1])

        #return values
        i_x = 0;i_y = 0 

        s1_x = 0;s1_y = 0
        s2_x = 0;s2_y = 0

        s1_x = p1_x - p0_x  
        s1_y = p1_y - p0_y
        s2_x = p3_x - p2_x
        s2_y = p3_y - p2_y

        s = 0;t = 0

        try: 
            s = (-s1_y * (p0_x - p2_x) + s1_x * (p0_y - p2_y)) / (-s2_x * s1_y + s1_x * s2_y)
            t = ( s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x)) / (-s2_x * s1_y + s1_x * s2_y)
        except:
            #print('something went wrong in 2D line intersection')
            return 0
        
        if (s >= 0 and s <= 1 and t >= 0 and t <= 1):
            # Collision detected
            i_x = p0_x + (t * s1_x);
            i_y = p0_y + (t * s1_y);
            return (i_x,i_y)  

        return 0; # No collision

    @property
    def normal(self):
        """ untested - normal for a 2d vector """ 
        invLength = 1.0/math.sqrt(self.x*self.x + self.y*self.y)
        self.x *= invLength
        self.y *= invLength
            
    def calc_normal(self):
        """
        Quoting from http://www.opengl.org/wiki/Calculating_a_Surface_Normal
        A surface normal for a triangle can be calculated by taking the vector cross product of two edges of that triangle. 
        The order of the vertices used in the calculation will affect the direction of the normal 
        (in or out of the face w.r.t. winding).
        So for a triangle p1, p2, p3, 
        if the vector U = p2 - p1 and the vector V = p3 - p1 then the normal N = U x V and can be calculated by:
        Nx = UyVz - UzVy
        Ny = UzVx - UxVz
        Nz = UxVy - UyVx
        """
        """
        The cross product of two sides of the triangle equals the surface normal. So, if VV = P2P2 - P1P1 and WW = P3P3 - P1P1, and NN is the surface normal, then:
        Nx=(Vy∗Wz)−(Vz∗Wy)Nx=(Vy∗Wz)−(Vz∗Wy)
        Ny=(Vz∗Wx)−(Vx∗Wz)Ny=(Vz∗Wx)−(Vx∗Wz)
        Nz=(Vx∗Wy)−(Vy∗Wx)Nz=(Vx∗Wy)−(Vy∗Wx)
        If AA is the new vector whose components add up to 1, then:
        Ax=Nx/(|Nx|+|Ny|+|Nz|)Ax=Nx/(|Nx|+|Ny|+|Nz|)
        Ay=Ny/(|Nx|+|Ny|+|Nz|)Ay=Ny/(|Nx|+|Ny|+|Nz|)
        Az=Nz/(|Nx|+|Ny|+|Nz|)
        """

        out = []

        #Nx = UyVz - UzVy
        #Ny = UzVx - UxVz
        #Nz = UxVy - UyVx

        return out

###############################################
class vec3(object):    

    def __init__(self,x=0,y=0,z=0):
        self.x=x;self.y=y;self.z=z  
        self.mu = math_util() 

    def __repr__(self):
        return '(%s, %s, %s)' % (self.x, self.y, self.z)

    def __abs__(self):
        return type(self)(abs(self.x), abs(self.y), abs(self.z))

    def __add__(self, other):
        if isinstance(other, np.ndarray):
            return type(self)(self.x+other[0], self.y+other[1], self.z+other[2])
        if isinstance(other, float) or isinstance(other, int):
            return type(self)(self.x+other, self.y+other, self.z+other)
        if isinstance(other, vec3):                    
            return type(self)(self.x+other.x, self.y+other.y, self.z+other.z)
        if isinstance(other, tuple) or isinstance(other, list):
            return type(self)(self.x+other[0], self.y+other[1], self.z+other[2])  

    def __sub__(self, other):
        if isinstance(other, np.ndarray):
            return type(self)(self.x-other[0], self.y-other[1], self.z-other[2])
        if isinstance(other, float) or isinstance(other, int):
            return type(self)(self.x-other, self.y-other, self.z-other)
        if isinstance(other, vec3):                    
            return type(self)(self.x-other.x, self.y-other.y, self.z-other.z)
        if isinstance(other, tuple) or isinstance(other, list):
            return type(self)(self.x-other[0], self.y-other[1], self.z-other[2])  

    def __mul__(self, other):
        if isinstance(other, np.ndarray):
            return type(self)(self.x*other[0], self.y*other[1], self.z*other[2])
        if isinstance(other, float) or isinstance(other, int):
            return type(self)(self.x*other, self.y*other, self.z*other)
        if isinstance(other, vec3):
            return type(self)(self.x*other.x, self.y*other.y, self.z*other.z)
        if isinstance(other, tuple) or isinstance(other, list):
            return type(self)(self.x*other[0], self.y*other[1], self.z*other[2])          

    def __truediv__(self, other):
        #untested - normalized? 
        #return type(self)(self.x / other.x, self.y / other.y)
        self.x = other.x/other.length()
        self.y = other.y/other.length()
        self.z = other.y/other.length()
        return type(self)(self.x, self.y, self.z)

    def __float__(self):
        return type(self)(float(self.x), float(self.y), float(self.z) )

    def __getitem__(self, index):
        if index==0:
            return self.x
        if index==1:
            return self.y
        if index==2:
            return self.z

    def __setitem__(self, key, item):
        if key==0:
            self.x = item
        if key==1:
            self.y = item
        if key==2:
            self.z = item

    @property
    def as_np(self):
        """  vec3 as np.array """ 
        return self.copy(vtype='numpy')

    def insert(self, iterable):
        """ convert an np.array, tuple or list  to vec3  
            does not check size, so just assume 3 items (x,y,z)
        """

        if isinstance(iterable, np.ndarray):
            self.x = iterable[0]
            self.y = iterable[1]            
            self.z = iterable[2]

        if isinstance(iterable, list) or isinstance(iterable, tuple):
            self.x = iterable[0]
            self.y = iterable[1]            
            self.z = iterable[2]
        return self 

 
    def copy(self, vtype=None):
        if vtype is None:
            return type(self)(self.x,self.y,self.z)
        if vtype is 'numpy':
            return np.array( (self.x,self.y,self.z) )
        if vtype is 'tuple':
            return ( (self.x,self.y,self.z) )

    @property
    def length(self):
        """ output - scalar representing the length of this vector """
        return math.sqrt(self.x*self.x + self.y*self.y + self.z*self.z )

    def dot (self, v2):
         """ input  - vector
             output - scalar  
         """
         return self.x*v2[0] + self.y*v2[1] + self.z*v2[2]

    def cross (self, other):
        """ input  - vector
            output - vector 
        """
        x = self.y*other.z - self.z*other.y
        y = self.z*other.x - self.x*other.z
        z = self.x*other.y - self.y*other.x
        return type(self)(x,y,z)

    @property
    def dtr(self):
        """ degree to radian """
        return ( self.mu.drt( self.x ),
                 self.mu.dtr( self.y ),
                 self.mu.dtr( self.z )                 
               )

    @property
    def rtd(self):
        """ radian to degree """
        return ( self.mu.rtd( self.x ),
                 self.mu.rtd( self.y ),
                 self.mu.rtd( self.z )                 
               )  
    @property
    def normal(self) :
        """ normal of this vector """
        try:
           invLength = 1.0/math.sqrt(self.x*self.x + self.y*self.y + self.z*self.z)
           return type(self)(self.x *invLength, self.y * invLength,  self.z * invLength)
        except:
           print('normal vec3: divide by zero error.') 
           return [0,0,1]

    @property
    def np_normal(self):
        """ unit vector of the vector using numpy """
        return self.as_np / np.linalg.norm(self.as_np)

    def look_at(self):
        """
        https://stackoverflow.com/questions/1251828/calculate-rotations-to-look-at-a-3d-point
        
        ###################################

        rotx = Math.atan2( y, z )
        roty = Math.atan2( x * Math.cos(rotx), z )
        rotz = Math.atan2( Math.cos(rotx), Math.sin(rotx) * Math.sin(roty) )
        ###################################
         rotx = Math.atan2( y, z );
         if (z >= 0) {
            roty = -Math.atan2( x * Math.cos(rotx), z );
         }else{
            roty = Math.atan2( x * Math.cos(rotx), -z );
         }        
        ###################################

        About X: -atan2(y, z)
        About Y: atan2(x, sqrt(y*y + z*z))
        About Z: 0 

        """
        pass

    def np_angle_between(self, v1, v2):
        """ 
           UNTESTED 

           Returns the angle in radians between vectors 'v1' and 'v2'::

                >>> angle_between((1, 0, 0), (0, 1, 0))
                1.5707963267948966
                >>> angle_between((1, 0, 0), (1, 0, 0))
                0.0
                >>> angle_between((1, 0, 0), (-1, 0, 0))
                3.141592653589793
        """
        
        if isinstance(v1, tuple):
            v1 = self.insert(v1)  #wrong?
        if isinstance(v2, tuple):
            v2 = self.insert(v2)  #wrong?

        v1_u = v1.np_normal;v2_u = v2.np_normal

        #print('### ', v1_u , v2_u )

        return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

    def angle_between(self, other):
        """ 
            result is in radians
            derived from law of cosines 
            range from 0 (colinear) to PI radians (180 deg)  

            FAILS WHEN 0 or 180 ?
        """
        o = math.acos( self.dot(other) / (self.length*other.length) ) 
        return o 


    def vector_mean(self, p_vecs, mode='vec3'):
        """ 
            UNTESTED 
         
            take a list of vectors and return a new vector that is an average

        """
        

        tmp_x = 0;tmp_y = 0;tmp_z = 0;

        count = 0
        a = 0 
        for v in p_vecs: 
            if mode=='angle':            
                a = a + self.np_angle_between( upvec, v )
            if mode=='vec3':

                if isinstance(v, vec3):
                    tmp_x = tmp_x+v.x
                    tmp_y = tmp_x+v.y
                    tmp_z = tmp_x+v.z 
                if isinstance(v, tuple) or isinstance(v, list):                
                    tmp_x = tmp_x+v[0]
                    tmp_y = tmp_x+v[1]
                    tmp_z = tmp_x+v[2] 

            count += 1 
        
        if mode=='angle':
            return a/count 
        if mode=='vec3':
            return type(self)(tmp_x/count, tmp_y/count, tmp_z/count)


        return None 


                


###############################################
class matrix33(object):
    """ 3X3 matrix from pure python 
        limited support to interface to numpy arrays

        the patern "return type(self) is nice to return copies of itself,
        but beware that this structure is not compatible for passing mutable types.
        Only primitive types work, in this case, 9 floats  
    """
    def __init__(self, a=1,b=0,c=0,
                       d=0,e=1,f=0,
                       g=0,h=0,i=1):
        self.m = [a,b,c,d,e,f,g,h,i]
        self.mu = math_util()

    def __getitem__(self, index):
        return self.m[index]

    def __setitem__(self, key, item):
        self.m[key] = item

    def __repr__(self):
        return '(%s, %s, %s, %s, %s, %s, %s, %s, %s )'%(
                self.m[0], self.m[1], self.m[2],  self.m[3],  self.m[4],  
                self.m[5],  self.m[6],  self.m[7], self.m[8] )

    def test_index(self):
        """ fill matrix with incrementing numbers to see indices 
            for testing transpose and other things 
        """
        tmp = type(self)()
        for i in range(9):
            tmp.m[i] = i
        return tmp

    @property
    def identity(self):
        """ using the "standard" 9 element, Row major, 3X3 rotation matrix """
        return type(self)()

    @property
    def determinant(self):
        """
            https://www.mathsisfun.com/algebra/matrix-determinant.html

                  a b c  |  0 1 2  |        
             A =  d e f  |  3 4 5  |       
                  g h i  |  6 7 8  |     
           
            |A| = a(ei − fh) − b(di − fg) + c(dh − eg)
            |A| = 0(48 − 57) − 1(38 − 56) + 2(37 − 46)

        """  
        a = self.copy() 
        o = a[0]* ((a[4]*a[8])-(a[5]*a[7])) - a[1]*((a[3]*a[8])-(a[5]*a[6])) +  a[2]*((a[3]*a[7])-(a[4]*a[6]))
        return o
      

    @property
    def np_inverse(self):
        """ seems to work, but I am not convinced """
        a = self.copy(mtype='numpy')
        b = inv(a)
        c = matrix33()
        c.insert(b)
        #print('## inverse \n\n', self  , ' \n\n', c , ' \n\n',  self*c )
        return b


    def serialize(self, inarray):
        """ serialize this array into a list 
            if you want it as a numpy array use self.copy(mtype='numpy')
        """

        out = []
        for i in self.m:
            out.append(i)
        return out


    def insert(self, iterable):
        """ load the first 9 things we find into this matrix 
            accepts numpy.ndarray, list, and tuple 
        """
       
        #numpy ND array 
        if isinstance(iterable, np.ndarray):
            out = [];idx=0
            for row in iterable:
                for col in row:
                    self.m[idx] = col
                    idx+=1
        
        #serialized simple array
        if isinstance(iterable, list) or isinstance(iterable, tuple):
            for idx,i in enumerate(iterable):
                if idx <= 9:
                     self.m[idx] = iterable[idx] 

    def copy(self, mtype=None):
        """ create a copy of this matrix - can be same type or a numpy.ndarray """
        if not mtype:
            return type(self)(
                self.m[0] , self.m[1] , self.m[2] ,
                self.m[3] , self.m[4] , self.m[5] ,
                self.m[6] , self.m[7] , self.m[8]
            )
        
        if mtype == 'numpy':
            return np.array(([self.m[0] , self.m[1] , self.m[2]] ,
                             [self.m[3] , self.m[4] , self.m[5]] ,
                             [self.m[6] , self.m[7] , self.m[8]])
                           )

    @property
    def transpose(self):
        """
        # standard indicies  |  #transposed indicies
        1  2  3              |   1 4 7
        4  5  6              |   2 5 8
        7  8  9              |   3 6 9 
       
        """

        return type(self)(
            self.m[0], self.m[3], self.m[6],
            self.m[1], self.m[4], self.m[7],
            self.m[2], self.m[5], self.m[8] 
        )

    @property
    def inverse(self):
        """
            NOT DONE 
            multiply a matrix by its inverse and we end up with the Identity Matrix.

            https://www.wikihow.com/Find-the-Inverse-of-a-3x3-Matrix
            
            https://www.mathsisfun.com/algebra/matrix-inverse.html
        a = self.copy()
        o = self.identity 

        det = a.determinant
        
        if det == 0:
            # determinant is 0, matrix has no inverse
            return None 
        else:     
            print( '### ### ' , a )

            b = a.transpose

            print ("TRANSPOSE IS  ", a.transpose )
       

            return type(self)(
                o.m[0], o.m[1], o.m[2],
                o.m[3], o.m[4], o.m[5],
                o.m[6], o.m[7], o.m[8] 
            )
        """
        pass


    def __add__(self, other):
        return type(self)(
            self.m[0]+n[0], self.m[1]+n[1], self.m[2]+n[2],
            self.m[3]+n[3], self.m[4]+n[4], self.m[5]+n[5],
            self.m[6]+n[6], self.m[7]+n[7], self.m[8]+n[8]     
        )
    
    def __sub__(self, other):
        return type(self)(
            self.m[0]-n[0], self.m[1]-n[1], self.m[2]-n[2],
            self.m[3]-n[3], self.m[4]-n[4], self.m[5]-n[5],
            self.m[6]-n[6], self.m[7]-n[7], self.m[8]-n[8]     
        )

    # def __truediv__(self, other):
    #     # matrices don't divide! there is no concept of dividing by a matrix.
    #     # multiply by an inverse, which achieves the same thing.

    def __mul__(self, n):
        """ multiply two 3X3 matricies together 

            the order that you mutliply will call a different object!!!
            make sure you do "this * other", NOT "other * this"

        """

        if isinstance(n, vec3):
            outx = self.m[0] * n.x + self.m[3] * n.y + self.m[6] * n.z 
            outy = self.m[1] * n.x + self.m[4] * n.y + self.m[7] * n.z 
            outz = self.m[2] * n.x + self.m[5] * n.y + self.m[8] * n.z 
            return  (outx, outy, outz)

        if isinstance(n, tuple) or isinstance(n, list) or isinstance(n, np.ndarray):
            outx = self.m[0] * n[0] + self.m[3] * n[1] + self.m[6] * n[2] 
            outy = self.m[1] * n[0] + self.m[4] * n[1] + self.m[7] * n[2] 
            outz = self.m[2] * n[0] + self.m[5] * n[1] + self.m[8] * n[2] 
            return  (outx, outy, outz)

        if type(n) == type(self):
            return type(self)(
                    self.m[0]*n[0]  + self.m[1]*n[3]  + self.m[2]*n[6],
                    self.m[0]*n[1]  + self.m[1]*n[4]  + self.m[2]*n[7],
                    self.m[0]*n[2]  + self.m[1]*n[5]  + self.m[2]*n[8],
                    self.m[3]*n[0]  + self.m[4]*n[3]  + self.m[5]*n[6],
                    self.m[3]*n[1]  + self.m[4]*n[4]  + self.m[5]*n[7],
                    self.m[3]*n[2]  + self.m[4]*n[5]  + self.m[5]*n[8],
                    self.m[6]*n[0]  + self.m[7]*n[3]  + self.m[8]*n[6],
                    self.m[6]*n[1]  + self.m[7]*n[4]  + self.m[8]*n[7],
                    self.m[6]*n[2]  + self.m[7]*n[5]  + self.m[8]*n[8]   
                   )


    def batch_mult_pts(self, pts):
        """ iterate a list of points and multiply them by this matrix """

        tmp_buffer = []
        out = None
        for pvec in pts:  
            tmp_buffer.append( self * pvec )
        return tmp_buffer


    def rotate_pts_3d(self, points, xrot, yrot, zrot):
        """
           previously named rotate_mat3 
           using the "standard" 9 element, Row major, 3X3 rotation matrix used by Maya
                
           [0  1  2]      xx xy xz 
           [3  4  5]      yx yy yz 
           [6  7  8]      zx zy zz 
           ------------------------------
           rotate Y matrix     
           |cos(y)  0      -sin(y) | 
           |0       1       0      | 
           |sin(y)  0       cos(y) | 
           ------------------------------
           rotate Z  matrix 
           |  cos(z)  sin(z)  0    | 
           | -sin(z)  cos(z)  0    |
           |  0       0       1    |
           ------------------------------
           rotate X matrix  
           | 1       0         0   |    
           | 0    cos(x)   sin(x)  |  
           | 0    -sin(x)   cos(x) |  
           ------------------------------            
 
        """
        dtr = self.mu.dtr

        # build rotationY (see diagram above) 
        y_matrix =  self.identity
        y_matrix[0]  =  math.cos(dtr( yrot ))
        y_matrix[2]  = -math.sin(dtr( yrot ))
        y_matrix[6]  =  math.sin(dtr( yrot ))
        y_matrix[8]  =  math.cos(dtr( yrot ))

        ####                
        # build rotationZ (see diagram above) 
        z_matrix    =  self.identity
        z_matrix[0] =  math.cos(dtr( zrot ))
        z_matrix[1] =  math.sin(dtr( zrot ))
        z_matrix[3] = -math.sin(dtr( zrot ))
        z_matrix[4] =  math.cos(dtr( zrot ))
        tmp_matr =  y_matrix * z_matrix 

        ####
        # build rotationX (see diagram above) 
        x_matrix =  self.identity
        x_matrix[4]  =   math.cos(dtr( xrot )) 
        x_matrix[5]  =   math.sin(dtr( xrot )) 
        x_matrix[7]  =  -math.sin(dtr( xrot ))
        x_matrix[8]  =   math.cos(dtr( xrot ))
        rotation_33 = x_matrix * tmp_matr 

        ############ 
        return rotation_33.batch_mult_pts(points)

###############################################
class matrix44(object):
    """ 4X4 matrix from pure python 
        limited support to interface to numpy arrays 

        the patern "return type(self) is nice to retrurn copies of itself,
        but beware that this structure is not compatible for passing mutable types.
        Only primitive types work, in this case floats  
    """    
    def __init__(self, a=1,b=0,c=0,d=0,
                       e=0,f=1,g=0,h=0,
                       i=0,j=0,k=1,l=0,
                       m=0,n=0,o=0,p=1):

        self.mu = math_util()
        self.m = [a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p]

    def __getitem__(self, index):
        return self.m[index]

    def __setitem__(self, key, item):
        self.m[key] = item

    def __repr__(self):
        return '(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s )'%(
                self.m[0], self.m[1], self.m[2],  self.m[3],  self.m[4],  self.m[5],  self.m[6],  self.m[7], 
                self.m[8], self.m[9], self.m[10], self.m[11], self.m[12], self.m[13], self.m[14], self.m[15])

    def __add__(self, n):
        return type(self)(
                self.m[0] +n[0] , self.m[1]+n[1]  , self.m[2]+n[2]  , self.m[3]+n[3]  , 
                self.m[4] +n[4] , self.m[5]+n[5]  , self.m[6]+n[6]  , self.m[7]+n[7]  ,
                self.m[8] +n[8] , self.m[9]+n[9]  , self.m[10]+n[10], self.m[11]+n[11],
                self.m[12]+n[12], self.m[13]+n[13], self.m[14]+n[14], self.m[15]+n[15]
               )
   
    
    def __sub__(self, n):
        return type(self)(
                self.m[0] -n[0] , self.m[1]-n[1]  , self.m[2]-n[2]  , self.m[3]-n[3]  , 
                self.m[4] -n[4] , self.m[5]-n[5]  , self.m[6]-n[6]  , self.m[7]-n[7]  ,
                self.m[8] -n[8] , self.m[9]-n[9]  , self.m[10]-n[10], self.m[11]-n[11],
                self.m[12]-n[12], self.m[13]-n[13], self.m[14]-n[14], self.m[15]-n[15]
               )
   
    def __mul__(self, n):
        """multiply two 4X4 matricies together """

        if isinstance(n, vec3) or isinstance(n, np.ndarray):
            outx = self.m[0] * n.x + self.m[4] * n.y + self.m[8]  * n.z + self.m[12]
            outy = self.m[1] * n.x + self.m[5] * n.y + self.m[9]  * n.z + self.m[13]
            outz = self.m[2] * n.x + self.m[6] * n.y + self.m[10] * n.z + self.m[14]
            return  (outx, outy, outz)

        if isinstance(n, tuple) or isinstance(n, list):
            outx = self.m[0] * n[0] + self.m[4] * n[1] + self.m[8]  * n[2] + self.m[12]
            outy = self.m[1] * n[0] + self.m[5] * n[1] + self.m[9]  * n[2] + self.m[13]
            outz = self.m[2] * n[0] + self.m[6] * n[1] + self.m[10] * n[2] + self.m[14]
            return  (outx, outy, outz)

        if type(n) == type(self):
            return type(self)(
                    self.m[0]*n[0]  + self.m[1]*n[4]  + self.m[2]*n[8]   + self.m[3]*n[12] ,
                    self.m[0]*n[1]  + self.m[1]*n[5]  + self.m[2]*n[9]   + self.m[3]*n[13] ,
                    self.m[0]*n[2]  + self.m[1]*n[6]  + self.m[2]*n[10]  + self.m[3]*n[14] ,
                    self.m[0]*n[3]  + self.m[1]*n[7]  + self.m[2]*n[11]  + self.m[3]*n[15] ,
                    self.m[4]*n[0]  + self.m[5]*n[4]  + self.m[6]*n[8]   + self.m[7]*n[12] ,
                    self.m[4]*n[1]  + self.m[5]*n[5]  + self.m[6]*n[9]   + self.m[7]*n[13] ,
                    self.m[4]*n[2]  + self.m[5]*n[6]  + self.m[6]*n[10]  + self.m[7]*n[14] ,
                    self.m[4]*n[3]  + self.m[5]*n[7]  + self.m[6]*n[11]  + self.m[7]*n[15] ,
                    self.m[8]*n[0]  + self.m[9]*n[4]  + self.m[10]*n[8]  + self.m[11]*n[12],
                    self.m[8]*n[1]  + self.m[9]*n[5]  + self.m[10]*n[9]  + self.m[11]*n[13],
                    self.m[8]*n[2]  + self.m[9]*n[6]  + self.m[10]*n[10] + self.m[11]*n[14],
                    self.m[8]*n[3]  + self.m[9]*n[7]  + self.m[10]*n[11] + self.m[11]*n[15],
                    self.m[12]*n[0] + self.m[13]*n[4] + self.m[14]*n[8]  + self.m[15]*n[12],
                    self.m[12]*n[1] + self.m[13]*n[5] + self.m[14]*n[9]  + self.m[15]*n[13],
                    self.m[12]*n[2] + self.m[13]*n[6] + self.m[14]*n[10] + self.m[15]*n[14],
                    self.m[12]*n[3] + self.m[13]*n[7] + self.m[14]*n[11] + self.m[15]*n[15]
                   )

    @property
    def identity(self):
        """ using the "standard" 16 element, Row major, 4X4 rotation matrix """
        #return  [ 1,0,0,0 ,0,1,0,0 ,0,0,1,0 ,0,0,0,1]
        return type(self)()


    @property
    def np_inverse(self):
        pass

    #@property
    def test_index(self):
        """ fill matrix with incermenting number to see indices """
        tmp = type(self)()
        for i in range(16):
            tmp[i] = i
        return tmp


    def serialize(self, inarray):
        """ serialize this array into a list """

        out = []
        for i in self.m:
            out.append(i)
        return out

    def insert(self, iterable):
        """ UNTESTED - load the first 16 things we find into this matrix """
       
         #numpy ND array 
        if isinstance(iterable, np.ndarray):
            out = [];idx=0
            for row in iterable:
                for col in row:
                    self.m[idx] = col
                    idx+=1

        #serialized simple array
        if isinstance(iterable, list) or isinstance(iterable, tuple):
            for idx,i in enumerate(iterable):
                if idx <= 16:
                     self.m[idx] = iterable[idx] 

    #@property
    def copy(self, mtype=None):
        """ UNTESTED """
        if not mtype:        
            return type(self)(
                self.m[0] , self.m[1] , self.m[2] , self.m[3] ,
                self.m[4] , self.m[5] , self.m[6] , self.m[7] ,
                self.m[8] , self.m[9] , self.m[10], self.m[11],
                self.m[12], self.m[13], self.m[14], self.m[15]
            )

        if mtype == 'numpy':
            return np.array((
                [self.m[0] , self.m[1] , self.m[2] , self.m[3]  ],
                [self.m[4] , self.m[5] , self.m[6] , self.m[7]  ],
                [self.m[8] , self.m[9] , self.m[10], self.m[11] ],
                [self.m[12], self.m[13], self.m[14], self.m[15] ]
            ))               

    @property
    def transpose(self):
        return type(self)(
            self.m[0], self.m[4], self.m[8] , self.m[12],
            self.m[1], self.m[5], self.m[9] , self.m[13],
            self.m[2], self.m[6], self.m[10], self.m[14],
            self.m[3], self.m[7], self.m[11], self.m[15]
        )


    @property
    def determinant(self):
        """
            Laplace expansion method 

            https://www.mathsisfun.com/algebra/matrix-determinant.html

            plus  a times the determinant of the matrix that is not in a's row or column,
            minus b times the determinant of the matrix that is not in b's row or column,
            plus  c times the determinant of the matrix that is not in c's row or column,
            minus d times the determinant of the matrix that is not in d's row or column,

                       |               | A*        -     B*      +        C*    -         D*|
            a  b  c  d | 0   1  2   3  | a -->     |  <--b-->    |     <--c-->  |      <--d |  
            e  f  g  h | 4   5  6   7  | |  f g h  |  e  |  g  h |  e  f  |  h  | e  f  g | |
            i  j  k  l | 8   9  10  11 |    j k l  |  i     k  l |  i  j     l  | i  j  k   |
            m  n  o  p | 12  13 14  15 |    n o p  |  m     o  p |  m  n     p  | m  n  o   | 

        """  


        def det33(a):
            #same as matrix33 determinant method 
            o = a[0]* ((a[4]*a[8])-(a[5]*a[7])) - a[1]*((a[3]*a[8])-(a[5]*a[6])) +  a[2]*((a[3]*a[7])-(a[4]*a[6]))
            return o

        a = self.copy() 
        
        a33 = [a[5],a[6],a[7],a[9],a[10],a[11],a[13],a[14],a[15] ]
        b33 = [a[4],a[6],a[7],a[8],a[10],a[11],a[12],a[14],a[15] ]
        c33 = [a[4],a[5],a[7],a[8],a[9] ,a[11],a[12],a[13],a[15] ]
        d33 = [a[4],a[5],a[6],a[8],a[9] ,a[10],a[12],a[13],a[14] ]
        o = (a[0] * det33(a33)) - (a[1]*det33(b33)) + (a[2]*det33(c33)) - (a[3]*det33(d33))
        return o

    def batch_mult_pts(self, pts):
        """ sub component of matrix rotate function 
            multiply a 3X3 matrix by a list of points 
        """

        #debug - make work with other types, like self.insert() 

        tmp_buffer = []
        out = None
        for pvec in pts:  
            tmp_buffer.append( self * pvec )
        return tmp_buffer


    def rotate_pts_3d(self, points, xrot, yrot, zrot):
        """
           -  previously named rotate_mat4
           
           using the "standard" 16 element, Row major, 4X4 rotation matrix used by Maya
         
           [0  1  2  3]     [   XVEC    0 ]
           [4  5  6  7]     [   YVEC    0 ]
           [8  9  10 11]    [   ZVEC    0 ]
           [12 13 14 15]    [ 0  0  0   0 ]

           [0  1  2  3]      xx xy xz 0
           [4  5  6  7]      yx yy yz 0
           [8  9  10 11]     zx zy zz 0
           [12 13 14 15]     0  0  0  0
           ------------------------------
           rotate Y matrix     
           |  cos(y)  0      -sin(y)  0 | 
           |  0       1       0       0 | 
           |  sin(y)  0       cos(y)  0 | 
           |  0       0       0       1 | 
           ------------------------------
           rotate Z  matrix 
           |  cos(z)  sin(z)  0       0 | 
           | -sin(z)  cos(z)  0       0 |
           |  0       0       1       0 |
           |  0       0       0       1 |
           ------------------------------
           rotate X matrix  
           |  1       0       0       0 |  
           |  0       cos(x)  sin(x)  0 |  
           |  0      -sin(x)  cos(x)  0 |  
           |  0       0       0       1 | 
        """
        dtr = self.mu.dtr

        ####
        #build rotationY (see diagram above) 
        y_matrix     =  self.identity
        y_matrix[0]  =  math.cos(dtr( yrot ))
        y_matrix[2]  = -math.sin(dtr( yrot ))
        y_matrix[8]  =  math.sin(dtr( yrot ))
        y_matrix[10] =  math.cos(dtr( yrot ))

        ####                
        #build rotationZ (see diagram above) 
        z_matrix    =  self.identity
        z_matrix[0] =  math.cos(dtr( zrot ))
        z_matrix[1] =  math.sin(dtr( zrot ))
        z_matrix[4] = -math.sin(dtr( zrot ))
        z_matrix[5] =  math.cos(dtr( zrot ))
        tmp_matr = y_matrix * z_matrix 

        ####
        #build rotationX (see diagram above) 
        x_matrix     =  self.identity
        x_matrix[5]  =   math.cos(dtr( xrot )) 
        x_matrix[6]  =   math.sin(dtr( xrot )) 
        x_matrix[9]  =  -math.sin(dtr( xrot ))
        x_matrix[10] =   math.cos(dtr( xrot ))
        rotation_44 = x_matrix * tmp_matr

        ############ 
        return rotation_44.batch_mult_pts(points)


    def buildPerspProjMat(self, fov, aspect, znear, zfar):
        """
            UNTESTED 
            
            http://stackoverflow.com/questions/8633034/basic-render-3d-perspective-projection-onto-2d-screen-with-camera-without-openg

            use homogenous transformations and coordinates. 

            take a point in space and:
                Position it relative to the camera using the model matrix.
                Project it either orthographically or in perspective using the projection matrix.
                Apply the viewport transformation to place it on the screen.
        """

        PI_OVER_360 = 0.00872664625
        xymax = 100
        ymax = znear * math.tan(fov * PI_OVER_360);

        ymin = -xymax;
        xmin = -xymax;

        width = xymax - xmin;
        height = xymax - ymin;

        depth =   zfar - znear;
        q     = -(zfar + znear) / depth;
        qn    = -2 * (zfar * znear) / depth;

        w = 2 * znear / width;
        w = w / aspect;
        h = 2 * znear / height;

        m = self.identity

        m[0]  = w; m[4]  = 0; m[8]  = 0 ; m[12] = 0
        m[1]  = 0; m[5]  = h; m[9]  = 0 ; m[13] = 0
        m[2]  = 0; m[6]  = 0; m[10] = q ; m[14] = -1
        m[3]  = 0; m[7]  = 0; m[11] = qn;
        #print(m)
        return m
