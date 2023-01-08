#!/usr/local/bin/python3


import itertools 
import math
import os

from gnelscript import NUMPY_IS_LOADED, SCIPY_IS_LOADED


if NUMPY_IS_LOADED:
    # print(' ## debug - loading numpy module. ')
    import numpy as np             #for testing - remove later 
    from numpy.linalg import inv

if SCIPY_IS_LOADED:    
    from scipy.linalg import expm

else:
    print(' ## debug - numpy module disabled. ')




DEG_TO_RAD = 0.0174532925 # degree = radian * (180 / PI) # PI = 3.14159265
RAD_TO_DEG = 57.29577951  # radian = degree * (PI/180)


###############################################
class math_util(object):    
    """ general math library - may contain some of the same functions 
        as the vector and matrix objects, but those will be implemented 
        to operate relative to themselves, these will work more generally

        for example self.dot_vec3() will be vec3.dot , etc 
    """

    ##-------------------------------------------##

    # # iterate through rows of X
    # for i in range(len(X)):
    #    # iterate through columns of Y
    #    for j in range(len(Y[0])):
    #        # iterate through rows of Y
    #        for k in range(len(Y)):
    #            result[i][j] += X[i][k] * Y[k][j]

    ##-------------------------------------------##

    # def matmult(a,b):
    #     zip_b = zip(*b)
    #     # uncomment next line if python 3 : 
    #     # zip_b = list(zip_b)
    #     return [[sum(ele_a*ele_b for ele_a, ele_b in zip(row_a, col_b)) 
    #              for col_b in zip_b] for row_a in a]

    ##-------------------------------------------##
    def calc_line_length(self, x1, y1, x2, y2):
        """ distance between two points 
            DEBUG -- make work in 2d and 3d!
            DEBUG -- merge with -  self.calc_line_length(pt1[0], pt1[1], pt2[0], pt2[1] )

        """
        return math.sqrt( ((x1-x2)**2)+ ((y1-y2)**2) )

    ##-------------------------------------------##        
    def dot_vec3 (self, v1, v2):
         """ scalar - mag but not direction """
         return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2] 

    ##-------------------------------------------##
    def cross_vec3 (self, v1, v2):
        return (v1[1]*v2[2]-v1[2]*v2[1], v1[2]*v2[0]-v1[0]*v2[2], v1[0]*v2[1]-v1[1]*v2[0])

    ##-------------------------------------------##
    def dtr (self, deg ):
       return deg * DEG_TO_RAD

    ##-------------------------------------------##
    def rtd (self, rad ):
       return rad * RAD_TO_DEG
    
    ##-------------------------------------------##
    def dtr_vec3(self, invec):
        return ( self.dtr( invec[0] ),
                 self.dtr( invec[1] ),
                 self.dtr( invec[2] ),
               )
    ##-------------------------------------------##
    def rtd_vec3(self, invec):
        return ( self.rtd( invec[0] ),
                 self.rtd( invec[1] ),
                 self.rtd( invec[2] ),
               )
    ##-------------------------------------------##
    def vect_add(self, v1, v2):
        return [v1[0]+v2[0], v1[0]+v2[1], v1[2]+v2[2]]
    ##-------------------------------------------##
    def mult_scalar (self, scalar, v):
        return [v[0]*scalar, v[1]*scalar, v[2]*scalar ]
    ##-------------------------------------------##
    def mult_m33_vec3(self, m, v):
        """ multiplies a 3X3 matrix by a 3D vector - returns a vector tuple 
            NUMPY DOT is the same as multiplying 
        """
        
        outx = m[0] * v[0] + m[3] * v[1] + m[6] * v[2] 
        outy = m[1] * v[0] + m[4] * v[1] + m[7] * v[2] 
        outz = m[2] * v[0] + m[5] * v[1] + m[8] * v[2] 
          
        return  (outx, outy, outz)
    ##-------------------------------------------##
    def mult_m44_vec3(self, m, v):
        """ multiplies a 4X4 matrix by a 3D vector - returns a vector tuple """
        
        outx = m[0]*v[0] + m[4]*v[1] + m[8]  * v[2]+m[12]
        outy = m[1]*v[0] + m[5]*v[1] + m[9]  * v[2]+m[13]
        outz = m[2]*v[0] + m[6]*v[1] + m[10] * v[2]+m[14]
          
        return  (outx, outy, outz)
    ##-------------------------------------------##
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
    ##-------------------------------------------##
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

    ##-------------------------------------------##
    def normalize_vec3(self, in_vec):
        try:
           invLength = 1.0/math.sqrt(in_vec[0]*in_vec[0] + in_vec[1]*in_vec[1]+ in_vec[2]*in_vec[2])
           return (in_vec[0] *invLength, in_vec[1] * invLength,  in_vec[2] * invLength)
        except:
           print('normalize_vec3: divide by zero error.') 
           return [0,0,1]

    ##-------------------------------------------##
    def mul_square_matrices(self, m1, m2):
        
        """
        # From the brilliant Bryant Ewert's script 
        matrix_playground.mel 

        matrix $cMatrix[4][4] = << 0, 0, 0, 0;  0, 0, 0, 0;  0, 0, 0, 0;  0, 0, 0, 0 >>;
        for ( $i = 0; $i < 4; $i++ )
        {
            for ( $j = 0; $j < 4; $j++ )
            {
                for ( $k = 0; $k < 4; $k++ )
                {
                    $cMatrix[$i][$j] = $cMatrix[$i][$j] + ( $aMatrix[$i][$k] * $bMatrix[$k][$j] );
                }
            }
        }
        return $cMatrix;
        """
        pass

    ##-------------------------------------------##
    def matrix_adjoint(self, m44, size):
        
        """
        # From the brilliant Bryant Ewert's script 
        matrix_playground.mel 

        # Returns a matrix which is the Adjoint of $aMatrix.
        # The $size of the input matrix must be specified (3 or 4).

        """
        cMatrix = matrix44()
        detSize = size - 1 

        if size > 2:
            """
            # Cofactor of top-left 3×3 matrix
            $cMatrix[0][0] =  determinant( matrixCollapseRowColumn( $aMatrix, $size, 1, 1 ), $detSize );
            $cMatrix[0][1] = -determinant( matrixCollapseRowColumn( $aMatrix, $size, 1, 2 ), $detSize );
            $cMatrix[0][2] =  determinant( matrixCollapseRowColumn( $aMatrix, $size, 1, 3 ), $detSize );

            $cMatrix[1][0] = -determinant( matrixCollapseRowColumn( $aMatrix, $size, 2, 1 ), $detSize );
            $cMatrix[1][1] =  determinant( matrixCollapseRowColumn( $aMatrix, $size, 2, 2 ), $detSize );
            $cMatrix[1][2] = -determinant( matrixCollapseRowColumn( $aMatrix, $size, 2, 3 ), $detSize );

            $cMatrix[2][0] =  determinant( matrixCollapseRowColumn( $aMatrix, $size, 3, 1 ), $detSize );
            $cMatrix[2][1] = -determinant( matrixCollapseRowColumn( $aMatrix, $size, 3, 2 ), $detSize );
            $cMatrix[2][2] =  determinant( matrixCollapseRowColumn( $aMatrix, $size, 3, 3 ), $detSize );
            """
            pass

        if size > 3:
            """
            # Cofactor of 4th column
            $cMatrix[0][3] = -determinant( matrixCollapseRowColumn( $aMatrix, $size, 1, 4 ), $detSize );

            $cMatrix[1][3] =  determinant( matrixCollapseRowColumn( $aMatrix, $size, 2, 4 ), $detSize );

            $cMatrix[2][3] = -determinant( matrixCollapseRowColumn( $aMatrix, $size, 3, 4 ), $detSize );

            # Cofactor of 4th row
            $cMatrix[3][0] = -determinant( matrixCollapseRowColumn( $aMatrix, $size, 4, 1 ), $detSize );
            $cMatrix[3][1] =  determinant( matrixCollapseRowColumn( $aMatrix, $size, 4, 2 ), $detSize );
            $cMatrix[3][2] = -determinant( matrixCollapseRowColumn( $aMatrix, $size, 4, 3 ), $detSize );
            $cMatrix[3][3] =  determinant( matrixCollapseRowColumn( $aMatrix, $size, 4, 4 ), $detSize );
            """ 
            pass 

        # Adjoint is TRANSPOSE of matrix containing cofactors
        #cMatrix = transpose( $cMatrix, $size );
        return cMatrix

    ##-------------------------------------------##

    def matrix_invert(self, aMatrix, size):
        
        """
        # From the brilliant Bryant Ewert's script 
        matrix_playground.mel 

        # Returns a matrix which is the Inverse of $aMatrix.
        # The $size of the input matrix must be specified (3 or 4).

        """

        iMatrix  = matrix44()
    
        #determinant = determinant( aMatrix, size )
    
        if  determinant != 0.0: 
            #iMatrix = ( 1 / determinant( aMatrix, size ) ) * adjoint( aMatrix, size )
            pass

        return iMatrix;

    ##-------------------------------------------##
    def matrix_rotate(self, aMatrix, size):
        
        """
        # From the brilliant Bryant Ewert's script 
        matrix_playground.mel 

        # Applies a Rotation Transformation to the specified Matrix ("A", "B" or "C").
        # The rotation value and unit is derived from the current UI settings.
        # Note: The rotation matrices used may seem Transposed to those typically
        #    documented, but they are correct for this implementation within Maya,
        #    specifically in regard to the Acquire and Apply functions (above).

        """

    
        sAxis = ( "X", "Y", "Z" )
        kRADIAN = 57.295779513082320876798154814105;
        aMatrix = matrix44()
        rMatrix = matrix44() 
   
        """ 
        if ( $unit == 1 )       // must convert to degrees
            $rotate = $fRotate / $kRADIAN;
            
        float $cos = `cos $rotate`;
        float $sin = `sin $rotate`;
            
        switch ( $axis )
        {
            case 1:     // X axis
                $rMatrix = << 1, 0, 0, 0;  0, $cos, $sin, 0;  0, -$sin, $cos, 0;  0, 0, 0, 1 >>;
                break;
                
            case 2:     // Y axis
                $rMatrix = << $cos, 0, -$sin, 0;  0, 1, 0, 0;  $sin, 0, $cos, 0;  0, 0, 0, 1 >>;
                break;
                
            case 3:     // Z axis
                $rMatrix = << $cos, $sin, 0, 0;  -$sin, $cos, 0, 0;  0, 0, 1, 0;  0, 0, 0, 1 >>;
                break;
        }
        
        $aMatrix = getMatrix( $which );
        $rMatrix = multiplyMatrix( $aMatrix, $rMatrix );
        
        populateMatrix( "C", $rMatrix, ( "Rotate Matrix A " + $fRotate + ( $unit == 1 ? " deg" : " rad" ) + " on " + $sAxis[$axis-1] ) );
        """

###############################################
class vec2(object):    

    def __init__(self,x=0,y=0):
        self.mu = math_util()
        self.x = x
        self.y = y    
        
    def __repr__(self):
        return '(%s, %s)' % (self.x, self.y)

    def __abs__(self):
        return type(self)(abs(self.x), abs(self.y))

    def __add__(self, other):
        return type(self)(self.x + other.x, self.y + other.y)

    def __sub__(self, other):
        return type(self)(self.x - other.x, self.y - other.y)

    def __mul__(self, other):
        
        #if isinstance(other, matrix22):
        #    matrix22

        print('## debug vec2 mult type other is ' , type(other) )

        if isinstance(other, tuple) or isinstance(other, list):
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

    def distance_to(self, other):
        """ """
        val = math.hypot((self.x - other.x), (self.y - other.y))
        return val

    def project_pt(self, A, B, offset):
        """ project a point along a vector """

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
        return type(self)(ptX, ptY)

    def rayintersectseg(p, edge):
        ''' UNTESTED 
            FROM https://rosettacode.org/wiki/Ray-casting_algorithm#Python
             takes a point p=Pt() and an edge of two endpoints a,b=Pt() of a line segment returns boolean
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
        """ untested - normalize a 2d vector """ 
        invLength = 1.0/math.sqrt(self.x*self.x + self.y*self.y)
        self.x *= invLength
        self.y *= invLength

##-------------------------------------------##
##-------------------------------------------##

class vec3(object):    

    def __init__(self,x=0,y=0,z=0):
        #this is sloppy - check the first item and assume we are initializing with a tuple xyz 
        if type(x) is tuple:
            self.x=x[0];self.y=x[1];self.z=x[2]         
        else:    
            self.x=x;self.y=y;self.z=z  
        self.mu = math_util() 

    def __repr__(self):
        return '(%s, %s, %s)' % (self.x, self.y, self.z)

    def __abs__(self):
        return type(self)(abs(self.x), abs(self.y), abs(self.z))

    def __add__(self, other):
        if NUMPY_IS_LOADED:
            if isinstance(other, np.ndarray):
                return type(self)(self.x+other[0], self.y+other[1], self.z+other[2])
        if isinstance(other, float) or isinstance(other, int):
            return type(self)(self.x+other, self.y+other, self.z+other)
        if isinstance(other, vec3):                    
            return type(self)(self.x+other.x, self.y+other.y, self.z+other.z)
        if isinstance(other, tuple) or isinstance(other, list):
            return type(self)(self.x+other[0], self.y+other[1], self.z+other[2])  

    def __sub__(self, other):
        if NUMPY_IS_LOADED:
            if isinstance(other, np.ndarray):
                return type(self)(self.x-other[0], self.y-other[1], self.z-other[2])
        if isinstance(other, float) or isinstance(other, int):
            return type(self)(self.x-other, self.y-other, self.z-other)
        if isinstance(other, vec3):                    
            return type(self)(self.x-other.x, self.y-other.y, self.z-other.z)
        if isinstance(other, tuple) or isinstance(other, list):
            return type(self)(self.x-other[0], self.y-other[1], self.z-other[2])  

    def __mul__(self, other):
        if NUMPY_IS_LOADED:
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
    def aspt(self):
        return (self.x, self.y, self.z)


    #def is_parallel(self, vec):

    def is_orthogonal(self, vec):
        """UNTESTED
        #https://www.learndatasci.com/glossary/orthogonal-and-orthonormal-vectors/
        """

        if self.dot(vec)==0:
            return True 
        else:
            return False     

    def to_euler(self):
        """ 
            UNFINISHED 
            experimental - inspired by this:
            https://stackoverflow.com/questions/21622956/how-to-convert-direction-vector-to-euler-angles
        
        """
        
        """
        D = (XD,YD,ZD) # direction 
        U = (XU,YU,ZU) # up
        H = (XD,YD,0)  # heading 
        
        angle_H = math.atan2(YD,XD)
        
        ZD=sin(angle_P)
        
        angle_P=asin(ZD)

        W0 = ( -YD, XD, 0 )
        U0 = W0 × D

        cos(angle_B) = Dot(U0,U) / abs(U0) / abs(U)
        sin(angle_B) = Dot(W0,U) / abs(W0) / abs(U)
        angle_B = atan2( Dot(W0,U) / abs(W0), Dot(U0,U) / abs(U0) )

        """
        pass  



    def project_pt(self, offset):
        """ project a point along a vector """
        
        # nX = B.x - A.x
        # nY = B.y - A.y
        # nZ = B.z - A.z 

        # distX = pow( (A.x - B.x ) , 2.0 ) 
        # distY = pow( (A.y - B.y ) , 2.0 ) 
        # distZ = pow( (A.z - B.z ) , 2.0 ) 
        # vecLength = math.sqrt(distX + distY + distZ )

        # normalized vector  
        # calcX = nX / self.length
        # calcY = nY / self.length
        # calcZ = nZ / self.length

        normal = self.normal

        # project point along vector with offset (can use negative too)
        ptX = self.x + (normal.x * offset)
        ptY = self.y + (normal.y * offset)
        ptZ = self.z + (normal.z * offset)
        return type(self)(ptX, ptY, ptZ)

    def project_vec3(self):
        """  UNFINISHED  

            from David Gould's book
            Complete Maya Programing II 
            page 32 
        

        """

        pass
    
    ##----------------------------------------

    def between(self, pt2):
        """ given 2 points in 3D, create a 3D vector 
            representing the offset between them 

            doesnt get much easier than this, just subtract 
             
            usage: 
                v1 = vec3(0,0,0)
                v2 = vec3(1,1,1)
                print(v2.between(v1) )

        """
   
        #if isinstance(pt1, tuple):
        #    pt1 = vec3( pt1[0], pt1[1], pt1[2]) 


        if isinstance(pt2, tuple):
            pt2 = vec3( pt2[0], pt2[1], pt2[2])

        return pt2 - self

    ##----------------------------------------
    def orthogonal_vec_from_pt(self, vecpt, unitvec, pt ):
        if NUMPY_IS_LOADED:        
            return (vecpt-pt) - ( np.dot((vecpt-pt), unitvec) ) * unitvec
        return None 
    
    ##----------------------------------------

    if NUMPY_IS_LOADED:
        @property
        def as_np(self):
            """  get this vec3 as an np.array """ 
            return self.copy(vtype='numpy')

    ##----------------------------------------

    def insert(self, iterable):
        """ convert an np.array, tuple or list  to vec3  
            does not check size, so just assume 3 items (x,y,z)
        """

        if isinstance(iterable, list) or isinstance(iterable, tuple):
            self.x = iterable[0]
            self.y = iterable[1]            
            self.z = iterable[2]

        if NUMPY_IS_LOADED:
            if isinstance(iterable, np.ndarray):
                self.x = iterable[0]
                self.y = iterable[1]            
                self.z = iterable[2]

        return self 

     ##----------------------------------------
    def copy(self, vtype=None):
        if vtype == None:
            return type(self)(self.x,self.y,self.z)

        if vtype == 'tuple':
            return ( (self.x,self.y,self.z) )

        if NUMPY_IS_LOADED:    
            if vtype == 'numpy':
                return np.array( (self.x,self.y,self.z) )

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

    ##------------------------------------- 
    @property
    def np_normal(self):
        """ unit vector of the vector using numpy """
        if NUMPY_IS_LOADED:
            return self.as_np / np.linalg.norm(self.as_np)
        else:
            pass 

    ##------------------------------------- 
    def lookat(self, pt):
        """
        UNFINISHED/NOT WORKING 

        https://stackoverflow.com/questions/1251828/calculate-rotations-to-look-at-a-3d-point
        
        https://math.stackexchange.com/questions/3139155/finding-a-vector-that-points-towards-a-coordinate

        ##--- 

        rotx = Math.atan2( y, z )
        roty = Math.atan2( x * Math.cos(rotx), z )
        rotz = Math.atan2( Math.cos(rotx), Math.sin(rotx) * Math.sin(roty) )

        About X: -atan2(y, z)
        About Y: atan2(x, sqrt(y*y + z*z))
        About Z: 0 

        """

        x = pt[0]
        y = pt[1]
        z = pt[2]

        rotx = math.atan2( y, z )
        if z >= 0:
            roty = -math.atan2( x * math.cos(rotx), z );
        else:
            roty = math.atan2( x * math.cos(rotx), -z );
        
        return [rotx,roty]


    ##------------------------------------- 
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

        if NUMPY_IS_LOADED:        
            if isinstance(v1, tuple):
                v1 = self.insert(v1)  #wrong?
            if isinstance(v2, tuple):
                v2 = self.insert(v2)  #wrong?

            v1_u = v1.np_normal;v2_u = v2.np_normal

            #print('### ', v1_u , v2_u )

            return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

    ##------------------------------------- 
    def angle_between(self, other):
        """ 
            result is in radians
            derived from law of cosines 
            range from 0 (colinear) to PI radians (180 deg)  

            FAILS WHEN 0 or 180 ?
        """
        try:
            o = math.acos( self.dot(other) / (self.length*other.length) ) 
            return o 
        except:
            return 0

    ##------------------------------------- 
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
    
    ##------------------------------------- 
    def ray_tri_intersect(self, orig, dir, v0, v1, v2):
        """ taken from C code at 
            https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/ray-triangle-intersection-geometric-solution.html 
        """

        t = 0 

        #DEBUG  a but hacky 
        if type(dir) is not vec3:
            dir = vec3(dir)
        if type(orig) is not vec3:
            orig = vec3(orig)

        # compute the plane's normal
        v0v1 = vec3()
        v0v2 = vec3()
        
        v0v1 = v1-v0
        v0v2 = v2-v0
        
        # no need to normalize
        N =vec3() 
        N = v0v1.cross(v0v2) 

        #// Step 1: check if the ray and plane are parallel.
        n_dot_ray_dir = N.dot(dir)

        # check almost 0
        #if abs(n_dot_ray_dir) < kEpsilon: 
        if n_dot_ray_dir <.001:
            # they are parallel, so they don't intersect! 
            return False 

        # compute (d) parameter using equation 2
        d = -N.dot(v0)
    
        # compute (t) - distance from the ray origin to the intersection point
        t = -(N.dot(orig) + d) / n_dot_ray_dir
        
        # check if the triangle is behind the ray
        # If is greater than 0, the triangle is "visible" to that ray
        if t < 0:
            return False 
     
        # compute the intersection point  
        P = vec3() 
        P = orig+(dir*t)
        
        #print("## distance %s point "%t, P)
        ##-----

        # inside-outside test - vector perpendicular to triangle's plane
        C = vec3() 
     
        # edge 0
        edge0 = vec3()
        edge0 = v1 - v0 
        vp0 =vec3()
        vp0 = P - v0
        C = edge0.cross(vp0)
        if N.dot(C) < 0:
            # P is on the right side
            return False; 
     
        # edge 1
        edge1 = vec3()
        edge1 = v2 - v1 
        vp1 = vec3() 
        vp1 = P - v1
        C = edge1.cross(vp1)
        if N.dot(C) < 0:  
            # P is on the right side
            return False; 
     
        # edge 2
        edge2 = vec3()
        edge2 = v0 - v2 
        vp2 = vec3()        
        vp2 = P - v2
        C = edge2.cross(vp2)
        if N.dot(C) < 0:
            # P is on the right side;
            return False; 
        
        # this ray hits the triangle
        return P

    

    ##------------------------------------- 

    def poly_intersect(self, ray, poly):
        """ DEBUG UNTESTED - 

            FROM :  https://stackoverflow.com/questions/312328/what-is-the-fastest-way-to-find-the-point-of-intersection-between-a-ray-and-a-po
            SEE ALSO: https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/ray-triangle-intersection-geometric-solution.html

        """

        point_is_inside = False
        
        e0=poly[0]
        e1=poly[1]
        e2=poly[2]

        #e0=e0.normal
        #e1=e1.normal
        #e2=e2.normal

        r1 = ray[0]
        r2 = ray[1]
        
        #r1=r1.normal
        #r2=r2.normal

        #point plane_normal = crossProduct(vectorSub(Poly.P[1], Poly.P[0]), vectorSub(Poly.P[2], Poly.P[0]))
        tmp = e1-e0
        plane_normal = tmp.cross(e2-e0)
        
        print("# poly_intersect normal ", plane_normal , ' length ', plane_normal.length )
        #keith is experimenting here - try normalizing?
        #plane_normal = plane_normal.normal


        #float denominator = dotProduct(vectorSub(Ray.R2, Poly.P[0]), plane_normal)
        tmp = r2-e0
        denominator = tmp.dot(plane_normal)

        ## ray is parallel to the polygon (compare to face normal to get this)
        if denominator == 0:
            print("ray is parallel to the polygon")
            return False 

        #float ray_scalar = dotProduct(vectorSub(Poly.P[0], Ray.R1), plane_normal)
        tmp = e0-r1
        ray_scalar = tmp.dot(plane_normal)

        #Answer = vectorAdd(Ray.R1, scalarMult(ray_scalar, Ray.R2))
        answer = r1 + (r2*ray_scalar)
  
        ##---- 

        # verify that the point falls inside the polygon
        #point test_line = vectorSub(Answer, Poly.P[0])
        test_line = answer-e0

        #point test_axis = crossProduct(plane_normal, test_line)
        test_axis = plane_normal.cross(test_line)

        #point test_point = vectorSub(Poly.P[1], Answer)
        test_point = e1-answer

        #bool prev_point_ahead = (dotProduct(test_line, test_point) > 0)
        prev_point_ahead = test_line.dot(test_point)

        #bool prev_point_above = (dotProduct(test_axis, test_point) > 0)
        prev_point_above = test_axis.dot(test_point)

        this_point_ahead = False
        this_point_above = False
        

        index = 2
        while index < len(poly):
            test_point = poly[index]-answer
            this_point_ahead = test_line.dot(test_point)>0
            if prev_point_ahead or this_point_ahead:
                this_point_above = test_axis.dot(test_point)>0
                if prev_point_above != this_point_above:
                    point_is_inside = not point_is_inside

            prev_point_ahead = this_point_ahead
            prev_point_above = this_point_above
            index+=1
 
 

        return [point_is_inside, answer, plane_normal ]



    ##------------------------------------- 


    ##------------------------------------- 


    """
    # https://stackoverflow.com/questions/2549708/intersections-of-3d-polygons-in-python

    if NUMPY_IS_LOADED:
        def np_intersection(self, ray):

           def cmp_floats(a,b, atol=1e-12):
               return abs(a-b) < atol

           # Returns a intersection point with a ray and the polygon. 

           n = self.normal()

           #Ray is parallel to the polygon
           if cmp_floats( np.dot( np.array(ray.direction), n ), 0. ):
               return None

           t = 1/(np.dot(np.array(ray.direction),n)) * ( np.dot(n,np.array(self.pts[0])) - np.dot(n,np.array(ray.position)) )
           
           #Intersection point is behind the ray
           if t < 0.0:
               return None

           #Calculate intersection point
           point = np.array(ray.position) + t*np.array(ray.direction)
           
           #Check if intersection point is really in the polygon or only on the (infinite) plane
           if self.on_surface(point):
               return [list(point)]

           return None

    def is_planar(self):
        #the determinant of the vectors (volume) must always be 0
        x_i = np.array(self.pts[i])
        x_i1 = np.array(self.pts[i+1])
        x_i2 = np.array(self.pts[i+2])
        det = np.linalg.det([x_0-x_i, x_0-x_i1, x_0-x_i2])
        assert cmp_floats( det, 0.0 ), "Points must be in a plane to create a Polygon"


    def on_surface(self, point):
        # Returns True if the point is on the polygon's surface and false otherwise. 
        n = len(self.pts)
        anglesum = 0
        p = np.array(point)

        for i in range(n):
            v1 = np.array(self.pts[i]) - p
            v2 = np.array(self.pts[(i+1)%n]) - p

            m1 = magnitude(v1)
            m2 = magnitude(v2)

            if cmp_floats( m1*m2 , 0. ):
                return True #point is one of the nodes
            else:
                # angle(normal, vector)
                costheta = np.dot(v1,v2)/(m1*m2)
            anglesum = anglesum + np.arccos(costheta)
        return cmp_floats( anglesum , 2*np.pi )

    """
            
    ##----------------------




##-------------------------------------------##
##-------------------------------------------##

class vec4(object):
    """ untested - 
        homogeneous coordinate experiment  
    """

    def __init__(self,x=0,y=0,z=0,w=1):
        self.x=x
        self.y=y
        self.z=z
        self.w=w  
        
        #self.mu = math_util() 

    def __mul__(self, other):
        # https://www.tomdalling.com/blog/modern-opengl/explaining-homogenous-coordinates-and-projective-geometry/
        
        if isinstance(other, float) or isinstance(other, int):
            #return type(self)(self.x*other, self.y*other, self.z*other)
            
            new_x = self.x * other 
            new_y = self.y * other 
            new_z = self.z * other 
            new_w = self.w * other 
            return type(self)( (new_x/new_w), (new_y/new_w), (new_z/new_w), new_w )


    def __repr__(self):
        return '(%s, %s, %s, %s)' % (self.x, self.y, self.z, self.w)

    def __getitem__(self, index):
        if index==0:
            return self.x
        if index==1:
            return self.y
        if index==2:
            return self.z
        if index==3:
            return self.w

    def __setitem__(self, key, item):
        if key==0:
            self.x = item
        if key==1:
            self.y = item
        if key==2:
            self.z = item
        if key==3:
            self.w = item

    def insert(self, iterable):
        """ convert an np.array, tuple or list  to vec3  
            does not check size, so just assume 3 items (x,y,z)
        """
        if isinstance(iterable, vec3):
            self.from_vec3(iterable) 

        if isinstance(iterable, list) or isinstance(iterable, tuple):
            self.x = iterable[0]
            self.y = iterable[1]            
            self.z = iterable[2]

            if len(iterable)==3: 
                self.w = 1  
            if len(iterable)==4:
                self.w = iterable[3]                 

        if NUMPY_IS_LOADED:
            if isinstance(iterable, np.ndarray):
                self.x = iterable[0]
                self.y = iterable[1]            
                self.z = iterable[2]
                if len(iterable)==3: 
                    self.w = 1  
                if len(iterable)==4:
                    self.w = iterable[3] 
        return self 

    def to_vec3(self):
        """ vec4 to vec3 - 
            divide the weight across the other 3 axis 
        """
        w = self.w
        return vec3( (self.x/w), (self.y/w), (self.z/w) )

    def from_vec3(self, vec3):
        """ vec4 from a vec3, simply add a weight of 1 """
        self.x = vec3.x        
        self.y = vec3.y  
        self.z = vec3.z
        self.w = 1  

#[xyzw]∗⎡⎣⎢⎢⎢m00m10m20m30m01m11m21m31m02m12m22m32m03m13m23m33⎤⎦⎥⎥⎥
#x′=x∗m00+y∗m10+z∗m20+w∗m30y′=x∗m01+y∗m11+z∗m21+w∗m31z′=x∗m02+y∗m12+z∗m22+w∗m32w′=x∗m03+y∗m13+z∗m23+w∗m33


##-------------------------------------------##
##-------------------------------------------##

class matrix22(object):
    """ 2D matrix experiment """

    def __init__(self, a=1, b=0, c=0, d=1):
        self.m = [a,b,c,d]
        self.mu = math_util()
    
    def __getitem__(self, index):
        return self.m[index]

    def __setitem__(self, key, item):
        self.m[key] = item

    def __repr__(self):
        return '(%s, %s, %s, %s)'%(self.m[0], self.m[1], self.m[2],  self.m[3])

    @property
    def identity(self):
        return type(self)()

    @property
    def transpose(self):
        """
        UNTESTED 
        # standard indicies  |  #transposed indicies
        1  2                 |   1 3 
        3  4                 |   2 4 
       
        """

        return type(self)(
            self.m[0], self.m[1],
            self.m[2], self.m[3]
        )


    def __add__(self, other):
        """ UNTESTED """
        return type(self)(
            self.m[0]+n[0], self.m[2]+n[1],
            self.m[1]+n[2], self.m[3]+n[3]
        )
    
    def __sub__(self, other):
        """ UNTESTED """        
        return type(self)(
            self.m[0]-n[0], self.m[1]-n[2],
            self.m[2]-n[1], self.m[3]-n[3]
        )

    # def __truediv__(self, other):
    #     # matrices don't divide! there is no concept of dividing by a matrix.
    #     # multiply by an inverse, which achieves the same thing.

    def __mul__(self, n):
        """ multiply two 2X2 matricies together 

            the order that you mutliply will call a different object!!!
            make sure you do "this * other", NOT "other * this"

        """

        if isinstance(n, vec2):
            outx = self.m[0]*n.x + self.m[2]*n.y   
            outy = self.m[1]*n.x + self.m[3]*n.y  
            return  (outx, outy)

        if isinstance(n, tuple) or isinstance(n, list):
            outx = self.m[0]*n[0] + self.m[2]*n[1]   
            outy = self.m[1]*n[0] + self.m[3]*n[1]  
            return  (outx, outy)

        if type(n) == type(self):
            return type(self)(
                    self.m[0]*n[0] + self.m[1]*n[2], 
                    self.m[0]*n[1] + self.m[1]*n[3],
                    self.m[2]*n[0] + self.m[3]*n[2], 
                    self.m[2]*n[1] + self.m[3]*n[3]                    
                   )

    def from_euler(self, rot):

        dtr = self.mu.dtr

        self.m[0]  =  math.cos(dtr( rot ))
        self.m[1]  = -math.sin(dtr( rot ))
        self.m[2]  =  math.sin(dtr( rot ))
        self.m[3]  =  math.cos(dtr( rot ))


    def batch_mult_pts(self, pts):
        """ iterate a list of points and multiply them by this matrix """

        tmp_buffer = []
        out = None
        for pt in pts:  
            tmp_buffer.append( self * pt )
        return tmp_buffer


##-------------------------------------------##
##-------------------------------------------##

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

    def np_inverse(self, mtype='numpy'):
        """ seems to work, but I am not convinced """
        if NUMPY_IS_LOADED:        
            a = self.copy(mtype='numpy')
            b = inv(a)
            if mtype=='numpy':
                return b 
            if mtype=='m33':    
                c = matrix33()
                c.insert(b)
                return c
            #print('## inverse \n\n', self  , ' \n\n', c , ' \n\n',  self*c )

    def serialize(self, inarray):
        """ serialize this array into a list 
            if you want it as a numpy array use self.copy(mtype='numpy')
        """

        out = []
        for i in self.m:
            out.append(i)
        return out


    def from_np(self, np):
        """ convert a numpy matrix to a matrix33 object
        """
        # return type(self)(
        #     np[0][0] , np[1][0] , np[2][0] ,
        #     np[0][1] , np[1][1] , np[2][1] ,
        #     np[0][2] , np[1][2] , np[2][2]
        # )

        return type(self)(
            np[0][0] , np[0][1] , np[0][2] ,
            np[1][0] , np[1][1] , np[1][2] ,
            np[2][0] , np[2][1] , np[2][2]
        )


    def insert(self, iterable):
        """ load the first 9 things we find into this matrix 
            accepts numpy.ndarray, list, and tuple 
        """
       
        if isinstance(iterable, matrix33):
            self.m = iterable.m
        
        #serialized simple array
        if isinstance(iterable, list) or isinstance(iterable, tuple):
            for idx,i in enumerate(iterable):
                if idx <= 9:
                     self.m[idx] = iterable[idx] 

        if NUMPY_IS_LOADED:
            #numpy ND array 
            if isinstance(iterable, np.ndarray):
                out = [];idx=0
                for row in iterable:
                    for col in row:
                        self.m[idx] = col
                        idx+=1

    def copy(self, mtype=None):
        """ create a copy of this matrix - can be same type or a numpy.ndarray """
        if not mtype:
            return type(self)(
                self.m[0] , self.m[1] , self.m[2] ,
                self.m[3] , self.m[4] , self.m[5] ,
                self.m[6] , self.m[7] , self.m[8]
            )

        if NUMPY_IS_LOADED:
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
            outx = self.m[0]*n.x + self.m[3]*n.y + self.m[6]*n.z 
            outy = self.m[1]*n.x + self.m[4]*n.y + self.m[7]*n.z 
            outz = self.m[2]*n.x + self.m[5]*n.y + self.m[8]*n.z 
            return  (outx, outy, outz)

        if isinstance(n, tuple) or isinstance(n, list) or isinstance(n, np.ndarray):
            outx = self.m[0]*n[0] + self.m[3]*n[1] + self.m[6]*n[2] 
            outy = self.m[1]*n[0] + self.m[4]*n[1] + self.m[7]*n[2] 
            outz = self.m[2]*n[0] + self.m[5]*n[1] + self.m[8]*n[2] 
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

    def from_vec3(self, vec3, angle):
        """ UNTESTED build a 3X3 matrix that will rotate around a vector 
  
            vec3 = the vector for an axis of rotation
            angle - angle in degrees to rotate 

        """
        
        theta = self.mu.dtr(angle)
        #axis = vec3.normal
        axis = vec3.as_np
        #tmpm33 = self.identity

        if NUMPY_IS_LOADED and SCIPY_IS_LOADED:
            tmpm33 = expm(np.cross(np.eye(3), axis / np.linalg.norm(axis) * theta)) 
            return self.from_np(tmpm33)
        else:
            return None 


    def align_two_vec3(self, a, b):
        """  UNFINISHED !! 

             https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
       
              ⎛cosθ −sinθ  0 ⎞
          G = ⎜sinθ  cosθ  0 ⎟ 
              ⎝0      0    1 ⎠

                --------------------------------------------------------------
               Given our unit vectors, we note that cosθ=A⋅B, and sinθ=||A×B||
      
              ⎛ A.B    -||A×B||    0 ⎞
          G = ⎜||A×B||   A.B       0 ⎟ 
              ⎝0          0        1 ⎠


        """
        theta = self.mu.dtr(angle)
        #axis = vec3.normal
        axis = vec3.as_np
        #tmpm33 = self.identity
        
        tmpm33 = expm(np.cross(np.eye(3), axis / np.linalg.norm(axis) * theta)) 

        return self.from_np(tmpm33)

    def from_euler(self, xrot, yrot, zrot):
        """
            derived from the rotate_pts_3d function 
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
 
        self.insert(rotation_33)

    def rotate_pts_3d(self, points, xrot, yrot, zrot):
        """
          The "standard" 9 element, Row major, 3X3 rotation matrix used by Maya
             

           ⎡0  1  2⎤      xx xy xz 
           ⎢3  4  5⎥      yx yy yz 
           ⎣6  7  8⎦      zx zy zz 
           ------------------------------
                Rotate Y matrix     
           ⎡cos(y)  0      -sin(y) ⎤ 
           ⎢0       1       0      ⎥ 
           ⎣sin(y)  0       cos(y) ⎦ 
           ------------------------------
                Rotate Z  matrix 
           ⎡  cos(z)  sin(z)  0    ⎤ 
           ⎢ -sin(z)  cos(z)  0    ⎥
           ⎣  0       0       1    ⎦
           ------------------------------
                Rotate X matrix  
           ⎡ 1       0         0   ⎤    
           ⎢ 0    cos(x)   sin(x)  ⎥  
           ⎣ 0   -sin(x)   cos(x)  ⎦  
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


##-------------------------------------------##
##-------------------------------------------##

class matrix44(object):
    """ 4X4 matrix from pure python 
        limited support to interface to numpy arrays 

        the patern "return type(self) is nice to return copies of itself,
        but beware that this structure is not compatible for passing mutable types.
        Only primitive types work, in this case floats  


       -------------------------------------------------------------------------
       standard affine transformation matrix.

       ⎡m00  m01 m02 0⎤
       ⎢m10  m11 m12 0⎥
       ⎢m20  m21 m22 0⎥
       ⎣Tx   Ty  Tz  1⎦
       -------------------------------------------------------------------------

    """    

    def __init__(self, a=1,b=0,c=0,d=0,
                       e=0,f=1,g=0,h=0,
                       i=0,j=0,k=1,l=0,
                       m=0,n=0,o=0,p=1):

        self.mu = math_util()
        self.m = [a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p]

    def load_file(self, filename, transpose=False):
       if os.path.lexists(filename):
            f = open( filename,"r", encoding='utf-8')
            ct = 0
            for lin in f.readlines():
                c = lin.replace('\n','').split(' ')
                if ct==0:
                    self.m[0]=float(c[0])
                    self.m[1]=float(c[1])
                    self.m[2]=float(c[2])
                    self.m[3]=float(c[3])                                                        
                if ct==1:
                    self.m[4]=float(c[0])
                    self.m[5]=float(c[1])
                    self.m[6]=float(c[2])
                    self.m[7]=float(c[3])  
                if ct==2:
                    self.m[8]=float(c[0])
                    self.m[9]=float(c[1])
                    self.m[10]=float(c[2])
                    self.m[11]=float(c[3])  
                if ct==3:
                    self.m[12]=float(c[0])
                    self.m[13]=float(c[1])
                    self.m[14]=float(c[2])
                    self.m[15]=float(c[3])  

                ct += 1

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
        """multiply this 4X4 by another 4X4 matrix or a vector3, vector4 """

        if isinstance(n, vec4):
            #untested
            #https://www.scratchapixel.com/lessons/3d-basic-rendering/perspective-and-orthographic
            #    -projection-matrix/building-basic-perspective-projection-matrix

            # row major * vec4 
            # outx = n.x*self.m[0]   + n.y*self.m[1]  +  n.z*self.m[2]   +  n.w *self.m[3] 
            # outy = n.x*self.m[4]   + n.y*self.m[5]  +  n.z*self.m[6]   +  n.w *self.m[7] 
            # outz = n.x*self.m[8]   + n.y*self.m[9]  +  n.z*self.m[10]  +  n.w *self.m[11]                         
            # outw = n.x*self.m[12]  + n.y*self.m[13] +  n.z*self.m[14]  +  n.w *self.m[15] 
            
            # column major * vec4 
            outx = n.x*self.m[0]  + n.y*self.m[4]  +  n.z*self.m[8]   +  n.w *self.m[12] 
            outy = n.x*self.m[1]  + n.y*self.m[5]  +  n.z*self.m[9]   +  n.w *self.m[13] 
            outz = n.x*self.m[2]  + n.y*self.m[6]  +  n.z*self.m[10]  +  n.w *self.m[14]                         
            outw = n.x*self.m[3]  + n.y*self.m[7]  +  n.z*self.m[11]  +  n.w *self.m[15] 

            return  (outx, outy, outz, outw)

   
        if isinstance(n, vec3): # or isinstance(n, np.ndarray):
            # column major -                      why add the last 12,13,14 ? (affine?)            
            # outx = self.m[0] * n.x + self.m[4] * n.y + self.m[8]  * n.z     + self.m[12]
            # outy = self.m[1] * n.x + self.m[5] * n.y + self.m[9]  * n.z     + self.m[13]
            # outz = self.m[2] * n.x + self.m[6] * n.y + self.m[10] * n.z     + self.m[14]

            # column major  , without elements 12,13,14 
            outx = self.m[0] * n.x + self.m[4] * n.y + self.m[8]  * n.z 
            outy = self.m[1] * n.x + self.m[5] * n.y + self.m[9]  * n.z 
            outz = self.m[2] * n.x + self.m[6] * n.y + self.m[10] * n.z 
            return  (outx, outy, outz)


        if isinstance(n, tuple) or isinstance(n, list):
            
            # what is the purspose of adding 12,13,14 ?
            outx = self.m[0] * n[0] + self.m[4] * n[1] + self.m[8]  * n[2]     + self.m[12]
            outy = self.m[1] * n[0] + self.m[5] * n[1] + self.m[9]  * n[2]     + self.m[13]
            outz = self.m[2] * n[0] + self.m[6] * n[1] + self.m[10] * n[2]     + self.m[14]

            # column major, same as first, without 12,13,14               
            #outx = self.m[0] * n[0] + self.m[4] * n[1] + self.m[8]  * n[2] 
            #outy = self.m[1] * n[0] + self.m[5] * n[1] + self.m[9]  * n[2] 
            #outz = self.m[2] * n[0] + self.m[6] * n[1] + self.m[10] * n[2] 

            # row major - same as first, without 12,13,14    
            # outx = self.m[0] * n[0] + self.m[1] * n[1] + self.m[2]  * n[2] 
            # outy = self.m[4] * n[0] + self.m[5] * n[1] + self.m[6]  * n[2] 
            # outz = self.m[8] * n[0] + self.m[9] * n[1] + self.m[10] * n[2] 
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

    def np_inverse(self, mtype='numpy'):
        """ untested """
        if NUMPY_IS_LOADED:
            a = self.copy(mtype='numpy')
            b = inv(a)
            if mtype=='numpy':
                return b 
            if mtype=='m44':    
                c = matrix44()
                c.insert(b)
                return c

    #@property
    def test_index(self):
        """ fill matrix with incrementing number to see indices """
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

        if isinstance(iterable, matrix44):
            self.m = iterable.m

        #serialized simple array
        if isinstance(iterable, list) or isinstance(iterable, tuple):
            for idx,i in enumerate(iterable):
                if idx <= 16:
                     self.m[idx] = iterable[idx] 

        if NUMPY_IS_LOADED:
             #numpy ND array 
            if isinstance(iterable, np.ndarray):
                out = [];idx=0
                for row in iterable:
                    for col in row:
                        self.m[idx] = col
                        idx+=1

    #@property
    def copy(self, mtype=None):
        """ UNTESTED """
        if mtype == None:        
            return type(self)(
                self.m[0] , self.m[1] , self.m[2] , self.m[3] ,
                self.m[4] , self.m[5] , self.m[6] , self.m[7] ,
                self.m[8] , self.m[9] , self.m[10], self.m[11],
                self.m[12], self.m[13], self.m[14], self.m[15]
            )

        if NUMPY_IS_LOADED:
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
            multiply a 4X4 matrix by a list of points 
        """
        #debug - make work with other types??

        tmp_buffer = []
        out = None
        for pvec in pts:  
            #total experiment for perpective/homogeneous coordinates  
            #pvec = vec4(pt[0], pt[1], pt[2], 1)
            #print(pvec)
            tmp_buffer.append( self * pvec )
 
        return tmp_buffer

    def from_euler(self, xrot, yrot, zrot):
        """
            derived from the rotate_pts_3d function 
            go read that doc for explanation  
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
 
        self.insert(rotation_44)

    def rotate_pts_3d(self, points, xrot, yrot, zrot):
        """
           The "standard" 16 element, Row major, 4X4 rotation matrix 

           ⎡0   1   2   3  ⎤    ⎡   XVEC    0 ⎤
           ⎢4   5   6   7  ⎥    ⎢   YVEC    0 ⎥
           ⎢8   9   10  11 ⎥    ⎢   ZVEC    0 ⎥
           ⎣12  13  14  15 ⎦    ⎣ 0  0  0   0 ⎦

           ⎡0   1   2   3 ⎤     ⎡xx  xy  xz  0⎤
           ⎢4   5   6   7 ⎥     ⎢yx  yy  yz  0⎥
           ⎢8   9   10  11⎥     ⎢zx  zy  zz  0⎥
           ⎣12  13  14  15⎦     ⎣0   0   0   0⎦
           ------------------------------
           rotate Y matrix     
           ⎡  cos(y)  0      -sin(y)  0 ⎤ 
           ⎢  0       1       0       0 ⎥ 
           ⎢  sin(y)  0       cos(y)  0 ⎥ 
           ⎣  0       0       0       1 ⎦
           ------------------------------
           rotate Z  matrix 
           ⎡  cos(z)  sin(z)  0       0 ⎤ 
           ⎢ -sin(z)  cos(z)  0       0 ⎥
           ⎢  0       0       1       0 ⎥
           ⎣  0       0       0       1 ⎦
           ------------------------------
           rotate X matrix  
           ⎡  1       0       0       0 ⎤  
           ⎢  0       cos(x)  sin(x)  0 ⎥  
           ⎢  0      -sin(x)  cos(x)  0 ⎥  
           ⎣  0       0       0       1 ⎦ 
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
            
            transformation matrix that changes the W element of each vertex. 
            After the the camera matrix is applied to each vertex, 
            but before the projection matrix is applied, 
            the Z element of each vertex represents the distance away from the camera. 
            the larger Z is, the more the vertex should be scaled down


            http://stackoverflow.com/questions/8633034/basic-render-3d-perspective-projection-onto-2d-screen-with-camera-without-openg

            use homogenous transformations and coordinates. 

            take a point in space and:
                Position it relative to the camera using the model matrix.
                Project it either orthographically or in perspective using the projection matrix.
                Apply the viewport transformation to place it on the screen.




        """

        PI_OVER_360 = 0.00872664625
        #xymax = 100
        xymax = znear * math.tan(fov * PI_OVER_360);

        ymin = -xymax
        xmin = -xymax

        width = xymax - xmin
        height = xymax - ymin

        depth =  zfar - znear
        
        #avoid those pesky divide by zero errors 
        if depth==0:
            depth=.1
        if width==0:
            width=1
        if height==0:
            height=1                        

        q     = -(zfar + znear) / depth
        qn    = -2 * (zfar * znear) / depth

        w = 2 * znear / width
        w = w / aspect
        h = 2 * znear / height

        m = self.identity

        # m[0]  = w; m[4]  = 0; m[8]  = 0 ; m[12] = 0
        # m[1]  = 0; m[5]  = h; m[9]  = 0 ; m[13] = 0
        # m[2]  = 0; m[6]  = 0; m[10] = q ; m[14] = qn 
        # m[3]  = 0; m[7]  = 0; m[11] =-1 ; m[15] = 0

        m[0]=w ; m[1]=0  ; m[2]=0  ; m[3]=0;
        m[4]=0 ; m[5]=h  ; m[6]=0  ; m[7]=0;
        m[8]=0 ; m[9]=0  ; m[10]=q ; m[11]=-1;
        m[12]=0; m[13]=0 ; m[14]=qn; m[15]=0;
        
        #print(m)
        return m

##-------------------------------------------##
##-------------------------------------------##

class quaternion(object):
    """ UNTESTED
        first stab at quaternion 
        mostly converted from C - 
        taken from this: https://github.com/mycmessia/3D-Math-Primer/blob/master/3D%20Math/Quaternion.cpp 

    """

    def __init__(self,w=1,x=0,y=0,z=0):
        self.w=w 
        self.x=x
        self.y=y
        self.z=z

    #def inverse(self)
    #def difference(self)
    #def to_euler(self)

    def __repr__(self):
        return '(%s, %s, %s, %s)' % (self.w, self.x, self.y, self.z)

    def __getitem__(self, index):
        if index==0:
            return self.w
        if index==1:
            return self.x
        if index==2:
            return self.y
        if index==3:
            return self.z

    def __setitem__(self, key, item):
        if key==0:
            self.w = item
        if key==1:
            self.x = item
        if key==2:
            self.y = item
        if key==3:
            self.z = item

    @property
    def identity(self):
        return type(self)(1,0,0,0)

    def set_identity(self):
        self.w=1 
        self.x=0
        self.y=0
        self.z=0

    def set_rotx (self, theta):
        theta_over2 = theta * .5
        self.w = math.cos(theta_over2);
        self.x = math.sin(theta_over2);
        #self.y = 0
        #self.z = 0

    def set_roty (self, theta):
        theta_over2 = theta * .5
        self.w = math.cos(theta_over2)
        #self.x = 0
        self.y = math.sin(theta_over2)
        #self.z = 0

    def set_rotz (self, theta):
        theta_over2 = theta * .5
        self.w = math.cos(theta_over2)
        #self.x = 0
        #self.y = 0
        self.z = math.sin(theta_over2)

    def from_euler(self, h, p, b, trans_type='obj2inertial'):
        sp=0;sb=0;sh=0
        cp=0;cb=0;ch=0
        
        sp = math.sin(p*.5) 
        sb = math.sin(b*.5) 
        sh = math.sin(h*.5) 

        cp = math.cos(p*.5) 
        cb = math.cos(b*.5) 
        ch = math.cos(h*.5) 
       
        if (trans_type == 'obj2inertial'):
            self.w = ch * cp * cb + sh * sp * sb
            self.x = ch * sp * cb + sh * cp * sb
            self.y = -ch * sp * sb + sh * cp * cb
            self.z = -sh * sp * cb + ch * cp * sb

        elif (trans_type == 'inertial2ob'):
            self.w = ch * cp * cb + sh * sp * sb
            self.x = -ch * sp * cb - sh * cp * sb
            self.y = ch * sp * sb - sh * cp * cb
            self.z = sh * sp * cb - ch * cp * sb
        
        else:
            print( "Invalid trans_type!" ) 

    def mag(self):
        mag = float( math.sqrt(  self.w*self.w + 
                                 self.x*self.x + 
                                 self.y*self.y + 
                                 self.z*self.z) ) 
        return mag 

    def normalize(self):
        mag = self.mag() 
        if mag > 0:
            oneOverMag = float( 1.0 / mag)
            self.w *= oneOverMag
            self.x *= oneOverMag
            self.y *= oneOverMag
            self.z *= oneOverMag
        else:
            self.set_identity()

    def dot_product(self, q): 
        #DEBUG - HOW IS THIS RIGHT? 
        return a.x * b.x + a.y * b.y + a.z * b.z + a.z * a.z;   

    def conjugate(self, q):
        result = type(self)()
        
        result.w =  q.w
        result.x = -q.x
        result.y = -q.y
        result.z = -q.z
        return result


    def __mul__(self, a):
        result = type(self)()
    
        w = self.w
        x = self.x
        y = self.y
        z = self.z

        result.w = w * a.w - x * a.x - y * a.y - z * a.z
        result.x = w * a.x + x * a.w + z * a.y + y * a.z
        result.y = w * a.y + y * a.w + x * a.z + z * a.x
        result.z = w * a.z + z * a.w + y * a.x + x * a.y
        
        return result;


    def from_m33(self, m33):
        """ 
           m11  m12 m13 
           m21  m22 m23 
           m31  m32 m33 
        """
        
        four_w_sq_min1 = m33[0] + m33[4] + m33[8]
        four_x_sq_min1 = m33[0] - m33[4] - m33[8]
        four_y_sq_min1 = m33[4] - m33[0] - m33[8]
        four_z_sq_min1 = m33[8] - m33[0] - m33[4]
        
        maxIndex = 0
        max = four_w_sq_min1
        
        if four_x_sq_min1 > max:
            max = four_x_sq_min1
            maxIndex = 1
        
        if (four_y_sq_min1 > max):
            max = four_y_sq_min1
            maxIndex = 2
        
        if (four_z_sq_min1 > max):
            max = four_z_sq_min1
            maxIndex = 3
        
        max = math.sqrt (max + 1.0) * 0.5
        mult = 0.25 / max
        
        if maxIndex==0:
            self.w = max;
            self.x = (m33[5] - m33[7]) * mult;
            self.y = (m33[6] - m33[2]) * mult;
            self.z = (m33[1] - m33[3]) * mult;

        if maxIndex==1:
            self.x = max;
            self.w = (m33[5] - m33[7]) * mult;
            self.y = (m33[1] + m33[3]) * mult;
            self.z = (m33[6] + m33[2]) * mult;

        if maxIndex==2:
            self.y = max;
            self.w = (m33[6] - m33[2]) * mult;
            self.x = (m33[1] + m33[3]) * mult;
            self.z = (m33[5] + m33[7]) * mult;

        if maxIndex==3:
            self.z = max;
            self.w = (m33[1] - m33[3]) * mult;
            self.x = (m33[6] + m33[2]) * mult;
            self.y = (m33[5] + m33[7]) * mult;
            self.z = (m33[5] + m33[7]) * mult;
 


    def to_m33(self, trans_type='inertial2obj'):
 
            mo = matrix33() 
            q = self 

            if (trans_type == 'inertial2obj'):
 
                mo[0] = 1.0 - 2.0 * (q.y * q.y + q.z * q.z)
                mo[1] = 2.0 * (q.x * q.y + q.w * q.z)
                mo[2] = 2.0 * (q.x * q.z + q.w * q.y)
                
                mo[3] = 2.0 * (q.x * q.y - q.w * q.z)
                mo[4] = 1.0 - 2.0 * (q.x * q.x + q.z * q.z)
                mo[5] = 2.0 * (q.y * q.z + q.w * q.x)
                
                mo[6] = 2.0 * (q.x * q.z + q.w * q.y)
                mo[7] = 2.0 * (q.y * q.z - q.w * q.x)
                mo[8] = 1.0 - 2.0 * (q.x * q.x + q.y * q.y)
 
                return mo 

            elif (trans_type == 'obj2inertial'):
 
                mo[0] = 1.0 - 2.0 * (q.y * q.y + q.z * q.z)
                mo[1] = 2.0 * (q.x * q.y - q.w * q.z)
                mo[2] = 2.0 * (q.x * q.z + q.w * q.y)
                
                mo[3] = 2.0 * (q.x * q.y + q.w * q.z)
                mo[4] = 1.0 - 2.0 * (q.x * q.x + q.z * q.z)
                mo[5] = 2.0 * (q.y * q.z - q.w * q.x)

                mo[6] = 2.0 * (q.x * q.z + q.w * q.y)
                mo[7] = 2.0 * (q.y * q.z + q.w * q.x)
                mo[8] = 1.0 - 2.0 * (q.x * q.x + q.y * q.y)
 
                return mo 

            else:
                print("Invalid trans_type!")
   
    def set_rot_zxis(self, axis, theta):
        #assert((vectorMag(axis) - 1.0f) < 0.01f);
        thetaOver2 = theta * .5
        sinThetaOver2 = math.sin(thetaOver2)
       
        w = math.cos(thetaOver2)
        x = axis.x * sinThetaOver2
        y = axis.y * sinThetaOver2
        z = axis.z * sinThetaOver2

    def get_rot_angle(self):
        thetaOver2 = math.acos(self.w)
        return thetaOver2 * 2.0
 
    def get_rot_axis(self):
        sin_theta_over2Sq = 1.0 - self.w * self.w
        one_over_sin_theta = 1.0 / math.sqrt(sin_theta_over2Sq)
        
        nx = self.x * one_over_sin_theta
        ny = self.y * one_over_sin_theta
        nz = self.z * one_over_sin_theta
        
        return vec3(nx, ny, nz)
 

    def slerp (self, q0, q1, t):
 
        if (t <= 0):
            return q0
        if (t >= 1):
            return q1

        cosOmega = self.dot_product(q0, q1)
        
        q1w = q1.w
        q1x = q1.x
        q1y = q1.y
        q1z = q1.z
        
        # avoid getting different results
        if (cosOmega < 0):
            q1w = -q1w
            q1x = -q1x
            q1y = -q1y
            q1z = -q1z
            cosOmega -= cosOmega
        
        k0 = 0
        k1 = 0
        # avoid sth over 0 happening
        if cosOmega > 0.999:
            k0 = 1.0 - t
            k1 = t
        else:
        
            sinOmega = math.sqrt(1.0 - cosOmega * cosOmega)
            omega = math.atan2 (sinOmega, cosOmega)
            oneOverSinOmega = 1.0/sinOmega
            
            k0 = math.sin ((1.0- t) * omega) * oneOverSinOmega
            k1 = math.sin (t*omega) * oneOverSinOmega
        
        ####
        result = type(self)()
        
        result.x = k0 * q0.x + k1 * q1x;
        result.y = k0 * q0.y + k1 * q1y;
        result.z = k0 * q0.z + k1 * q1z;
        result.w = k0 * q0.w + k1 * q1w;
        
        return result


##-------------------------------------------##
##-------------------------------------------##

class spherical(object):
    """ UNTESTED -   polar and spherical coordinates 

    # https://stackoverflow.com/questions/4116658/faster-numpy-cartesian-to-spherical-coordinate-conversion

    def cart2sph(x,y,z):
        XsqPlusYsq = x**2 + y**2
        r = m.sqrt(XsqPlusYsq + z**2)               # r
        elev = m.atan2(z,m.sqrt(XsqPlusYsq))     # theta
        az = m.atan2(y,x)                           # phi
        return r, elev, az

    def cart2sphA(pts):
        return np.array([cart2sph(x,y,z) for x,y,z in pts])

    def appendSpherical(xyz):
        np.hstack((xyz, cart2sphA(xyz)))

    """

    def __init__(self,r=0,t=0,p=0):
        """ r = radius 
            t = theta 
            p = phi  (not needed for polar)
        """
        #self.mu = math_util() 
        self.r=r
        self.t=t
        self.p=p

    def __repr__(self):
        return '(%s, %s, %s)' % (self.r, self.p, self.t)

    def __getitem__(self, index):
        if index==0:
            return self.r
        if index==1:
            return self.t
        if index==2:
            return self.p

    def __setitem__(self, key, item):
        if key==0:
            self.r = item
        if key==1:
            self.t = item
        if key==2:
            self.p = item

    def to_cartesian(self, deg_rad='rad'):
        """ from David Gould's book
            Complete Maya Programing II 
            page 14 
        """
        x = self.r * math.sin(self.p) * math.cos(self.t)
        y = self.r * math.sin(self.p) * math.sin(self.t)
        z = self.r * math.cos(self.p)    
        return (x,y,z)
 
    def from_cartesian(self, vec3):
        """ UNTESTED 
            from David Gould's book
            Complete Maya Programing II 
            page 14 

            r = length( x y z )
            p = tan-1 (length(x y), z)) 
            t = tan-1 (y,x)
        """        
        if isinstance(vec3,tuple):
            vec3 = vec3(vec3[0],vec3[1],vec3[2]) 

        r = vec3.length  
        p = math.atan( -1*math.sqrt(vec3[0]*vec3[0]+vec3[1]*vec3[1])*vec3[2]  ) 
        t = math.atan( -1*(vec3[1]*vec3[0]))
        print('### r %s p %s t %s '%(r,p,t))
        self.r=r
        self.t=t
        self.p=p
        #return type(self)(r,p,t)

    def polar_to_cartesian(self):
        """ UNTESTED
            polar coordinates are spherical without the phi element 
            from David Gould's book
            Complete Maya Programing II 
            page 13         
        """ 
        pass

    def cartesian_to_polar(self, vec3):
        """ UNTESTED
            polar coordinates are spherical without the phi element 
            from David Gould's book
            Complete Maya Programing II 
            page 13         
        """ 
        if isinstance(vec3,tuple):
            vec3 = vec3(vec3[0],vec3[1],0) #Z is ignored  

        r = math.sqrt(vec3[0]*vec3[0] + vec3[1]*vec3[1])  
        t = math.atan( vec3[1]*vec3[0] ) 

        print('### cartesian to polar  r %s t %s '%(r,t))
        self.r=r
        self.t=t
        self.p=0


