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
# ***** END GPL LICENCE BLOCK *****


""" 
DEBUG NOT DONE - JUST COPIED FROM OBJECT3D
"""


import math 


from gnelscript.pygfx.point_ops_2d import polygon_operator_2d
from gnelscript.pygfx.math_ops import vec3 



class object25d(polygon_operator_2d):
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


    ##-------------------------------------------## 
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

    ##-------------------------------------------##
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

    ##-------------------------------------------##
    def prim_triangle(self, pos=(0,0,0), rot=(0,0,0), size=1):
        pts =  [(-size,0,0), (0,size,0), (size,0,0) ]
        poly = [(1,2,3)]

        #not a good solution - this doesnt take into account existing geom 
        self.points.extend(pts)
        self.polygons.append( (1,2,3) )


    # def prim_circle(self,  center=(0,0,0), dia=1):
    #     pass


    ##-------------------------------------------##  

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