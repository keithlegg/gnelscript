
""" 
DEBUG NOT DONE - JUST COPIED FROM OBJECT3D
"""

import math 


from gnelscript.pygfx.point_ops_2d import *
from gnelscript.pygfx.obj2d import *



class object2d(polygon_operator_2d):
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

        self.rot             = [0,0]
        self.pos             = [0,0]
        self.scale           = [1,1]

    def reset(self):
        self.rot          = [0,0]
        self.pos          = [0,0]
        self.scale        = [1,1]

    ##-------------------------------------------## 
    def show(self):
        data = []
 
        lines = 0
        tris = 0
        quads = 0
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
        data.append('  position      : %s %s ' %( self.pos[0], self.pos[1]     ))
        data.append('  rotation      : %s %s ' %( self.rot[0], self.rot[1]     ))
        data.append('  scale         : %s %s ' %( self.scale[0], self.scale[1] ))
        data.append(' --------------------------- ' )        
        #data.append('  num face normals   : %s' %  self.numfacnrml )   
        #data.append('  num verts          : %s' %  self.numpts     )
        #data.append('  num polygons       : %s' %  self.numply     )
        data.append(' --------------------------- ' )         
        data.append('  num triangles      : %s' %  tris  )
        data.append('  num quads          : %s' %  quads )  
        data.append('  num lines          : %s' %  lines ) 
        data.append('  num other          : %s' %  other ) 
        data.append('############################\n')

        for d in data:
            print(d) 

    def prim_square(self,  pos=(0,0,0), rot=(0,0,0), size=1):
        pts = [] 

        #keep Z but leave it zero
        pts.append((-size, -size))  
        pts.append((-size,  size))  
        pts.append(( size,  size)) 
        pts.append(( size, -size)) 

        #not a good solution - this doesnt take into account existing geom 
        self.points.extend(pts)
        self.polygons.append( (1,2,3,4) )


    def prim_triangle(self, pos=(0,0), rot=(0,0), size=1):
        pts =  [(-size,0,0), (0,size,0), (size,0,0) ]
        poly = [(1,2,3)]

        #not a good solution - this doesnt take into account existing geom 
        self.points.extend(pts)
        self.polygons.append( (1,2,3) )


    # def prim_circle(self,  center=(0,0,0), dia=1):
    #     pass






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


##-------------------------------------------##

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

##-------------------------------------------##

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