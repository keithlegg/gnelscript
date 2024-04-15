
""" 
general tools for grids, coordinates and a little GIS even 

"""


import os, sys, math



from gnelscript.pygfx.point_ops import *
from gnelscript.pygfx.obj3d import *

from gnelscript.pygfx.point_ops_2d import *
from gnelscript.pygfx.dag_ops import *

#not done 
#from gnelscript.pygfx.obj2d import *


class cell(node_base):
    """ what is idx VS coord ?? 
        IDX (index) is an index to help locate the cell on a 3d grid and optionally name the cell 
        the coord is the location in space , which has no need to coresepond to the index itself, but may if its a cubic structure 
        
    """

    def __init__(self, name,
                       width, height, 
                       idx_x , idx_y , idx_z, 
                       coordx, coordy, coordz 
                ):
        super().__init__()     
        self.name = name

        self.width = width
        self.height = height

        self.coord_x = coordx
        self.coord_y = coordy
        self.coord_z = coordz

        self.cell_x = idx_x
        self.cell_y = idx_y
        self.cell_z = idx_z

        self.points = [] 
        self.polys  = [] 

    # @property
    # def bbox(self):

    # @property
    # def centroid(self):

    @property
    def boundary_pts(self):
        pts = []
        halfwidth = self.width/2 
        halfheight = self.height/2 

        pts.append( (self.coord_x-halfwidth , self.coord_y-halfheight) )  
        pts.append( (self.coord_x-halfwidth , self.coord_y+halfheight) ) 
        pts.append( (self.coord_x+halfwidth , self.coord_y+halfheight) ) 
        pts.append( (self.coord_x+halfwidth , self.coord_y-halfheight) ) 
        #make periodic
        pts.append( (self.coord_x-halfwidth , self.coord_y-halfheight) ) 

        shifted=[]
        for pt in pts:
            shifted.append( (pt[0]+halfwidth, pt[1]+halfheight) ) 
        return shifted

##------------------------------##

class tessellator(data_graph):

    def __init__(self):
        super().__init__()         

        #dimensions of 2d grid
        self.width  = 0
        self.height = 0
        self.depth = 0

        #number of cells in 2d grid  
        self.x_squares = 0
        self.y_squares = 0
        self.z_squares = 0

        self.minx = 0
        self.miny = 0
        self.maxx = 0        
        self.maxy = 0
    
    def show(self):
        for c in self.nodes:
            print('DEBUG type cell ', type(c))

            print("# name %s %s %s "%(c.name, c.coord_x+(c.width/2), c.coord_y+(c.height/2) ) )

    ##----------- 
    def _set_extents(self, bbox):
        """ set global extents for generating data 
            based on PIL coordinate which is [left, top, right, bottom] 
            [minx, miny, maxx, maxy]  
        """
        
        self.minx = bbox[0]
        self.miny = bbox[1]         
        self.maxx = bbox[2]
        self.maxy = bbox[3]

    ##----------- 

    def new_cell_2d(self, name, 
                          width, height,
                          idx_x, idx_y, idx_z, 
                          coordx, coordy, coordz,
                          attrs=None):

        c = cell(name,
                 width, height,
                 idx_x, idx_y, idx_z, 
                 coordx, coordy, 0)

        c.set_position([coordx, coordy, coordz])
        #c.set_rotation ([, , ])
        c.addattr('centroid', [coordx,coordy])
        c.addattr('idx_x', idx_x)
        c.addattr('idx_y', idx_y)
        c.addattr('idx_z', idx_z)
        c.addattr('width', width)
        c.addattr('height', height)
    
        ##calc the edges (midpoints of edges) for a sqaure polygon 
        o = object3d() 
        o.prim_rect(pos=(coordx,coordy,coordz), sizex=width, sizey=height)
        s=o.points 
         
        ## DEBUG midpoints of opposing edges can also derive an angle per cell  
        e1mid = o.locate_pt_along3d(s[0], s[1], 1)
        e2mid = o.locate_pt_along3d(s[1], s[2], 1)
        e3mid = o.locate_pt_along3d(s[2], s[3], 1)
        e4mid = o.locate_pt_along3d(s[3], s[0], 1)
        c.addattr('e1mid', e1mid)
        c.addattr('e2mid', e2mid)
        c.addattr('e3mid', e3mid)
        c.addattr('e4mid', e4mid)

        self.add(c)

    ##----------- 
    def build_2d_cells(self, numx, numy, scale=1):
        """   
        top edge verts match bottom edge
        right edge verts match left edge 
        
        -------------           
        | 0_0 | 0_1 |
        -------------
        | 1_0 | 1_1 |
        ------------- 
        """
   
        self.x_squares = numx
        self.y_squares = numy

        res_x = self.maxx - self.minx
        res_y = self.maxy - self.miny

        res_x = res_x*scale
        res_y = res_y*scale

        spacing_x = res_x/numx 
        spacing_y = res_y/numy

        print( "\n\ntesselation size: resx %s resy %s space x %s space y %s "% (res_x, res_y, spacing_x, spacing_y ) )

        ##---------------------
        for idx_x in range(self.x_squares):
            for idx_y in range(self.y_squares):
                coordx = self.minx+spacing_x*idx_x
                coordy = self.miny+spacing_y*idx_y

                self.new_cell_2d( "cell_%s_%s"%(idx_x,idx_y), 
                                  spacing_x,
                                  spacing_y,
                                  idx_x, idx_y, 0,
                                  coordx, coordy, 0
                                )


    ##----------- 
    def build_3d_cells(self):
        """   
            DEBUG - NOT DONE 
        """
        for x in range(self.x_squares):
            for y in range(self.y_squares):
                for z in range(self.z_squares):
                    print("%s %s %s"%(x,y,z))

    ##----------- 
    def from_square_outline(self, array=None):
        """   
            DEBUG - NOT DONE  - builds a grid from X,Y outline 

            input  - [a,b,c,d], [e,f,g,h], [i,j,k,l], [m,n,o,p] 
            
            output -  a b c d
                    p . . . . e
                    o . . . . f
                    n . . . . g
                    m . . . . h
                      l k j i 

            started out as a simple cube generator 
            I got tyhe idea to build a grid based on square nested arrays from another function 

            ARGS - 
                array (must be square) grid of 3d points  

        """
        output = [] 

        # DEBUG - add a check here 
        is_square=True 

        ##tested with a 2D grid of 3D points (a polygon face) 
        #if is_square:
        #    for x in array:
        #        for y in x_cells:
        #            print(xy_cells) 

        width = len(array) 
        height = len(array[0]) 

        # matrix transpose 
        tmp  = array[0] 
        tmp2 = array[1] 
        array[0] = tmp2
        array[1] = tmp
        array[0].reverse()
        array[1].reverse()

        #only works with 4 sided  (N per side)
        for y in range(height):
            for x in range(width):
                output.append( array[x][y] ) 

        return output 
        
    ##----------- 
    def arrange_cell(self):
        """ """
        pass 

    ##----------- 
    def load_tesselation(self, filename):
        #DEBUG - need to make a proper wrapper for cells VS DAG nodes 

        #for c in self.nodes:
        #    print("# name %s %s %s "%(c.name, c.coord_x+(c.width/2), c.coord_y+(c.height/2) ) ) 

        self.load_graph_file(filename) 
        
        for c in self.nodes:
            print(c)

    ##----------- 
    def save_tesselation(self, filename):
        #DEBUG - need to make a proper wrapper for cells VS DAG nodes

        #for c in self.nodes:
        #    print("# name %s %s %s "%(c.name, c.coord_x+(c.width/2), c.coord_y+(c.height/2) ) ) 

        self.save_graph_file(filename)           





