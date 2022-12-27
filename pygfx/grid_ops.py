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




#class teselator_cell(object2d):
class cell(node_base):
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


class teselator(data_graph):

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
            print("# name %s %s %s "%(c.name, c.coord_x+(c.width/2), c.coord_y+(c.height/2) ) )

    def _set_extents(self, bbox):
        """ set global extents for generating data 
            based on PIL coordinate which is [left, top, right, bottom] 
            [minx, miny, maxx, maxy]  
        """
        
        self.minx = bbox[0]
        self.miny = bbox[1]         
        self.maxx = bbox[2]
        self.maxy = bbox[3]

    def build_3d_cells(self):
        """   
        DEBUG - NOT DONE 
        """

        for x in range(self.x_squares):
            for y in range(self.y_squares):
                for z in range(self.z_squares):
                    print("%s %s %s"%(x,y,z))


    def new_cell_2d(self, name, 
                          width, height,
                          idx_x, idx_y, idx_z, 
                          coordx, coordy, coordz):

        c = cell(name,
                 width, height,
                 idx_x, idx_y, idx_z, 
                 coordx, coordy, 0)

        self.add(c)


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
        
        print('#######################')
        print( "resx %s resy %s spacex %s spacey %s "% (res_x, res_y, spacing_x, spacing_y ) )

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




    def arrange_cell(self):
        """ """
        pass 




#############################################################

def arc_to_degree(NS, degrees, minutes, seconds, EW):
    """ arc minutes to decimal degree ( example n50d0'02"e ) """
    
    outdegrees = 0.0

    if NS =='n':
        outdegrees = degrees
        outdegrees = outdegrees + (minutes*.0166667) #1/60
        outdegrees = outdegrees + (seconds*.0166667*.0166667) #1/60
    if NS =='s':
        outdegrees = 180.0
        outdegrees = outdegrees + degrees
        outdegrees = outdegrees + (minutes*.0166667) #1/60
        outdegrees = outdegrees + (seconds*.0166667*.0166667) #1/60
    if EW =='w' and NS =='s':
        outdegrees = outdegrees * -1
    if EW =='e' and NS =='n':
        outdegrees = outdegrees * -1
  
    return outdegrees


