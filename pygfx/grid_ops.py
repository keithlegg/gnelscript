""" 
general tools for grids, coordinates and a little GIS even 

"""


import os, sys, math


from gnelscript.pygfx.point_ops import *
from gnelscript.pygfx.obj3d import *

from gnelscript.pygfx.point_ops_2d import *


#not done 
#from gnelscript.pygfx.obj2d import *




#class teselator_cell(object2d):
class teselator_cell(object):
    def __init__(self):
        self.points = [] 





class teselator(object):

    def __init__(self):
        self.cells = []
        self.width  = 10
        self.height = 10
        self.num_per_edge = 2


    def define_cell(self):
        """   

        top edge verts match bottom edge
        right edge verts match left edge 
        
        ----------------           
        | a1 | a2 | a3 |
        ---------------|
        | b1 | b2 | b3 | 
        ---------------| 
        | c1 | c2 | c3 | 
        ----------------         
        """
        pass 

    def arrange_cell(self):
        """ """
        pass 




# def import WKT()


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


