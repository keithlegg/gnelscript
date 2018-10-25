


import os, sys, math

from pycv.point_ops import PointOperator




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


