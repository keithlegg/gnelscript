


#from pygfx.point_ops import *

from pygfx.math_ops import  vec2



"""
a = vec2(2,5)
b = vec2(1,4)

print( a.distance_to(b) ) 
"""

#################


"""
a = vec2(1,1)
b = vec2(2,2)

print( a.project_pt(a, b, 0) ) 

"""


#################
a = vec2()

s1 = vec2( 5, 5)
e1 = vec2(-5,-5)

s2 = vec2( 6, 3)
e2 = vec2(-3,-3)

print( a.intersect(s1,e1,s2,e2) )

