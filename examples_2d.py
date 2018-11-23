


#from pygfx.point_ops import *

from pygfx.raster_ops import *

from pygfx.math_ops import  vec2



"""
a = vec2(2,5)
b = vec2(1,4)

print( a.distance_to(b) ) 
"""

###################################################

"""
a = vec2(1,1)
b = vec2(2,2)

print( a.project_pt(a, b, 0) ) 

"""


###################################################

a = vec2()

s1 = vec2( 5, 5)
e1 = vec2(-5,-5)

s2 = vec2( 6, 3)
e2 = vec2(-3,-3)

print( a.intersect(s1,e1,s2,e2) )


###################################################

def bloody_simple_2drender( imagename, pts=None, vecs=None, gridsize=50):
    """ draw some points and lines on a grid """

    fb = PixelOp()   
    fb.create_buffer(800, 800)
    fb.graticule(gridsize)

    #pt_size = 3

    #if line_geom is not None:
    # connect the dots , scaled 

    if vecs is not None:
        for v in vecs:
            fb.draw_vector_2d(   v,  scale=gridsize)

    if pts is not None:
        for p in pts:
            fb.draw_point_2d( p , scale=gridsize )

    fb.save(imagename)




# pts = [(1,1),(2,2),(3,3)]
# bloody_simple_2drender('2d_render.png', pts=pts, vecs=pts)


###################################################


def project_point_along_2Dvector():
    
    # 2d vector 
    a = vec2(  1  ,1   )
    b = vec2( 1.01 ,1.01 )
    com = vec2() 
    print( com.project_pt(a, b, 1) )



###################################################
def render_2d_vector(v1, gridsize=50):
    """draw a vector and tell us the angle of it in degrees
       
       vector    : 2 2D tuples e.g. (1.5,1), (0,0)
       gridsize  : specify pixels per linear unit

    """

    fb = PixelOp()   
    fb.create_buffer(800, 800)
    fb.graticule(gridsize)
    fb.draw_vector_2d(   v1,  scale=gridsize)
    fb.save('vec.png')


# render_2d_vector( (0,5)  ) 

