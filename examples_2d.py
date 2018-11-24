

from pygfx.obj3d import  *
from pygfx.point_ops import *

from pygfx.raster_ops import *
from pygfx.math_ops import  vec2





###################################################
def render_2d_vector(v1, gridsize=50):
    """draw a vector and tell us the angle of it in degrees
       
       vector    : 2 2D tuples e.g. (1.5,1), (0,0)
       gridsize  : specify pixels per linear unit

    """

    fb = PixelOp()   
    fb.create_buffer(800, 800)
    fb.graticule(gridsize)
    fb.render_vector_2d(   v1,  scale=gridsize)
    fb.save('vec.png')


# render_2d_vector( (0,5)  ) 

###################################################



###################################################

def bloody_simple_2drender( imagename, pts=None, vecs=None, lines=None, obj=None, gridsize=50):
    """ draw some points and lines on a grid 

        ARGS:

            pts      - list of points to render 
            vecs     - list of vectors to render 
            lines    - list of lines to render 
            obj      - list of 3D obect models to render
            gridsize - parameter to set a grid units to pixels ratio

        TO RUN:

            pts = [(1,1),(2,2),(3,3)]
            lines = [ [ (1,1), (1,2), (2,1)], [ (6,1), (1,6), (5,-1)] ]
            bloody_simple_2drender('2d_render.png', pts=pts, vecs=pts, lines=lines )

    """

    fb = PixelOp()   
    fb.create_buffer(800, 800)
    fb.graticule(gridsize)

    #pt_size = 3
    pointcolor = ()
    linecolot = () 

    if obj is not None:
        for o in obj:
            # look! I wrote a renderer in 5 lines!! 
            for ply in o.polygons:
                poly = []
                for fid in ply:
                    poly.append( o.points[fid-1] )
                fb.render_line_2d( poly,  scale=gridsize)   

    if lines is not None:
        for line in lines:        
            fb.render_line_2d( line,  scale=gridsize)

    if vecs is not None:
        for v in vecs:
            fb.render_vector_2d(   v,  scale=gridsize)

    if pts is not None:
        for p in pts:
            fb.render_point_2d( p , scale=gridsize )

    fb.save(imagename)



## ---------------------------------------


def example_BSR():
    """ BSR = bloody simple renderer """
    pts = [(1,1),(2,2),(3,3)]
    lines = [ [ (1,1), (1,2), (2,1)], [ (6,1), (1,6), (5,-1)] ]

    bloody_simple_2drender('2d_render.png', pts=pts, vecs=pts, lines=lines )



## ---------------------------------------


def load_obj_render_BSR(objfile):
    """ BSR = bloody simple renderer """
    #load a 3d model and render it in 2D
    obj = object3d() 
    obj.load(objfile)
    
    obj2 = object3d() 
    obj2.load('objects/monkey.obj')

    bloody_simple_2drender('2d_render.png', obj=[obj,obj2], gridsize=100)


load_obj_render_BSR('original_sin.obj') 



###################################################
###################################################




def distance_between_2vecs():

    a = vec2(1,1)
    b = vec2(2,2)

    c = a-b 
  
    print( '## this number %s ' % a.distance_to(b) ) 
    print( '.. should be equal to this %s '% c.length ) 


    vecs = [ a, b, c]

    # pts = [(1,1),(2,2),(3,3)]
    bloody_simple_2drender('2d_render.png', pts=None, vecs=vecs)
 

# distance_between_2vecs() 

###################################################

"""
a = vec2(1,1)
b = vec2(2,2)

print( a.project_pt(a, b, 0) ) 

"""


def test_2d_intersect():
    a = vec2()

    s1 = vec2( 5, 5)
    e1 = vec2(-5,-5)

    s2 = vec2( 6, 3)
    e2 = vec2(-3,-3)


    vecs = [ a, b, c]
    # pts = [(1,1),(2,2),(3,3)]
    bloody_simple_2drender('2d_render.png', pts=None, vecs=vecs)

    print( a.intersect(s1,e1,s2,e2) )





###################################################


def project_point_along_2Dvector():
    
    # 2d vector 
    a = vec2(  1  ,1   )
    b = vec2( 1.01 ,1.01 )
    com = vec2() 
    print( com.project_pt(a, b, 1) )





