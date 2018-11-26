

from pygfx.obj3d import  *
from pygfx.point_ops import *

from pygfx.raster_ops import *
from pygfx.math_ops import  *

mu = math_util() 



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

def bloody_simple_2drender( imagename, pts=None, vecs=None, lines=None, obj=None, gridsize=50, fb=None):
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

    if fb is None:
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
            fb.render_point_2d( p , scale=gridsize, dotsize=10)

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
    
    #obj.scale_pts( (.1,.2,.1) )
    #obj.rotate_pts( (.1,.1,.1) )

    #obj2 = object3d() 
    #obj2.load('objects/monkey.obj')

    bloody_simple_2drender('2d_render.png', obj=[obj], gridsize=100)


#load_obj_render_BSR('original_sin.obj') 



###################################################

"""
def partial_bresenham( x1, y1, x2, y2):
    # only works when the first line is numerically lower than the second 
    # I suppose you can presort them and get around this problem ?
    dy = y2-y1
    dx = x2-x1
    d = 2*dy - dx
    x = x1
    y = y1
    
    print(' distances x %s y %s  '%(dx,dy) )
    print(' d %s start x %s y %s '%(d,x,y) )
 
    while x <= x2:
        print( 'Draw pixel at (%s,%s) \n'%(x,y))

        x+=1
        if d<0 :
            d += dy + dy
            # print('D < 0 IN LOOP %s\n'% d )

        else:
            d += 2*(dy-dx)
            # print('D !<0 IN LOOP %s\n'% d )

            y+=1
"""

        


def animate_bresenham( x1, y1, x2, y2):
    """ 
        make an animation of the bresenham algorithm 

        TO RUN:
            animate_bresenham(-3,3,5,6)

    """

    dy = y2-y1
    dx = x2-x1
    d = 2*dy - dx
    x = x1
    y = y1
    
    pixels_per_unit = 50 

    #print(' distances x %s y %s  '%(dx,dy) )
    #print(' d %s start x %s y %s '%(d,x,y) )
 
    fb = PixelOp()   
    fb.create_buffer(800, 800)
    fb.graticule(pixels_per_unit)

    fr_cnt = 0 
    pt_bufr = [] 
    while x <= x2:
   
        
        bloody_simple_2drender('unit_circle_%s.png'%fr_cnt, pts=[(x1,y1), (x2,y2)], gridsize=pixels_per_unit, fb=fb)

        pt_bufr.append( (x,y) )
        bloody_simple_2drender('unit_circle_%s.png'%fr_cnt, pts=pt_bufr, gridsize=pixels_per_unit, fb=fb)

        x+=1
        if d<0 :
            d += dy + dy
            # print('D < 0 IN LOOP %s\n'% d )

        else:
            d += 2*(dy-dx)
            # print('D !<0 IN LOOP %s\n'% d )

            y+=1
        
        fr_cnt += 1 


        



###################################################





def unit_circle_anim( ):
    """ animate the right triangle(s) sweeping the full 360 degrees of a unit circle 
    """

    pixels_per_unit = 200 

    fr_cnt = 0
    for theta in range(1,360,1):
        
        fb = PixelOp()   
        fb.create_buffer(800, 800)
        fb.graticule(pixels_per_unit)

        # print('THETA IS ', theta )

        x = math.cos(mu.dtr( theta) )    
        y = math.sin(mu.dtr( theta) ) 

        hypot  = vec3(x,y,0)
        adaj   = vec3(x,0,0)
        oppos  = vec3(0,y,0)
      
        #form the 3 vectors of the right triangle 
        obj = object3d() 
        obj.one_vec_to_obj(hypot)
        obj.one_vec_to_obj(adaj)
        obj.one_vec_to_obj(oppos, adaj)

        #put a cube at the current theta angle 
        obj.prim_cube(pos=(x, y, 0), size=.05,linecolor=(255,0,0),rot=(0,0,0),pivot='world')


        #calculate the points between 1 and theta to form a circle 
        dots = []
        for dot in range(1,theta,1):
            xx = math.cos(mu.dtr( dot) )    
            yy = math.sin(mu.dtr( dot) )             
            dots.append( (xx,yy) ) 
        bloody_simple_2drender('unit_circle_%s.png'%fr_cnt, pts=dots, gridsize=pixels_per_unit, fb=fb)


        #draw the OBJ file forming a right triangle 
        bloody_simple_2drender('unit_circle_%s.png'%fr_cnt, obj=[obj], gridsize=pixels_per_unit, fb=fb)
        fr_cnt += 1


# unit_circle_anim()


###################################################



###################################################



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





