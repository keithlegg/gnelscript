
from gnelscript.pygfx.math_ops import  *
from gnelscript.pygfx.raster_ops import *
from gnelscript.pygfx.point_ops import *
from gnelscript.pygfx.point_ops_2d import *
from gnelscript.pygfx.obj3d import  *


from gnelscript.examples.render import bloody_simple_2drender


mu = math_util() 



def model_geom_from_scratch(): 
    """ build a new polygon object in memory from points 
        then insert it into an object and export  
    """ 

    geom  = [[],[]]
    geom2 = [[],[]]

    obj = object25d()

    #add new geom and auto increment the ids
    polys = [(1,2,3), (2,3,4) ]
    pts = [(1,1,1),(0,1,1),(-1,-1,1),(2,-2,1)]
    geom = obj.insert_polygons(polys, pts, geom=geom) 

 
    polys = [(1,2,3,4) ]
    pts = [(4,-4.3,-3),(1.5,-2.5,-2.1),(-2,2,-4),(4,-4.2,1)]
    geom2 = obj.insert_polygons(polys, pts, geom=geom2) 

    # use insert to add geom to object 
    obj.insert(geom) 
    obj.insert(geom2) 
 
    # see what we have done, or not done 
    obj.show() 

    obj.save("3d_obj/foo.obj")

#model_geom_from_scratch()


##-------------------------------------------## 
#def make_checker_map( , ):

##-------------------------------------------## 
def render_texture_map(input_img, gridsize=50):
    
    fb = raster_op() 
    #fb.load(input_img) 

    #fb.set_pix()


    print() 

    pass


##-------------------------------------------##
def render_2d_vector(v1, gridsize=50):
    """draw a vector and tell us the angle of it in degrees
       
       vector    : 2 2D tuples e.g. (1.5,1), (0,0)
       gridsize  : specify pixels per linear unit

    """

    fb = pixel_op()   
    fb.create_buffer(800, 800)
    fb.graticule(gridsize)
    fb.render_vector_2d(   v1,  scale=gridsize)
    fb.save('vec.png')


# render_2d_vector( (0,5)  ) 







##-------------------------------------------##


def example_BSR():
    """ BSR = bloody simple renderer """
    pts = [(1,1),(2,2),(3,3)]
    lines = [ [ (1,1), (1,2), (2,1)], [ (6,1), (1,6), (5,-1)] ]

    bloody_simple_2drender('2d_render.png', pts=pts, vecs=pts, lines=lines )



##-------------------------------------------##

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


##-------------------------------------------##
def test_matrix22(gridsize=50):
    """ first test of the 2X2 matrix object
        demonstrates how to:
            multiply 2 matrices together 
            multiply a vector by a matrix  
    """

    v1 = vec2(3,0)
    v2 = vec2(0,3)

    #rotate 45 degrees 
    m22 = matrix22()
    m22.from_euler(45)

    # make a second matrix, also 45 degrees, should give us 90 total 
    m22_2 = matrix22()
    m22_2.from_euler(45)
    m22 =  m22_2 * m22

    # mutliply a vector by the matrix 
    v3 =  m22 * v2 

    fb = pixel_op()   
    fb.create_buffer(800, 800)
    fb.graticule(gridsize)
    
    pts = [ (0,0), (0,1), (2,1), (0,2) ]
    #bloody_simple_2drender('2d_rotation.png', pts=pts, gridsize=50, pfb=fb)

    vecs = [v2,v3]
    bloody_simple_2drender('2d_rotation.png', vecs=vecs, gridsize=50, pfb=fb)

    #rotate the points by matrix multiplication 
    pts = m22.batch_mult_pts(pts) 
    bloody_simple_2drender('2d_rotation.png', pts=pts, gridsize=50, pfb=fb)
    fb.save('2d_rotation.png')


#test_matrix22() 

##-------------------------------------------##
def test_2d_object(gridsize=50):
    """ first test of new 2D object """

    obj = object25d() 

    #you can load a 3D object, Z axis gets ignored
    obj.load('objects/sphere.obj')

    #obj.prim_square()
    #obj.prim_triangle()

    #obj.save('2d_square.obj')

    #rotate 45 degrees 
    m22 = matrix22()
    m22.from_euler(45)

    fb = pixel_op()   
    fb.create_buffer(800, 800)
    fb.graticule(gridsize)

    #rotate the points by matrix multiplication 
    obj.points = m22.batch_mult_pts( obj.points ) 

    #saving a 2d object from 3D flattens it on Z axis. 
    #utterly mangles the topology 
    obj.save('2d_rotated.obj')

    #bloody_simple_2drender('2d_rotation.png', pts=pts, gridsize=200, pfb=fb)
    bloody_simple_2drender('2d_rotation.png', obj=[obj], gridsize=200, pfb=fb)

    fb.save('2d_rotation.png')


# test_2d_object() 


##-------------------------------------------##

def draw_fractal_tree():
    tree = []

    tn = tnode()
    tn.set(0, -20) 
    tree.append(tn)

    fractal(tree, 0,10) 

    pts = tree_to_lines(tree) 

    bloody_simple_2drender('frac_tree.png', lines=pts, gridsize=1, gratic=False)


# draw_fractal_tree() 

##-------------------------------------------##

def project_point_along_2Dvector():
    """ visualize the output of vector function 
        not sure this even works ? debug  
    """ 
    
    # 2d vector 
    a = vec2(  1 , 1  )
    b = vec2( -1 , -1 )
    com = vec2() 

    #fb = pixel_op()   
    #fb.create_buffer(800, 800)
    #fb.graticule(pixels_per_unit)

    vecs = [a,b]
    pts = [com.project_pt(a, b, 2)]

    bloody_simple_2drender('2d_render.png', vecs=vecs, pts=pts, gridsize=40)



##-------------------------------------------##

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

        
##-------------------------------------------##

def animate_bresenham( x1, y1, x2, y2):
    """ 
        make an animation of the bresenham algorithm 

        this is not a full implementation of bresenham, but it works 

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
 
    fb = pixel_op()   
    fb.create_buffer(800, 800)
    fb.graticule(pixels_per_unit)

    fr_cnt = 0 
    pt_bufr = [] 
    while x <= x2:
   
        #draw the start and end point on each framebuffer 
        bloody_simple_2drender('unit_circle_%s.png'%fr_cnt, pts=[(x1,y1), (x2,y2)], gridsize=pixels_per_unit, pfb=fb)

        #calc and draw each step in between 
        pt_bufr.append( (x,y) )
        bloody_simple_2drender('unit_circle_%s.png'%fr_cnt, pts=pt_bufr, gridsize=pixels_per_unit, pfb=fb)

        x+=1
        if d<0 :
            d += dy + dy
            # print('D < 0 IN LOOP %s\n'% d )

        else:
            d += 2*(dy-dx)
            # print('D !<0 IN LOOP %s\n'% d )

            y+=1
        
        fr_cnt += 1 

       



##-------------------------------------------##

def unit_circle_anim( ):
    """ animate the right triangle(s) sweeping the full 360 degrees of a unit circle 
    """

    pixels_per_unit = 200 

    fr_cnt = 0
    for theta in range(1,360,1):
        
        fb = pixel_op()   
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
        bloody_simple_2drender('unit_circle_%s.png'%fr_cnt, pts=dots, gridsize=pixels_per_unit, pfb=fb)


        #draw the OBJ file forming a right triangle 
        bloody_simple_2drender('unit_circle_%s.png'%fr_cnt, obj=[obj], gridsize=pixels_per_unit, pfb=fb)
        fb.save( 'unit_circle_%s.png'%fr_cnt )
        fr_cnt += 1


# unit_circle_anim()


##-------------------------------------------##


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

##-------------------------------------------##

"""
a = vec2(1,1)
b = vec2(2,2)

print( a.project_pt(a, b, 0) ) 

"""


def test_2d_intersect():

    pixels_per_unit = 50 

    fb = pixel_op()   
    fb.create_buffer(800, 800)
    fb.graticule(pixels_per_unit)

    com = vec2()  #container for commands 

    s1 = vec2( 5, 5)
    e1 = vec2(-5,-5)

    s2 = vec2( 2, 3)
    e2 = vec2( 2,-3)
      
    #debug - make auto convert from vec2 so we dont have to re enter these 
    #lines = [ ( s1,e1 ), (s2,e2 ) ] 

    lines = [ ( s1,e1 ),  (s2,e2  ) ] 
    bloody_simple_2drender('XXX', lines=lines, pfb=fb)
    
    pt = com.intersect(s1,e1,s2,e2)
    bloody_simple_2drender('XXX', pts=[pt],  pfb=fb)

    fb.save('intersect_2d.png')

    print( 'vectors intersect at point ', pt )





##-------------------------------------------##






