
import sys 


from pygfx.render import *
from pygfx.point_ops import *
from pygfx.math_ops import  *
from pygfx.obj3d import  *


from pygfx.raytrace import  *








#######################################################


def make_gif_from_frames():
    # ref - https://engineering.giphy.com/how-to-make-gifs-with-ffmpeg/

    ##convert frames to video 
    ## ffmpeg -r 30 -i my3d_%d.png -s hd480 -vcodec libx264 output.mp4

    ##convert frames to video preserve aspect ratio 
    ## ffmpeg -r 30 -i unit_circle_%d.png -vf scale=w=320:h=240:force_original_aspect_ratio=decrease  -vcodec libx264 output.mp4


    ## convert video to GIF 
    # ffmpeg -i output.mp4 -f gif output.gif

    ## convert video to GIF with speed change
    # ffmpeg -i output.mp4 -f gif -filter:v "setpts=3*PTS" output.gif

    pass





def basic_animation():
    obj = object3d() 
    obj.prim_cube() 

    #obj.load('objects')

    ropr = simple_render()
    ropr.anim(objs=[obj], init_rots=(0,0,0), linethick=5, numframes=5, scale=150 )
#    anim(self, objs, init_rots=(0,0,0), linethick=5, numframes=5, scale=150):

#basic_animation()

#######################################################

# run the raytracer 

def raytrace():
    
    rtrace = raytracer() 
    rtrace.save_image( rtrace.main() )
 



#######################################################



def build_perspective_matrix():
    #debug - NOT WORKING!  Work In Progress 

    obj = object3d()
    
    #obj.prim_cube()
    obj.load('objects/monkey.obj')

    obj.xform_pts( (0,0,-10) )
    #obj.scale_pts((3,3,30))
    #obj.rotate_pts((30,30,30))
    ropr = simple_render()
    
    persp_m44 = matrix44()
    #                          fov, aspect, znear, zfar)
    persp_m44 = persp_m44.buildPerspProjMat( 300, 1, 100, 1000)
    
    print( persp_m44)

    # apply_matrix_pts( pts, m33=None, m44=None):
    obj.points = obj.apply_matrix_pts(obj.points, m44=persp_m44) 

    ropr.render_obj((100,0,255), 0, 0, 0, 1, 150, object3d=obj)
    ropr.save_image('simple_render.png')

    # render_matrix_obj (self,  m33, m44, thick, scale, filename, object3d =None):
    #ropr.render_matrix_obj( m33=None, m44=persp_m44, thick=1, scale=200, filename='out.png', object3d=obj)


# build_perspective_matrix() 


#####################################################
def pass_matrix_to_render():
    """ use a 3X3 or 4X4 matrix to adjust a render 
        attempt to "visualize" a matrix 
    """

    obj = object3d()
    obj.prim_cube(pos=(0,0,0), rot=(0,0,0), linecolor=(255,0,0))
    ropr = simple_render()
    m44 = matrix44()
    m44.from_euler(45,45,0)
    ropr.render_matrix_obj( None, m44, 3, 100, 'custom_render.png' , obj      )


#######################################################

def texmapping_test(fnum=1):
    
    obj = object3d()
    
    #obj.load('objects/monkey.obj')
    obj.load('objects/sphere.obj')
    #obj.load('objects/cube.obj')

    #obj.prim_quad(axis='z',  pos=(0,0,0), rot=(0,0,0))
    obj.triangulate() 
    obj.rotate_pts( (200, -20, 0) )


    # obj.rotate_pts( (180,0,0) )

    # load the texture to map to polygons 
    img_op = PixelOp()   
    
    #img_op.load('tex.png') 
    img_op.load('uvmap_sm.jpg') 

    #------------
    #you can do PIL operations on the images at anytime!
    """
    from PIL import ImageEnhance
    enhancer = ImageEnhance.Brightness(img_op.fb)
    img_op.fb = enhancer.enhance(3) #.show("Sharpness %f" % factor)
    """
    #------------

    render_scale = 200
    lightpos = (0,-1,-3)

    ropr = simple_render()

 
    #ropr.COLOR_MODE = 'lighted'
    ropr.COLOR_MODE = 'lightedshaded'
    
    ropr.SHOW_FACE_CENTER = False
    ropr.SHOW_EDGES       = False     

    ropr.scanline(obj, render_scale, lightpos=lightpos, texmap=img_op ) 
    ropr.save_image('simple_render_%s.png'%fnum)


# texmapping_test(30)

#ANIMATE IT 
#for a in range(20):
#    texmapping_test(a) 

#######################################################

def lighting_test( lightpos, fnum=1):
    """ run the scanline render with a lighting model 
        lighting model shades the polygons based on angle to a point
        the point becomes a cheap lighting simulation 

        #to run : 

           lighting_test( (0,10, 10), 1 )

    """
    save_lighting_data_as_objects = True 

    obj = object3d()
    obj.load('objects/monkey.obj')
    obj.rotate_pts( (180,0,0) )

    #obj.load('objects/sphere.obj')
    #obj.prim_quad(axis='z',  pos=(0,0,0), rot=(0,0,0)) 
    #obj.points = obj.rotate_pts((-10,180,180),pts=obj.points)
    
    obj.triangulate() 

    ropr = simple_render()

    render_linecolor = (255,0,255)
    render_scale = 200 

    #################
    # some render properties you can tweak 

    ## ropr.SHOW_EDGES = False
    ropr.SHOW_FACE_CENTER = False
    ## ropr.COLOR_MODE = 'flat'
    ropr.SHOW_EDGES = False 
    #ropr.COLOR_MODE = 'normal'
    ropr.COLOR_MODE = 'lighted'

    #################
    ## scanline render 
    ropr.scanline(obj, render_scale, lightpos=lightpos ) 
    ropr.save_image('simple_render_%s.png'%fnum)

    #################    
    """ this is cool! 
        we can take the data generated during the render process 
        and dump it back out to a 3D object to verify it is working 
        how we think it should be. 

    """
    if save_lighting_data_as_objects:

        #keep a copy of the triangulated mesh for reference 
        obj.save('triangulated_mesh.obj')

        # do a dump of the lighting data vectors as an OBJ file 
        obj2 = object3d() 
        obj2.prim_cube(pos=lightpos,size=.05,linecolor=(255,0,0),rot=(0,0,0),pivot='world')
        
        LV = ropr.lighting_vectors

        # print( "## debug size of lighting vectors is %s "%len( LV ) ) 

        # lighting_vectors format is:
        # fcntr,   nrml,   vec_to_light,   angle,  poly (polygon_indices)     

        #----
        for v in  LV:    
            # vector from face center to light
            obj2.one_vec_to_obj( v[2] , pos=v[0] )  
        obj2.save('light_to_face_vectors.obj')
        obj2.flush() 

        #----
        for v in  LV:  
            # unit length face normal, from world origin   
            obj2.one_vec_to_obj( v[1] , pos=v[0])  
        obj2.save('face_normal_vectors.obj')
        obj2.flush() 

        #----
        # experiment to extract the parts of the model that are facing the light 
        # it works!!
        threshold = 75  #angle to determine which face to keep or not 

        points = []
        polys  = []
        for v in  LV:  
            if v[3]<threshold:  
                polys.append(  v[4] ) 

        obj2.insert_polygons(polys, obj.points) 
        obj2.save('visible_faces.obj')
        obj2.flush() 


#lighting_test( (10,-3, 0) )

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
    
def animate_light_in_spherical_coords():
    """ generate some 3d positions in a spherical coordinates 
        and call the renderer in a loop with a rotating light 
        slow, but it works 
    """

    mu = math_util() 
    obj = object3d()

    fnum = 0
    for theta in range(-180,180,30):
        print('## theta ', theta )
        for phi in range(-180,180,30):        
            sp = spherical(1.5, mu.dtr(theta), mu.dtr(phi) ) 
            pt=  sp.to_cartesian() 
           
            lighting_test(pt,fnum)
            fnum+=1 


# animate_light_in_spherical_coords()
#######################################################


def three_renderers():
    """ example of the 3 main ways to render  
            - single object 
            - multi object (single in a loop)
            - scanline 
     """

    obj = object3d()
    obj.load('objects/sphere2.obj')
    #obj.load('original_sin.obj')

    #obj.rotate_pts((20,170,170))
    obj.triangulate() 

    ropr = simple_render()

    render_linecolor = (255,0,255)
    render_scale = 200 

    ####

    ## # some render properties you can tweak 
    ## ropr.SHOW_EDGES = False
    ## ropr.SHOW_FACE_CENTER = False
    ## ropr.COLOR_MODE = 'normal'
    ## ropr.COLOR_MODE = 'flat'
    ## ropr.SHOW_EDGES = True 

    ####

    # render single object 
    ropr.render_obj((100,0,255), 0, 0, 0, 1, 150, object3d=obj)
    ropr.save_image('simple_render.png')
    
    ####

    ##  render multiple objects
    ## obj2 = object3d()
    ## obj2.prim_quad()
    ## ropr.render_objects.append(obj) 
    ## ropr.render_objects.append(obj2) 
    ## #                GS (  color,  rx, ry, rz, linethick, scale)
    ## ropr.render_multiobj( render_linecolor, 45, 45, 45, 4, render_scale) 

    ####
    ## scanline render 
    #ropr.scanline(obj, render_scale) 
    #ropr.save_image('simple_render.png')




#three_renderers()



########################################################



#obj = object3d() 
#obj.load('original_sin.obj')
#bloody_simple_2drender('2d_render.png', obj=obj, gridsize=200)
