
import sys 

from PIL import Image, ImageOps


from gnelscript.pygfx.render import *
from gnelscript.pygfx.point_ops import *
from gnelscript.pygfx.math_ops import  *
from gnelscript.pygfx.obj3d import  *

from gnelscript import NUMPY_IS_LOADED

if NUMPY_IS_LOADED:
    from gnelscript.tools.raytrace import  *








#######################################################


def make_gif_from_frames():
    # ref - https://engineering.giphy.com/how-to-make-gifs-with-ffmpeg/

    ##convert frames to video 
    ## ffmpeg -r 30 -i my3d_%d.png -s hd480 -vcodec libx264 output.mp4

    ##convert frames to video preserve aspect ratio 
    ## ffmpeg -r 30 -i unit_circle_%d.png -vf scale=w=640:h=480:force_original_aspect_ratio=decrease  -vcodec libx264 output.mp4


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


#######################################################

if NUMPY_IS_LOADED:
    # run the raytracer 
    def raytrace():
        
        rtrace = raytracer() 
        rtrace.save_image( rtrace.main() )


#######################################################
def vector_render(objfile):
    obj = object3d()
    
    obj.load(objfile)

    #obj.prim_triangle( "z", (0,0,0), (0,0,0)  )

    #PIL Images have origin in top left - so flip 
    #obj.rotate_pts( (0, 0, 180) )


    #need to go back and check all functions in gnelscript
    #make work in 2d too 

    #GEOM types 
    #BBOX types 
    #PTGROUP types (zero index, option to shift)

    
    #XFROM_PTS IS BROKEN 
    #transform half of picture width (in pixels) to get to center
    #obj.xform_pts( (-200, -200, 0) )
    
    #TRS_POINTS SEEMS TO WORK 
    #obj.trs_points( ply, translate=(0,0,0), rotate=(0,0,0), scale=(gs,gs,gs) ))
    obj.points = obj.trs_points( obj.points, translate=(-3,-3,0) )


    #obj.triangulate() 

    ropr = simple_render()

    render_linecolor = (255,0,255)
    render_scale = 100 

    ##--------------------------------------##

    ## # some render properties you can tweak 
    ## ropr.SHOW_EDGES         = False
    ## ropr.SHOW_FACE_CENTER   = False
    ## ropr.COLOR_MODE         = 'normal'
    ropr.COLOR_MODE         = 'flat'
    ropr.SHOW_EDGES         = True 
    ropr.SHOW_VTXS             = True 
    ropr.USE_PERSPECTIVE       = False


    ##--------------------------------------##

    if ropr.USE_PERSPECTIVE:
        ropr.SHOW_VTXS = False
        persp_m44 = matrix44()    
        persp_m44 = persp_m44.buildPerspProjMat( 60, 1, 1, 10)
        obj.points = obj.apply_matrix_pts(obj.points, m44=persp_m44)

    # render single object 
    ropr.render_obj((100,0,255), 1, 1, 1, 1, render_scale, object3d=obj)
    ropr.save_image('simple_render_.png' )

    ropr.vec_fr_buffer.show()
    ropr.vec_fr_buffer.scale_pts((.01,.01,.01))
    ropr.vec_fr_buffer.save("3d_obj/vec_buffer.obj")
    

#######################################################


def three_renderers(objfile, fnum):
    """ PIL USES TOP LEFT AS (0,0)
        BECAUSE OF THIS IMAGES APPEAR "UPSIDE DOWN"
        YOU CAN SIMPLY HORIZONTAL FLIP AT THE END TO FIX THIS 

        example of the 3 main ways to render  
            - single object 
            - multi object (single in a loop)
            - scanline 
     """

    obj = object3d()
    obj.load(objfile)

    #obj.xform_pts( (1, 1, 1) )

    obj.rotate_pts( (45, fnum*10, 0) )
    #obj.triangulate() 

    ropr = simple_render()

    render_linecolor = (255,0,255)
    render_scale = 500 

    ##--------------------------------------##

    ## # some render properties you can tweak 
    ## ropr.SHOW_EDGES         = False
    ## ropr.SHOW_FACE_CENTER   = False
    ## ropr.COLOR_MODE         = 'normal'
    ## ropr.COLOR_MODE         = 'flat'
    ## ropr.SHOW_EDGES         = True 
    ropr.SHOW_VTXS             = True 
    ropr.USE_PERSPECTIVE       = False
    
    ##--------------------------------------##

    if ropr.USE_PERSPECTIVE:
        ropr.SHOW_VTXS = False
        persp_m44 = matrix44()    
        persp_m44 = persp_m44.buildPerspProjMat( 60, 1, 1, 10)
        obj.points = obj.apply_matrix_pts(obj.points, m44=persp_m44)

    # render single object 
    ropr.render_obj((100,0,255), 1, 1, 1, 1, render_scale, object3d=obj)
    ropr.save_image('simple_render_%s.png'%fnum)
    
    ##--------------------------------------##

    ##  render multiple objects
    ## obj2 = object3d()
    ## obj2.prim_quad()
    ## ropr.render_objects.append(obj) 
    ## ropr.render_objects.append(obj2) 
    ## #                GS (  color,  rx, ry, rz, linethick, scale)
    ## ropr.render_multiobj( render_linecolor, 45, 45, 45, 4, render_scale) 

    ##--------------------------------------##

    ## ropr.render_normal_lines  = True 
    ## ropr.SHOW_VEC_HITS        = False  #faces cover this up!
    ## ropr.SHOW_FACES           = True
    #ropr.SHOW_NORMALS         = True  
    ## ropr.SHOW_SCREEN_CLIP     = False # broken 
    ## ropr.DO_SCREEN_CLIP       = False # broken 

    #ropr.COLOR_MODE         = 'normal'
    #ropr.COLOR_MODE         = 'flat'

    ## scanline render 
    po = pixel_op()
    po.load("images/in/refer.jpg")
    ropr.scanline(obj, render_scale, texmap=po) 
    ropr.save_image('simple_scanline_%s.png'%fnum)

    ##--------------------------------------##

    
def animate_three_renders():
    for x in range(100):
        three_renderers(x)


#animate_three_renders() 

#######################################################



def build_perspective_matrix(near, far, num):
    #debug - NOT WORKING!  Work In Progress 
 
    render_scale = 20

    obj = object3d()
    
    #obj.prim_cube()
    obj.load('objects/monkey.obj')

    #obj.xform_pts( (0,0,0) )
    obj.rotate_pts( (0,num*10,0) )

    #obj.scale_pts((3,3,30))
    #obj.rotate_pts((30,30,30))
    ropr = simple_render()

    #ropr.COLOR_MODE = 'lighted'
    #ropr.COLOR_MODE = 'lightedshaded'
    #ropr.SHOW_FACE_CENTER = False
    #ropr.SHOW_EDGES       = False     

    persp_m44 = matrix44()
    #                          fov, aspect, znear, zfar)
    #persp_m44 = persp_m44.buildPerspProjMat( 60, 1, near, far)
    persp_m44 = persp_m44.buildPerspProjMat( 60, 1, 1, 10)

    # apply_matrix_pts( pts, m33=None, m44=None):
    #obj.points = obj.apply_matrix_pts(obj.points, m44=persp_m44) 

    ropr.render_obj((100,0,255), 0, 0, 0, 1, render_scale, object3d=obj)
    
    print('rendering frame %s'%num)
    ropr.save_image('anim/persp_render_%s.png'%num)

    # render_matrix_obj (self,  m33, m44, thick, scale, filename, object3d =None):
    #ropr.render_matrix_obj( m33=None, m44=persp_m44, thick=1, scale=200, filename='out.png', object3d=obj)



def animate_persp():
   step = 1 

   ct=0
   for n in range(2,20,step):
       #for f in range(2,50,step):    
       build_perspective_matrix(1,n,ct) 
       ct+=1


#animate_persp() 

#build_perspective_matrix( 1, 10, 1 ) 





#######################################################

def texmapping_test(objfile, texfile, pathout, fnum=1):
    
    obj = object3d()
    
    #obj.load('objects/teapot.obj') #really slow!
    
    #obj.load('objects/cube.obj')    
    
    obj.load(objfile)
    #obj.load('objects/sphere.obj')
    #obj.load('objects/cube.obj')

    #obj.prim_quad(axis='z',  pos=(0,0,0), rot=(0,0,0))
    obj.triangulate() 
    
    obj.rotate_pts( (45, -45, 0) ) #for monkey ,teapot
    
    #obj.rotate_pts( (45, fnum, 45) ) 
    # obj.scale_pts((.5,fnum,.5))
    obj.points = obj.xform_pts( (0, 0, 1),  pts=obj.points ) 

    #obj.scale_pts( (math.sin( fnum/40 ), math.sin( fnum/40 ), math.sin( fnum/40 )) ) #for teapot 


    # obj.rotate_pts( (180,0,0) )

    # load the texture to map to polygons 
    img_op = pixel_op()   
    
    #img_op.load('tex.png') 
    img_op.load(texfile) 


    #------------
    #you can do PIL operations on the images at anytime!
    """
    from PIL import ImageEnhance
    enhancer = ImageEnhance.Brightness(img_op.fb)
    img_op.fb = enhancer.enhance(3) #.show("Sharpness %f" % factor)
    """
    #------------

    render_scale = 10
    lightpos = (0, 1 ,3)

    ropr = simple_render()
 
    #ropr.COLOR_MODE = 'lighted'
    ropr.COLOR_MODE = 'lightedshaded'
    ropr.SHOW_FACE_CENTER = False
    ropr.SHOW_EDGES       = False     

    ropr.USE_PERSPECTIVE  = True 

    if ropr.USE_PERSPECTIVE:
        ropr.SHOW_VTXS = False
        persp_m44 = matrix44()    
        #                                       fov, aspect, znear, zfar):
        persp_m44 = persp_m44.buildPerspProjMat( fnum*2,  1,     1,     10)
        obj.points = obj.apply_matrix_pts(obj.points, m44=persp_m44)

    ropr.scanline(obj, render_scale, lightpos=lightpos, texmap=img_op ) 
    ropr.save_image('%s/simple_render_%s.png'%(pathout,fnum))



# single frame render 
#texmapping_test(30)


def animate_scanline():
    for a in range(11,20,1):
        texmapping_test(a) 


#animate_scanline() 

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

    
    obj.load('3d_obj/teapot.obj')
    image = pixel_op()
    image.load("images/uvmap_sm.jpg")



    #obj.rotate_pts( (180,0,0) )
    obj.scale_pts( (.5,.5,.5) )

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
    ropr.scanline(obj, render_scale, lightpos=lightpos, texmap=image.fb ) 
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
        obj2.save('light_to_face_vectors_%s.obj'%fnum)
        obj2.flush() 

        #----
        for v in  LV:  
            # unit length face normal, from world origin   
            obj2.one_vec_to_obj( v[1] , pos=v[0])  
        obj2.save('face_normal_vectors_%s.obj'%fnum)
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
        obj2.save('visible_faces_%s.obj'%fnum)
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


########################################################



#obj = object3d() 
#obj.load('original_sin.obj')
#bloody_simple_2drender('2d_render.png', obj=obj, gridsize=200)
