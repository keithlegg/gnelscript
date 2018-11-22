
from pygfx.render import *
from pygfx.point_ops import *
from pygfx.math_ops import  *
from pygfx.obj3d import  *


from pygfx.raytrace import  *



#######################################################

# run the raytracer 

def raytrace():
    
    rtrace = raytracer() 
    rtrace.save_image( rtrace.main() )
 



#######################################################



def build_perspective_matrix():
    #debug - NOT WORKING!  Work In Progress 

    obj = object3d()
    obj.prim_cube()
    #obj.scale_pts((3,3,30))
    obj.rotate_pts((30,30,30))
    ropr = simple_render()
    #                          fov, aspect, znear, zfar)
    #mx = m44.buildPerspProjMat( 200, 1, 1, 100)
    ropr.render_obj((100,0,255), 0, 0, 0, 1, 150, object3d=obj)
    ropr.save_image('simple_render.png')




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


# lighting_test( (10,-3, 0) )

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

#######################################################


def three_renderers():
    """ example of the 3 main ways to render  
            - single object 
            - multi object (single in a loop)
            - scanline 
     """

    obj = object3d()
    #obj.load('objects/sphere2.obj')
    obj.load('extrudez.obj')

    obj.rotate_pts((20,170,170))
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

