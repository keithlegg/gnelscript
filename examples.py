

from pygfx.render import *
from pygfx.point_ops import *
from pygfx.math_ops import  *
from pygfx.obj3d import  *



#######################################################
# test of new vec4 type 

"""
 
v3 = vec3(41,32,13)
v4 = vec4()

m44 = matrix44(1,3,4,5, 6,7,23,5, 6,3,23,3 ,5,6,7,8 )

v4.from_vec3(v3)

print( v4 )
 
"""

# print( m44.np_inverse )
# print(m44*v4)

#######################################################
# test of new quaternion type 

""" 
q1 = quaternion()
m33 = matrix33()
m33.from_euler(45,0,45)
q1.from_m33(m33)
m9 = q1.to_m33() 
obj = object3d()
obj.prim_cube()
ropr = simple_render()
ropr.render_matrix_obj( m9 , None ,     1,   100, 'custom_render.png' , obj      )
""" 

#######################################################

"""
obj = object3d()
obj.load_obj('objects/sphere.obj')

fid = 23
print( obj.get_face_data(fid ) ) #reindex=True 
print( obj.get_face_edges(fid ) ) #DEBUG - add reindex 
print( obj.get_face_normal(fid ) )
print( obj.get_face_centroid(fid ) )
"""





def loft_test():
    obj = object3d()    
    obj.prim_circle() 




def test_rotate_points():
    obj = object3d()
    obj.load_obj('objects/monkey.obj')
    #pts = [(2,2,2), (4,4,4), (8,8,8)]
    pts2 = obj.rotate_pts((45,45,45) )
    #print(pts2)
    obj.save_obj('foo.obj')



def modify_a_subselect():
    """ UNFINSIHED ! """

    obj = object3d()
    obj.load_obj('objects/sphere.obj')
    geom = obj.sub_select_geom( slice=[1,200]  , reindex=True )
    newpts = obj.rotate_pts((45,45,45), points=geom[1])

    obj2 = object3d() 
    obj2.insert_polygons(geom[0], newpts  )      
    obj2.save_obj('sphere_modify.obj')

modify_a_subselect()



def modify_part_of_an_object():
    """ UNFINSIHED ! """

    obj = object3d()
    obj.load_obj('objects/sphere.obj')

    geom = obj.sub_select_geom( slice=(10,50), reindex=True )
    newpts = obj.rotate_pts((45,45,45), points=geom[1])

    obj2 = object3d() 
    obj2.insert_polygons(geom[0], newpts  )      
    obj2.save_obj('sphere_modify.obj')



################################################


"""
    IDEAS 
def dagops_l_system_extrude():
def time_based_simulation():
def time_based_simulation_bresenham():
def time_based_simulation_fft():
def time_based_simulation_quaternion_slerp():
def visualize_edges_as_little_arrows(obj1, obj2, slice): 
def copy_obj_rotate_to_each_face(obj1, obj2, slice): 
"""

def extrude_single_face(fid): 
    """ UNFINISHED! """
    obj = object3d()
    obj.load_obj('objects/sphere.obj')

    print( obj.get_face_data(fid ) ) #reindex=True 
    print( obj.get_face_edges(fid ) ) #DEBUG - add reindex 
    print( obj.get_face_normal(fid ) )
    print( obj.get_face_centroid(fid ) )
     
#extrude_face(23)



def extrude_single_edge(fid): 
    """ UNFINISHED! """
    obj = object3d()
    obj.load_obj('objects/sphere.obj')

    print( obj.get_face_data(fid ) ) #reindex=True 
    print( obj.get_face_edges(fid ) ) #DEBUG - add reindex 
    print( obj.get_face_normal(fid ) )
    print( obj.get_face_centroid(fid ) )



def triangulate_test():
    """ UNFINISHED - 
        test stacks of operations, 
        repeating this ame op over , etc 
    """
    obj = object3d()
    obj.prim_circle() 
    obj.radial_triangulate_obj( offset=None)#as_new_obj=False
    #obj.radial_triangulate_obj()
    obj.save_obj('triangulated.obj')



def multi_face_triangulate_offset():
    """ broken - DEBUG """
    
    obj = object3d()
    obj.load_obj('objects/sphere.obj')

    nrmls = []
    for i,p in enumerate(obj.points):
        nrmls.append( (i, obj.get_face_normal(i)*3) )
    for n in nrmls[1:4]:
        print(n[0])
        obj.radial_triangulate_face(n[0], offset=n[1] )


    obj.save_obj("durrian.obj")



def circle_with_cube_all_pts():
    """ BROKEN - FIX THIS 
        make a circle with a rotated cube at each point 
    """

    obj = object3d()
    obj.prim_circle(axis='z', pos=(0,0,0), spokes=42) 
    ctr = obj.get_face_centroid(0)
    obj.triangulate(force=True)
    pts = obj.get_face_pts(0) 
    ct = 0
    for pt in pts:
        tmp = object3d()
        tmp.prim_cube(size=.05, pos=pt, rot=(ct,ct,ct), pivot='world')
        ct += 10
        obj.insert(tmp)  
    obj.save_obj("cubey.obj")

#######################################################
#######################################################
#######################################################
#######################################################
#these are all tested-ish 


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



def pass_matrix_to_render():
    """ use a 3X3 or 4X4 matrix to adjust a render 
        attempt to "visualize" a matrix 
    """

    obj = object3d()
    obj.prim_cube()
    ropr = simple_render()
    m44 = matrix44()
    m44.from_euler(45,45,0)
    ropr.render_matrix_obj( None, m44, 3, 100, 'custom_render.png' , obj      )








def test_copysop():
    obj = object3d() 
    obj.prim_circle(axis='y') 
    #copy_sop( slice=None, ids=None, reindex=False, offset=(0,0,0), num=1):
    obj.copy_sop(ids=[0], offset=(0,2,0), num=4)
    obj.save_obj('stax.obj')


def slice_extract_and_makenew():

    """ load two models, extract parts of them, weld them into a new model 
        fekkin awesome mate!  
    """
    
    obj = object3d() 
    obj.load_obj('objects/sphere2.obj')
    geom = obj.sub_select_geom( slice=(0,51) , ids=[100,120,105,53,55,73], reindex=True)


    obj3 = object3d() 
    obj3.load_obj('objects/monkey.obj')
    geom2 = obj3.sub_select_geom( slice=(30,100) , ids=[101,105,148], reindex=True)


    obj2 = object3d() 
    # weld two models together 
    obj2.insert_polygons(geom[0], geom[1]  ) 
    obj2.insert_polygons(geom2[0], geom2[1]  )

    obj2.save_obj('new.obj')



def object_primitives():
    """ demo various built in primitive objects """

    obj = object3d() 

    position = (0,0,0)
    rotation = (0,0,0)
    size = 1 

    do_flush = False

    obj.prim_line( pos=position, rot=rotation, size=size)
    obj.save_obj("new_line.obj")
    if do_flush:
        obj.flush()

    obj.prim_triangle( pos=position, rot=rotation, size=size)
    obj.save_obj("new_triangle.obj")
    if do_flush:
        obj.flush()

    obj.prim_quad( pos=position, rot=rotation, size=size)
    obj.save_obj("new_quad.obj")
    if do_flush:
        obj.flush()

    obj.prim_circle( pos=position, rot=rotation, size=size)
    obj.save_obj("new_circle.obj")
    if do_flush:
        obj.flush()

    obj.prim_sphere( pos=position, rot=rotation, size=size)
    obj.save_obj("new_sphere.obj")
    if do_flush:
        obj.flush()

    obj.prim_cone( pos=position, rot=rotation, size=size)
    obj.save_obj("new_cone.obj")
    if do_flush:
        obj.flush()

    obj.prim_sphere( pos=position, rot=rotation, size=size)
    obj.save_obj("new_sphere.obj")
    if do_flush:
        obj.flush()

    obj.prim_locator( pos=position, rot=rotation, size=size)
    obj.save_obj("new_locator.obj")
    if do_flush:
        obj.flush()

    obj.prim_locator_xyz( pos=position, rot=rotation, size=size)
    obj.save_obj("new_locator_xyz.obj")
    if do_flush:
        obj.flush()

    #obj.prim_arrow( pos=position, rot=rotation, size=size)
    #obj.save_obj("new_arrow.obj")
    #if do_flush:
    #    obj.flush()


def three_renderers():
    """ example of the 3 main ways to render  
            - single object 
            - multi object (single in a loop)
            - scanline 
     """

    obj = object3d()
    obj.load_obj('objects/sphere2.obj')
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
    #ropr.render_obj((100,0,255), 0, 0, 0, 1, 150, object3d=obj)

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
    ropr.scanline(obj, render_scale) 

    ropr.save_image('simple_render.png')




def angle_between_vectors():
    v1 = vec3(1, 0, 0)
    v2 = vec3(0, 1, 0)
    v3 = vec3()
    mu = math_util() 
    print( mu.rtd(v2.angle_between(v1))        ) 
    print( mu.rtd(v3.np_angle_between(v1, v2)) )




def model_from_scratch(): 
    """ build a new polygon object from points """ 

    obj = object3d()

    # pass one - you can do as many as you like
    pts = [(1,1,1),(0,1,1),(-1,-1,1),(2,-2,1)]
    polys = [(1,2,3,4) ]
    obj.insert_polygons(polys, pts)

    # pass two - add as many times as you want 
    pts = [(0,-3,-1),(2,-2,1),(3,-1,1)]
    polys = [(3,2,1) ]
    #reindex True  -  assume this is a new object to merge with old (default) 
    #reindex False -  assume this is a continuation of past geometry 
    obj.insert_polygons(polys, pts)

    obj.save_obj("my_new_object.obj")





def load_obj_build_another_from_it(objectpath):
    """ load an object, 
        turn its normals into another object, 
        render and save image and new object 

        load_obj_build_another_from_it('objects/sphere.obj')
        
    """

    obj = object3d()
    obj.load_obj(objectpath)

    obj2 = object3d()

    obj2.save_obj("edges.obj")

    ropr = simple_render()
    ropr.render_obj((100,0,255), 0, 0, 0, 1, 150, object3d=obj2)
    ropr.save_image('simply_render.png')




def load_obj_build_another_from_normals(objectpath):
    """ load an object, 
        turn its normals into another object, 
        render and save image and new object 

        load_obj_build_another_from_normals('objects/sphere.obj')

    """

    obj = object3d()
    obj.load_obj(objectpath)

    obj2 = object3d()
    for i in range(obj.numpts):
        edges  = obj.get_face_edges(i)  
        normal = obj.get_face_normal(i)
        pos    = obj.get_face_centroid(i) 

        obj2.vectorlist_to_obj(edges[1])
        obj2.vectorlist_to_obj( [normal], pos)

    obj2.save_obj("edges.obj")

    ropr = simple_render()
    ropr.render_obj((100,0,255), 0, 0, 0, 1, 150, object3d=obj2)
    ropr.save_image('simply_render.png')




