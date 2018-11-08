

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
obj.load('objects/sphere.obj')

fid = 23
print( obj.get_face_data(fid ) ) #reindex=True 
print( obj.get_face_edges(fid ) ) #DEBUG - add reindex 
print( obj.get_face_normal(fid ) )
print( obj.get_face_centroid(fid ) )
"""


def test_extrude():
    obj = object3d()
    obj.load('objects/sphere.obj')
    #obj.load('objects/kube.obj')
    obj.extrude_face(3)

    obj.save('extrudez.obj')

test_extrude()


def loft_test():
    obj = object3d()    
    obj.prim_circle() 


def test_copysop():
    """ copy SOP is a subselect, copy and transform 
        optional loop and increment 
    """
    obj = object3d() 
    #obj.prim_circle(axis='y') 
    obj.load('objects/sphere.obj')
    #copy_sop( slice=None, ids=None, reindex=False, offset=(0,0,0), num=1):
    obj.copy_sop(slice=(1,50), offset=(0,2,0), num=5)
    obj.save('stax.obj')


def grab_all_pts():
    obj = object3d()
    #obj.load('objects/sphere.obj')
    obj.load('objects/kube.obj')

    all_pts = obj.points
    plyidx1 =  obj.get_pt_ids([0,1,2])  
    plyidx2 =  obj.get_pt_ids([3,5])

    obj2 = object3d() 
    obj2.points = all_pts

    #obj2.insert_polygons( gr1, all_pts   ) 
    #gr3 = obj2.xform_pts((2,2,2), gr2 )

    obj2.polygons.extend(plyidx1)
    obj2.polygons.extend(plyidx2)

    obj2.save('kube_modify.obj')






def modify_a_subselect():
    """ UNFINSIHED ! """

    obj = object3d()
    #obj.load('objects/sphere.obj')
    obj.load('kube.obj')

    geom = obj.sub_select_geom( slice=[1,5]  , reindex=True )
    newpts = obj.rotate_pts((45,45,45), points=geom[1])

    print(geom)
    geom2 = obj.sub_select_geom( slice=[5,6]  , reindex=True )
    #newpts2 = obj.rotate_pts((-45,0,45), points=geom2[1])

    obj2 = object3d() 
    #obj2.insert_polygons(geom[0], newpts  )      
    obj2.insert_polygons(geom2[0], geom2[1], asnewgeom=False  ) 
    # obj2.insert_polygons(geom2[0], newpts2  , asnewgeom=False) 
    obj2.save('kube_modify.obj')












def test_rotate_points():
    obj = object3d()
    obj.load('objects/monkey.obj')
    #pts = [(2,2,2), (4,4,4), (8,8,8)]
    pts2 = obj.rotate_pts((45,45,45) )
    #print(pts2)
    obj.save('foo.obj')



def modify_part_of_an_object():
    """ UNFINSIHED ! """

    obj = object3d()
    obj.load('objects/sphere.obj')

    geom = obj.sub_select_geom( slice=(10,50), reindex=True )
    newpts = obj.rotate_pts((45,45,45), points=geom[1])

    obj2 = object3d() 
    obj2.insert_polygons(geom[0], newpts  )      
    obj2.save('sphere_modify.obj')



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
    obj.load('objects/sphere.obj')

    print( obj.get_face_data(fid ) ) #reindex=True 
    print( obj.get_face_edges(fid ) ) #DEBUG - add reindex 
    print( obj.get_face_normal(fid ) )
    print( obj.get_face_centroid(fid ) )
     
#extrude_face(23)



def extrude_single_edge(fid): 
    """ UNFINISHED! """
    obj = object3d()
    obj.load('objects/sphere.obj')

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
    obj.save('triangulated.obj')



def multi_face_triangulate_offset():
    """ broken - DEBUG """
    
    obj = object3d()
    obj.load('objects/sphere.obj')

    nrmls = []
    for i,p in enumerate(obj.points):
        nrmls.append( (i, obj.get_face_normal(i)*3) )
    for n in nrmls[1:4]:
        print(n[0])
        obj.radial_triangulate_face(n[0], offset=n[1] )


    obj.save("durrian.obj")



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
    obj.save("cubey.obj")

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











def slice_extract_and_makenew():
    """ load two models, extract parts of them using the subselect tool 
        subselect grabs polys and points at the same time, with option to reindex
         
        reindex effectively makes a new object 

        weld them into a new model 
        fekkin awesome mate!  
    """
    
    obj = object3d() 
    obj.load('objects/sphere2.obj')
    geom = obj.sub_select_geom( slice=(0,51) , ids=[100,120,105,53,55,73], reindex=True)


    obj3 = object3d() 
    obj3.load('objects/monkey.obj')
    geom2 = obj3.sub_select_geom( slice=(30,100) , ids=[101,105,148], reindex=True)

    obj2 = object3d() 
    # weld two models together 
    obj2.insert_polygons(geom[0], geom[1]  ) 
    obj2.insert_polygons(geom2[0], geom2[1]  )

    obj2.save('new.obj')



def object_primitives():
    """ demo various built in primitive objects """

    obj = object3d() 

    position = (0,0,0)
    rotation = (0,0,0)
    size = 1 

    do_flush = False

    obj.prim_line( pos=position, rot=rotation, size=size)
    obj.save("new_line.obj")
    if do_flush:
        obj.flush()

    obj.prim_triangle( pos=position, rot=rotation, size=size)
    obj.save("new_triangle.obj")
    if do_flush:
        obj.flush()

    obj.prim_quad( pos=position, rot=rotation, size=size)
    obj.save("new_quad.obj")
    if do_flush:
        obj.flush()

    obj.prim_circle( pos=position, rot=rotation, size=size)
    obj.save("new_circle.obj")
    if do_flush:
        obj.flush()

    obj.prim_sphere( pos=position, rot=rotation, size=size)
    obj.save("new_sphere.obj")
    if do_flush:
        obj.flush()

    obj.prim_cone( pos=position, rot=rotation, size=size)
    obj.save("new_cone.obj")
    if do_flush:
        obj.flush()

    obj.prim_sphere( pos=position, rot=rotation, size=size)
    obj.save("new_sphere.obj")
    if do_flush:
        obj.flush()

    obj.prim_locator( pos=position, rot=rotation, size=size)
    obj.save("new_locator.obj")
    if do_flush:
        obj.flush()

    obj.prim_locator_xyz( pos=position, rot=rotation, size=size)
    obj.save("new_locator_xyz.obj")
    if do_flush:
        obj.flush()

    #obj.prim_arrow( pos=position, rot=rotation, size=size)
    #obj.save("new_arrow.obj")
    #if do_flush:
    #    obj.flush()


def three_renderers():
    """ example of the 3 main ways to render  
            - single object 
            - multi object (single in a loop)
            - scanline 
     """

    obj = object3d()
    obj.load('objects/sphere2.obj')
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




def load_build_another_from_normals(objectpath):
    """ load an object, 
        turn its normals into another line object, 
        render and save image and new object 

        load_build_another_from_normals('objects/sphere.obj')
    """

    obj = object3d()
    obj.load(objectpath)

    obj2 = object3d()
    for i in range(obj.numpts):
        edges  = obj.get_face_edges(i)  
        normal = obj.get_face_normal(i)
        pos    = obj.get_face_centroid(i) 

        obj2.vectorlist_to_obj(edges[1])
        obj2.vectorlist_to_obj( [normal], pos)

    obj2.save("edges.obj")

    ropr = simple_render()
    ropr.render_obj((100,0,255), 0, 0, 0, 1, 150, object3d=obj2)
    ropr.save_image('simply_render.png')




def model_from_scratch(): 
    """ build a new polygon object from points """ 

    obj = object3d()

    #add new geom and auto increment the ids
    pts = [(1,1,1),(0,1,1),(-1,-1,1),(2,-2,1)]
    polys = [(1,2,3,4) ]
    obj.insert_polygons(polys, pts) 

    #add new geom and auto increment the ids
    pts = [(0,-3,-1),(2,-2,1),(3,-1,1)]
    polys = [(3,2,1) ]
    obj.insert_polygons(polys, pts)

    #add polys without new points 
    obj.insert_polygons( [(1,2,3,4,5,6,7),(1,7,2)], None, asnewgeom=False)

    obj.save("my_new_object.obj")



def extract_by_hack():
    """ *slightly* higher level than raw geom 
        use the lookup util to get the pt ids by face index 

    """
    obj = object3d()
    obj.load('objects/kube.obj')

    #mirror all points with out thinking about it 
    all_pts = obj.points

    #extract two chunks of poly ids 
    poly1 =  obj.get_pt_ids([0,1,2])  
    poly2 =  obj.get_pt_ids([3,5])

    #make a new object and dump data into it
    obj2 = object3d() 
    obj2.points = all_pts
    obj2.polygons.extend(poly1)
    obj2.polygons.extend(poly2)

    obj2.save('kube_modify.obj')

