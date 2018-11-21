

from pygfx.render import *
from pygfx.point_ops import *
from pygfx.math_ops import  *
from pygfx.obj3d import  *






#######################################################
#######################################################

"""

    THIS IS A PLACE TO WORK ON NEW TOOlS, EXPEIMENTS ETC

    If you want to see the "goodies", go look in unit_tests.py 

"""


#######################################################
#######################################################




"""
    IDEAS / TODO list 

def dagops_l_system_extrude():

def time_based_simulation():
def time_based_simulation_bresenham():
def time_based_simulation_fft():
def time_based_simulation_quaternion_slerp():

def visualize_edges_as_little_arrows(obj1, obj2, slice): 
def copy_obj_rotate_to_each_face(obj1, obj2, slice): 
"""




#######################################################

# run the raytracer 

"""
from pygfx.raytrace import  *
rtrace = raytracer() 
rtrace.save_image( rtrace.main() )
"""

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

#test of 3D arrow primitive model 


def vec_to_euler():

    obj = object3d() 
    obj.prim_arrow( axis='y')
    #obj.save('arrow.obj')


def make_normal_all_faces():
    """ the zero indexing is a big problem  
        all functions should use zero indexing, but OBJ face indices are 1 indexed
        not sure how to handle it, but need to get it sorted 
    """
    obj = object3d() 
    obj.load('objects/sphere.obj')

    vectors = [] 

    for fid in range(len(obj.polygons)):
        nrml = obj.get_face_normal(fid) 
        cntr = obj.get_face_centroid(fid)
        vectors.append( (nrml, cntr) )

    obj2 = object3d()    
    obj2.vectorlist_to_obj(vectors)
    obj2.save('normals.obj') 


## vector_between_to_obj(self, r3_1, r3_2):

#######################################################
def select_polygons_spatially():
    """ UNFINISHED -  
        half working but it looks like the wrong faces are being selected
    """
    obj = object3d() 
    obj.load('objects/monkey.obj')

    from_pt = (-1, 0, -1)

    # get the IDS of polygons near a point in space 
    fids = obj.select_by_location('polygons', from_pt, 1.3 ) 

    # use sub select to grab the geometry from the face IDs 
    geom = obj.sub_select_geom( ids=fids , reindex=True)

    # dump that into a new file 
    obj2 = object3d() 
    
    obj2.prim_cube( size=.05, pos=from_pt, pivot='world' )  
    #tmp.prim_cube(linecolor=(255,0,0), size=.05, pos=pt, rot=(ct,ct,ct), pivot='world')

    obj2.insert_polygons(geom[0], geom[1]  ) 
    obj2.save('polygons_near.obj')

    # #dump the vectors in the work buffer to a new object 
    # obj3 = object3d() 
    # #obj3.vectorlist_to_obj( obj.vec_buffer )
    # # r3, pos=None, arrowhead=False 
    # for v_pt in obj.vec_buffer:
    #     obj3.one_vec_to_obj( v_pt, pos=from_pt ) 
    # obj3.save('dist_vectors.obj')



#######################################################


def extrude_single_edge(fid): 
    """ UNFINISHED! """
    obj = object3d()
    obj.load('objects/sphere.obj')

    print( obj.get_face_geom(fid ) ) #reindex=True 
    print( obj.get_face_edges(fid ) ) #DEBUG - add reindex 
    print( obj.get_face_normal(fid ) )
    print( obj.get_face_centroid(fid ) )



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




#######################################################





def loft_test():
    """ UNFINISHED ! """
    obj = object3d()    
    obj.prim_circle() 


#######################################################


def modify_a_subselect():
    """ UNFINSIHED ! """

    obj = object3d()
    #obj.load('objects/sphere.obj')
    obj.load('kube.obj')

    geom = obj.sub_select_geom( slice=[1,5]  , reindex=True )

    # rotate_pts( rot, pts=None, ptgrp=None):
    newpts = obj.rotate_pts(rot=(45,45,45), pts=geom[1])

    print(geom)
    geom2 = obj.sub_select_geom( slice=[5,6]  , reindex=True )
    #newpts2 = obj.rotate_pts((-45,0,45), points=geom2[1])

    obj2 = object3d() 
    #obj2.insert_polygons(geom[0], newpts  )      
    obj2.insert_polygons(geom2[0], geom2[1], asnew_shell=False  ) 
    # obj2.insert_polygons(geom2[0], newpts2  , asnew_shell=False) 
    obj2.save('kube_modify.obj')



#######################################################



def modify_part_of_an_object():
    """ UNFINSIHED ! """

    obj = object3d()
    obj.load('objects/sphere.obj')

    geom = obj.sub_select_geom( slice=(10,50), reindex=True )
    newpts = obj.rotate_pts((45,45,45), points=geom[1])

    obj2 = object3d() 
    obj2.insert_polygons(geom[0], newpts  )      
    obj2.save('sphere_modify.obj')




#######################################################

def triangulate_test():
    """ UNFINISHED - 
        test stacks of operations, 
        repeating this ame op over , etc 
    """
    obj = object3d()
    obj.prim_circle(axis='y', pos=(0,0,0), spokes=124) 
    obj.radial_triangulate_obj( offset=None)#as_new_obj=False
    #obj.radial_triangulate_obj()
    obj.save('triangulated.obj')



#######################################################

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

    obj.save("durian_fruit.obj")







