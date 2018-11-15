

from pygfx.render import *
from pygfx.point_ops import *
from pygfx.math_ops import  *
from pygfx.obj3d import  *





#######################################################

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


#####################################################





"""
dia = .1

obj = object3d()
axis = 'y'
pos = (.75,0,0)

obj.prim_circle(axis=axis, pos=pos, dia=dia)

tiplen = dia*2

if axis=='x':
    oset = (-tiplen,0,0)
if axis=='y':
    oset = (0,-tiplen,0)            
if axis=='z':
    oset = (0,0,-tiplen) 
obj.radial_triangulate_face(1, offset=oset )


obj2 = object3d()
obj2.prim_circle(axis=axis, pos=pos, spokes=4 , dia=.02)
obj2.extrude_face(1, distance=.75)
#obj.insert(obj2)


obj.save("kone.obj")
"""


"""
mu = math_util() 
obj = object3d()
sp = spherical(1.73,-.95,-.78) 

#cart = vec3(1,1,1)
#sp.from_cartesian( cart )
#print('### spherical ', sp)

print( '##### to cartesian ', sp.to_cartesian() ) 

#obj.prim_cube(pos=pt, size=.1, linecolor=(255,0,0), rot=(0,0,0), pivot='world')

#obj.save('ball_of_cubes.obj') 
"""



def lighting_test():
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
    ropr.COLOR_MODE = 'lighted'
    ## ropr.SHOW_EDGES = True 


    ## scanline render 
    ropr.scanline(obj, render_scale) 
    ropr.save_image('simple_render.png')



lighting_test() 



def offset_between_2vecs():
    obj = object3d()
    com = vec3() #container for commands
    
    a = vec3(2, 0, 2)  
    b = vec3(3, 0, 3)  

    c = b - a
    print("## offset is ",  c) 




def extrude_single_edge(fid): 
    """ UNFINISHED! """
    obj = object3d()
    obj.load('objects/sphere.obj')

    print( obj.get_face_geom(fid ) ) #reindex=True 
    print( obj.get_face_edges(fid ) ) #DEBUG - add reindex 
    print( obj.get_face_normal(fid ) )
    print( obj.get_face_centroid(fid ) )



def build_orthogonal_vector():
    obj = object3d()
    com = vec3() #container for commands

    # the point we are "looking" from 
    pt1 = vec3(-1,1,1)
    obj.prim_cube(pos=pt1,size=.05,linecolor=(255,0,0),rot=(0,0,0),pivot='world')

    # the point of the line origin
    pt2 = vec3(-2,.3,.9)
    obj.prim_cube(pos=pt2,size=.1,linecolor=(255,0,0),rot=(0,0,0),pivot='world')    
    
    # the line, needs to be normalized for the math to work  
    display_unitvec = vec3(0,5,0)
    unitvec = display_unitvec.normal
    
    #render it as full size, not unit length 
    obj.one_vec_to_obj( display_unitvec , pos=pt2) 
    #make a negative version as well, to really get the idea of the size 
    display_unitvec = display_unitvec * -1
    obj.one_vec_to_obj( display_unitvec , pos=pt2) 

    d= com.orthogonal_vec_from_pt(pt2, unitvec, pt1)
    obj.one_vec_to_obj( d , pos=pt1)  

    obj.save('perpendicular.obj')







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

def prim_arrow(axis): 
    
    obj = object3d() 
    #obj.prim_circle(axis=axis, spokes=4 , dia=.02)
    #obj.extrude_face(1, distance=.75)
    
    dist = 3
    if axis=='x':
        #posi=[dist,0,0]
        obj.prim_cone(axis=axis, pos=[2,0,0], dia=.1 )

    if axis=='y':
        #posi = [0,dist,0]
        obj.prim_cone(axis=axis, pos=[0,2,0], dia=.1 )    
    
    if axis=='z':
        #posi = [0,0,dist]
        obj.prim_cone(axis=axis, pos=[0,0,2], dia=.1 ) 


    obj.save("arrow.obj")





#######################################################








def loft_test():
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



################################################

def extrude_single_face(fid): 
    """ UNFINISHED! """
    obj = object3d()
    obj.load('objects/sphere.obj')

    print( obj.get_face_geom(fid ) ) #reindex=True 
    print( obj.get_face_edges(fid ) ) #DEBUG - add reindex 
    print( obj.get_face_normal(fid ) )
    print( obj.get_face_centroid(fid ) )
    

#extrude_single_face(20)






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




def circle_with_cube_all_pts():
    """ BROKEN - FIX THIS 
        make a circle with a rotated cube at each point 
    """

    obj = object3d()
    obj.prim_circle(axis='z', pos=(0,0,0), spokes=22) 
    ctr = obj.get_face_centroid(0)
    obj.triangulate(force=True)
    pts = obj.get_face_pts(0) 
    ct = 0
    for pt in pts:
        tmp = object3d()
        tmp.prim_cube(linecolor=(255,0,0), size=.05, pos=pt, rot=(ct,ct,ct), pivot='world')
        ct += 10
        obj.insert(tmp)  
    obj.save("cube_pts.obj")



def rotate_around_vec():
    obj = object3d()

    #add first vector that will be rotated Y axis
    obj.one_vec_to_obj( (0,1,0) ) 
    ptgrp = obj.get_pt_grp()    

    # construct a matrix to transform them 
    rotated_m33 = matrix33()
    m = rotated_m33.from_vec3( vec3(1,0,0) , -45) 
    #rotated_m33.from_euler(90,0,0)
    
    # apply the matrix to the points in the model 
    rotated_points = obj.apply_matrix_ptgrp(ptgrp, m33=m) 
    obj.insert_pt_grp(rotated_points)

    #add second vector to compare X axis *2 
    obj.one_vec_to_obj( (2,0,0) )     
    obj.save('vectorDaCleaner.obj')

