

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



 


def speherical_to_point():
    mu = math_util() 

    sp = spherical(1.5, mu.dtr(10), mu.dtr(80) ) 
    pt=  sp.to_cartesian() 
    
    obj = object3d()
    obj.prim_cube(pos=pt, size=.2, linecolor=(255,0,0), rot=(0,0,0), pivot='world')
    obj.save('ballz.obj') 




def rotate_matrix_to_vec():
    obj = object3d()

    obj.one_vec_to_obj( (0,55,0) ) 
    #obj.rotate_pts( (45,0,0))

    # get the points for this object 
    ptgrp = obj.get_pt_grp()    

    # contruct a matrix to transform them 
    rotated_m33 = matrix33()

    m = rotated_m33.from_vec3( vec3(1,0,0) , 45) 
    
    

    """
    # apply the matrix to the points in the model 
    rotated_points = obj.apply_matrix_ptgrp(ptgrp, m33=rotated_m33) 
    obj.insert_pt_grp(rotated_points)
    obj.one_vec_to_obj( (0,1,0) )     
    obj.save('vectorDaCleaner.obj')
    """

#rotate_matrix_to_vec() 



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


#API  
#get_face_geom



# insert_polygons

# get_geom_edges

# get_face_edges

# sub_select_geom

# get_face_edges2

# extrude_face











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

    #be cautious of large number of polys. It gets slow real quick!
    obj.copy_sop(slice=(1,10), offset=(0,2,0), num=5, distance=.75)
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
    obj2.insert_polygons(geom2[0], geom2[1], asnew_shell=False  ) 
    # obj2.insert_polygons(geom2[0], newpts2  , asnew_shell=False) 
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

    print( obj.get_face_geom(fid ) ) #reindex=True 
    print( obj.get_face_edges(fid ) ) #DEBUG - add reindex 
    print( obj.get_face_normal(fid ) )
    print( obj.get_face_centroid(fid ) )
    

#extrude_single_face(20)


def extrude_single_edge(fid): 
    """ UNFINISHED! """
    obj = object3d()
    obj.load('objects/sphere.obj')

    print( obj.get_face_geom(fid ) ) #reindex=True 
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

    obj.save("durian_fruit.obj")




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



