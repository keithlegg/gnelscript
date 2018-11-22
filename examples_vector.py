

from pygfx.render import *
from pygfx.point_ops import *
from pygfx.math_ops import  *
from pygfx.obj3d import  *






#####################################################

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



#####################################################
    
def visualize_cross_product():
    obj = object3d()
    a = vec3(2,2,2)
    b = vec3(-4,1,10)

    obj.one_vec_to_obj(a)
    obj.one_vec_to_obj(b)
    obj.one_vec_to_obj( a.cross(b).normal )

    obj.save('cross_product.obj')





#####################################################

def offset_between_2vecs():
    """ create a vector representing offset between two points """

    obj = object3d()
    com = vec3() #container for commands
    
    a = vec3(2, 0, 2)  
    b = vec3(3, 0, 3)  

    c = a.between(b)

    print("## offset is ",  c) 

    
#####################################################

def angle_between_vectors():
    v1 = vec3(1, 0, 0)
    v2 = vec3(0, 1, 0)
    v3 = vec3()
    mu = math_util() 
    print( mu.rtd(v2.angle_between(v1))        ) 
    print( mu.rtd(v3.np_angle_between(v1, v2)) )




#####################################################
def rotate_around_vec():
    """ test of function to generate a matrix 
        that will rotate around a vector 

        see also matrix.from_euler()

    """

    obj = object3d()

    #add first vector that will be rotated Y axis
    obj.one_vec_to_obj( (0,1,0) ) 
    ptgrp = obj.get_pt_grp()    

    # construct a matrix to transform them 
    rotated_m33 = matrix33()
    m = rotated_m33.from_vec3( vec3(1,0,0) , -45) 
    
    
    # apply the matrix to the points in the model 
    rotated_points = obj.apply_matrix_ptgrp(ptgrp, m33=m) 
    obj.insert_pt_grp(rotated_points)

    #add second vector to compare X axis *2 
    obj.one_vec_to_obj( (2,0,0) )     
    obj.save('vectorDaCleaner.obj')


#####################################################

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


#####################################################
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

        #obj2.vectorlist_to_obj(edges[1])
        #obj2.vectorlist_to_obj( [normal], pos)

        obj2.vectorlist_to_obj(edges[1])

    obj2.save("all_face_normals.obj")

    ropr = simple_render()
    ropr.render_obj((100,0,255), 0, 0, 0, 1, 150, object3d=obj2)
    ropr.save_image('pretty_render.png')

#####################################################

def model_geom_from_scratch_calc_normals(): 

    obj   = object3d() # container for 3D object 
    geom  = [[],[]]    # container for some random geometry (optional)

    # add new geom  
    polys = [(1,2,3,4) ]  #if you reverse the indexing, the normal inverts   
    pts   = [ (1,1,0), (1,-1,0), (-1,-1,0), (-1,1,0) ]

    # if you pass a geom container it will populate it 
    # instead of putting the geometry into the object  
    geom = obj.insert_polygons(polys, pts, geom=geom) 
    # use insert to add geom to object 
    obj.insert(geom) 

    # get the data from face ID with helper functions 
    # normal    = obj.get_face_normal(0)
    # centroid  = obj.get_face_centroid(0) 

    # ... or calculate them yourself.  
    normal   = obj.calc_tripoly_normal(pts[0:3], True)
    centroid = obj.centroid_pts(pts[0:3]) 

    # see what we have done, or not done 
    # obj.show() 
    obj.save("new_geom.obj")

    #######################
    obj2 = object3d()
    obj2.vectorlist_to_obj([normal.normal]) #, pos=centroid)
    obj2.save("new_normal.obj")

#model_geom_from_scratch_calc_normals() 