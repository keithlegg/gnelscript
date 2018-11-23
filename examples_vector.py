


from pygfx.math_ops import  *
from pygfx.point_ops import *
from pygfx.obj3d import  *
from pygfx.render import *






mu = math_util() 

#####################################################

def angle_between_vectors():
    """ print the angle between two vectors """

    v1 = vec3(1, 0, 0)
    v2 = vec3(0, 1, 0)

    v3 = vec3()
    print( mu.rtd(v2.angle_between(v1))        ) 
    print( mu.rtd(v3.np_angle_between(v1, v2)) )



#####################################################

def unit_circle_viewer( ):
    """ UNFINSIHED 
        playground for sin, cos, tan , etc 
        good for looking at polar coords as well 
    """

    radius = 1 
    theta = 63.9 

    #x = math.cos(xcrd)    
    #y = math.sin(ycrd) 

    obj = object3d()

    a = vec3(1,0,0)
    #b = vec3(-4,1,10)

    obj.one_vec_to_obj(a)
    #obj.calc_circle(dia=radius, axis='z', periodic=True, spokes=23)
    obj.prim_circle(dia=radius, axis='z', spokes=23)

    obj.save('original_sin.obj')


#####################################################


def right_triangle_viewer( xcrd, ycrd ):
    """ UNFINSIHED 
        playground for sin, cos, tan , etc 
        good for looking at polar coords as well 

        good reference:
            http://tutorial.math.lamar.edu/pdf/Trig_Cheat_Sheet.pdf

    """

    radius = 1 

    obj = object3d()

    def make_right_triangle(theta):
        """ put a cube at the point the triangle meets the unit circle """
        x = math.cos(mu.dtr( theta) )    
        y = math.sin(mu.dtr( theta) ) 

        hypot  = vec3(x,y,0)
        adaj   = vec3(x,0,0)
        oppos  = vec3(0,y,0)

        print('## ---------- ')
        print('## theta is             %s'%theta)
        print('## length of hypotenuse %s'%hypot.length)
        print('## length of adjacent   %s'%adaj.length)      
        print('## length of opposite   %s'%oppos.length)
        print('## a^2 + b^           = %s'% str(adaj.length**2 + oppos.length**2)  )

        #print('## angle between  h a   %s'%hypot.angle_between(adaj)  )
        #print('## angle between  a o   %s'%adaj.angle_between(oppos)  )      
        #print('## angle between  o h   %s'%oppos.angle_between(hypot) )

        obj.one_vec_to_obj(hypot)
        obj.one_vec_to_obj(adaj)
        obj.one_vec_to_obj(oppos, adaj)

        obj.prim_cube(pos=(x, y, 0), size=.05,linecolor=(255,0,0),rot=(0,0,0),pivot='world')

    ## -------- 
    for i in range(0,360,79):
        make_right_triangle(i)


    #obj.calc_circle(dia=radius, axis='z', periodic=True, spokes=23)
    obj.prim_circle(dia=radius, axis='z', spokes=23)

    obj.save('original_sin.obj')



# right_triangle_viewer(2,2)


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
    
    a = vec3(-2, 0, 2)  
    b = vec3(3, 0, 3)  

    c = a.between(b)

    print("## offset is ",  c) 

    obj.one_vec_to_obj(a)
    obj.one_vec_to_obj(b)
    obj.one_vec_to_obj( c, a ) #move the vector to the end of another - forms a triangle 

    obj.save('offset_between.obj')



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


