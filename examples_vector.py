


from pygfx.math_ops import  *
from pygfx.point_ops import *
from pygfx.obj3d import  *
from pygfx.render import *



PYCORE_OBJ_IN = 'objects/sphere.obj'
PYCORE_OBJ_OUT = 'PYCORE.obj'


mu = math_util() 




#####################################################
def rotate_around_vec():
    """ uses numpy!

        test of function to generate a matrix 
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
    obj.save(PYCORE_OBJ_OUT)

# rotate_around_vec()


#####################################################

def render_m33_as_vec3s(m33, transpose=False, vlist=None):
    """ create renderable line geom to visualize a matrix 
        vlist is an optional list of vectors to append 
        to aid in visualizing 
    """

    # a=1,b=0,c=0,
    # d=0,e=1,f=0,
    # g=0,h=0,i=1):

    if transpose:
        v1 = vec3(m33[0], m33[3], m33[6])
        v2 = vec3(m33[1], m33[4], m33[7])
        v3 = vec3(m33[2], m33[5], m33[8])  
    else:
        v1 = vec3(m33[0], m33[1], m33[2])
        v2 = vec3(m33[3], m33[4], m33[5])
        v3 = vec3(m33[6], m33[7], m33[8])
    
    obj = object3d()
    rendervecs = [v1,v2,v3]
    rendervecs.extend(vlist)

    obj.vectorlist_to_obj( rendervecs )
    obj.save(PYCORE_OBJ_OUT)       



def numpy_m33_fromvec():
    my_m33 = matrix33()
    axis = vec3(1,1,1)
    m = my_m33.from_vec3( axis , 90 )

    render_m33_as_vec3s(m, vlist=[axis])


#####################################################

def make_right_triangle(theta, obj=None):
    """ UNTESTED 
        calulate three 3D vectors to form a right triangle based on a theta angle 
        
        if object is passed in : 
            bake those vectors into a 3d line geometry 
            insert a 3D cube at the point the triangle meets the unit circle 

    """
    x = math.cos(mu.dtr( theta) )    
    y = math.sin(mu.dtr( theta) ) 

    axis = 'y'

    if axis == 'x':
        hypot  = vec3(x,y,0)
        adaj   = vec3(x,0,0)
        oppos  = vec3(0,y,0)
        #obj.prim_cube(pos=(x, y, 0), size=.05,linecolor=(255,0,0),rot=(0,0,0),pivot='world')

    if axis == 'y':
        hypot  = vec3(x,0,y)
        adaj   = vec3(x,0,0)
        oppos  = vec3(0,0,y)
        #obj.prim_cube(pos=(x, y, 0), size=.05,linecolor=(255,0,0),rot=(0,0,0),pivot='world')

    if axis == 'z':
        hypot  = vec3(0,x,y)
        adaj   = vec3(0,x,0)
        oppos  = vec3(0,0,y)
        #obj.prim_cube(pos=(x, y, 0), size=.05,linecolor=(255,0,0),rot=(0,0,0),pivot='world')

    # print('## ---------- ')
    # print('## theta is             %s'%theta)
    # print('## length of hypotenuse %s'%hypot.length)
    # print('## length of adjacent   %s'%adaj.length)      
    # print('## length of opposite   %s'%oppos.length)
    # print('## a^2 + b^           = %s'% str(adaj.length**2 + oppos.length**2)  )
    
    # print('## ---------- ')
    #print('## angle between  h a   %s'%hypot.angle_between(adaj)  )
    #print('## angle between  a o   %s'%adaj.angle_between(oppos)  )      
    #print('## angle between  o h   %s'%oppos.angle_between(hypot) )
    
    if obj is None:
        return [oppos, adaj, hypot]

    if obj is not None:
        obj.one_vec_to_obj(hypot)
        obj.one_vec_to_obj(adaj)
        obj.one_vec_to_obj(oppos, adaj)


#####################################################
def homogeneous_test():
    camera = vec3( 0,0,-5 )
 
    #pt_4d = vec4(15,21,0,3)
    
    pt    = vec4(15,21,3)
    pt_4d = vec4(15,21,0,3)


    m44= matrix44()
    
    #m44.from_euler(45,45,45)
    
    persp = m44.buildPerspProjMat(49,1,1,1.1)

    #def (self, fov, aspect, znear, zfar):
    
    print( pt_4d )
    print( persp * pt_4d) 
    print( persp * pt) 




# homogeneous_test() 


"""
def homogeneous_test():
    camera = vec3( 0,0,-5 )
 
    fov = 25
    
    vecs = make_right_triangle(fov)
    #view frustrum boundaries 
    vf_tl = vec4()
    vf_tl.insert(vecs[0])
    vf_tr = vec4()
    vf_tr.insert(vecs[1])

    obj = object3d() 
    obj.one_vec_to_obj(vf_tl)    
    obj.one_vec_to_obj(vf_tr)
    obj.save('frustrum.obj')
"""

#####################################################
"""
def homogeneous_test_2d():
    camera = vec3(0,0,-5)
 
    fov = 90

    #view frustrum boundaries 
    vf_tl = vec4() 
    vf_tr = vec4() 
    #vf_br = vec4() 
    #vf_bl = vec4() 
"""

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

    obj.save(PYCORE_OBJ_OUT)


# unit_circle_viewer()

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

    ## -------- 
    for i in range(0,360,79):
        make_right_triangle(i, obj)

    #obj.calc_circle(dia=radius, axis='z', periodic=True, spokes=23)
    obj.prim_circle(dia=radius, axis='z', spokes=23)
    obj.save(PYCORE_OBJ_OUT)

#right_triangle_viewer(2,2)


#####################################################
    
def visualize_cross_product():
    obj = object3d()
    a = vec3(-1,0,1)
    b = vec3(1,0,1)

    obj.one_vec_to_obj(a)
    obj.one_vec_to_obj(b)

    obj.one_vec_to_obj( a.cross(b).normal )
    #obj.one_vec_to_obj( a.cross(b) )

    obj.save(PYCORE_OBJ_OUT)


#visualize_cross_product()

#####################################################

def offset_between_2vecs():
    """ create a vector representing offset between two points """

    obj = object3d()
    com = vec3() #container for commands
    
    a = vec3(-2, 4, 2)  
    b = vec3(3, 0, 3)  

    c = a.between(b)

    print("## offset is ",  c) 

    obj.one_vec_to_obj(a)
    obj.one_vec_to_obj(b)
    
    #move the vector to the end of another - forms a triangle
    #it becomes a lot more illustrative if we connect the two endpoints 
    obj.one_vec_to_obj( c, a )  

    obj.save(PYCORE_OBJ_OUT)


#offset_between_2vecs() 



#####################################################

def make_normal_all_faces():
    """ the zero indexing is a big problem  
        all functions should use zero indexing, but OBJ face indices are 1 indexed
        not sure how to handle it, but need to get it sorted 
    """
    obj = object3d() 
    obj.load('objects/sphere.obj')

    vectors = [] 

    # iterate all faces, calculate a normal for it 
    # transfrom that normal vector to the 3D center of the face 
    for fid in range(len(obj.polygons)):
        nrml = obj.get_face_normal(fid) 
        cntr = obj.get_face_centroid(fid)
        vectors.append( (nrml, cntr) )

    obj2 = object3d()    
    obj2.vectorlist_to_obj(vectors)
    obj2.save(PYCORE_OBJ_OUT) 

# make_normal_all_faces() 

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

    obj2.save(PYCORE_OBJ_OUT)

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
    obj.save(PYCORE_OBJ_OUT)

    #######################
    obj2 = object3d()
    obj2.vectorlist_to_obj([normal.normal]) #, pos=centroid)
    obj2.save("new_normal.obj")

#model_geom_from_scratch_calc_normals() 





#####################################################

def build_orthogonal_vector():
    """ treats the "line" as an infinite vector 

    """

    obj = object3d()
    com = vec3() #container for commands

    # the point we are "looking" from 
    pt1 = vec3(-1,1,-4)
    obj.prim_cube(pos=pt1,size=.05,linecolor=(255,0,0),rot=(0,0,0),pivot='world')

    # the point of the line origin
    pt2 = vec3(10,-5, 17)
    obj.prim_cube(pos=pt2,size=.1,linecolor=(255,0,0),rot=(0,0,0),pivot='world')    
    
    # the line, needs to be normalized for the math to work  
    display_unitvec = vec3(0,1,0)
    unitvec = display_unitvec.normal
    
    #render it as full size, not unit length 
    obj.one_vec_to_obj( display_unitvec , pos=pt2) 
    #make a negative version as well, to really get the idea of the size 
    display_unitvec = display_unitvec * -1
    obj.one_vec_to_obj( display_unitvec , pos=pt2) 

    d = com.orthogonal_vec_from_pt(pt2, unitvec, pt1)

    #obj.one_vec_to_obj( d*-1 )   
    obj.one_vec_to_obj( d , pt1 )   

    obj.save(PYCORE_OBJ_OUT)


# build_orthogonal_vector() 
