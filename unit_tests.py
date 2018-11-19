

from pygfx.render import *
from pygfx.point_ops import *
from pygfx.math_ops import  *
from pygfx.obj3d import  *


#######################################################
"""
     "UNIT TESTS" are just examples that got put into long term storage.

     Not a unit test in what you would think of as a test suite. 
     Not yet, anyway. 

"""
#######################################################



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



#lighting_test( (10,-3, 0) )


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
def test_copysop():
    """ copy SOP is a subselect, copy and transform the result
        ala Houdini 
        optional loop and increment 
    """
    obj = object3d() 
    obj.load('objects/sphere.obj')

    #be cautious of large number of polys. It gets slow real quick!
    obj.copy_sop(slice=(1,10), offset=(0,2,0), num=5, distance=.75)
    obj.save('stax.obj')

#####################################################
def test_rotate_points():
    """ simple example to use one of the 3 standard transform tools 
        
        the 3 tools are:
            xfrom_pts , rotate_pts, scale_pts 

        rotate and scale use matrices internally
        transform uses simple addition 

        the transform tools can be used on 
             a point group
             a list of points
             the whole object (as a pointgroup)
    """
    obj = object3d()
    obj.load('objects/monkey.obj')
    #pts = [(2,2,2), (4,4,4), (8,8,8)]
    pts2 = obj.rotate_pts((45,45,45) )
    #print(pts2)
    obj.save('foo.obj')

#####################################################
def modify_part_of_an_object():
    """ extract a slice of polygons
        spatially trasform them
        insert result into a new object  
    """

    obj = object3d()
    obj.load('objects/sphere.obj')

    geom = obj.sub_select_geom( slice=(10,50), reindex=True )
    newpts = obj.rotate_pts((45,45,45), pts=geom[1])

    obj2 = object3d() 
    obj2.insert_polygons(geom[0], newpts  )      
    obj2.save('sphere_modify.obj')

#####################################################
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

#####################################################
def multi_face_triangulate_offset():
    """ broken - DEBUG """
    
    obj = object3d()
    obj.load('objects/sphere.obj')


    nrmls = []
    for i,p in enumerate(obj.polygons):
        nrmls.append( (i, obj.get_face_normal(i)*3) )

    print('## NUM NORMALS ', len(nrmls) )

    #for n in nrmls:
    #    obj.radial_triangulate_face(n[0], offset=n[1] , as_new_obj= True )

    obj.save("durian_fruit.obj")

multi_face_triangulate_offset() 

#####################################################
def circle_with_cube_all_pts():
    """ make a circle with a rotated cube at each point """

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

#####################################################
def spherical_to_point():
    """ test of spherical coordinates to a cartesian point 
        done in a nested loop to make a sphere
    """

    mu = math_util() 
    obj = object3d()

    for theta in range(-180,180,20):
        print('## theta ', theta )
        for phi in range(-180,180,20):        
            sp = spherical(1.5, mu.dtr(theta), mu.dtr(phi) ) 
            pt=  sp.to_cartesian() 
            obj.prim_cube(pos=pt, size=.1, linecolor=(255,0,0), rot=(0,0,0), pivot='world')

    obj.save('ball_of_cubes.obj') 

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

#####################################################
def object_primitives():
    """ demo various built in primitive objects """

    obj = object3d() 

    position = (0,0,0)
    rotation = (0,0,0)
    size = 1 
    axis = 'y'

    do_flush = True

    obj.prim_line( axis=axis, pos=position, rot=rotation, size=size)
    obj.save("new_line.obj")
    if do_flush:
        obj.flush()

    obj.prim_triangle( axis=axis, pos=position, rot=rotation, size=size)
    obj.save("new_triangle.obj")
    if do_flush:
        obj.flush()

    obj.prim_quad( axis=axis, pos=position, rot=rotation, size=size)
    obj.save("new_quad.obj")
    if do_flush:
        obj.flush()

    obj.prim_circle( axis=axis, pos=position, dia=size) #rot=rotation
    obj.save("new_circle.obj")
    if do_flush:
        obj.flush()

    obj.prim_sphere(  pos=position, rot=rotation, size=size)
    obj.save("new_sphere.obj")
    if do_flush:
        obj.flush()

    obj.prim_locator(  pos=position, rot=rotation, size=size)
    obj.save("new_locator.obj")
    if do_flush:
        obj.flush()

    obj.prim_locator_xyz(  pos=position, rot=rotation, size=size)
    obj.save("new_locator_xyz.obj")
    if do_flush:
        obj.flush()

    obj.prim_cone( axis=axis, pos=position, dia=size) #rot=rotation
    obj.save("new_cone.obj")
    if do_flush:
        obj.flush()
#object_primitives() 
#####################################################
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

#####################################################
def angle_between_vectors():
    v1 = vec3(1, 0, 0)
    v2 = vec3(0, 1, 0)
    v3 = vec3()
    mu = math_util() 
    print( mu.rtd(v2.angle_between(v1))        ) 
    print( mu.rtd(v3.np_angle_between(v1, v2)) )

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
def test_subsel_point_transform(): 
    """ example of translate, rotate, scale of a point group 
        translate tools work with "ptgroups", or raw points 
    """
    obj = object3d()
    obj.load('objects/monkey.obj')

    ptgrp = obj.get_pt_grp( slice=(1,300) )
    xformed_pts = obj.scale_pts(2.5, ptgrp=ptgrp)   

    ptgrp = obj.get_pt_grp( slice=(1,300) )
    xformed_pts = obj.rotate_pts((0,30,0),  ptgrp=ptgrp)

    ptgrp = obj.get_pt_grp( slice=(1,100) )
    xformed_pts = obj.xform_pts( (0,2,0),  ptgrp=ptgrp) 

    obj.save('ptgrp.obj')

#####################################################
def test_point_transform(): 
    """ example of translate, rotate, scale of raw points 
        translate tools work with "ptgroups", or raw points
    """

    obj = object3d()
    obj.load('objects/monkey.obj')

    obj.points = obj.scale_pts(1.5, pts=obj.points )   

    obj.points = obj.rotate_pts((0,30,0), pts=obj.points ) 

    obj.points = obj.xform_pts( (0,2,0),  pts=obj.points ) 

    obj.save('ptgrp.obj')

#####################################################
def slice_extract_and_makenew():
    """ load two models, extract parts of them using the subselect tool 
        subselect grabs polys and points at the same time, with option to reindex
         
        reindex effectively makes a new object 

        weld them into a new model 
        fekkin awesome mate!  
    """
    
    obj = object3d() 
    obj.load('objects/sphere2.obj')
    geom = obj.sub_select_geom( slice=(1,50), ids=[53,55,73], reindex=True)

    obj3 = object3d() 
    obj3.load('objects/monkey.obj')
    geom2 = obj3.sub_select_geom( slice=(3,100) , ids=[101,105,148], reindex=True)

    obj2 = object3d() 
    ## weld two models together 
    obj2.insert_polygons(geom[0], geom[1]  ) 
    obj2.insert_polygons(geom2[0], geom2[1]  )

    obj2.save('new.obj')

#####################################################
def model_geom_from_scratch(): 
    """ build a new polygon object in memory from points 
        then insert it into an object and export  
    """ 

    geom  = [[],[]]
    geom2 = [[],[]]

    obj = object3d()

    #add new geom and auto increment the ids
    polys = [(1,2,3), (2,3,4) ]
    pts = [(1,1,1),(0,1,1),(-1,-1,1),(2,-2,1)]
    geom = obj.insert_polygons(polys, pts, geom=geom) 

 
    polys = [(1,2,3,4) ]
    pts = [(4,-4.3,-3),(1.5,-2.5,-2.1),(-2,2,-4),(4,-4.2,1)]
    geom2 = obj.insert_polygons(polys, pts, geom=geom2) 

    # use insert to add geom to object 
    obj.insert(geom) 
    obj.insert(geom2) 
 
    # see what we have done, or not done 
    obj.show() 

    obj.save("my_new_object.obj")

#####################################################
def extract_by_copy_hack():
    """ *slightly* higher level than raw geom 
        use the face lookup util to get the geom by face index 
    """

    obj = object3d()
    obj.load('objects/kube.obj')

    # duplicate all points with out thinking about it 
    all_pts = obj.points

    # extract two chunks of poly ids 
    polygr1 =  obj.get_pt_ids([0,1,2])  
    polygr2 =  obj.get_pt_ids([3,5])

    # make a new object and dump data into it
    obj2 = object3d() 
    obj2.points = all_pts  #move all points over - DEBUG add a clean func to remove unused

    for ply in polygr1:
        obj2.polygons.extend(ply)
    for ply in polygr2:
        obj2.polygons.extend(ply)        


    obj2.save('kube_modify.obj')

#####################################################
def model_obj_from_scratch(): 
    """ build a new polygon object from points directly into an object """ 

    obj = object3d()

    #add new geom and auto increment the ids
    polys = [(1,2,3,4) ]
    pts = [(1,1,1),(0,1,1),(-1,-1,1),(2,-2,1)]
    obj.insert_polygons([], pts) 

    #add new geom and auto increment the ids
    pts = [(0,-3,-1),(2,-2,1),(3,-1,1)]
    obj.insert_polygons([], pts)

    #add polys without new points into same "shell"
    obj.insert_polygons( [(1,2,3,4,5,6,7),(1,7,2)], None, asnew_shell=False)
    
    #add new polygon in a new "shell" 
    obj.insert_polygons( [(1,2,3,4)], [(3,3,3), (3,-4,5), (-4,-2.5,3.1), (6.2,-2.7,8)], asnew_shell=True)

    obj.save("my_new_object.obj")

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

    obj2.save("edges.obj")

    ropr = simple_render()
    ropr.render_obj((100,0,255), 0, 0, 0, 1, 150, object3d=obj2)
    ropr.save_image('simply_render.png')

#####################################################
def face_extrude():
    """ brure force test of face extrude 
        extrudes all faces in a polygon object 
        also will display a cheapo test "progress bar"
        because it can be slow 
    """

    obj = object3d()
    obj.load('objects/monkey.obj')
   
    tenths = int(obj.numply/10)
    ct = 1
    for i in range(1,len(obj.polygons) ):   
        if i % tenths == 0:
            print('%%%s0 processed.'%ct) 
            ct+=1
        obj.extrude_face(i, 10/i)

    obj.save('extrudez.obj')

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

def offset_between_2vecs():
    """ create a vector representing offset between two points """

    obj = object3d()
    com = vec3() #container for commands
    
    a = vec3(2, 0, 2)  
    b = vec3(3, 0, 3)  

    c = a.between(b)

    print("## offset is ",  c) 


#####################################################

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
