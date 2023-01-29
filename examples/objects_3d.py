

from gnelscript.pygfx.render import *
from gnelscript.pygfx.point_ops import *
from gnelscript.pygfx.math_ops import  *
from gnelscript.pygfx.obj3d import  *






## obj = object3d()
## obj.load('objects/cube.obj')
## #pids = obj.indexer( span=(0,20), nth=4 ) 
## geom = obj.sub_select_geom( span=(0,'n') )
## 
## for g in geom:
##     print(g)

        
PYCORE_OBJ_IN = 'objects/sphere.obj'
PYCORE_OBJ_OUT = 'PYCORE.obj'



def do_obj_intersection(filename):
    """ 
        UNTESTED 
        test intersection of a single vector with a multiple triangles  
    """

    # project_pt
    # locate_pt_along3d
    pop = object3d()

    flipray = True 

    pop.load(filename)
    pop.triangulate()

    if flipray:
        ray = (vec3(0,0,-1), vec3(0, 0, 1))
    else: 
        ray = (vec3(0,0,1), vec3(0, 0, -1))

    hits = pop.ray_hit( ray[0],  ray[1])

    n = object3d()

    #ray origin
    n.prim_locator(ray[0], size=.1)
    #ray vector 
    n.one_vec_to_obj(ray[1], ray[0], arrowhead=True)

    fids = []
    if hits:
          #hit location
          for h in hits:
              fids.append(h[0])
              n.prim_locator((h[1]), size=.1)
    
    #geom = pop.sub_select_geom(span=None, ids=fids, reindex=False)
    geom = pop.get_face_geom( fids, reindex=False )
    print(geom)
    n.insert(geom)
    n.save('obj_intersect.obj')



##-------------------------------------------------------## 

def test_copysop():
    """ copy SOP is a subselect, copy and transform the result
        ala Houdini 
        optional loop and increment 
    """
    obj = object3d() 
    obj.load(PYCORE_OBJ_IN)

    #be cautious of large number of polys. It gets slow real quick!
    obj.copy_sop(span=(1,20), offset=(0,.5,0), num=2, distance=.5)
    obj.save(PYCORE_OBJ_OUT)

##-------------------------------------------------------## 
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
    obj.save(PYCORE_OBJ_OUT)

##-------------------------------------------------------## 
def modify_part_of_an_object():
    """ extract a span of polygons
        spatially trasform them
        insert result into a new object  
    """

    obj = object3d()
    obj.load('objects/sphere.obj')

    geom = obj.sub_select_geom( span=(10,50), reindex=True )
    newpts = obj.rotate_pts((45,45,45), pts=geom[1])

    obj2 = object3d() 
    obj2.insert_polygons(geom[0], newpts  )      
    obj2.save(PYCORE_OBJ_OUT)

##-------------------------------------------------------## 
def triangulate_test():
    """ UNFINISHED - 
        test stacks of operations, 
        repeating this ame op over , etc 
    """
    obj = object3d()
    obj.prim_circle() 
    obj.radial_triangulate_obj( offset=None)#as_new_obj=False
    #obj.radial_triangulate_obj()
    obj.save(PYCORE_OBJ_OUT)

##-------------------------------------------------------## 
def multi_face_triangulate_offset():
    """ broken - DEBUG """
    
    obj = object3d()
    obj.load('objects/sphere.obj')


    nrmls = []
    for i,p in enumerate(obj.polygons):
        nrmls.append( (i, obj.get_face_normal(i)*3) )

    print('## NUM NORMALS ', len(nrmls) )

    # something is broken here 
    #for i,n in enumerate(nrmls):
    #    print(i)
    #    obj.radial_triangulate_face(i, offset=n[1] , as_new_obj= True )

    obj.save(PYCORE_OBJ_OUT)


##-------------------------------------------------------## 
def circle_with_cube_all_pts():
    """ make a circle with a rotated cube at each point """

    obj = object3d()
    obj.prim_circle(axis='z', pos=(0,0,0), spokes=42) 
    obj.triangulate(force=True)
    pts = obj.get_face_pts(0) 
    ct = 0
    for pt in pts:
        tmp = object3d()
        tmp.prim_cube(size=.05, pos=pt, rot=(ct,ct,ct), pivot='world')
        ct += 10
        obj.insert(tmp)  
    obj.save(PYCORE_OBJ_OUT)


##-------------------------------------------------------## 
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

    obj.save(PYCORE_OBJ_OUT) 



##-------------------------------------------------------## 
def object_primitives():
    """ demo various built in primitive objects """

    obj = object3d() 

    position = (0,0,0)
    rotation = (0,0,0)
    size = 1 
    axis = 'y'

    do_flush = True

    obj.prim_line( axis=axis, pos=position, rot=rotation, size=size)
    obj.save(PYCORE_OBJ_OUT)
    if do_flush:
        obj._flush()

    obj.prim_triangle( axis=axis, pos=position, rot=rotation, size=size)
    obj.save(PYCORE_OBJ_OUT)
    if do_flush:
        obj._flush()

    obj.prim_quad( axis=axis, pos=position, rot=rotation, size=size)
    obj.save(PYCORE_OBJ_OUT)
    if do_flush:
        obj._flush()

    obj.prim_circle( axis=axis, pos=position, dia=size) #rot=rotation
    obj.save(PYCORE_OBJ_OUT)
    if do_flush:
        obj._flush()

    obj.prim_sphere(  pos=position, rot=rotation, size=size)
    obj.save(PYCORE_OBJ_OUT)
    if do_flush:
        obj._flush()

    obj.prim_locator(  pos=position, rot=rotation, size=size)
    obj.save(PYCORE_OBJ_OUT)
    if do_flush:
        obj._flush()

    obj.prim_locator_xyz(  pos=position, rot=rotation, size=size)
    obj.save(PYCORE_OBJ_OUT)
    if do_flush:
        obj._flush()

    obj.prim_cone( axis=axis, pos=position, dia=size) #rot=rotation
    obj.save(PYCORE_OBJ_OUT)
    if do_flush:
        obj._flush()




##-------------------------------------------------------## 
def test_point_transform(): 
    """ example of translate, rotate, scale of raw points 
        translate tools work with "ptgroups", or raw points
    """

    obj = object3d()
    obj.load('objects/monkey.obj')

    obj.points = obj.scale_pts(1.5      , pts=obj.points )   

    obj.points = obj.rotate_pts((0,30,0), pts=obj.points ) 

    obj.points = obj.xform_pts( (0,2,0),  pts=obj.points ) 

    obj.save(PYCORE_OBJ_OUT)



##-------------------------------------------------------## 
def face_extrude():
    """ brute force test of face extrude 
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

    obj.save(PYCORE_OBJ_OUT)


##-------------------------------------------------------## 

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
    obj.save(PYCORE_OBJ_OUT)
