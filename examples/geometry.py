

from gnelscript.pygfx.math_ops import  *
from gnelscript.pygfx.point_ops import *
from gnelscript.pygfx.obj3d import  *


""" 

    examples of building geometry from scratch - lower level stuff 

"""





##-------------------------------------------------------## 

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

    obj.save("scratch.obj")
    

##-------------------------------------------------------## 
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

    obj.save(PYCORE_OBJ_OUT)


##-------------------------------------------------------## 
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
    centroid = obj.centroid(pts[0:3]) 

    # see what we have done, or not done 
    # obj.show() 
    obj.save(PYCORE_OBJ_OUT)

    #######################
    obj2 = object3d()
    obj2.vectorlist_to_obj([normal.normal]) #, pos=centroid)
    obj2.save(PYCORE_OBJ_OUT)

#model_geom_from_scratch_calc_normals() 



