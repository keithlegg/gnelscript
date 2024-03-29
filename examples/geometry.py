

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
    geom = obj.insert_polygons(plyids=polys, points=pts, geom_obj=geom) 

    print('#########')
    print(geom)

    polys = [(1,2,3,4) ]
    pts = [(4,-4.3,-3),(1.5,-2.5,-2.1),(-2,2,-4),(4,-4.2,1)]
    geom2 = obj.insert_polygons(plyids=polys, points=pts, geom_obj=geom2) 

    # use insert to add geom to object 
    obj.insert(geom) 
    obj.insert(geom2) 
 
    # see what we have done, or not done 
    obj.show() 

    obj.save("scratch.obj")
    

##-------------------------------------------------------## 
def model_obj_from_scratch(fileout): 
    """ build a new polygon object from points directly into an object """ 

    obj = object3d()

    #add new geom and auto increment the ids
    polys = [(1,2,3,4) ]

    pts = [(1,1,1),(0,1,1),(-1,-1,1),(2,-2,1)]
    obj.points.extend(pts) 

    #add new geom and auto increment the ids
    pts = [(0,-3,-1),(2,-2,1),(3,-1,1)]
    obj.points.extend(pts) 

    #add polys without new points into same "shell"
    obj.insert_polygons( plyids=[(1,2,3,4,5,6,7),(1,7,2)], ans=False)
    
    #add new polygon in a new "shell" 
    obj.insert_polygons( plyids=[(1,2,3,4)], points=[(3,3,3), (3,-4,5), (-4,-2.5,3.1), (6.2,-2.7,8)], ans=True)

    obj.save(fileout+'.obj')


##-------------------------------------------------------## 
def model_geom_from_scratch_calc_normals(fileout): 

    obj   = object3d() # container for 3D object 
    geom  = [[],[]]    # container for some random geometry (optional)

    # add new geom  
    polys = [(1,2,3,4) ]  #if you reverse the indexing, the normal inverts   
    pts   = [ (1,1,0), (1,-1,0), (-1,-1,0), (-1,1,0) ]

    # if you pass a geom container it will populate it 
    # instead of putting the geometry into the object  
    geom = obj.insert_polygons(polys, pts) 
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
    obj.save(fileout+'.obj')

    #######################
    obj2 = object3d()
    obj2.vectorlist_to_obj([normal.normal]) #, pos=centroid)
    obj2.save(fileout+'_vect.obj')

#model_geom_from_scratch_calc_normals() 



