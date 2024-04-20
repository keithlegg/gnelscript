

from gnelscript.pygfx.math_ops import  *
from gnelscript.pygfx.point_ops import *
from gnelscript.pygfx.obj3d import  *


""" 

    examples of building geometry from scratch - lower level stuff 

"""

def lathe(path, name):
    """ needs to have the same num U and V to work  

        usage:
            # simplest possible lathe example (must be square , num pts == num revolutions)
            pts = [(.1,.1,0),(1,1,0),(2,2,0),(3,3,0)]
            obj.lathe(pts, 4)


            # using bezier curve function 
            num = 23
            start = (1 ,  0, 0)
            ctrl1 = (.5,  0, 0)
            ctrl2 = ( 0, .5, 0)
            end   = (0 ,  1, 0)
            curve = obj.cubic_bezier(num, start, ctrl1, ctrl2, end)
            obj.lathe(curve, num)
    
    """ 
    
    obj = object3d()

    #------ 

    num = 10

    start  = (2 ,  0, 0)
    ctrl1  = (.5,  .5, 0)
    ctrl2  = (3, .5, 0)
    end    = (1 ,  1, 0)
    curve1 = [  start, ctrl1, ctrl2, end ]

    curvepts = obj.cubic_bezier(num, start, ctrl1, ctrl2, end )

    start  = (1   ,  1   , 0)
    ctrl1  = (2.5 ,  1.0   , 0)
    ctrl2  = (1   ,  1.5 , 0)
    end    = (2   ,  2   , 0)
    curve2 = [ start, ctrl1, ctrl2, end ]
  
    tmp      = obj.cubic_bezier(num, start, ctrl1, ctrl2, end )
   
    curvepts.extend(tmp)

    #obj.draw_splines( num, [curve1, curve2], drawctrls=True, drawhulls=True) 

    #obj.lathe2(curvepts, num*2)
    
    obj.lathe(curvepts, num*2)

    #for i in range(120):
    #    obj.extrude_face(i, -.4)

    #scal command kills COLOR - DOH !
    #obj.scale_pts( (1,2,1) )

    #------ 
    
    # num = 20
    # #cross_section = obj.calc_circle(pos=(1,0,0), rot=(0,180,0), dia=.5, axis='z', periodic=True, spokes=num-1)
    # cross_section = obj.calc_circle(pos=(1,0,0), dia=.5, axis='z', periodic=True, spokes=num-1)
    # cross_section.reverse()
    # obj.lathe(cross_section, num)


    obj.save('%s/%s'%(path,name))

##-------------------------------------------------------## 

def matrix_extrusion():
    #rotate a 3d object to a series of points to "extrude" a line

    obj = object3d()

    seedpts = obj.calc_circle(dia=3,spokes=10, start=0, end=180)
    obj.points = seedpts 
    obj.save('vec_connect.obj')

    #obj.vec_connect_pts(pts=seedpts, draw_obj='rect_2d')
    obj.vec_connect_pts(pts=seedpts, draw_obj='rect_2d', axis='x', width=.05)



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



