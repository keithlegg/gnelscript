

from pygfx.render import *
from pygfx.point_ops import *
from pygfx.math_ops import  *
from pygfx.obj3d import  *



#######################################################


"""
obj = object3d()
obj.load_obj('objects/sphere2.obj')
#obj.prim_quad(axis='z')
#obj.prim_cube()
obj.triangulate() 
#obj.calc_face_normals() 
#obj.scale_pts((.2,.2,.2)) 
#obj.rotate_pts((30,30,30)) 
obj.show() 
rscale = 150
ropr = simple_render()
ropr.render_objects.append(obj) 
"""


"""
#ropr.render_multiobj( (100,0,255), 0, 0, 0, 1, rscale )
ropr.SHOW_EDGES = False
#ropr.SHOW_FACE_CENTER = False
#ropr.COLOR_MODE = 'normal'
ropr.COLOR_MODE = 'zdepth'
#ropr.SHOW_EDGES = False 
ropr.scanline(obj, rscale) 
#ropr.render_obj((100,0,255), 0, 0, 0, 1, rscale, object3d=obj)
ropr.save_image('simple_render.png')
"""


#######################################################

""" 


obj = polygon_operator()
obj.vector2_to_object3d(wirez[0],wirez[1])


obj.vector2_to_object3d(wirez2[0],wirez2[1])

obj.prim_cube(pos=(1.2,6.2,-2))
obj.prim_locator(pos=(3,3,3))
obj.prim_locator_xyz(pos=(5,5,5)) 

#dump the wire vectors to an OBJ file 
obj.save_obj('wire_viz.obj', as_lines=True)
""" 

 
#######################################################

"""
obj = object3d()


obj.load_obj('objects/sphere2.obj')
#obj.prim_cube()

obj.triangulate() 
obj.calc_face_normals() 
#obj.scale_pts((.2,.2,.2)) 

rscale = 150

#obj.rotate_pts((30,30,30)) 



ropr = simple_render()
#ropr.render_obj((100,0,255), 0, 0, 0, 1, rscale, object3d=obj)
ropr.SHOW_FACE_CENTER = False
ropr.COLOR_MODE = 'normal'
#ropr.SHOW_EDGES = False 

#ropr.save_image('simple_render.png')
ropr.scanline(obj, rscale)  

"""




#######################################################

""" 
v1 = vec3(1, 0, 0)
v2 = vec3(0, 1, 0)
v3 = vec3()
mu = math_util() 
print( mu.rtd(v2.angle_between(v1))        ) 
print( mu.rtd(v3.np_angle_between(v1, v2)) )
""" 

#######################################################

""" 
obj = object3d()
obj.prim_cube()
#obj.scale_pts((3,3,30))
obj.rotate_pts((30,30,30))
ropr = simple_render()
#                          fov, aspect, znear, zfar)
#mx = m44.buildPerspProjMat( 200, 1, 1, 100)
ropr.render_obj((100,0,255), 0, 0, 0, 1, 150, object3d=obj)
ropr.save_image('simple_render.png')
""" 

#######################################################

""" 
obj = object3d()
obj.prim_cube()
obj.rotate_pts((30,30,30))
obj.scale_pts((.1,.1,.1))

ropr = simple_render()
ropr.anim([obj], linethick=1, numframes=100, scale=100)
"""
#######################################################

"""
# build polygons from scratch points 

obj = object3d()

pts = [(1,1,1),(0,1,1),(-1,-1,1),(2,-2,1)]
polys = [(1,2,3,4), (2,4,1)]
obj._insert_poly_idxs(polys)
obj._insert_points(pts)

pts = [(0,-3,-1),(2,-2,1),(3,-1,1)]
polys = [(1,2,3)]
obj._insert_poly_idxs(polys)
obj._insert_points(pts)

obj.save_obj("radial_trichop.obj")
"""

#######################################################


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

#######################################################
# examples of direct matrix manipulation 


"""
obj = object3d()
obj.prim_cube()
ropr = simple_render()
m33 = matrix33(.4,.4,0,.5,-.5,0,0,0,1)
ropr.render_matrix_obj( m33 , None ,     1,   200, 'custom_render.png' , obj      )
"""


""" 
obj = object3d()
obj.prim_cube()
ropr = simple_render()
m44 = matrix44()
m44.from_euler(0,45,0)
ropr.render_matrix_obj( None , m44 ,     1,   100, 'custom_render.png' , obj      )
""" 

"""

obj = object3d()
#obj.prim_triangle(axis='z', pos=(0,0,0)) 
#obj.prim_circle(axis='z', pos=(0,0,0), spokes=42) 

obj.load_obj('objects/sphere.obj')


obj2 = object3d()
for i in range(obj.numpts):

    edges  = obj.get_face_edges(i)  
    normal = obj.get_face_normal(i)
    pos    = obj.get_face_centroid(i) 

    obj2.vectorlist_to_obj(edges[1])
    obj2.vectorlist_to_obj( [normal], pos)

obj2.save_obj("edges.obj")


ropr = simple_render()
ropr.render_obj((100,0,255), 0, 0, 0, 1, 150, object3d=obj2)
ropr.save_image('simply_render.png')
"""

#print( obj.calc_circle() )


#######################################################

#x = matrix33(-5,5,9,10,0,1,66,77,88)
#x = matrix33(5,1,5,4,4,22,6,0,0)
#iv = x.np_inverse






#######################################################





"""
#render(color, rx, ry, rz, thick, scale, framebuffer=None, object3d =None,):
ropr = render_3d(800,800)
ropr.render((0,255,0), 45, 45, 45, 3, 100) 
ropr.save_image('my3d.png')
"""


#######################################################



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
        tmp.prim_cube(size=.05, pos=pt, rot=(ct,ct,ct))
        ct += 10
        obj.insert(tmp)  
    obj.save_obj("cubey.obj")




def load_obj_build_another_fromit(objectpath):
    """ load an object, 
        turn its normals into another object, 
        render and save image and new object 
    """

    obj = object3d()
    obj.load_obj(objectpath)

    obj2 = object3d()
    for i in range(obj.numpts):
        edges  = obj.get_face_edges(i)  
        normal = obj.get_face_normal(i)
        pos    = obj.get_face_centroid(i) 

        obj2.vectorlist_to_obj(edges[1])
        obj2.vectorlist_to_obj( [normal], pos)

    obj2.save_obj("edges.obj")

    ropr = simple_render()
    ropr.render_obj((100,0,255), 0, 0, 0, 1, 150, object3d=obj2)
    ropr.save_image('simply_render.png')



load_obj_build_another_fromit('objects/sphere.obj')

