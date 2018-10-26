

from pygfx.render import *
from pygfx.point_ops import *




######################################################################



obj = object3d()


obj.load_obj('objects/sphere2.obj')
#obj.prim_cube()

obj.triangulate() 
#obj.calc_face_normals() 
#obj.scale_pts((.2,.2,.2)) 
#obj.rotate_pts((30,30,30)) 


obj.show() 

rscale = 150



ropr = simple_render()
#ropr.render_obj((100,0,255), 0, 0, 0, 1, rscale, object3d=obj)
ropr.SHOW_FACE_CENTER = False
#ropr.COLOR_MODE = 'normal'
#ropr.SHOW_EDGES = False 

#ropr.save_image('simple_render.png')
ropr.scanline(obj, rscale)  



######################################################################

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

 
#############

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




#############

"""
v1 = Vec3d(0,.9,0)
v2 = Vec3d(0, 1,0)
v3 = Vec3d()
print( v3.math.rtd(v2.angle_between(v1))        ) 
print( v3.math.rtd(v3.np_angle_between(v1, v2)) )
"""

#############

#mx = my.buildPerspProjMat( .0001, 1.0, 1, 100)
#ropr.render_matrix_obj( None , mx ,     1,   10000, 'custom_render.png' , obj      )



##############################
# examples of direct matrix manipulation and acroname code  
##############################

"""
obj = polygon_operator()
obj.prim_cube()
ropr = render_3d()
mu = math_util()
ac = acroname() 

rx = 20;ry = 10;rz = 0

dtr = mu.dtr


ropr.render_matrix_obj( matrix33 , None ,     1,   100, 'custom_render.png' , obj      )
 
""" 

#x = matrix_33(-5,5,9,10,0,1,66,77,88)
#x = matrix_33(5,1,5,4,4,22,6,0,0)
#iv = x.np_inverse










######################################################################





"""
#render(color, rx, ry, rz, thick, scale, framebuffer=None, object3d =None,):
ropr = render_3d(800,800)
ropr.render((0,255,0), 45, 45, 45, 3, 100) 
ropr.save_image('my3d.png')
"""


#######################################################



