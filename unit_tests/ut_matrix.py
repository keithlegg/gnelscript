

from gnelscript import NUMPY_IS_LOADED


from gnelscript.pygfx.render import *
from gnelscript.pygfx.point_ops import *
from gnelscript.pygfx.math_ops import  *
from gnelscript.pygfx.obj3d import  *




if NUMPY_IS_LOADED:

    def ut_transpose(m33):
        
        print(m33.transpose() )
        print(np.transpose(m33) )


    def ut_matrix_rotate(theta):
        #theta =(45,45,0)

        m33 = matrix33()

        rotmat1=m33.np_rotate( vec3(1,0,0), angle_deg=theta )
        rotmat2=m33.from_euler(theta)
        rotmat3=m33.np_euler_rotation_matrix(theta)

        print("####################################################")

        print('## output of .np_rotate()                \n', rotmat1)
        print('## output of .from_euler()               \n', rotmat2)
        print('## output of .np_euler_rotation_matrix() \n', rotmat2)
    







