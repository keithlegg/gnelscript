
from PIL import Image, ImageOps

from gnelscript.pygfx.math_ops import  *
from gnelscript.pygfx.raster_ops import *
from gnelscript.pygfx.point_ops_2d import *

#from pygfx.point_ops import *
#from pygfx.obj3d import  *


mu = math_util() 



COMMON_IN  = 'foobar'
COMMON_OUT = 'rasterfoo'
COMMON_EXT = 'png'

##----------------------------------------------------






##----------------------------------------------------
def vectorize( inputfile, outputfile=None):
    """
        Experiment to to an "Adobe streamline" type tool 
        Far from working, not even started. 

        load image 
        posterize 
        sort pixels into those bands 
        trace the chape of those bands and convert to vectors 
    """

    simg = Image.open( inputfile )
  
    dimg = ImageOps.posterize(simg, bits=1)
    dimg.save("%s%d.%s" % (COMMON_OUT,1,COMMON_EXT) )



##----------------------------------------------------
def rasterize(in_img, out_img):
    """ load a bitmap, do stuff and make another from it """
    
    # in_img 
    # out_img
    pass


##----------------------------------------------------

def generate_image(hres, yres, outputfile=None):
    
    color  = (0,0,255)
    color2 = (0,255,0)

    #fb = RasterObj() 
 
    fb = PixelOp()    
    

    #fb.load(input_img)
    fb.create_buffer(hres,yres)
    fb.fill_color(color)
  
    for p_x in range(hres):
        for p_y in range(hres):
            # print('## pixel is %s %s '%(p_x, p_y))

            if p_x < fb.center[0]:
                fb.set_pix( (p_x, p_y) , color)
            elif p_y <fb.center[1]:
                fb.set_pix( (p_x, p_y) , color2)

    if outputfile is None:
        return fb
    else:    
        fb.save(outputfile)
   


#fb = generate_image(115,115, "neu.bmp")








##----------------------------------------------------

