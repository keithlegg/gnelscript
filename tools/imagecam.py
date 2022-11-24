from __future__ import print_function
import binascii
import struct
from PIL import Image
import numpy as np
import scipy
import scipy.misc
import scipy.cluster

###

import sys

from PIL import Image, ImageFilter, ImageEnhance, ImageOps

from gnelscript.pygfx.math_ops import  *
from gnelscript.pygfx.raster_ops import *
from gnelscript.pygfx.point_ops_2d import *

#from pygfx.point_ops import *
#from pygfx.obj3d import  *


mu = math_util() 

potrace_command = "WHERE_YOU_PUT/potrace"



COMMON_IN  = 'foobar'
COMMON_OUT = 'rasterfoo'
COMMON_EXT = 'png'


##----------------------------------------------------


def firstpass( iters, chops, inputfile, outputfolder, outputfile ):
    """
        try to break a bitmap up into the most basic shapes of color regions 

    """
    import subprocess
     
    # max_zones = 5 
    # luminance = .5 

    simg = Image.open( inputfile )

    blur = 2.0 
    contrast = 1.2 
    bright = 1.1 

    do_posterize = False 
    do_potrace = False 
    
    tsize = 10 #suppress "turd" size - speckles 

    
    posterfile = "%s/%s_%d.%s"%(outputfolder, outputfile,1,"bmp")


    width = int(simg.width/chops)
    height = int(simg.height/chops) 

    #simg = simg.filter(ImageFilter.GaussianBlur(radius = blur))
    ##contrast pass 
    #cont_en = ImageEnhance.Contrast(simg)
    #simg = cont_en.enhance(2)

    for i in range(iters):
        
        ### blur pass
        simg = simg.filter( ImageFilter.GaussianBlur(radius=blur) )
        
        ### contrast pass 
        cont_en = ImageEnhance.Contrast(simg)
        simg = cont_en.enhance(contrast)
        
        ### bright pass  
        #bright_en = ImageEnhance.Brightness(simg)
        #simg = bright_en.enhance(.8)
        
        simg.save( "%s/%s_%d.%s"%(outputfolder, outputfile,i,"bmp") )

    if do_posterize:
        simg = ImageOps.posterize(simg, bits=1)
        simg.save( posterfile)

    if do_potrace:
        
        #"-i" -invert 

        command = [potrace_command, posterfile, "-b", "svg", "-W", str(width), "-H", str(height), "-t", str(tsize)]
        
        print(command)

        subprocess.run(command)



def secondpass(inputimage, outputimage):
    #from stack overflow  - attempt to get most common colors in an image 

    NUM_CLUSTERS = 5

    im = Image.open(inputimage )

    #im = im.resize((150, 150))      # optional, to reduce time

    ar = np.asarray(im)
    shape = ar.shape
    ar = ar.reshape(scipy.product(shape[:2]), shape[2]).astype(float)

    codes, dist = scipy.cluster.vq.kmeans(ar, NUM_CLUSTERS)
    print('cluster centres:\n', codes)

    vecs, dist = scipy.cluster.vq.vq(ar, codes)         # assign codes

    # counts, bins = scipy.histogram(vecs, len(codes))    # count occurrences
    # index_max = scipy.argmax(counts)                    # find most frequent
    # peak = codes[index_max]
    # colour = binascii.hexlify(bytearray(int(c) for c in peak)).decode('ascii')
    # print('most frequent is %s (#%s)' % (peak, colour))

    import imageio
    
    c = ar.copy()
    for i, code in enumerate(codes):
        c[scipy.r_[scipy.where(vecs==i)],:] = code
    imageio.imwrite('%s_commonbands.png'%outputimage, c.reshape(*shape).astype(np.uint8))


def thirdpass( inputfile, outputfolder  ):
    """ break an already posterized image into seperate iamges X colors """

    import subprocess

    simg = PixelOp()
    simg.load( inputfile )

    width  = int( simg.res_x )
    height = int( simg.res_y ) 
    
    tsize = 10 #suppress "turd" size - speckles 
    chops = 100

    vwidth = int(width/chops)
    vheight = int(height/chops) 

    ## colors = [ 
    ##        ["black" ,  (  0  , 0  , 0   ) ],
    ##        ["red"   ,  (  255, 0  , 0   ) ],
    ##        ["green" ,  (  0  , 255, 0   ) ],
    ##        ["blue"  ,  (  0  , 0  , 255 ) ],
    ##        ["white" ,  (  255, 255, 255 ) ],
    ##        ]

    colors= [ ["a", (0,0,0)],
               ##["b", (7,3,2)],
               ##["c", (94,96,82)],
               ##["d", (98,79,7)],
               ["e", (255,255,255)],
             ]

    for c in colors:
        simg.extract_by_color( outputfolder, c[0], c[1] )

        command = [potrace_command, "%s/%s.bmp"%(outputfolder, c[0] ), "-b", "dxf", "-W", str(vwidth), "-H", str(vheight), "-t", str(tsize)]
        
        print(command)

        subprocess.run(command)




##----------------------------------------------------


#firstpass(10, "images/stop.webp", "poster", "postergirl")
#firstpass(10, "images/d.jpg",     "poster", "postergirl")
#firstpass(2, "images/refer.jpg",  "poster", "postergirl")


#firstpass(10, 250, "images/me2.jpg", "poster", "postergirl")

# ( /usr/local/opt/python@3.10/bin/python3.10 ./imagecam.py ) 

#secondpass("images/me4.png", "poster" )

##


simg = PixelOp()
color = [(255,255,255)]

#print( simg.closest_color(color, (128,128,128))) 


#print( simg.color_distance((0,100,80), (0,100,100))) 




































##----------------------------------------------------



# Text for the file header, the parameter is the name of the module, ex "LOGO".
header = """(module %(name)s
  (layer F.Cu)
  (at 0 0)
  (fp_text reference "VAL***" (at 0 10) (layer F.SilkS) hide
    (effects (font (thickness 0.3)))
  )
  (fp_text reference "%(name)s" (at 0 5) (layer F.SilkS) hide
    (effects (font (thickness 0.3)))
  )
"""


# Text for the file footer, the only parameter is the name of the module
footer = """)
"""


# Places a pixel with (x, y) at the upper-left corner
# Size is in units of mm
def make_pixel(x, y, px_size):
    return """  (fp_poly
    (pts
      (xy %(0)s %(1)s)
      (xy %(2)s %(1)s)
      (xy %(2)s %(3)s)
      (xy %(0)s %(3)s)
      (xy %(0)s %(1)s)
    )
    (layer F.SilkS)
    (width 0.01)
  )
""" % {"0": x, "1": y, "2": x + px_size, "3": y + px_size}


def conv_image_to_module(image, module_name, scale_factor):

    """ Returns the text for the module, and the size: (x, y) in mm. """
    w, h = image.size
    #print "Original image dimensions: {0} x {1}".format(w, h)
    #print "Writing module file to \"{0}\"".format(output_filename)
    module = header % {"name": module_name}
    for y in range(h):
        for x in range(w):
            #print image.getpixel((x,y))
            if image.getpixel((x, y)) == 0:
                module += make_pixel(scale_factor * x, scale_factor * y, scale_factor)
    module += footer % {"name": module_name}
    return module, (scale_factor * w, scale_factor * h)






##----------------------------------------------------

def vectorizer(inputfile, outputfile=None):
    """ example of PIL edge detect """
    output = "poster"

    # Opening the image (R prefixed to string
    # in order to deal with '\' in paths)
    image = Image.open(inputfile)

    # Converting the image to grayscale, as edge detection 
    # requires input image to be of mode = Grayscale (L)
    image = image.convert("L")

    # Detecting Edges on the Image using the argument ImageFilter.FIND_EDGES
    image = image.filter(ImageFilter.FIND_EDGES)

    # Saving the Image Under the name Edge_Sample.png
    image.save("%s%d.%s" % ("%s/%s"%(output,"edge_detect"),1,COMMON_EXT) )






##----------------------------------------------------

def vectorizer2(inputfile, outputfile=None):
    """ based on kicad img convert """
    output = "poster"

    input_image = "images/dot.png"                              #sys.argv[1]
    module_name = "foobar"  #sys.argv[3]
    scale_factor = int(2)

    print("Reading image from \"%s\"" % input_image)

    # If we .convert("1") it makes a binary image with a 127 threshold
    # Keeping it at greyscale allows us to make the threshold configurable
    image = Image.open(input_image).convert("L")
    # The above conversion will assume a black background for removing transparency.
    # If we want to use a different background color, use this:
    #white = Image.new("RGB", image.size, (255,255,255))
    #r,g,b,a = image.split()
    #image = Image.composite(image, white, a)

    # If you want to invert the image:
    #image = ImageOps.invert(image)

    # If you want to do non-127 thresholding, change the 127 below
    image = image.point(lambda i: 0 if i < 127 else 255)

    module, size = conv_image_to_module(image, module_name, scale_factor)
    print ("Output image size: %f x %f mm" % (size[0], size[1]) )

    fid = open(outputfile, "w")
    fid.write(module)
    fid.close()




    


