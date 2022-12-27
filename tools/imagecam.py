#from __future__ import print_function
""" 

requires 
   Python:
       scipy 
       numpy
       geojson 

   Bin:   
      potrace 

"""

##-- 
import imageio

import numpy as np
import scipy
import scipy.misc
import scipy.cluster


##-- 

import binascii
import struct

import sys
import subprocess

from PIL import Image, ImageFilter, ImageEnhance, ImageOps

from gnelscript.pygfx.math_ops import  *
from gnelscript.pygfx.raster_ops import *

from gnelscript.pygfx.point_ops_2d import *
from gnelscript.pygfx.point_ops import *

from gnelscript.pygfx.obj3d import  *
from gnelscript.pygfx.obj25d import  *


from gnelscript.pygfx.kicad_ops import *
from gnelscript.pygfx.gis_vector_ops import *



from gnelscript.examples.selection import * 
#from gnelscript.examples.milling import *
#from gnelscript.examples.raster import *
#from gnelscript.examples.render import *
#from gnelscript.examples.vector import *
#from gnelscript.examples.wip import *

mu = math_util() 


#you need to install and confifgure this tool 
potrace_command = "potrace"



COMMON_IN  = 'foobar'
COMMON_OUT = 'rasterfoo'
COMMON_EXT = 'png'



##----------------------------------------------------
##   /usr/local/opt/python@3.10/bin/python3.10 ./imagecam.py  
#firstpass_bw(10, 1.5, 1.5, 1, 250, "images/in/art.jpg", "images/out", "output")
## (iteration , blur , contrast, bright, scaling(divs) , in, out )
#firstpass(10, 0, 1, 1, 250, "images/in/oil.png", "images/out", "output")

#firstpass(10, 1.5, 1.2, 1, 250, "images/in/er.jpg", "images/out", "output")



##----------------------------------------------------
##   /usr/local/opt/python@3.10/bin/python3.10 ./imagecam.py  
#secondpass("images/out/er.bmp", "images/out" , 16, False)


##----------------------------------------------------
#set the RGB values from last tool and run this 
#thirdpass( "images/out/commonbands.png", "images/out", "dxf" )
#thirdpass( "images/out/commonbands.png",  "images/out" , "geojson")

#thirdpass( "images/out/commonbands.png",  "images/out" , "dxf")


##----------------------------------------------------

def firstpass_bw( iters, blur, contrast, bright, chops, inputfile, outputfolder, outputfile ):
    """
        UNTESTED 
        try to break a BW bitmap up into the most basic shapes of color regions 

    """
    simg = Image.open( inputfile )

    simg = simg.convert('L') # convert 8 bit

    do_posterize = False 
    do_potrace = False 

    width = int(simg.width/chops)
    height = int(simg.height/chops) 


    for i in range(iters):
        
        ### blur pass
        simg = simg.filter( ImageFilter.GaussianBlur(radius=blur) )
        
        ### contrast pass 
        cont_en = ImageEnhance.Contrast(simg)
        simg = cont_en.enhance(contrast)
        
        ### bright pass  
        bright_en = ImageEnhance.Brightness(simg)
        simg = bright_en.enhance(bright)

        copy = simg.convert("RGBA")        
        copy.save( "%s/%s_%d.%s"%(outputfolder, outputfile,i,"bmp") )


##----------------------------------------------------

def firstpass( iters, blur, contrast, bright, chops, inputfile, outputfolder, outputfile ):
    """
        filter out noise of an image but keep the main shapes and color regions 

    """
    simg = Image.open( inputfile )

    width = int(simg.width/chops)
    height = int(simg.height/chops) 

    for i in range(iters):
        
        ### blur pass
        simg = simg.filter( ImageFilter.GaussianBlur(radius=blur) )
        
        ### contrast pass 
        cont_en = ImageEnhance.Contrast(simg)
        simg = cont_en.enhance(contrast)
        
        ### bright pass  
        bright_en = ImageEnhance.Brightness(simg)
        simg = bright_en.enhance(bright)
        
        simg.save( "%s/%s_%d.%s"%(outputfolder, outputfile,i,"bmp") )

    # posterfile = "%s/%s_%d.%s"%(outputfolder, outputfile,1,"bmp")
    # if do_posterize:
    #     simg = ImageOps.posterize(simg, bits=1)
    #     simg.save( posterfile)



##----------------------------------------------------

def secondpass(inputimage, outputpath, numbands, fast=False):
    """
    #from stack overflow  - get most common colors in an image 
    
    """


    im = Image.open(inputimage )

    if fast:
        im = im.resize((300, 300))      # optional, to reduce time

    ar = np.asarray(im)
    shape = ar.shape
    ar = ar.reshape(scipy.product(shape[:2]), shape[2]).astype(float)

    codes, dist = scipy.cluster.vq.kmeans(ar, numbands)

    colorfile = '%s/commonbands.txt'%outputpath
    print("writing color bands %s"%colorfile)

    f = open(colorfile, "a")

    #write a file - make a list of RGB ints to read in later 
    for i,clr in enumerate(codes):
        f.write("%s %s %s %s \n"%(i, int(clr[0]), int(clr[1]), int(clr[2])) )
    f.close()

    print(codes)

    vecs, dist = scipy.cluster.vq.vq(ar, codes)         # assign codes

    # counts, bins = scipy.histogram(vecs, len(codes))    # count occurrences
    # index_max = scipy.argmax(counts)                    # find most frequent
    # peak = codes[index_max]
    # colour = binascii.hexlify(bytearray(int(c) for c in peak)).decode('ascii')
    # print('most frequent is %s (#%s)' % (peak, colour))

    c = ar.copy()
    for i, code in enumerate(codes):
        c[scipy.r_[scipy.where(vecs==i)],:] = code

    print("writing file ", '%s/commonbands.png'%outputpath) 
    imageio.imwrite('%s/commonbands.png'%outputpath, c.reshape(*shape).astype(np.uint8))


##----------------------------------------------------

def thirdpass( inputfile, outputfolder, fileformat, bmpinvert=False, po_invert=False, fastmode=False  ):
    """ break an already posterized image into seperate images X colors """

    simg = pixel_op()
    simg.load( inputfile )

    width  = int( simg.res_x )
    height = int( simg.res_y ) 
    
    tsize = 10 #suppress "turd" size - speckles 
    chops = 100

    vwidth = int(width/chops)
    vheight = int(height/chops) 

    
    #fileformat = "dxf" 
    #fileformat = "geojson" 

    ##populate this from the output of second pass to get the five best tasty colors I know 
    ## colors= [  ["a", (78,27,40)],
    ##            ["b", (163,91,94)],
    ##            ["c", (14,4,12)],
    ##            ["d", (254,246,241)],
    ##            ["e", (243,162,158)],
    ##          ]

    colorfile = '%s/commonbands.txt'%outputfolder

    print("loading rgb commonbands file %s"%colorfile )

    colors = [] 

    f = open(colorfile, "r")
    for x in f:
        l = x.split(' ')
        colors.append( [l[0], (int(l[1]), int(l[2]), int(l[3]) ) ] ) 

    for c in colors:
        
        #extract_by_color( path, name, color, slowmode, exactcolor, invert, framebuffer)
        simg.extract_by_color( outputfolder, c[0], c[1], False, False, bmpinvert, False )

        if po_invert:
            command = [potrace_command, "%s/%s.bmp"%(outputfolder, c[0] ), "-i", "-b", fileformat, "-W", str(vwidth), "-H", str(vheight), "-t", str(tsize)]        
        else:    
            command = [potrace_command, "%s/%s.bmp"%(outputfolder, c[0] ), "-b", fileformat, "-W", str(vwidth), "-H", str(vheight), "-t", str(tsize)]
        #print(command)
        subprocess.run(command)


##----------------------------------------------------


def geojson_to_ngc(folder, fnames, onefile=False):
    kiparser = generic_ngc()

    ##-- 

    #indexer(ids=None, span=None, unique=True, nth=None)
    ids = kiparser.indexer(span=[1,2])

    ##-- 

    passnum = 1

    sort = True
    doscale=False 

    #kiparser.load_geojson('images/out/%s.json'%passnum, 0, getfids=None, getids=None)

    if onefile:
        for file in fnames:
            kiparser.load_geojson('%s/%s'%(folder,file), 0, getfids=None, getids=None)

        if sort:
            #DEBUG BROKEN  - see below 
            pass 

        ##--
        #scale  
        if doscale: 
            gs = kiparser.global_scale
            xformed = []
            for ply in kiparser.gr_polys:
                xformed.append(kiparser.trs_points( ply, translate=(0,0,0), rotate=(0,0,0), scale=(gs,gs,gs) ))
            kiparser.gr_polys = xformed

        kiparser.calculate_paths()
        #kiparser.save_line_obj('3d_obj/foo.obj')
        kiparser.export_ngc("images/out/pass%s.ngc"%passnum)

    #export each in a loop
    else:                     
        for file in fnames:
            fspl = file.split('.')
    
            kiparser.load_geojson('%s/%s'%(folder,file), 0, getfids=None, getids=None)
            if sort:
                kiparser.index_grsort()
                kiparser.export_sorted_centroids(fspl[0], folder)
                kiparser.export_sorted_extents(fspl[0], folder)
                kiparser.show_buffers()

                kiparser.make_grid(fspl[0], folder, 5, 5)

            ##--
            #scale  
            if doscale: 
                gs = kiparser.global_scale
                xformed = []
                for ply in kiparser.gr_polys:
                    xformed.append(kiparser.trs_points( ply, translate=(0,0,0), rotate=(0,0,0), scale=(gs,gs,gs) ))
                kiparser.gr_polys = xformed
            kiparser.calculate_paths()
            #kiparser.save_line_obj('3d_obj/foo.obj')
            kiparser.export_ngc("%s/%s.ngc"%(folder, fspl[0]) )

    ##--
    #bbox = kiparser.calc_bbox_pt(2, (5,5))
    #pts = kiparser.cvt_2d_to_3d(kiparser.extents_fr_bbox(bbox, periodic=True))
    #debug - need to solve the clean_pts_str debacle?
    #kiparser.gr_polys.append(pts)

    ##--
    #bbox = kiparser.calc_bbox_pt(1.75, (-3,3))
    #pts = kiparser.cvt_2d_to_3d(kiparser.extents_fr_bbox(bbox, periodic=True))
    #debug - need to solve the clean_pts_str debacle?    
    #kiparser.gr_polys.append(pts)

    ##--
    #print(pts)
    #kiparser.grply_inspect()
    #kiparser.cvt_grpoly_obj3d()
    #kiparser.save("3d_obj/foo.obj")

    ##-- 

    # ADD BETTER TOOLS TO XFROM POINTS 

    # ADD BETTER TOOLS OMIT POLYGONS  

    # TODO ADD A SPATIAL POLYGON SORT TO REDUCE SPINDLE TRAVEL ON RETRACTS 





























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




    


