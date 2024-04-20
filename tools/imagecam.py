#from __future__ import print_function
""" 


extract data from images and turn it into vectors.


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
from gnelscript.pygfx.obj2d import  *


from gnelscript.pygfx.kicad_ops import *
from gnelscript.pygfx.vector_ops import *
from gnelscript.pygfx.milling_ops import *



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
        
        if bright:
            ### bright pass  
            bright_en = ImageEnhance.Brightness(simg)
            simg = bright_en.enhance(bright)

        copy = simg.convert("RGBA")        
        copy.save( "%s/%s_%d.%s"%(outputfolder, outputfile,i,"bmp") )


##----------------------------------------------------
def firstpass( iters, blur, contrast, bright, chops, inputfile, outputfolder, outputfile, rotate=None ):
    """
        filter out noise of an image but keep the main shapes and color regions 
        
        todo:
            add max size on long edge option 
    """

    simg = Image.open( inputfile )

    width = int(simg.width/chops)
    height = int(simg.height/chops) 

    if rotate:
        simg = simg.rotate(rotate, expand=True)

    for i in range(iters):
        
        ### blur pass
        simg = simg.filter( ImageFilter.GaussianBlur(radius=blur) )
        
        ### contrast pass 
        cont_en = ImageEnhance.Contrast(simg)
        simg = cont_en.enhance(contrast)

        if bright:        
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

    TODO:
        warn if numbands is too high? does it matter?
        warn if image is too large? 
        resize if too large?

    
    """
    print("calculating %s common colors. may take some time, especially on large images"%numbands)

    im = Image.open(inputimage )

    if fast:
        im = im.resize((fast, fast))      # optional, to reduce time

    ar = np.asarray(im)
    shape = ar.shape
    ##DEBUG - will crash here if num bands exceedes what is in image data 
    ar = ar.reshape(scipy.product(shape[:2]), shape[2]).astype(float)

    #sort colors  WRONG 
    #np.sort(ar)

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
def thirdpass( inputfile, outputfolder, fileformat, bmpinvert=False, po_invert=False, fastmode=False, blurrad=None  ):
    """ break an already posterized image into seperate images X colors 

        
        DEBUG NEED TO ADD THESE:
        
        #tsize = 10    - suppress "turd" size - speckles 
        #chops = 100   - scale 
 
        inputfile      -
        outputfolder   -
        blurrad        - bitmap blur (prior to tracing vectors)
        fileformat     - format of vctor output (dxf,json) 
                         if no format - dont run potrace - comic mode only uses BMP files 
        bmpinvert      -

        po_invert      - not working?
        fastmode       -

    """

    simg = pixel_op()
    simg.load( inputfile )

    width  = int( simg.res_x )
    height = int( simg.res_y ) 
    
    tsize = 10 #suppress "turd" size - speckles 
    chops = 100

    vwidth = int(width/chops)
    vheight = int(height/chops) 

    ## format of input - five best tasty colors I know 
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

        #if no format - dont run potrace - comic mode only uses BMP files 
        if fileformat:
            if po_invert:
                command = [potrace_command, "%s/%s.bmp"%(outputfolder, c[0] ), "-i", "-b", fileformat, "-W", str(vwidth), "-H", str(vheight), "-t", str(tsize)]        
            else:    
                command = [potrace_command, "%s/%s.bmp"%(outputfolder, c[0] ), "-b", fileformat, "-W", str(vwidth), "-H", str(vheight), "-t", str(tsize)]
            #print(command)
            subprocess.run(command)
        else:
            print('thirdpass: POTRACE DISABLED')

##----------------------------------------------------
def comic_pass( infolder, names, outfolder, colorbg_file, lineweight, bluramt, brightamt, blacklines=True ):
    """ trace edges and merge with color background 
    """

    edge_cases = [] 

    rop = raster_op()

    # get size we are working with - assume they are all same size 
    infile = '%s/%s.bmp'%(infolder,names[0])
    image = Image.open(infile)    
    res_x = image.size[0]
    res_y = image.size[1] 

    ##-- 

    ## iterate list of BW files to edge detect and composite back onto color image 
    for fname in names:
        infile = '%s/%s.bmp'%(infolder,fname)

        image = Image.open(infile)
     
        # Converting the image to grayscale, as edge detection
        # requires input image to be of mode = Grayscale (L)
        image = image.convert("L")
     
        # Detecting Edges on the Image using the argument ImageFilter.FIND_EDGES
        image = image.filter(ImageFilter.FIND_EDGES)
      
        #invert to make black edges 
        inverted_image = ImageOps.invert(image)

        #image.save('%s/%s.bmp'%(outfolder,fname))
        edgefile='%s/%s_edge.bmp'%(outfolder,fname)
        edge_cases.append(edgefile)
        image.save(edgefile)        
        #invedgefile='%s/%s_inv_edge.bmp'%(outfolder,fname)
        #inverted_image.save(invedgefile)


    ##----
    # merge all the edge files into one
    # THIS IS NOT RIGHT , BUT IT IS INTERESTNG - IT ACCUMULATES LINE WEIGHT 

    if blacklines:
        alledges = Image.new('RGBA', [res_x,res_y], (255,)*4)
    else:
        alledges = Image.new('RGBA', [res_x,res_y], (0,)*4)

    allwhite = Image.new('RGBA', [res_x,res_y], (255,)*4) # white mask to composite
    allblack = Image.new('RGBA', [res_x,res_y], (0,)*4)   # black mask to color lines where mask shows through 

    
    weighted = True

    if weighted:
        for i,e in enumerate(edge_cases):
            tmp = Image.open(e)         
            tmp = tmp.convert("RGBA")
            alledges =  Image.blend(tmp, alledges, .5) 
            #alledges = Image.combine((tmp, alledges), lambda pixels: max(pixels))
            #r, g, b, a = tmp.split()
            #alledges = Image.merge("RGB", (rr, gg, bb))
        alledges = alledges.convert("L")
        alledges = ImageOps.invert(alledges)


    #attempt to get even line weight (all egdes combined in one)
    #we get a near similar effect below if we skip blur and overdrive contrast
    else:
        ro = raster_op()
        ro.create_buffer(res_x,res_y)
        ro.bitmode = 'RGB'
        alledges = ro.empty_buffer( (255,255,255) )
        outpix = alledges.load() 
        for i,e in enumerate(edge_cases):
            image = Image.open(edge_cases[i])
            inpix = image.load() 
            for x in range(res_x):
                for y in range(res_y):
                    #edges are inverted so white == black
                    if inpix[x,y]>240:
                        outpix[x,y]=0 

        alledges = alledges.convert("L")
        alledges = ImageOps.invert(alledges)

    ########### 
    ## process lines   
 
    if lineweight:
        ## blur lines pass
        if bluramt:
            alledges = alledges.filter( ImageFilter.GaussianBlur(radius=bluramt) )
        ## contrast lines 
        cont_en = ImageEnhance.Contrast(alledges)
        alledges = cont_en.enhance(lineweight)


    alledges.save("%s/alledges.bmp"%(outfolder))

    ## bright pass  
    #bright_en = ImageEnhance.Brightness(simg)
    #simg = bright_en.enhance(bright)

    #cheapo dilation
    #alledges = alledges.filter(ImageFilter.MaxFilter(3))

    #cheapo erode
    #alledges = alledges.filter(ImageFilter.MinFilter(1))

    ## end process lines
    ########### 

    #if blacklines:
    #    alledges = alledges.convert("L")
    #    alledges = ImageOps.invert(alledges)
 

    ##----
    # now composite them onto the color background 
    background = '%s/%s'%(infolder,colorbg_file)
    colorbg = Image.open(background) 

    if brightamt:
        ## bright pass color background  
        bright_en = ImageEnhance.Brightness(colorbg)
        colorbg = bright_en.enhance(brightamt)

    final = Image.new('RGBA', [res_x,res_y]) 
    final = Image.composite(colorbg, allblack, alledges) 
    final.save("%s/cartoon.bmp"%(outfolder))

    ## FIRST STAB - composite using white edges - black BG - it works, but only for each layer at a time ##
    # background = '%s/%s'%(infolder,colorbg_file)
    # colorbg = Image.open(background) 
    # final = Image.new('RGBA', colorbg.size) 
    # for i,e in enumerate(edge_cases):
    #     edge_mask = Image.open(e) 
    #     final =  Image.composite(allblack, colorbg, edge_mask) 
    #     final.save("%s/foofu_%s.bmp"%(outfolder,i))

##----------------------------------------------------
def redraw(  infolder, colorbg_file, lineweight, bluramt, brightamt, blacklines=True ):
    """ trace edges 
        and merge with external lines , generated ny hand or from comic_pass() 

        todo:
            run a second iteration:
               1> take commonbands and run through firstpass
               2> composite those regions with the BW masks generated in thirdpass()  

    """

    edge_cases = [] 

    # get size we are working with - assume they are all same size 
    aefile = '%s/alledges.bmp'%(infolder)
    alledges = Image.open(aefile)    
    res_x = alledges.size[0]
    res_y = alledges.size[1] 

    ########### 
    ## process lines   
    
    if lineweight:
        ## blur lines pass
        alledges = alledges.filter( ImageFilter.GaussianBlur(radius=bluramt) )
        ## contrast lines 
        cont_en = ImageEnhance.Contrast(alledges)
        alledges = cont_en.enhance(lineweight)

    #allwhite = Image.new('RGBA', [res_x,res_y], (255,)*4) # white mask to composite
    allblack = Image.new('RGBA', [res_x,res_y], (0,)*4)   # black mask to color lines where mask shows through 

    ## bright pass  
    #bright_en = ImageEnhance.Brightness(simg)
    #simg = bright_en.enhance(bright)

    #cheapo dilation
    #alledges = alledges.filter(ImageFilter.MaxFilter(3))

    #cheapo erode
    #alledges = alledges.filter(ImageFilter.MinFilter(1))

    ## end process lines
    ########### 

    ##----
    # now composite them onto the color background 
    background = '%s/%s'%(infolder,colorbg_file)
    colorbg = Image.open(background) 

    if brightamt:
        ## bright pass color background  
        bright_en = ImageEnhance.Brightness(colorbg)
        colorbg = bright_en.enhance(brightamt)

    final = Image.new('RGBA', [res_x,res_y]) 
    final = Image.composite(colorbg, allblack, alledges) 
    final.save("%s/cartoon.bmp"%(infolder))


##----------------------------------------------------
##----------------------------------------------------
##----------------------------------------------------
##----------------------------------------------------


def pcb_to_tesselated(path, inpcb):
    pcb = pcbfile()
    vflo = vectorflow()


    # zval, path, infname):
    pcb.load_kicadpcb(0, path, inpcb)
    o = object3d() 

    pts = o.scale_pts(.05, pts=pcb.gr_polys[0])
    
    vflo.gr_polys = [pts] 
    vflo._sort() 
    vflo.gl_move_center()
    
    seed = vflo.gr_sort[0][4]
    #(scale, numx, numy, numrots, seedpts, folder, injson, outjson):
    tesselate_json(4, 10, 10, 3, seed, path, 'warped.json', 'tesslated.json')

    #vflo.export_geojson_polygon(path, 'kicad')


##----------------------------------

def tesselate_json(scale, numx, numy, numrots, seedpts, folder, injson, outjson):
    """scan 2D polygons, look for 4 sided only, bisect the edges  
    """
    do_export = False  

    vflo = vectorflow()

    #builds a DAG node for each centroid of all polys - (QUAD polys only)
    if False:
        vflo.load_geojson( '%s/%s'%(folder, injson) )
        vflo.cvt_grsort_todag()

    #or build a grid of nodes 
    if True:
        vflo.tesl.minx = -1 
        vflo.tesl.miny = -1 
        vflo.tesl.maxx = 1 
        vflo.tesl.maxy = 1 
        vflo.tesl.build_2d_cells( numx, numy, scale=scale)
    
    if False:
        ## example of how to add a single cell    
        ## DEBUG - add a function to do this and procederally set the attrs from data passed in 
        vflo.tesl.new_cell_2d('bob', 
                              2, 2,
                              0, 0, 0, 
                              0, 0, 0)

    #procedurally build some basic geom in the cells     
    #vflo.build_tesselation_sample()

    #pts = [(-.1,.3,0), (.3,1,0), (.1,-.3,0), (-.1,-.3,0)]
    vflo.build_tesselation_test2(numrots, seedpts)

    vflo.gr_polys = [] 
    vflo.gr_sort = [] 

    vflo.cvt_tessl_grsort()

    #flatten the polygons into some bad topology 
    points=[]
    for ply in vflo.gr_sort:
        points.extend(ply[4])
    

    ##
    #convert geom to shapely so we can run buffer on it 
    v2 = vectorflow()
    v2.gr_polys.append(points)
    v2._sort()
    fc = v2.cvt_grsort_shapely()

    #buffer with 0 distance can "repair" topology
    bply = v2.shapely_buffer(fc[0],0)

    #clear any existing geom 
    v2.gr_polys=[]
    v2.gr_sort=[]

    v2.gr_polys.append(v2.cvt_2d_to_3d(bply))
    v2._sort()
    

    #build 3D geom from the points
    if True:
        v2.cvt_grpoly_obj3d()
        v2.vec_connect_pts(pts=v2.cvt_2d_to_3d(bply), axis='x', draw_obj='rect_2d')
        #v2.vec_connect_pts(pts=v2.cvt_2d_to_3d(bply), axis='z', draw_obj='arrow')

    #buffered (fixed) geom    
    if True:
        v2.export_geojson_lines(  folder, 'fixed')
        v2.export_geojson_polygon(  folder, 'fixed')

    #original broken geom
    if False:
       vflo.export_geojson_lines(  folder, outjson)
       vflo.export_geojson_polygon(  folder, outjson)



##----------------------------------------------------
def tesselation_examples(scale, numx, numy, numrots, seedpts, folder, injson, outjson):
    """scan 2D polygons, look for 4 sided only, bisect the edges  
    """
    do_export = False  

    vflo = vectorflow()

    #builds a DAG node for each centroid of all polys - (QUAD polys only)
    if False:
        vflo.load_geojson( '%s/%s'%(folder, injson) )
        vflo.cvt_grsort_todag()

    #or build a grid of nodes 
    if True:
        vflo.tesl.minx = -1 
        vflo.tesl.miny = -1 
        vflo.tesl.maxx = 1 
        vflo.tesl.maxy = 1 
        vflo.tesl.build_2d_cells( numx, numy, scale=scale)
    
    if False:
        ## example of how to add a single cell    
        ## DEBUG - add a function to do this and procederally set the attrs from data passed in 
        vflo.tesl.new_cell_2d('bob', 
                              2, 2,
                              0, 0, 0, 
                              0, 0, 0)

    #procedurally build some basic geom in the cells     
    #vflo.build_tesselation_sample()

    #pts = [(-.1,.3,0), (.3,1,0), (.1,-.3,0), (-.1,-.3,0)]
    vflo.build_tesselation_test2(numrots, seedpts)

    vflo.gr_polys = [] 
    vflo.gr_sort = [] 

    vflo.cvt_tessl_grsort()

    #flatten the polygons into some bad topology 
    points=[]
    for ply in vflo.gr_sort:
        points.extend(ply[4])
    
    ####

    v2 = vectorflow()
    v2.gr_polys.append(points)
    v2._sort()
    fc = v2.cvt_grsort_shapely()

    #buffer with 0 distance can "repair" topology
    bply = v2.shapely_buffer(fc[0],0)

    v2.gr_polys=[]
    v2.gr_sort=[]
    v2.gr_polys.append(v2.cvt_2d_to_3d(bply))
    v2._sort()

    #buffered (fixed) geom    
    if True:
        v2.export_geojson_lines(  folder, 'fixed')
        v2.export_geojson_polygon(  folder, 'fixed')

    #original broken geom
    if False:
       vflo.export_geojson_lines(  folder, outjson)
       vflo.export_geojson_polygon(  folder, outjson)


##----------------------------------------------------


def make_rect(sx, sy, path, name): 
    vflo = vectorflow()
    vflo.prim_rect(axis='z', sizex=sx, sizey=sy, sizez=1, periodic=True)
    vflo.cvt_obj3d_grpoly()
    vflo.export_geojson_polygon(path, name)
    #vflo.export_geojson_lines(path, name)
    
    # export_ngc(  rh, ch, cdpi, cmax, filename, do3d=False, do_retracts=True, do_laser=False, laserpwm=400,  do_gpio=False):
    vflo.export_ngc(1, 0, .1   ,    2, '%s/%s.ngc'%(path, name), do_laser=False, do3d=False, do_retracts=False, do_gpio=0)  


#make_rect(GLOBAL_PROJ, 'quad')


##----------------------------------------------------
def move_json(path, name): 
    vflo = vectorflow()
    vflo.load_geojson('%s/%s.json'%(path,name))
    #vflo.gl_translate(0,-.3,0)
    #vflo.gl_scale([.9,.9,.9])
    vflo.export_geojson_polygon(path, '%s.json'%name)


#move_json(GLOBAL_PROJ, 'spacewarp')

##----------------------------------------------------

def scale_fit_test(path, infile, outfile):
    vflo = vectorflow()
    vflo.load_geojson(path+'/'+infile)

    sc = vflo.gl_scale_to_fit(4,4)
    vflo.gl_scale( [sc[0], sc[1], 1] )

    vflo.gr_sort = vflo.filter_by_bbox(.2,.2)

    #vflo.gl_move_extents_corner(which='bl')
    vflo.export_extents_ngc(path, outfile)
    vflo.export_ngc(1, 0, .1, 2, '%s/%s.ngc'%(path, outfile), do_laser=False, do3d=False, do_retracts=False, do_gpio=0)  


##----------------------------------------------------
def geojson_to_ngc(folder, fnames, onefile=False):
    kiparser = vectorflow()

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

        #if sort:
            #DEBUG BROKEN  - see below 

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
                kiparser.export_sorted_centroids(folder, fspl[0])
                kiparser.export_sorted_extents(folder, fspl[0])
                #kiparser.show_buffers()
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
    #bbox = kiparser.calc_2d_bbox_pts(2, (5,5))
    #pts = kiparser.cvt_2d_to_3d(kiparser.extents_fr_bbox(bbox, periodic=True))
    #debug - need to solve the clean_pts_str debacle?
    #kiparser.gr_polys.append(pts)

    ##--
    #bbox = kiparser.calc_2d_bbox_pts(1.75, (-3,3))
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


##------------------------------------------------------------




class streamline(object):
    """ you can call each pass manually or use this tool to try and do it all at once 

        TODO:
            make this a "script engine" for processing images 
            each art project will get a script to tell it how to work 

    """

    def __init__(self):
        self.vflo = vectorflow()

        self.p1_settings = [7, 0, 1.1, 1.1, -99]
        self.p2_settings = []
        self.p3_settings = []
        
        self.inputimg = ''
        self.input_folder = ''
        self.output_folder = ''

    ###########

    def load_json(self, jsonfile):
        """WIP 
           test to work with some json data 
           we want this to become scriptable later on 
        """

        #def load_geojson(self, inputfile, zaxis, getfids=None, getids=None):
        self.vflo.load_geojson(jsonfile)
          
                
        ## global scale the data  
        # self.vflo.gl_scale()

        ## global centroid  
        # print(self.vflo.gl_centroid())


        #this will dump the geometry into an object3d instance 
        self.vflo.cvt_grpoly_obj3d()

        self.vflo.show()
        #print(len(self.vflo.points)) 



        # kip.load_geojson('%s/images/out/%s.json'%(GLOBAL_PROJ, 'vselect') )
        # kip.show_buffers()
        # kip.show_setup()
        # kip.export_global_extents('%s/images/'%GLOBAL_PROJ, 'foo_extents.ngc', ngc=True)
        # kip.export_global_extents('%s/images/'%GLOBAL_PROJ, 'foo_extents.json', ngc=False)
        # kip.export_sorted_centroids('%s/images/'%GLOBAL_PROJ, 'kentroidz.json')
        # kip.export_ngc(1, 0, .1, 2, '%s/images/out/%s.ngc'%(GLOBAL_PROJ, 'cowtri.ngc') , do3d=True)
        # print(kip.gr_sort[0])
        # kip.export_ngc(1, 0, .1, 2, '%s/images/out/%s'%(GLOBAL_PROJ, 'wargames.ngc') , do3d=False)
    
    ###########

    #def load_script(self, scriptpath):
    
    """
    def sample_script(self, scriptpath):
        out = []

        out.append('## sample worfarts script ')
        out.append(' ')
        out.append('## SETUP PARAMETERS  ')        
        out.append('set path $global /Users/klegg/serv/gnolmec/images')
        out.append(' ')

        out.append('## IMAGE $input ')                
        out.append('load image $input $global predator.jpg')
        out.append('frame image $input 10px')
        out.append('scale image $input ".8')
        out.append('rotate image $input -45')
        out.append('save image $input . color_input.jpg')

        f = open(scriptpath, "w")
        for line in out:
            f.write('%s\n'%line)
        f.close()

    def run(self, scriptfile):

        self.inputimg = infile
        path = Path(infile)
        self.input_folder =  path.parent.absolute()
        self.output_folder = outfolder
    
    """

    def oldrun(self, whichiter, infile, outfolder, outname, numbands=4, invbmp=False, invpo=False, dopass='all'):
        """whichiter
           infile
           outfolder
           outname
           numbands
           invbmp
           invpo
           dopass
        """

        self.inputimg = infile
        path = Path(infile)
        self.input_folder =  path.parent.absolute()
        self.output_folder = outfolder
        
        if dopass==1:                
            self._pass1(infile, outfolder, outname) 

        if dopass==2: 
            print('using iteration %s/%s_%s as input'%(outfolder, outname, whichiter) )   
            self._pass2( numbands, outname, outfolder, outname, whichiter) 

        if dopass==3:                
            #infilename, outfolder, mode='geojson',bmpinvert=False, po_invert=False
            self._pass3('%s/commonbands.png'%outfolder, outfolder, mode='geojson',bmpinvert=invbmp, po_invert=invpo)

        if dopass=='all':                
            self._pass1(infile, outfolder, outname) 
            self._pass2( numbands, outname, outfolder, outname, whichiter) 
            #infilename, outfolder, mode='geojson',bmpinvert=False, po_invert=False
            self._pass3('%s/commonbands.png'%outfolder, outfolder, mode='geojson',bmpinvert=invbmp, po_invert=invpo)

        if dopass=='skipfirst':
            #use the source file and dont try to cleanup first  
            print(infile, outfolder, outname) 

            self._pass2( numbands, infile, outfolder, outname, whichiter) 
            
            #infilename, outfolder, mode='geojson',bmpinvert=False, po_invert=False
            #self._pass3('%s/commonbands.png'%outfolder, outfolder, mode='geojson',bmpinvert=invbmp, po_invert=invpo)

    def _pass1(self, infile, outfolder, outname, dobw=False):
        """
            0 iterations   - number of times to repeat  
            1 blur         - pixels per step  
            2 contrast     - contrast per step 
            3 bright       - brightness per step 
            4 scaling(divs)  - UNIMPLEMENTED
            5 in            - full path to image 
            6 out           - output folder 
            7 out           - output filename 
        """
        args = self.p1_settings
        
        print('\n # running first pass on %s '%infile)

        if dobw:
            firstpass_bw(args[0], args[1], args[2], args[3], args[4], infile, outfolder, outname)
        else:
            firstpass(args[0], args[1], args[2], args[3], args[4], infile, outfolder, outname)

    def _pass2(self, numbands, infilename, outfolder, outname, which=None):
        """which is the iteration from first pass to use
           if none specified, determine how many there are and take the middle  
        """

        args = self.p2_settings
        
        #figure out which image - take the middle if none specified
        found = []
        
        usefile = ''

        # look in folder and count the source images
        # sloppy method to take the middle file
        count = 0
        for file in os.listdir(self.output_folder):
            if outname in file:
                count=count+1
                found.append(file)
        if which is None:        
            half = int(len(found)/2)
            for file in found:
                if str(half) in file:
                    usefile=file
        else:
            for file in found:
                if str(which) in file:
                    usefile=file           

        infile = '%s/%s'%(outfolder,usefile)
        print('# running secondpass on %s '%infile, outfolder)

        secondpass(infile, outfolder , numbands, fast=False )

    def _pass3(self, infilename, outfolder, mode='geojson',bmpinvert=False, po_invert=False):
        """
           0 inputfile
           1 outputfolder
           2 fileformat
           3 bmpinvert=False
           4 po_invert=False
        """

        print('#running third pass on %s '%(infilename))

        #set the RGB values from last tool and run this 

        #thirdpass( infilename, outfolder, mode, bmpinvert=True )
        thirdpass( infilename, outfolder, mode, bmpinvert=bmpinvert, po_invert=po_invert, fastmode=False )




##----------------------------------------------------




#obj_to_tesselation()

##----------------------------------------------------
def draw_test(folder, filename):
    vflo = vectorflow()
    vflo.prim_quad( sizex=7, sizey=4, axis='z') 
    vflo.cvt_obj3d_grpoly()
    vflo.export_geojson_polygon(folder, '%s.json'%filename)


##----------------------------------------------------
##----------------------------------------------------
def iso_flat_render(folder, rotation, objfile, outjson):
    #render by simply dropping the Z axis
    #Z SORT goes all to hell if any polygons happen to intersect
    vflo = vectorflow()
    vflo.load( '%s/3d_obj/%s'%(folder, objfile) )
    vflo.rotate_pts(rot=rotation)
    vflo.cvt_obj3d_grpoly()

    #vflo.export_ngc(1, 0, .1, 2, '%s/%s.ngc'%(folder, 'vvvxxbef'), do3d=False, do_retracts=False)
    #vflo.export_sorted_extents(folder,'vvvxxbef.json')
    #vflo.export_sorted_centroids(folder, 'centers.json')
    
    vflo.export_geojson_lines(  folder, outjson, periodic=False)
    vflo.export_geojson_polygon(folder, outjson)


##----------------------------------------------------

def vector_render_3dobj(folder, objname, outngc):
    #use the render object to make 2d path segments from a 3d object  
    vflo = vectorflow()
    
    #vflo.load('%s/3d_obj/spacetime.obj'%folder)
    pts = vflo.vec_render_obj( 35,35, 90, 1.0, '%s/3d_obj'%folder, objname)

    for pt in pts:
        vflo.gr_polys.append(vflo.cvt_2d_to_3d(pt))
    vflo._sort()
    
    #DEBUG this should be done by render  
    #vflo.move_center()

    #vflo.export_ngc(1, 0, .1, 2, '%s/%s.ngc'%(folder, outngc), do3d=False, do_retracts=False)
    vflo.export_geojson_lines(folder, outngc)
    #vflo.export_geojson_polygon(folder, outngc)

#vector_render_ngc()

##----------------------------------------------------
##----------------------------------------------------
def test_streamline(folder, bands, name, dopass='all', whichiter=2): 
    adbe = streamline()
    #adbe.sample_script('test.wrf')
    #adbe.load_json('%s/images/out/%s.json'%(folder, 'voronoi2'))

    #infile, outfolder, outname, numbands=4, dopass='all', invbmp=False, invpo=False
    adbe.oldrun(whichiter,
                '%s/images/in/%s'%(folder,name),
                '%s/images/out'%folder,
                'smithz',
                numbands=bands,
                dopass=dopass,
                invbmp = True, 
                invpo = False
                )

    def makegcode(num):
        vflo = vectorflow()
        for x in range(num-1):
            vflo.load_geojson('%s/images/out/%s.json'%(folder,x))
        vflo.export_ngc(1, 0, .1, 2, '%s/%s.ngc'%(folder, name), do3d=False, do_retracts=False)

    def makegcode_single(num, fid):
        vflo = vectorflow()
        vflo.load_geojson('%s/images/out/1.json'%folder,getfids=fid)
        vflo.export_ngc(1, 0, .1, 2, '%s/%s.ngc'%(folder, name), do3d=False, do_retracts=False)

    #makegcode(bands)




##----------------------------------------------------
def optimize_gcode(size, path, jsonfile, outname='processed'):
    """
        TODO 
            - make work with batch (path, name) in loop 
            - show how many polys reduced  
            - export json extents for each layer 
            - options for polygon, line, ngc 

        sort polygons by size from biggest to smallest 
        anything smaller than 1/4 of total gets sorted by quadtree 
        quadtree 
    """
    vflo = vectorflow()
    
    vflo.load_geojson(jsonfile) 

    filter_small_polys = True
    reorient           = False
    
    #vflo.load_geojson('/Users/klegg/serv/gnolmec/images/out/1.json') 
    #vflo.load_geojson('/Users/klegg/serv/gnolmec/images/out/2.json') 
           
    #vflo.polysize_info()
    
    if filter_small_polys:
        vflo.gr_sort = vflo.filter_by_bbox(size,size)
 
    #####################
    ## various common tranforms 
    
    #if model is negative corner 
    if reorient:
        vflo.gl_rotate(-90)

        #if model needs center shift 
        width,height = vflo.gl_width_height()
        vflo.gl_translate(width*.5,height/2,0)
        vflo.gl_scale([3,3,1])


    #if you need to center 
    #vflo.gl_move_center()


    ## NOTDONE ## #vflo.gl_move_extents_corner()

    #######################
    
    vflo.export_geojson_polygon(path, outname)
    vflo.export_geojson_lines(path, outname)

    #vflo.export_extents_ngc(path, outname)
    #vflo.export_ngc(1, 0, .1, 2, '%s/%s.ngc'%(path, outname), do_laser=True, do3d=False, do_retracts=False)




##----------------------------------------------------
def json_to_ngc(folder, name, scale, frame=1):
    vflo = vectorflow()
    vflo.load_geojson('%s/%s.json'%(folder,frame))
    
    vflo.gl_scale(scale)
    
    #my machine has long axis on (Y?)
    vflo.gl_rotate(-90)

    #vflo.filter_by_bbox
    
    vflo.export_extents_ngc()

    vflo.export_ngc(1, 0, .1, 2, '%s/%s.ngc'%(GLOBAL_PROJ, '%s%s'%(name,frame)), do_laser=True, do3d=False, do_retracts=False)



##----------------------------------------------------
##----------------------------------------------------
def json_to_laser(jsonfile, outpath, outname):
    vflo = vectorflow()
    vflo2 = vectorflow()
    
    vflo.load_geojson(jsonfile)
    vflo2.load_geojson(jsonfile)
    
    vflo2.rh = .01

    for ply in vflo.gr_sort:
        polypts = ply[4]
        pts = vflo.scanline_ngon(5, 20, polypts)
        vflo2.scanlines_to_segments(pts)

    vflo2.cvt_obj3d_grpoly()
    vflo2.export_ngc(1, 0, .1, 2, '%s/%s.ngc'%(outpath, outname) , do3d=False, do_retracts=False, do_laser=True)
    #vflo2.save('scanz.obj') 


#def batch_laser(folder, name):
#    for a in range(5):
#        json_to_laser('/Users/klegg/serv/gnolmec/images/out/%s.json'%a, folder, 'holez%s'%a)

##----------------------------------------------------
##----------------------------------------------------
def test_scanline(jsonfile, outpath, outname):
    vflo = vectorflow()
    vflo2 = vectorflow()
    
    vflo.load_geojson(jsonfile)
    vflo2.load_geojson(jsonfile)

    polypts = vflo.gr_sort[230][4]
    print(polypts)
    
    pts = vflo.scanline_ngon(5, 20, polypts)
    vflo2.scanlines_to_segments(pts)

    vflo2.cvt_obj3d_grpoly()
    vflo2.export_ngc(1, 0, .1, 2, '%s/%s.ngc'%(outpath, outname) , do3d=True, do_retracts=True)





##----------------------------------------------------
##----------------------------------------------------


def break_into_many_ngc(folder, infile, outfile):
    vflo = vectorflow()
    vflo.load_geojson(infile )
    numpolys =(len(vflo.gr_sort))
    vflo.flush()
    print("FOUND %s POLYGONS "%numpolys)

    for x in range(numpolys):
        vflo.flush()
        vflo.load_geojson(infile,getfids=x)
        vflo.scale_gcode([.5,.5,.5])
        vflo.export_ngc(1, 0, .1, 2, '%s/%s.ngc'%(folder, '%s_%s'%(outfile,x)), do3d=False, do_retracts=True, do_laser=True)


##----------------------------------------------------
##----------------------------------------------------
def test_3dprinting(folder, infile, outfile):
    """
       #DEBUG - very experimental 
    
    """

    gop = cnc_op()

    #gn_dir(gop)
    #gop.prim_sphere()
    #gop.prim_cube()
    
    #gop.load_geojson(infile, 0) 
    
    gop.rh = 10 

    gop.show_setup()

    origin = [100,100,0]
    gop.prim_quad( axis='z', pos=origin, rot=(0,0,0), size=10, periodic=True)
    
    #origin = [140,140,0]
    #gop.prim_quad( axis='z', pos=origin, rot=(0,0,0), size=10, periodic=True)

    gop.cvt_obj3d_grpoly()

    print(gop.gr_polys)

    #             (rh, ch, cdpi, cmax, filename                           , do3d=False)
    gop.export_ngc(1,   0,   .1,    2, '%s/%s.ngc'%(folder, outfile), do3d=True )

    #gop.export_3dprint( '%s/%s.gcode'%(folder, 'isofoo') )


##----------------------------------------------------
def test_milling(folder, infile, outfile):
    """
       #DEBUG - very experimental 
    
    """
        
    gop = cnc_op()

    #gn_dir(gop)
    #gop.prim_sphere()
    #gop.prim_cube()
    origin = [100,100,0]

    #gop.load_geojson(infile, 0)  

    for x in range(20):
        gop.prim_quad( axis='z', pos=origin, rot=(0,0,0), size=10)

    gop.cvt_obj3d_grpoly()

    #             (rh, ch, cdpi, cmax, filename                           , do3d=False)
    #gop.export_ngc(1,   0,   .1,    2, '%s/%s.ngc'%(GLOBAL_PROJ, 'isofoo'), do3d=True )

    gop.export_3dprint( '%s/%s.gcode'%(folder, outfile) )

    #gop.show_setup()



##----------------------------------------------------

def load_gcode(basepath, path, name):
    vflo = vectorflow()

    gop = cnc_op()
    gop.load_gcode('%s/%s/%s.ngc'%(basepath, path, name) )
    for ply in  gop.segments:
        #print(len(ply))
        vflo.gr_polys.append(ply) 
    vflo._sort()

    #if you want to prove the data is in true 3D, do this. 
    #vflo.gl_rotate([45,45,45]) 

    #vflo.export_geojson_polygon(basepath, 'textgrave')
    vflo.export_geojson_lines(basepath, '%s.json'%name)

    # OBJ not working very well 
    #vflo.cvt_grpoly_obj3d() 
    #vflo.save('%s/%s.obj'%(basepath, name) )



##----------------------------------------------------



##----------------------------------------------------
def kicad_to_json(path, infile, outfile):
    vflo = vectorflow()
    pcb = pcbfile()
    pcb.load_kicadpcb(0, path,infile )
    
    vflo.gr_polys = pcb.gr_polys
    vflo._sort()
    
    #move_center is for 3d obj buffer
    #vflo.move_center()
    
    #this is for json/ngc buffer
    vflo.gl_move_center() 

    
    #kicad uses top left origin - we need to flip it 
    vflo.gl_scale([1,-1,1])

    vflo.export_geojson_polygon(path, outfile)



##----------------------------------------------------

def test_plotter(path, infile, outfile):
    vec = vectorflow()
    vec.load_geojson(path+'/'+infile)
    
    vec.gl_scale(4)

    #name = infile.split('.')

    #DEBUG extents is not scaling!! grr 
    #print(vec.gl_extents())
    vec.export_extents_ngc(path , outfile, type='plot_ngc')

    #vec.export_geojson_polygon()

    #self,       rh,ch,cdpi,cmax, filename, do3d=False, do_retracts=True, do_laser=False, laserpwm=400):
    vec.export_ngc(0,0,0, 0, filename='%s/%s.ngc'%(path,outfile), do3d=False, do_retracts=True, do_laser=False, laserpwm=400 ,  do_gpio=0)




##----------------------------------------------------
##----------------------------------------------------


def recursive_spiral(basepath):
    """DEBUG not done - converting from external code  
    """ 
    vflo = vectorflow()
    def draw_spiral(x, y, length, direction):
        L = length
        c = 0
        while length>1 or c<3:
            if length>2:
                draw_spiral(x,y,length*0.255,160+direction)
            
            #turtle.up()
            #turtle.seth(direction)
            #turtle.goto(x,y)
            v = vec3(x,y,0)
            v = v* vec3(direction,0,0)

            vflo.prim_circle(axis='y', pos=v, spokes=4)
            #def prim_circle(self, axis, pos=(0,0,0), rot=(0,0,0), dia=1, spokes=9):

            #if length <= 2: 
            #    turtle.down()
            #turtle.fd(length)
            #x,y = turtle.xcor(), turtle.ycor()

            #print(length, direction,x,y)

            length *= 0.93
            direction += 20
            c += 1
    draw_spiral(1,-1,30,9)
    vflo.cvt_obj3d_grpoly()
    #vflo.export_geojson_polygon(basepath, 'spiral.json')
    vflo.export_geojson_lines(basepath, 'spiral.json') 






##----------------------------------------------------
##----------------------------------------------------
##----------------------------------------------------
"""
def vectorizer(inputfile, outputfile=None):
    # example of PIL edge detect 
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
    # based on kicad img convert 
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
"""




    


