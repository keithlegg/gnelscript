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
    #bbox = kiparser.calc_2d_bbox_pt(2, (5,5))
    #pts = kiparser.cvt_2d_to_3d(kiparser.extents_fr_bbox(bbox, periodic=True))
    #debug - need to solve the clean_pts_str debacle?
    #kiparser.gr_polys.append(pts)

    ##--
    #bbox = kiparser.calc_2d_bbox_pt(1.75, (-3,3))
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

    def oldrun (self, whichiter, infile, outfolder, outname, numbands=4, invbmp=False, invpo=False, dopass='all'):
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






def test_scanline2(jsonfile, outpath, outname):
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




    


