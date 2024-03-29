
import os, sys, math

from math import sqrt

from PIL import Image, ImageOps

from gnelscript.pygfx.math_ops import  NUMPY_IS_LOADED
from gnelscript.pygfx.point_ops_2d import pop2d


if NUMPY_IS_LOADED:
    import numpy as np  
else:
    print(' ## debug - numpy module disabled in point ops. ')



class raster_op(object):
    
    def __init__(self):
        self.ptgen = pop2d()

        self.debug_msg = False
        self.res_x = None
        self.res_y = None
        self.bitmode = 'RGBA' #PIL mode
        self.fb = None #main framebuffer  
 
    @property  
    def extents(self):
        return self.ptgen.calc_square_diag( (0,0), (self.size[0],self.size[1]) )  

    @property  
    def center(self):
        return ( int(self.res_x/2), int(self.res_y/2) ) 

    @property  
    def size(self):
        return ( self.fb.size )   

    def log_msg(self, *args ):
        if self.debug_msg:
            msg=''
            for value in args :
                msg+=( "%s " % str(value) )
            print(msg)

    def image_info(self):
        stats = {}
        stats['resolution']='W:'+str(self.res_x)+' H:'+str(self.res_y)
        print( stats )

    def save(self, name, noalpha=False):
        print("Saving file to: " + name)
        if noalpha:
            self.fb.convert('RGB').save(name)
        else:    
            self.fb.save(name)

    #def image_from_bytes(self, bytes_data):
    #    self.fb = Image.open( bytes_data )

    def load(self, name):
        self.fb = Image.open(name)
        self.res_x = self.fb.size[0]
        self.res_y = self.fb.size[1]
   
    def set_res(self, rx, ry):
        self.res_x = rx
        self.res_y = ry 

    def empty_buffer(self, color):
        """ make a new PIL Image that has the same properties of this one 
            fill with a color 
            
            usage:
                ro = raster_op()
                ro.create_buffer(800,600)
                ro.save("images/out/foobar.jpg", noalpha=True)
                po = ro.empty_buffer( (255,255,0) )
                po.convert('RGB').save("images/out/foobar2.jpg", noalpha=True)

        """
        
        image =  Image.new(self.bitmode, (self.res_x, self.res_y) )
        dpix = image.load() 
        for x in range(self.res_x):
            for y in range(self.res_y):
                dpix[ x, y ] = color
        return image  

    def create_buffer(self, rx, ry):
        """ create a brand new image in memory  

            usage:

                ro = raster_op()
                ro.create_buffer(800,600)
                ro.save("images/out/foobar.jpg", noalpha=True)

        """        
        self.res_x = rx
        self.res_y = ry        
        self.fb = Image.new(self.bitmode, (self.res_x, self.res_y) )
   
    def read_buffer(self, pilBuffer):
        #make sure you pass a PIL Image object 
        if type(pilBuffer) == Image.Image:
            self.fb = pilBuffer
            self.res_x = pilBuffer.size[0]
            self.res_y = pilBuffer.size[1] 
            #print("debug raster op buffer read ", self.fb.show() )
        else:
            print(type(pilBuffer))
            raise ValueError('object passed to read_buffer that was not a PIL Image') 

    def invert_bw(self):
        self.invert()
        if self.fb.mode != 'L':
            self.fb = self.fb.convert('L')

    def invert(self):
        if self.fb.mode == 'RGB':
            self.fb= ImageOps.invert(self.fb)
        if self.fb.mode == 'RGBA':
            self.fb= ImageOps.invert(self.fb)            
            #self.fb = self.fb.convert('RGBA')

    def cvt_1bit(self):    
        #img.point(lambda x: bool(x))
        self.fb = self.fb.convert('L') # convert 8 bit
        self.fb = self.fb.convert('1') # convert 1 bit

    def cvt_24bit_alpha(self):    
        self.fb = self.fb.convert("RGBA")
    
    ##---------------------------------------------##     
    def cvt_24bit(self):    
        self.fb = self.fb.convert("RGB")

    ##---------------------------------------------## 
    def get_pix(self, pt ):
        """ DEBUG - add a "wrap around feature  
            if a pixel is requested out of bounds - loop around to start 
        """
        get_pt = list(pt) 
        # print(get_pt)

        if pt[0] >= self.size[0]:
            new_x = pt[0] - self.size[0]
            get_pt = [new_x, pt[1]]
            #print( '## error U get pix ', pt, get_pt)
        if pt[1] >= self.size[1]:
            new_y = pt[1] - self.size[1]
            get_pt = [pt[0], new_y]
            #print( '## error V get pix ', pt, get_pt)

        return self.fb.getpixel(tuple(get_pt)) 
    
    ##---------------------------------------------## 
    def set_pix(self, pt , color ):
        dpix = self.fb.load()
        dpix[pt[0], pt[1]] = color
    
    ##---------------------------------------------##  
    def rotate_pil_raw(self, rotation):
        
        #rotate and expand - nothing else
        self.fb = self.fb.rotate(rotation, expand=1)        
        self.res_x = self.fb.size[0]
        self.res_y = self.fb.size[1] 
    
    ##---------------------------------------------## 
    def rotate_pil(self, rotation):
        
        #rotate, expand and composite white in the empty areas
        if self.fb.mode != 'RGBA':
            self.fb = self.fb.convert('RGBA')

        rot = self.fb.rotate(rotation, expand=1)        
        self.res_x = rot.size[0]
        self.res_y = rot.size[1] 

        fff = Image.new('RGBA', rot.size, (255,)*4) #white mask to composite
        self.fb = Image.composite(rot, fff, rot)

    ##---------------------------------------------## 
    def add_margins(self, size):

        old_size = self.fb.size
        new_size = (old_size[0]+size, old_size[1]+size)

        #new_im = Image.new(self.fb.mode, new_size)    #for black
        new_im = Image.new(self.fb.mode, new_size, (255,)*4) #for white 

        new_im.paste(self.fb, (new_size[0]-old_size[0])/2, (new_size[1]-old_size[1])/2  )
        new_im.show()

    ##---------------------------------------------## 
    def get_island(self, offset=None):
        """       
            this uses PIL.getbbox to exclude "empty" data from the edges of an image  
        """

        self.cvt_24bit()
        tmp_fb = ImageOps.invert(self.fb)
        inside_data = tmp_fb.getbbox()
        self.cvt_24bit_alpha()
        if not offset:
            return self.ptgen.extents_fr_bbox(inside_data)
        if offset:
            return self.ptgen.extents_fr_bbox(inside_data, offset)

    ##---------------------------------------------## 
    def crop_island(self, margin=None):
        """ crop out the image borders with no data in them 
            optional margin will buffer the image borders with white 
            nagative margins will trim the image edges
        """
    
        #you cant invert an image with alpha in PIL  
        self.cvt_24bit() #first we convert to RGB 
       
        tmp_fb = ImageOps.invert(self.fb)
        inside_data = tmp_fb.getbbox() #crops out black pixels around edges

        if not margin:
            self.fb = self.fb.crop( inside_data )
            self.cvt_24bit_alpha()#convert back to RGBA   
        
        if margin:
            if margin <0:
                inside_data=self.ptgen.bbox_buffer_2d(inside_data, margin)
                self.fb = self.fb.crop( inside_data )
                self.cvt_24bit_alpha()#convert back to RGBA   

            if margin >0:
                self.fb = self.fb.crop( inside_data )
                self.cvt_24bit_alpha()#convert back to RGBA 

                double = int(margin*2)
                bgimg = Image.new('RGBA', (self.fb.size[0]+double, self.fb.size[1]+double), (255,)*4) #white mask to composite
                img_w, img_h = self.fb.size
                bg_w, bg_h = bgimg.size
                bgimg.paste(self.fb, (margin, margin ) )
                ##
                self.fb = bgimg
                self.cvt_24bit_alpha()#convert back to RGBA  
                self.res_x = bg_w
                self.res_y = bg_h 

    ##---------------------------------------------## 
    def crop_pt(self, size, pt_coord):
        """
        # crop area from point

        usage: 
            ro = raster_op()
            ro.load("images/in/refer.jpg")
            out = ro.crop_pt( 400, ro.center) 
            out.save("images/out/foobar.png")       
        """
        xtntx = tuple(self.ptgen.calc_bbox( size, pt_coord) )
        if xtntx[0]<0 or xtntx[1]<0  or xtntx[2]>self.res_x  or xtntx[3]>self.res_y:
            print('# ERROR raster_ops.crop_pt - out of image bounds') 
 
        #self.fb = self.fb.crop( xtntx )
        return self.fb.crop( xtntx )

    ##---------------------------------------------##       
    def crop_corner(self, size, mode):
        """
        #crop the corners in a square 

        usage:
            ro = raster_op()
            ro.load("images/in/refer.jpg")
            out = ro.crop_corner( 100, 'bl') 
            out.save("images/out/foobar.png")        
        """
        if mode == 'bl':
            xtntx = (0, self.res_y - size, size, self.res_y)
        if mode == 'tl':
            xtntx = (0, 0, size, size)            
        if mode == 'tr':
            xtntx = (self.res_x - size, 0, self.res_x, size)            
        if mode == 'br':
            xtntx = (self.res_x - size, self.res_y - size, self.res_x, self.res_y)                  

        return self.fb.crop( xtntx ) 



        
class pixel_op (raster_op):
    """ 
        Pixel operator with raster goodies for drawing and sampling pixels 
        TODO:
           deal with cases where the sampling runs off the page 
    """

    ##---------------------------------------------##     
    def __init__(self):
        super().__init__()  
        self.filter = pixelFilter() 
        #super(pixel_op , self).__init__() 

    ##---------------------------------------------## 
    def graticule(self, divs=10, spacing=None, scale=1):
        """  make a graticule grid  
            start at center and go out from there based on spacing value
            spacing is in pixels 
  
            usage:
                ro = pixel_op()
                ro.create_buffer(800,800)
                ro.graticule(8,1)
                ro.save("images/out/foobar.png")

        """
        
        clr_backg = (0,50,90)
        clr_lines = (0,150,190)
        clr_dots  = (0,255,0)
        gridcolor = (145,120,180)

        #draw a dot in the center
        cen_x = self.center[0]  
        cen_y = self.center[1] 

        #flood fill back ground
        self.fill_color( clr_backg )

        #optional zoom 
        if spacing:
            spacing = spacing*scale
        res_x = self.res_x*scale
        res_y = self.res_y*scale

        #spacing puts lines at fixed spaces 
        if spacing:
            x = cen_x
            while(x<self.res_x):
                self.connect_the_dots( [(x, 0), (x, res_y)],
                                        clr_lines, 1 )
                x+=spacing
           
            x = cen_x
            while(x>0):
                self.connect_the_dots( [(x, 0), (x, res_y)],
                                        clr_lines, 1 )
                x-=spacing

            y = cen_y
            while(y>0):
                self.connect_the_dots( [(0, y), (res_x, y)],
                                        clr_lines, 1 )
                y-=spacing

            y = cen_y
            while(y<self.res_y):
                self.connect_the_dots( [(0, y), (res_x, y)],
                                        clr_lines, 1 )
                y+=spacing

        #divs divides up the image based on num breaks
        if divs:
            pass


        #draw lines from center across image 
        self.vert_line(self.center[0], gridcolor)
        self.horiz_line(self.center[1], gridcolor)

        #put a dot at the center
        self.draw_fill_circle(self.center[0],self.center[0], 2, (200,255,0) ) 

    ##---------------------------------------------## 
    def draw_cntr_line(self, points, color=(0,255,200), size=1, mag=1, framebuffer=None):
        """ DEBUG use offset feature of connect_the_dots  """
        
        if mag >1:
            tmp = []
            for pt in points:
                tmp.append( (pt[0]*mag, pt[1]*mag, pt[2]*mag ) )
            points = tmp 
        else:
            tmp = points

        if framebuffer==None:
            framebuffer = self.fb
        self.connect_the_dots( tmp, color, size, origin=(self.center[0] ,self.center[1]), framebuffer=framebuffer)

    ##---------------------------------------------## 
    def draw_cntr_pt(self, dot, size=1, origin=(0,0), color=(255,0,0), framebuffer=None):
        """ draw a point relative to center of image """

        sp = (self.center[0]+origin[0]) +  dot[0]
        ep = (self.center[1]+origin[1]) +  dot[1] #make y negative to flip "UP" -PIL uses top left origin

        #put a dot at the center
        self.draw_fill_circle(sp, ep, size, color ) 

    ##---------------------------------------------## 
    def render_line_2d(self, pts, invert_y=False, origin=(0,0), scale=10):
        """ this will draw a point in 2d from the center of the image 
            the only slightly special thing going on here is the scale factor 
            this allows for an image with variable grid size 
        """
        
        #make y negative to flip "UP" -PIL uses top left origin
        #-1 will flip , 1 will NOT flip 
        if invert_y:
            invert = -1
        else:
            invert = 1    
    
        transformed_pts = []

        for pt in pts:
            x = (self.center[0]+origin[0]) + (pt[0]*scale)
            y = (self.center[1]+origin[1]) + (pt[1]*scale) * invert 
            
            transformed_pts.append( (x,y) ) 
        #print('### ', pts )
        self.connect_the_dots( transformed_pts, (230, 200, 0 ), 1 )

    ##---------------------------------------------## 
    def render_point_2d(self, vec, invert_y=True, origin=(0,0), scale=10, dotsize=3):
        """ this will draw a point in 2d from the center of the image 
            the only slightly special thing going on here is the scale factor 
            this allows for an image with variable grid size 
        """
        
        #make y negative to flip "UP" -PIL uses top left origin
        #-1 will flip , 1 will NOT flip 
        if invert_y:
            invert = -1
        else:
            invert = 1    

        x = (self.center[0]+origin[0]) + (vec[0]*scale)
        y = (self.center[1]+origin[1]) + (vec[1]*scale) * invert 

        self.draw_fill_circle( x, y, dotsize, (230,0,0) ) 

    ##---------------------------------------------## 
    def render_vector_2d(self, vec, invert_y=True, origin=(0,0), scale=10):
        """ this will draw a vector in 2d from the center of the image 
            the only slightly special thing going on here is the scale factor 
            this allows for an image with variable grid size 
        """

        #make y negative to flip "UP" -PIL uses top left origin
        #-1 will flip , 1 will NOT flip 
        if invert_y:
            invert = -1
        else:
            invert = 1    

        ex = (self.center[0]+origin[0]) + (vec[0]*scale)
        ey = (self.center[1]+origin[1]) + (vec[1]*scale) * invert 

        #self.graticule(50) # - better to do this outside of this fucntion for more control over layering 

        self.connect_the_dots([ ((self.center[0]+origin[0]),(self.center[1]+origin[1])),
                                (ex,ey)], (0,230,0), 1 )

        #calc_angle = self.ptgen.old_calc_theta_vert( ( (self.center[0]+origin[0]), (self.center[1]+origin[1]) ),  (ex,ey) )
        # print("ANGLE OF VECTOR FROM VERTICAL (UP) %s" % calc_angle)

    ##---------------------------------------------## 
    def render_obj_2d(self, obj, invert_y=True, origin=(0,0), scale=10):
        """ this will draw an obj file in orthogrphic 2d from the center of the image 

            the only slightly special thing going on here is the scale factor 
            this allows for an image with variable grid size 
        """
        
        #make y negative to flip "UP" -PIL uses top left origin
        #-1 will flip , 1 will NOT flip 
        if invert_y:
            invert = -1
        else:
            invert = 1    
    
        transformed_pts = []

        #for pt in pts:
        #    x = (self.center[0]+origin[0]) + (pt[0]*scale)
        #    y = (self.center[1]+origin[1]) + (pt[1]*scale) * invert 
        #    
        #    transformed_pts.append( (x,y) ) 
        self.connect_the_dots( transformed_pts, (230, 200, 0 ), 1 )

    ##---------------------------------------------## 
    def normal_to_color(self, norml):
        out = [0,0,0]
        
        out[0]=int(norml[0]*255)
        out[1]=int(norml[1]*255)
        out[2]=int(norml[2]*255)

        if out[0]>255:
             out[0]=255
        if out[1]>255:
             out[1]=255
        if out[2]>255:
             out[2]=255             
        return tuple(out)

    ##---------------------------------------------## 
    def tint(self, color, com):
        """ i needed a way to procedurally tweak color 
          used for the glowing neon line effect to darken linear borders
        """
        amt = 120
        clamp_low = amt  #poor - make this better!
        clamp_high = amt #poor - make this better!

        tmp = 0
        tl = [0,0,0]

        #minus_red
        if com == 'drkr':
            t = color    
            if t[0]>clamp_low:
                tl[0]=t[0]-amt
            if t[1]>clamp_low:
                tl[1]=t[1]-amt
            if t[2]>clamp_low:
                tl[2]=t[2]-amt
                                                
            return ( tl[0], tl[1], tl[2] )

        #minus_red
        if com == 'mr':
            t = color    
            if t[0]>clamp_low:
                tmp=t[0]-amt
            return ( tmp, t[1], t[2] )

        #minus_green
        if com == 'mg':
            t = color    
            if t[1]>clamp_low:
                tmp =t[1]-amt
            return ( t[0], tmp, t[2] )

        #minus_blue
        if com == 'mb':
            t = color    
            if t[2]>clamp_low:
                tmp =t[2]-amt
            return ( t[0], t[1], tmp )
  
    ##---------------------------------------------## 
    def insert_numpy(self, data, brightness=255 ):
        for y,row in enumerate(data):
            for x,col in enumerate(row):
                color = ( int(col[0]*brightness), int(col[1]*brightness), int(col[2]*brightness))
                self.set_pix( (x, y) , color)

    ##---------------------------------------------## 
    def insert_image(self, px, py, foregroundfile, backgroundfile, outfile):
        """ Filename1 and 2 are input files; outfile is a path where results are saved (with extension).

        """

        img = Image.open(foregroundfile ,'r')
        img_w, img_h = img.size
        
        #background = Image.new('RGBA', (1024,512), (255, 255, 255, 255))
        bgimg = Image.open(backgroundfile ,'r')        
        bg_w, bg_h = bgimg.size
        bgimg.paste(img, (px, py ) )
        bgimg.save(outfile)

    ##---------------------------------------------## 
    def fill_color(self, color, framebuffer=None):    
        """ fills image with solid color 

        """
        
        if framebuffer:
            self.read_buffer(framebuffer)
        else:
            framebuffer= self.fb
        


        dpix = framebuffer.load() 
        for x in range(self.res_x):
            for y in range(self.res_y):
                dpix[ x, y ] = color

    ##---------------------------------------------## 
    def center_square(self, tl, br, color, framebuffer=None):  
        """ fills a centered square from the top left to bottom right corner """  
        
        if framebuffer==None:
            framebuffer = self.fb

        dpix = framebuffer.load() 
        for x in range(tl, self.res_x):
            for y in range(br, self.res_y):
                if x <self.res_x-tl and y <self.res_y-br:
                    dpix[ x, y ] = color
                #if y <self.res_y-br:
                #    dpix[ x, y ] = color

    ##---------------------------------------------##     
    def vert_line(self, xloc, color, framebuffer=None):    
        """ draw vertical line across entire image """
        
        if framebuffer:
            self.read_buffer(framebuffer)
        else:
            framebuffer= self.fb

        dpix = framebuffer.load() 
        for x in range(self.res_x):
            if x == xloc:
                for y in range(self.res_y):
                    dpix[ x, y ] = color

    ##---------------------------------------------##     
    def horiz_line(self, yloc, color, framebuffer=None):    
        """ draw horizontal line across entire image """
        
        if framebuffer:
            self.read_buffer(framebuffer)
        else:
            framebuffer= self.fb

        dpix = framebuffer.load() 
        for y in range(self.res_y):
            if y == yloc:
                for x in range(self.res_x):
                    dpix[ x, y ] = color

    ##---------------------------------------------## 
    def vert_line_thick(self, xloc, width, color, framebuffer=None):    
        """ draw horizontal line with thickness """
        
        if framebuffer:
            self.read_buffer(framebuffer)
        else:
            framebuffer= self.fb

        dpix = framebuffer.load() 
        for x in range(self.res_x):
            if x == xloc:
                for w in range(x, x+width):
                    for y in range(self.res_y):
                        dpix[ w, y ] = color

    ##---------------------------------------------## 
    def batch_draw_pixels(self, data, framebuffer=None):
        """ draw scanned data back into an image [ (value, (x,y)) .. ] """

        if framebuffer:
            self.read_buffer(framebuffer)
        else:
            framebuffer= self.fb
        
        dpix = framebuffer.load() 
        for px in data:
            #dpix[px[1][0], px[1][1]] = px[0]
            if px[0] ==1:
                dpix[px[1][0], px[1][1]] = red
            if px[0] ==0:                
                dpix[px[1][0], px[1][1]] = green
 
    ##---------------------------------------------## 
    def draw_fill_circle(self, x_orig, y_orig, dia, color, framebuffer=None):
        if framebuffer:
            self.read_buffer(framebuffer)
        else:
            framebuffer= self.fb

        for x in range(dia):
            self.draw_circle( x_orig, y_orig, x, color, framebuffer)
    
    ##---------------------------------------------##         
    def draw_circle(self, x_orig, y_orig, dia, color, framebuffer=None):
        plot_x = 0;plot_y = 0;
        if framebuffer:
            self.read_buffer(framebuffer)
        else:
            framebuffer= self.fb

        dpix = framebuffer.load()  
        
        if framebuffer.mode =='P':
            if color[0] or color[1] or color[2]:
                color = 128
            else:
                color = 0

        for i in self.ptgen.calc_circle_2d(x_orig, y_orig, dia):
            try:
                dpix[ i[0], i[1] ]= color
            except IndexError:
                pass

    ##---------------------------------------------##         
    def draw_points_batch(self, points, color, dia, origin=(0,0),framebuffer=None):
        """ debug - add check to make sure it doesnt go off edge of page 

            usage:
                ro = pixel_op()
                ro.create_buffer(800,800)
                ro.graticule(spacing=100, scale=1)
                pts = [(100,100),(-200,100), (-50,-100) ]
                ro.draw_points_batch(pts, (255,0,0), 5, origin=ro.center )
                ro.save("images/out/foobar.png")

        """

        if framebuffer:
            self.read_buffer(framebuffer)
        else:
            framebuffer= self.fb

        for pt in range(len(points)):
            p1 = list(points[pt])    

            #shift to another place before drawing
            #if origin[0]!=0 or origin[1]!=0:
            p1[0] = p1[0]+origin[0]
            p1[1] = p1[1]+origin[1]

            self.draw_fill_circle(p1[0], p1[1], dia, color, framebuffer)

    ##---------------------------------------------##     
    def connect_the_dots(self, points, color, thickness, origin=(0,0), framebuffer=None):
        """ debug - add check to make sure it doesnt go off edge of page 
                    clipping is a whole thing unto itself and maybe better put elsewhere? 

            usage:
                ro = pixel_op()
                ro.create_buffer(800,800)
                ro.graticule(spacing=100, scale=1)
                pts = [(1,1),(-20,3), (-20,-10) ]
                ro.connect_the_dots(pts, (255,0,0), 2, origin=ro.center )
                ro.save("images/out/foobar.png")

        """

        if framebuffer:
            self.read_buffer(framebuffer)
        else:
            framebuffer= self.fb
        ##
        #count = 0
        for pt in range(len(points)-1):
            p1 = list(points[pt])    
            p2 = list(points[pt+1])  

            #shift to another place before drawing
            if origin[0]!=0 or origin[1]!=0:
                p1[0] = p1[0]+origin[0]
                p1[1] = p1[1]+origin[1]
                p2[0] = p2[0]+origin[0]
                p2[1] = p2[1]+origin[1]

            #if count>0:
            color=color
            self.draw_line(tuple(p1), tuple(p2), color, thickness, framebuffer)
            #count += 1

    ##---------------------------------------------##     
    def draw_vector(self, vec, color, thickness=0, framebuffer=None):
        self.draw_line(vec[0] , vec[1] , color, thickness, framebuffer  )

    ##---------------------------------------------## 
    def draw_line(self, pt1, pt2, color, thickness=0, framebuffer=None):
        if framebuffer:
            self.read_buffer(framebuffer)
        else:
            framebuffer= self.fb

        pts = self.ptgen.calc_line(pt1[0], pt1[1], pt2[0], pt2[1])
        
        #attempt to make it work at different bit depths 
        if framebuffer.mode =='P':
            if color[0] or color[1] or color[2]:
                color = 0 #black
            else: 
                color = 128 #white 

        dpix = framebuffer.load()          
        for pt in pts:
            if not thickness:
                try:
                    dpix[ pt[0], pt[1] ] = color 
                except IndexError:
                    #print('clipping error, check yo coordinates dawg. ')
                    pass

            #really crappy way to add line thickness - makes a point a "plus sign"      
            if thickness:
                for pthick in range(0, thickness):
                    try:
                        dpix[ pt[0],        pt[1]        ] = color  
                        dpix[ pt[0],        pt[1]+pthick ] = color  
                        dpix[ pt[0],        pt[1]-pthick ] = color                                                         
                        dpix[ pt[0]+pthick, pt[1]        ] = color  
                        dpix[ pt[0]-pthick, pt[1]        ] = color
                    except IndexError:
                        #print('clipping error, check yo coordinates dawg. ')
                        pass

    ##---------------------------------------------## 
    def draw_pt_along_vector(self, pt1, pt2, num, color, dia=1, framebuffer=None):
        """ draw any number of points along a vector """

        if framebuffer:
            self.read_buffer(framebuffer)
        else:
            framebuffer= self.fb

        pts = self.ptgen.locate_pt_along( pt1[0], pt1[1], pt2[0], pt2[1], num )
        dpix = framebuffer.load()  
            
        for pt in range(len(pts)):
            self.draw_fill_circle( pts[pt][0], pts[pt][1], 5, color, framebuffer)
        

    ##---------------------------------------------## 
    def color_distance(self, rgb1, rgb2):
        """ number between 0 and 440 ish """

        if len(rgb1)==3:
            r,  g,   b = rgb1
        if len(rgb1)==4:
            r,  g,   b, a  = rgb1

        if len(rgb2)==3:
            cr, cg, cb = rgb2
        if len(rgb2)==4:
            cr, cg, cb, ca = rgb2

        return sqrt((r - cr)**2 + (g - cg)**2 + (b - cb)**2)
    
    ##---------------------------------------------## 
    def closest_color(self, colors, rgb):
        """
        https://stackoverflow.com/questions/54242194/python-find-the-closest-color-to-a-color-from-giving-list-of-colors

        colors = (
            (181, 230, 99),
            (23, 186, 241),
            (99, 23, 153),
            (231, 99, 29),
        )
        """

        r, g, b = rgb
        color_diffs = []
        for color in colors:
            cr, cg, cb = color
            color_diff = sqrt((r - cr)**2 + (g - cg)**2 + (b - cb)**2)
            color_diffs.append((color_diff, color))
        return min(color_diffs)[1]

    ##---------------------------------------------## 
    def extract_by_color(self, path, name, color, slowmode=False, exactcolor=False, invert=False, framebuffer=None):
        """ INVERT IS BROKEN - FIX IT DEBUG !
            self.empty_buffer() MAY ALSO BE BROKEN - THINK MORE ABOUT IT 

            return a new image that matches a color ignoring all others 

        """

        if framebuffer:
            self.read_buffer(framebuffer)
        else:
            framebuffer= self.fb

    
        new_fb = self.empty_buffer(  (0,0,0,255) ) 

        siz = framebuffer.size
        
        src_pix = framebuffer.load() 
        dpix    = new_fb.load() 

        for x in range(siz[0]):
            for y in range(siz[1]):

                ## reeeeaaallly slow because it iterates ALL pixels with math on each one 
                ## should poduce a better image because it will average colors into one 
                if slowmode: 
                    if self.color_distance(src_pix[ x, y ], color)<100: 
                         if exactcolor:
                             dpix[ x, y ]=color
                         else:  
                             dpix[ x, y ]=(255,255,255)
                    #else:
                    #    if invert: 
                    #         dpix[ x, y ]=(255,255,255)

                #fast but only catches colors that match exactly 
                else: 
                    if src_pix[ x, y ] == color: 
                         if exactcolor:
                             dpix[ x, y ]=color
                         else:  
                             dpix[ x, y ]=(255,255,255)
                    #else:
                    #    if invert: 
                    #         dpix[ x, y ]=(255,0,255)

        if invert:
            new_fb = new_fb.convert('RGB')
            print("INVERTING ", new_fb.mode)
            im_invert = ImageOps.invert(new_fb)
            im_invert.save("%s/%s.bmp"%(path,name) )

        else:              
            new_fb.save("%s/%s.bmp"%(path,name) ) 




class pycv(pixel_op):
    """ sketch of a computer vision system to locate fiducuals 
        moved into its own class from pixel_op
    """

    def __init__(self):
        super().__init__()  

    def pretty_much_yellow(self, pixel):
        if pixel[0]>250 and  pixel[1]>250 and  pixel[2]<15:
            return True
        return False

    def line_scan(self, pt1, pt2 , filterNoise=False, framebuffer=None):
        """
            filternoise is a tuple/ (#places to look forward/back , replace value) 

            scan a row of pixels along a line and return array of 1's and 0's 
            this is useful in two ways:
             - it checks averages pixels into black or white and stores them serialized
             - it also stores location of each pixel in XY catesian space
             
             [(PIXEL, ( XCOORD, YCOORD) ) .. ]
        """

        if framebuffer:
            self.read_buffer(framebuffer)
        else:
            framebuffer= self.fb


        x_1 = pt1[0]
        y_1 = pt1[1]
        
        x_2 = pt2[0]
        y_2 = pt2[1]
        
        pts = self.ptgen.calc_line( x_1, y_1, x_2, y_2 )
        output = []
        for pt in pts:
            pixel_bit = 0
            if self.scanner_darkly( framebuffer.getpixel(pt)):
                pixel_bit = 1
            output.append( (pixel_bit,pt) ) #( (color,coordinate),.. )

        if filterNoise:
            output = self.filter.filter_noise(output, filterNoise[0], filterNoise[1], True, False)
        return output

    def get_luminance(self, point, framebuffer=None):
        """ Returns the perceived luminance, from 0 - 1, of given point.
            
            Works with 'RGBA', 'RGB', 'L' and '1' image modes, and if an unexpected mode is 
            encountered raise an assertion error. For RGB uses ITU-R 601-2 luma transform.
            Alpha channels are ignored.
        """
        
        if framebuffer:
            self.read_buffer(framebuffer)
        else:
            framebuffer= self.fb

        mode = framebuffer.mode
        assert mode in ('RGBA', 'RGB', 'L', '1')
        
        # Running off the edge of the page shall be black
        try:
            color = framebuffer.getpixel(point)
        except IndexError:
            return 0
        
        if mode == 'RGBA':
            brightness = (0.299 * color[0] + 0.587 * color[1] + 0.114 * color[2]) # * (color[3] / 255) # ignore alpha
            brightness = brightness / 255
        elif mode == 'RGB':
            brightness = 0.299 * color[0] + 0.587 * color[1] + 0.114 * color[2]
            brightness = brightness / 255
        elif mode == 'L':
            brightness = color / 255
        elif mode == '1':
            brightness = color
        return brightness

    def scanner_darkly(self, pixel):
        """ quantize a 1/8/24 bit color pixel into a 1 bit boolean value """

        #1 bit  
        if isinstance( pixel, int ):
            if pixel==0:
                return True 
            else:
                return False   

        #24 or 32 bit                        
        else:
            #avg = sum([val for val in pixel]) / len(pixel)
            if pixel[0]<15 and pixel[1]<15 and pixel[2]<15:
                return True
        return False
            
    def line_scan_simple(self, pt1, pt2, calc_brightness=True, framebuffer=None):
        """ Returns a list of tuples of (color, (x, y)) for each pt from pt1 to pt2.
            
            use_brightness - If True, returns color as a scalar 0 - 1 representing brightness.
                             Otherwise, returns the color as represented by PIL.
        """
        
        if framebuffer:
            self.read_buffer(framebuffer)
        else:
            framebuffer= self.fb
        
        pts = self.ptgen.calc_line(pt1[0], pt1[1], pt2[0], pt2[1])
        output = []
        for pt in pts:
            pixel_bit = 0
            if calc_brightness:
                color = self.get_luminance(pt)
            else:
                color = framebuffer.getpixel(pt)
            output.append((color, pt))
        return output
    
    def line_scan_frames(self, pt1, pt2, filterNoise=False, framebuffer=None):
        """
            filternoise is a tuple/ (#places to look forward/back , replace value) 

            sort by contiguous blocks , put in list of lists (each sub list = a length)
            for example 11110010101110000 = [ [1111], [00], [1], [0], [1], [0], [111], [0000] ]

            True is a black pixel 
            stores tuple of tuples - (((True/False), (x,y)), ... )
        """
        
        if framebuffer:
            self.read_buffer(framebuffer)
        else:
            framebuffer= self.fb
            
        if filterNoise:    
            self.log_msg('# scanning with noise filter active')
            pixels = self.line_scan( pt1, pt2, filterNoise   )
        else:
            pixels = self.line_scan( pt1, pt2 )


        changecount = 0
        lastpix = None
        
        data_frames = [] #sort into groups of changes 
        buff = []
        for p in pixels:
            if p[0] == lastpix:
                 buff.append(p)
            if p[0] != lastpix:
                if buff!=[]:
                    data_frames.append(buff)
                buff = [] 
                buff.append(p)           
                changecount += 1
            lastpix = p[0]  

        #grab the last block if it is different  
        if buff != lastpix:
            data_frames.append(buff)
        return data_frames
   
    def circle_scan(self, x_orig, y_orig, dia, framebuffer=None):
        """ orignial scoring tool, looks in a circle for dark pixels from a center point """
        
        if framebuffer:
            self.read_buffer(framebuffer)
        else:
            framebuffer= self.fb

        pts = self.ptgen.calc_circle_2d(x_orig, y_orig, dia)
        pts.append((x_orig, y_orig))  #add the center too 

        is_checked = False
        for pt in pts:
            if ( self.scanner_darkly( framebuffer.getpixel(pt) ) ):
                is_checked = True
        return is_checked



class pixelFilter(object):
    """
       home for various image filters for serialized data, pixels, etc 
    """

    def mean_pix(self, listv):
        """ average the sampled pixel data [(Value, (X,Y))] """
        value = 0
        for v in listv:
            value+=v[0]
        return round( value/len(listv) )

    def filter_noise(self, scandata, filterSize, repval, shift_data=False, bookend=False):
        """
           TODO - whatever the first and last value is gets "carried"  
           
           (like the openCV healing feature - grapefruit example in book)

           basically averages a linear buffer of pixels
        """

        output = []
        lastpixel = 0

        total_cnt = 0
        size_data = len(scandata)-1 
        future_pix = -1 #the sample ahead pixel value (data is pre thresholded to 1 or 0) 

        if shift_data:
            shiftsize = round(filterSize/2)

        for b in scandata:
            filt = None 
            
            #lets look into the future of our data 
            sample_ahead = total_cnt+filterSize
            if sample_ahead<=size_data:
                future_pix = scandata[sample_ahead]
                fwindow =[]
                for w in range(total_cnt, sample_ahead):
                    fwindow.append(scandata[w])
       
                #sample behind now
                pwindow = []
                for w in range(0, -filterSize, -1):
                    pwindow.append( scandata[w])
  
                avgf = self.mean_pix(fwindow)
                avgp = self.mean_pix(pwindow)

                if (round(avgf+avgp/2)):    
                    if shift_data:
                        filt = (1, scandata[total_cnt+shiftsize][1] ) 
                    else:
                        filt = (1, b[1] )            
                else:
                    filt = (repval, b[1] )                   

                output.append(filt)
            if sample_ahead>size_data:
                future_pix = -1 #null value thats not 0 or 1 

            total_cnt+=1
        if bookend:
            return output[bookend[0]:-bookend[1]]
        else:    
            return output
