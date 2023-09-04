# ***** BEGIN GPL LICENSE BLOCK *****
#
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# ***** END GPL LICENCE BLOCK *****

from PIL import Image 


from gnelscript.pygfx.raster_ops import *
from gnelscript.pygfx.point_ops import *

from gnelscript.pygfx.math_ops import math_util as mu
from gnelscript.pygfx.math_ops import NUMPY_IS_LOADED

from gnelscript.pygfx.obj3d import object3d

#from pygfx.math_ops import matrix33, vec2, vec3  




if NUMPY_IS_LOADED:
    # print(' ## debug - loading numpy module in point ops. ')
    import numpy as np  
else:
    print(' ## debug - numpy module disabled in point ops. ')





"""
class render3d(object):  

    def __init__(self, resx=800, resy=600, framebuffer=None):
        self.pg   = polygon_operator()      #use this for the matrix rotation and more
        self.m33  = matrix33()       
        self.m44  = matrix44()
        self.pixop   = pixel_op()   #framebuffer and raster freinds
        
        self.rp   = []          # render path data  (lines)
        self.rpts = []          # render point data (points)

        self.tex_fb = None  #first stab at texture mapping 

        if framebuffer==None:
            self.res = [resx, resy]
            self.pixop.create_buffer(resx, resy)
        
        if framebuffer:
            self.pixop = framebuffer
            self.res = [framebuffer.res_x, framebuffer.res_y]


        ###############################################
        self.render_objects = [] #list of object3d's 

        # renderer properties 
        self.COLOR_MODE           = 'flat'    ####'flat', 'zdepth',  'normal'
        self.is_orthographic      = True
        self.SHOW_VTXS            = True 
        self.render_normal_lines  = True 
        self.SHOW_VEC_HITS        = False  #faces cover this up!
        self.SHOW_EDGES           = True
        self.SHOW_FACES           = True
        self.SHOW_NORMALS         = True  
        self.SHOW_FACE_CENTER     = True
        self.SHOW_SCREEN_CLIP     = False # broken 
        self.DO_SCREEN_CLIP       = False # broken 
 

    ## ## ## ## ## 
    def save_image(self,filename='output.png', noalpha=True):
        self.post_process()
        self.pixop.save(filename, noalpha=noalpha)     
        


    ## ## ## ## ##  
    def project_polygons(self, object3d,  rx, ry, rz, scale, res_x=None, res_y=None):
        
         # rotate and project 3D line geometry into 2D.
         # You can convert this data INTO 3d by adding an empty Z axis.
         # iterate through each point (vector) and mutliply it by a rotation 4X4 matrix

        
        if res_x==None:
            res_x = self.res[0]
        
        if res_y==None:
            res_y = self.res[1]

        center = (int(res_x/2), int(res_y/2) ) #we will need to know center of image when we project geometry into screen space
        
        ############################
        # various matrix ops
        pvtxs = self.m33.rotate_pts_3d( object3d.points, rx, ry, rz )  # tested and works

        #pvtxs = self.m44.rotate_pts_3d( object3d.points, rx, ry, rz )  # tested and works

        # rotate the points by a 3X3 matrix directly 
        #matrix = self.m33.identity
        #pvtxs =  matrix.batch_mult_pts( object3d.points )

        ############################

        lines_to_draw = []
        ##########################
        #print('render geom info  points:%s poly:%s'%(len(object3d.points), len(object3d.polygons) )  )
        ##########################



        #project rotated points into screen space  
        for ply in object3d.polygons:
            num_idx = len(ply) #walk array of indeces to vertices
            for pt in range(num_idx):

                if pt<num_idx-1:
                    #a line needs two points, simply walk through the list two at a time
                    idx  = int(ply[pt])-1 #index start of line in 2d
                    idx2 = int(ply[pt+1])-1 #index end of line in 2d
                    
                    # #start of line
                    x = pvtxs[idx][0] #first vtx - x component  
                    y = pvtxs[idx][1] #first vtx - y component 
                    z = pvtxs[idx][2] #first vtx - z component 

                    # #end of line
                    x2 = pvtxs[idx2][0] #second vtx - x component  
                    y2 = pvtxs[idx2][1] #second vtx - y component 
                    z2 = pvtxs[idx2][2] #second vtx - z component 

                    #attempt at perspective rendering . math, BLAH!
                    if self.is_orthographic==False:
                        #z = (pvtxs[idx][2] + pvtxs[idx2][2])/2
                        # here is my sad attempt at cheapo perspective:
                        #scale = scale -(z*20)   #terrible perspective illusion,  but a really cool effect
                        pass


                    #start of line to draw in 2d 
                    sx =  ((x*scale) +center[0])  
                    sy =  ((y*scale) +center[1])  
                    #end of line to draw in 2d
                    ex =  ((x2*scale)+center[0])  
                    ey =  ((y2*scale)+center[1])   

                    ############################
                    if self.DO_SCREEN_CLIP:
                        # clip lines to 2D screen size  
                        cxmin = center[0]-self.clip_x[0] 
                        cxmax = center[0]+self.clip_x[1] 
                        
                        cymin = center[1]-self.clip_y[0] 
                        cymax = center[1]+self.clip_y[1] 

                        if sx < cxmin:
                            sx = cxmin
                        if sx > cxmax:
                            sx = cxmax

                        if ex < cxmin:
                            ex = cxmin
                        if ex > cxmax:
                            ex = cxmax

                        if sy < cymin:
                            sy = cymin
                        if sy > cymax:
                            sy = cymax

                        if ey < cymin:
                            ey = cymin
                        if ey > cymax:
                            ey = cymax

                    lines_to_draw.append(  ( (sx,sy), (ex, ey) ) )
   
        return lines_to_draw
"""


##########################################################################################################
##########################################################################################################


class simple_render(object):  
    """         
        orthographic projection of 3D data. 
        Images can be rendered with PIL or dumped out as vectors. 

        The vectors are 2D in 3d space by adding 0 as the Z axis. 
        The "pixels" are coordinates so vector renders can be quite large.
        They need to be scaled and centered around the zero axis. 
        
    """

    def __init__(self, resx=800, resy=600, pixop=None):
        self.pg   = polygon_operator()      #use this for the matrix rotation and more
        self.m33  = matrix33()       
        self.m44  = matrix44()
        self.pixop   = pixel_op()   #framebuffer and raster freinds
        
        self.rp   = []          # render path data  (lines)
        self.rpts = []          # render point data (points)
        
        #store render output to a vector framebuffer - a 2D "3D" object 
        self.vec_fr_buffer = object3d()

        self.uvx = 0 #pointer for UV map X axis  
        self.uvy = 0 #pointer for UV map X axis
        self.ptgen = point_operator_2d()

        half_x = resx/2 
        half_y = resy/2 

        #centered to screen - so use HALF of res        
        self.clip_x = [0,half_x]
        self.clip_y = [0,half_y]

        if pixop==None:
            self.res = [resx, resy]
            self.pixop.create_buffer(resx, resy)
        
        if pixop:
            if type(pixop)==pixel_op:
                self.pixop = pixop
                self.res = [pixop.res_x, pixop.res_y]

        self.render_objects = []   # list of object3d's 
        self.lighting_vectors = [] # data about lighting - to visualize upstream from here

        ###############################################
        # renderer properties 
        self.SHOW_VTXS            = True 
        self.render_normal_lines  = True 
        self.SHOW_VEC_HITS        = False  #faces cover this up!
        self.SHOW_EDGES           = True
        self.SHOW_FACES           = True
        self.SHOW_NORMALS         = True  
        self.SHOW_FACE_CENTER     = True

        self.USE_PERSPECTIVE      = False # experimental -         
        self.SHOW_SCREEN_CLIP     = False # broken 
        self.DO_SCREEN_CLIP       = False # broken 
        
        self.STORE_VECTOR_RENDER  = True   


        self.COLOR_MODE           = 'flat'    ####'flat', 'zdepth',  'normal' 
        ###############################################

    ##--------------------------------## 
    def save_image(self,filename='output.png', noalpha=True):
        self.post_process()
        self.pixop.save(filename, noalpha=noalpha)     

    def post_process(self):
        if self.SHOW_SCREEN_CLIP:
            center = (int(self.pixop.res_x/2), int(self.pixop.res_y/2) )

            cxmin = center[0]-self.clip_x[0] 
            cxmax = center[0]+self.clip_x[1] 
            cymin = center[1]-self.clip_y[0] 
            cymax = center[1]+self.clip_y[1] 

            ## square = [ (cxmin,cymin),(cxmin,cymax),
            ##            (cxmin,cymax),(cxmax,cymax),
            ##            (cxmax,cymax),(cxmax,cymin),                       
            ##            (cxmax,cymin),(cxmin,cymin)                       
            ## ]


            self.pixop.connect_the_dots( [(cxmin,cymin),(cxmax,cymax)] , (0,255,0), 2  )   
    
        #PIL USES TOP LEFT FOR (0,0) SO WE NEED TO FLIP TO CORRECT 
        self.pixop.fb = self.pixop.fb.transpose(Image.FLIP_LEFT_RIGHT) 


    ##--------------------------------##  
    def project_points(self, object3d,  rx, ry, rz, scale, res_x=None, res_y=None):
        """
           project 3D point geometry into 2D
        """

        #I dont think these belong here  
        if res_x==None:
            res_x = self.res[0]
        if res_y==None:
            res_y = self.res[1]

        center = (int(res_x/2), int(res_y/2) ) #we will need to know center of image when we project geometry into screen space
        
        ##############################
        # various matrix ops 
        pvtxs = self.m33.rotate_pts_3d( object3d.points, rx, ry, rz )  #3X3 matrix - tested and works

        # rotate with a 4X4 matrix 
        #pvtxs = self.m44.rotate_pts_3d( object3d.points, rx, ry, rz )  #4X4 matrix - tested and works
        
        # rotate the points by a 3X3 matrix directly 
        #matrix = self.m33.identity
        #pvtxs =  matrix.batch_mult_pts( object3d.points )
        
        ##############################

        points_projected = []

        #crappy orthographic projection (no Z at all) 
        ## for p in pvtxs:
        ##     sx = (p[0]*scale) + center[0]  
        ##     sy = (p[1]*scale) + center[1]  
        ##     points_projected.append( (sx,sy) )
   
        for p in pvtxs:
            sx = p[0]  
            sy = p[1] 
            sz = p[2]   

            #dont think this works - debug 
            if self.USE_PERSPECTIVE:
                sx = sx/sz 
                sy = sy/sz

            points_projected.append( ((sx*scale)+center[0] , (sy*scale)+center[1] ) )

        return points_projected

    ##--------------------------------##   
    def project_polygons(self, object3d,  rx, ry, rz, scale, res_x=None, res_y=None):
        """
          rotate and project 3D line geometry into 2D.

          You can convert this data INTO 3d by adding an empty Z axis.
             - - - 
          iterate through each point (vector) and mutliply it by a rotation 4X4 matrix

        """
        if res_x==None:
            res_x = self.res[0]
        
        if res_y==None:
            res_y = self.res[1]

        center = (int(res_x/2), int(res_y/2) ) #we will need to know center of image when we project geometry into screen space
        
        ############################
        # various matrix ops
        #pvtxs = self.m33.rotate_pts_3d( object3d.points, rx, ry, rz )  # tested and works

        pvtxs = self.m44.rotate_pts_3d( object3d.points, rx, ry, rz )  # tested and works

        # rotate the points by a 3X3 matrix directly 
        #matrix = self.m33.identity
        #pvtxs =  matrix.batch_mult_pts( object3d.points )

        ############################

        lines_to_draw = []
        ##########################
        #print('render geom info  points:%s poly:%s'%(len(object3d.points), len(object3d.polygons) )  )
        ##########################

        #project rotated points into screen space  
        for ply in object3d.polygons:
            num_idx = len(ply) #walk array of indeces to vertices
            for pt in range(num_idx):

                if pt<num_idx-1:
                    #a line needs two points, simply walk through the list two at a time
                    idx  = int(ply[pt])-1 #index start of line in 2d
                    idx2 = int(ply[pt+1])-1 #index end of line in 2d
                    
                    if idx>len(pvtxs)-1:
                        raise ValueError("\n\nbad index '%s'- check your data "%ply[pt])
                    if idx2>len(pvtxs)-1:
                        raise ValueError("\n\nbad index '%s' - check your data"%ply[pt+1]) 

                    # # #start of line
                    # x = pvtxs[idx][0] #first vtx - x component  
                    # y = pvtxs[idx][1] #first vtx - y component 
                    # z = pvtxs[idx][2]/10 #first vtx - z component 
                    # # #end of line
                    # x2 = pvtxs[idx2][0] #second vtx - x component  
                    # y2 = pvtxs[idx2][1] #second vtx - y component 
                    # z2 = pvtxs[idx2][2]/10 #second vtx - z component 

                    # #start of line
                    x = pvtxs[idx][0] #first vtx - x component  
                    y = pvtxs[idx][1] #first vtx - y component 
                    z = pvtxs[idx][2]/10 #first vtx - z component 

                    # #end of line
                    x2 = pvtxs[idx2][0] #second vtx - x component  
                    y2 = pvtxs[idx2][1] #second vtx - y component 
                    z2 = pvtxs[idx2][2]/10 #second vtx - z component 

                    
                    #first shot at perspective - this is a good start 
                    if self.USE_PERSPECTIVE:
                        x = x/z 
                        y = y/z
                        x2 = x2/z2 
                        y2 = y2/z2 


                    ###################################################

                    # print('## Z IS ', z , z2 )

                    # orthographic - NO Z coordinate  
                    ## #start of line to draw in 2d 
                    ## sx =  ((x*scale) +center[0])  
                    ## sy =  ((y*scale) +center[1])  
                    ## #end of line to draw in 2d
                    ## ex =  ((x2*scale)+center[0])  
                    ## ey =  ((y2*scale)+center[1])   
  

                    #start of line to draw in 2d 
                    sx =  ((x*scale) +center[0])  
                    sy =  ((y*scale) +center[1])  
                    #end of line to draw in 2d
                    ex =  ((x2*scale)+center[0])  
                    ey =  ((y2*scale)+center[1])   

                    ############################
                    if self.DO_SCREEN_CLIP:
                        # clip lines to 2D screen size  
                        cxmin = center[0]-self.clip_x[0] 
                        cxmax = center[0]+self.clip_x[1] 
                        
                        cymin = center[1]-self.clip_y[0] 
                        cymax = center[1]+self.clip_y[1] 

                        if sx < cxmin:
                            sx = cxmin
                        if sx > cxmax:
                            sx = cxmax

                        if ex < cxmin:
                            ex = cxmin
                        if ex > cxmax:
                            ex = cxmax

                        if sy < cymin:
                            sy = cymin
                        if sy > cymax:
                            sy = cymax

                        if ey < cymin:
                            ey = cymin
                        if ey > cymax:
                            ey = cymax

                    lines_to_draw.append(  ( (sx,sy), (ex, ey) ) )
                
                # ###################################
                # # test to draw last line of poly 
                # if pt==num_idx-1:
                #     idx  = int(ply[pt])-1   # index start of line in 2d
                #     idx2 = int(ply[0])-1   # index end of line in 2d
                #     # start of line
                #     x = pvtxs[idx][0] #first vtx - x component  
                #     y = pvtxs[idx][1] #first vtx - y component 
                #     z = pvtxs[idx][2]/10 #first vtx - z component 
                #     # end of line
                #     x2 = pvtxs[idx2][0] #second vtx - x component  
                #     y2 = pvtxs[idx2][1] #second vtx - y component 
                #     z2 = pvtxs[idx2][2]/10 #second vtx - z component 
                #     #start of line to draw in 2d 
                #     sx =  ((x*scale) +center[0])  
                #     sy =  ((y*scale) +center[1])  
                #     #end of line to draw in 2d
                #     ex =  ((x2*scale)+center[0])  
                #     ey =  ((y2*scale)+center[1])   
                #     lines_to_draw.append(  ( (sx,sy), (ex, ey) ) )

        return lines_to_draw

    ##--------------------------------## 
    def render_obj (self, color, rx, ry, rz, thick, scale, framebuffer=None, object3d =None):
        """ 
            render a single object 
        """
       
        if not object3d.points:
            print('## errror - object has no point geometry.')
            return None

        #output image properties 
        res_x = self.res[0]
        res_y = self.res[1]
       
        ###########################
        # optional framebuffer passed in so we can render on top of other images
        if framebuffer:
            if type(framebuffer==pixel_op):
                rndr_bfr = framebuffer
        else:
            # default is to just make a new framebuffer from scratch
            rndr_bfr = self.pixop 
            #rndr_bfr.create_buffer(res_x, res_y)#make a new image in memory
            #rndr_bfr.fill_color( (25,20,25) ) #make bg all dark 

        ###########################
        # render line geometry  
        if self.SHOW_EDGES:
            self.rp = self.project_polygons(object3d, rx, ry, rz, scale, res_x, res_y)

            for i,l in enumerate(self.rp):
                
                #cool idea - needs more work to implement 
                #if object3d.linecolors:
                #    color = object3d.linecolors[i]
                
                rndr_bfr.connect_the_dots( l, color, int(thick/2) )  #points, color, thickness

                if self.STORE_VECTOR_RENDER:    
                    self.vec_fr_buffer.insert_line(l)


        ###########################
        # render point geometry 
        if self.SHOW_VTXS:
            self.rpts = self.project_points(object3d, rx, ry, rz, scale, res_x, res_y)
            rndr_bfr.draw_points_batch( self.rpts ,  (255,255,0) , int(thick)        )  #points, color, thickness
   
    ##--------------------------------## 
    def render_multiobj(self, color, rx, ry, rz, thick, scale):
        if not self.render_objects:
            print ('error no objects to render')
            return None        

        #rebuild a new framebuffer with each "render pass" 
        #if you comment these 3 lines out, it renders over the old one each pass 
        self.pixop = pixel_op()
        self.pixop.create_buffer(self.res[0], self.res[1])#make a new image in memory
        # self.pixop.fill_color( (25,20,25) ) #make bg all dark

        if self.render_objects:
            for obj in self.render_objects:
                self.render_obj(color, rx, ry, rz, thick, scale, framebuffer=self.pixop, object3d=obj) 

    ##--------------------------------## 
    def anim(self, objs, init_rots=(0,0,0), linethick=5, numframes=5, scale=150):
        """
            DEBUG - add interpolation for "keyframes"

            use imagemagik to make animated gif:
                convert loop 1 *.png spin.gif
        
        """

        outfolder   = 'anim'
        outfilename = 'my3d'

        output_type = 'png'

        step_degrees = 5
        
        #linethick = 5

        RX=init_rots[0];RY=init_rots[1];RZ=init_rots[2]

        ### single object example 
        # obj = polygon_operator() 
        # obj.prim_cube() 
        # for f in range(numframes):
        #     self.render_obj((255,0,0), RX, RY+(f*step_degrees), RZ, 3, 100,  object3d=obj ) 
        #     self.save( '%s/%s_%s.png'%(outfolder,outfilename,f) )

        self.render_objects = objs 
        
        # #obj2 = polygon_operator() 
        # #obj2.prim_locator() 
        # #self.render_objects.append(obj2) 
        # obj.prim_cube(linecolor=(0,255,0), size=.3,pos=(1,0,0)) 
        # self.render_objects.append(obj )

        for f in range(numframes):
            self.render_multiobj((255,0,0), RX, RY+(f*step_degrees), RZ, linethick, scale ) 
            self.save_image( '%s/%s_%s.%s'%(outfolder,outfilename,f,output_type) )

    ##--------------------------------## 
    def render_matrix_obj (self,  m33, m44, thick, scale, filename, object3d =None):
        """ high level wrapper to call render_custom_matrix """

        color = (0,0,255)
        res_x = self.res[0];res_y = self.res[1]
       
        # default is to just make a new framebuffer from scratch
        rndr_bfr = self.pixop 
        rndr_bfr.create_buffer(res_x, res_y)#make a new image in memory
        rndr_bfr.fill_color( (25,20,25) ) #make bg all dark 

        ###########################
        ## render line geometry  
        if self.SHOW_EDGES:
            render_data = self.render_custom_matrix(object3d, scale, m33, m44, res_x, res_y)
            for i,l in enumerate(render_data[1]):
                if object3d.linecolors:
                    color = object3d.linecolors[i]
                rndr_bfr.connect_the_dots( l, color, int(thick/2) )  #points, color, thickness
        ###########################
        ## render point geometry 
        if self.SHOW_VTXS:
            rndr_bfr.draw_points_batch( render_data[0] ,  (255,255,0) , int(thick)        )  #points, color, thickness

        rndr_bfr.save(filename) 

    ##--------------------------------##   
    def render_custom_matrix(self, object3d, scale, m33=None, m44=None, res_x=None, res_y=None):
        """ 
            Utility for viewing matricies. Pass in a 3X3, 4X4, ( point_ops matrix or numpy ndarray )
            This will render a single object using that matrix for rotation. 

            will render a 3X3 OR a 4X4 but not both at the same time.
        """
        if res_x==None:
            res_x = self.res[0]
        if res_y==None:
            res_y = self.res[1]

        center = (int(res_x/2), int(res_y/2) ) 

        ########################################################################## 
        ########################################################################## 
        if m33 is not None and m44 is not None:
            print('Only one size of matrix can be used at a time!! \n\n')
            return None 

        # rotate the points of the object by multiplying by a matrix 
        if m33 is not None:
            if NUMPY_IS_LOADED:
                if isinstance(m33, np.ndarray):
                    tmp = matrix33()
                    tmp.insert(m33)
                    pvtxs =  tmp.batch_mult_pts(object3d.points ) 
            if isinstance(m33, matrix33):                
                pvtxs =  m33.batch_mult_pts(object3d.points )          
           
        if m44 is not None:
            if NUMPY_IS_LOADED:            
                if isinstance(m44, np.ndarray):
                    tmp = matrix44()
                    tmp.insert(m44)
                    pvtxs =  tmp.batch_mult_pts(object3d.points ) 
            if isinstance(m44, matrix44):                
                pvtxs =  m44.batch_mult_pts(object3d.points )  

        ########################################################################## 
        ########################################################################## 
        points_to_draw  = []
        lines_to_draw   = []

        #project rotated lines into screen space  
        for ply in object3d.polygons:
            num_idx = len(ply) #walk array of indeces to vertices
            for pt in range(num_idx):

                if pt<num_idx-1:
                    # a line needs two points, simply walk through the list two at a time
                    idx  = int(ply[pt])-1   #index start of line in 2d
                    idx2 = int(ply[pt+1])-1 #index end of line in 2d
                    
                    # start of line in 3D
                    x = pvtxs[idx][0] #first vtx - x component  
                    y = pvtxs[idx][1] #first vtx - y component 
                    # Z is not being used for this ultra simple  projection 

                    # end of line in 3D
                    x2 = pvtxs[idx2][0] #first vtx - x component  
                    y2 = pvtxs[idx2][1] #first vtx - y component 
                    # Z is not being used for this ultra simple  projection 

                    # start of line to draw in 2D 
                    sx =  ((x*scale) +center[0])  
                    sy =  ((y*scale) +center[1])  
                    # end of line to draw in 2D
                    ex =  ((x2*scale)+center[0])  
                    ey =  ((y2*scale)+center[1])    
                    lines_to_draw.append( ( (sx,sy), (ex, ey) ) )
        
        #project rotated points to screen space 
        for p in pvtxs:
            sx = (p[0]*scale) + center[0]  
            sy = (p[1]*scale) + center[1]  
            points_to_draw.append( (sx,sy) )

        return [points_to_draw,lines_to_draw] #points and lines at same time 

    ##--------------------------------## 
    def sort_polys(self, obj):
        """ sorting logic is in the object3d class, 
            but this is where it gets called prior to rendering. 

            some things we do:
                check the number of faces,  
                triangulate if needed
                reindexing 
                zsorting 

        return [points[], polygons[] ]
        """
       
        polydata = []

        #obj.triangulate() 

        n_plys = len(obj.polygons)
        
        """
            #buggy! 

              File "XXXX/gnelscript/pygfx/point_ops.py", line 394, in z_sort
                tmp.sort(reverse=True)
            TypeError: '<' not supported between instances of 'list' and 'tuple'
       
        """

        # flip them becasue the perspective mangles the Z coords 
        if self.USE_PERSPECTIVE:
            obj.z_sort()  
        else:
            obj.z_sort(reverse=True)  

        polydata.append( obj.points       )
        polydata.append( obj.polygons     )
        polydata.append( obj.normals )

        print('## num points %s; num polygons %s; num normals %s '%(len(obj.points) , len(obj.polygons), len(obj.normals) ) )
        return polydata


    ##--------------------------------## 
    def paint_line(self, points, color_fb, framebuffer=None, xoffset=0, yoffset=0, bright=1):
        """ paint a line of pixels from an image  """

        if type(framebuffer)==pixel_op:
            dpix = framebuffer.fb.load() 
        else:
            raise ValueError("paint_line wrong object passed as framebuffer ")

        p1 = list(points[0])    
        p2 = list(points[1])  

        pts = self.ptgen.calc_line(p1[0], p1[1], p2[0], p2[1])
        
        pxct = 0

        for pt in pts:
            try:
                #mix luminance with color from texture
                pre_color = list(color_fb.get_pix( (pxct+xoffset,yoffset)) )
                
                pre_color[0] = int(pre_color[0] - bright) # R
                pre_color[1] = int(pre_color[1] - bright) # G
                pre_color[2] = int(pre_color[2] - bright) # B

                dpix[ pt[0], pt[1] ] = tuple(pre_color) 

                #or just use pure texture with full luminance 
                #dpix[ pt[0], pt[1] ] = color_fb.get_pix( (pxct+xoffset,yoffset)) 

            except:
                pass
            pxct += 1 

    ##--------------------------------## 
    def scanline(self, obj, scale=200, lightpos=(0,10,0) , texmap=None):
        """ 
            polydata = [ points[], polygons[], normals[] ]
        """
        vecmath = vec2()     # use for math operations

        polydata = self.sort_polys(obj) #pre-process polys before rendering 

        output = self.pixop

        output.fill_color((22,22,22))

        res_x = self.res[0]
        res_y = self.res[1]
        center = (self.res[0]/2, self.res[1]/2)
        
        light_intensity  = .6 

        if texmap is not None:
            print('loading texmap! ', texmap)
            self.tex_fb = texmap

        #mark the 0 point for convenience
        #output.horiz_line(res_y/2, (255,0,255) ) 
        #output.vert_line( res_x/2, (255,0,255) ) 

        ###############################################
         
        # polydata is [ points, polygons, normals ] 
        for ply in polydata[1]:
            num_idx = len(ply) # number of vertices per poly 
            drwply = []        # 3 points of triangle to draw
            
            self.uvx = 0
            self.uvy = 0

            #only look at 3 and four sided polys 
            if num_idx==3 or num_idx == 4:

                #  #pointer for UV map X axis  
                self.uvy += 1 #pointer for UV map X axis

                # DEBUG, if face is 4 sided - we need to triangulate it 

                for pt in range(num_idx):
                    idx = int(ply[pt])  
                    drwply.append( polydata[0][idx-1] )
               

                # ///////////////////////////////////// 
                #test of perspective with scanline - debug 
                if self.USE_PERSPECTIVE:
                    s1 = ( (drwply[0][0]/drwply[0][2]*scale)+center[0], (drwply[0][1]/drwply[0][2]*scale)+center[1] )
                    e1 = ( (drwply[1][0]/drwply[1][2]*scale)+center[0], (drwply[1][1]/drwply[1][2]*scale)+center[1] )
                    s2 = ( (drwply[1][0]/drwply[1][2]*scale)+center[0], (drwply[1][1]/drwply[1][2]*scale)+center[1] )
                    e2 = ( (drwply[2][0]/drwply[2][2]*scale)+center[0], (drwply[2][1]/drwply[2][2]*scale)+center[1] )
                    s3 = ( (drwply[2][0]/drwply[2][2]*scale)+center[0], (drwply[2][1]/drwply[2][2]*scale)+center[1] )
                    e3 = ( (drwply[0][0]/drwply[0][2]*scale)+center[0], (drwply[0][1]/drwply[0][2]*scale)+center[1] )
                else: 
                    # build up line data for three sides of triangle
                    s1 = ( (drwply[0][0]*scale)+center[0], (drwply[0][1]*scale)+center[1] )
                    e1 = ( (drwply[1][0]*scale)+center[0], (drwply[1][1]*scale)+center[1] )
                    #
                    s2 = ( (drwply[1][0]*scale)+center[0], (drwply[1][1]*scale)+center[1] )
                    e2 = ( (drwply[2][0]*scale)+center[0], (drwply[2][1]*scale)+center[1] )
                    #
                    s3 = ( (drwply[2][0]*scale)+center[0], (drwply[2][1]*scale)+center[1] )
                    e3 = ( (drwply[0][0]*scale)+center[0], (drwply[0][1]*scale)+center[1] )
      
                ##----------------------------------##      
                #edge lines
                l1 = [(s1[0], s1[1]), (e1[0], e1[1])]  
                l2 = [(s2[0], s2[1]), (e2[0], e2[1])]  
                l3 = [(s3[0], s3[1]), (e3[0], e3[1])]  

                ##----------------------------------##   
                #add some color to the polygons here  
                    
                if self.COLOR_MODE=='normal':  
                    # DEBUG 
                    # this whole thing seems broken
                    # it would be better to calc the normal while we have the face ! 
                    
                    """
                    # obj.calc_face_normals()

                    n2 = polydata[2][idx-1]
                    n2 = vec3(n2[0],n2[1],n2[2])
                    n2 = (  n2.angle_between( vec3(0,0,1)) ) 
                    angle  = mu().rtd(n2 )
                    
                    #facecolor = output.normal_to_color( (n2,n2,n2) )
                    if angle > 160:
                        facecolor = ( 255,0,0 )
                    else:    
                        facecolor = ( 0,0,255 )
                    """
                    facecolor = (128,0,0)

                if self.COLOR_MODE=='zdepth':                 
                    #COLOR BY Z DISTANCE FROM CAMERA
                    zgrad = int(drwply[0][2]*200)
                    facecolor = (zgrad,zgrad,zgrad )
                
                if self.COLOR_MODE=='flat':                 
                    facecolor = (100,100,100 )
           

                if self.COLOR_MODE=='lighted' or self.COLOR_MODE=='lightedshaded':  
                    """ light or dark depending on angle to a point (light) 
                         -- once that works add:
                              light intesity 
                              light color 
                              multiple lights 

                    """ 
                    if isinstance(lightpos, tuple):
                        # store position of "light" in a vec3 object 
                        # it is a position, not a vector but this allows the tools to be used on it 
                        lightpos = vec3(lightpos[0],lightpos[1],lightpos[2])   
                    
                    # get the center of face to move light vector to 
                    f_cntr = obj.centroid(drwply)
                    
                    # put face centroid point into a vec3 object 
                    fcntr = vec3(f_cntr[0],f_cntr[1],f_cntr[2])
                    
                    # build the face normal vector
                    nrml =  obj.calc_tripoly_normal( drwply, unitlen=False) 
                    
                    # build a vector between face center and light position  
                    vec_to_light = fcntr.between(lightpos) 
     
                    #--------
                    # calculate the angle between face and light vectors 
                    # treating the 3D light position as a vector
                    light_angle = (  nrml.angle_between( vec_to_light ) )
                    angle  = int(mu().rtd(light_angle ))
                    # this is a terrible way to do it, but it looks pretty good for a first try
                    light_pow = int(light_intensity*256)
                    facecolor = (light_pow-angle,light_pow-angle,light_pow-angle)

                    #--------
                    # store lighting vectors that were calulated so we can play with them later  
                    self.lighting_vectors.append( [fcntr, nrml, vec_to_light, angle, ply ] )  

                ##########
                # define the scanline geometry, iterate each horizontal line of image
                for hscan in range(1,res_x):
                    s_hvec = (-1*(res_x/2),hscan)
                    e_hvec = (       res_x,hscan)

                    # take the 3 edges of a triangle and determine if the horizontal scanline intersects any 
                    i = vecmath.intersect(s_hvec, e_hvec, s1, e1) #left  side of triangle
                    j = vecmath.intersect(s_hvec, e_hvec, s2, e2) #right side of triangle
                    k = vecmath.intersect(s_hvec, e_hvec, s3, e3) #top   side of triangle

                    if self.SHOW_VEC_HITS:
                        # debug tool - show the "hits" for the horizontal scanline
                        if i: 
                            output.draw_fill_circle( i[0], i[1], 1, (255,0,0) ) 
                        if j: 
                            output.draw_fill_circle( j[0], j[1], 1, (0,255,0) )  
                        if k: 
                            output.draw_fill_circle( k[0], k[1], 1, (0,0,255) )  

                    
                    ##--------------------------------## 
                    ## fill a polygon - test of texture mapping   
                    if self.SHOW_FACES:
                       
                        # get pixel color value from texture map 
                        # self.uvx and self.uvy are iterators that keep track of pixels drawn

                        # not a "true" UV coordinate, but an offset in image space
                        vtx_u = .0 
                        vtx_v = .0 

                        tex_u = int(self.tex_fb.size[0] * vtx_u) + self.uvx
                        tex_v = int(self.tex_fb.size[1] * vtx_v) + self.uvy

                        pix_clr=(0,0,0)
                        #if self.COLOR_MODE=='lighted': 
                        #    pix_clr = self.tex_fb.get_pix((tex_u, tex_v))                        

                        #add lighting into to pixel color 
                        #pix_clr = (pix_clr[0]-angle, pix_clr[1]-angle, pix_clr[2]-angle)

                        #print('## color x %s y %s is '%(self.uvx, self.uvy),pix_clr)
                        
                        ######################
                        # def paint_line( points, color_fb, framebuffer=None,):
                        if self.COLOR_MODE=='lightedshaded': 
                            
                            luminance = (256-angle)
                         
                            if i and j: 
                                drawlin = [i,j]
                                self.paint_line( drawlin, self.tex_fb , output, tex_u, tex_v, luminance) 
                                self.uvy += 1
                            if i and k:
                                drawlin = [i,k]
                                self.paint_line( drawlin, self.tex_fb , output, tex_u, tex_v, luminance)                             
                                self.uvy += 1                            
                            if j and k:
                                drawlin = [j,k]
                                self.paint_line( drawlin, self.tex_fb , output, tex_u, tex_v, luminance)                             
                                self.uvy += 1

                        ######################
                        if self.COLOR_MODE=='lighted': 
                            #ineffecient! why draw the whole vertical sweep ?
                            # should only go top to bottom of polygon 
                            if self.SHOW_FACES:
                                if i and j: 
                                    drawlin = [i,j]
                                    output.connect_the_dots( drawlin, facecolor, 1)
                                if i and k:
                                    drawlin = [i,k]
                                    output.connect_the_dots( drawlin, facecolor, 1)
                                if j and k:
                                    drawlin = [j,k]
                                    output.connect_the_dots( drawlin, facecolor, 1) 


                ##################                
                if self.SHOW_NORMALS:
                    #draw face normal 
                    pass

                if self.SHOW_FACE_CENTER:
                    #draw face normal 
                    cn =  obj.centroid(drwply) 
                    cntr_pt = ( (cn[0]*scale)+center[0], (cn[1]*scale)+center[1] ) 
                    output.draw_fill_circle( cntr_pt[0], cntr_pt[1], 2, (255,255,20) )

                if self.SHOW_EDGES:
                    #draw polygon edges 
                    output.connect_the_dots( l1, (0,255,0), 1) 
                    output.connect_the_dots( l2, (0,255,0), 1) 
                    output.connect_the_dots( l3, (0,255,0), 1)                 
        












