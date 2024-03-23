
import os 
import sys


from gnelscript.pygfx.math_ops import  *
from gnelscript.pygfx.point_ops import *
from gnelscript.pygfx.obj3d import  *
from gnelscript.pygfx.render import *

#from gnelscript.pygfx.vector_ops import *

# from gnelscript.pygfx.kicad_ops import * 
# from gnelscript.pygfx.milling_ops import * 



from gnelscript.examples.selection import * 
from gnelscript.examples.vector import *
from gnelscript.examples.wip import *
from gnelscript.examples.milling import *
#from gnelscript.examples.raster import *
#from gnelscript.examples.render import *



from gnelscript.tools.imagecam import * 
from gnelscript import GLOBAL_PROJ



from gnolinker import * 


mu = math_util() 

""" 

    TO RUN FROM COMMAND LINE 

    python3 pycore.py 3d_obj/sphere.obj runcommand 



    netstat -na | grep 2864


"""




##***********************************************************##

## print(sys.argv, len(sys.argv))

if __name__=="__main__":
    
    if len(sys.argv) <2:
        print("no arguments to command line pycore\n")
        exit()
    
    if len(sys.argv) > 2:
        if not os.path.exists(sys.argv[1]):
            print("pycore warn: "+sys.argv[1] +" does not exist \n"  )
            #exit() 

    PYCORE_OBJ_IN   = sys.argv[1]
    PYCORE_GEOMPATH = "3d_obj"
    PYCORE_OBJ_OUT  = "%s/%s"%(PYCORE_GEOMPATH, "PYCORE.obj")

    PYCORE_BMP_OUT  = "py_render.bmp"
    M44_DISK_FILE = "camera_matrix.olm"
    # print("# PYCORE %s --> %s "% (PYCORE_OBJ_IN, PYCORE_OBJ_OUT) )

    print("\n\n\n\n")
    print("### PYCORE INPUT %s"% PYCORE_OBJ_IN)


"""
print("## DEBUG PATHS ")  
print("## PYCORE_OBJ_IN   ", PYCORE_OBJ_IN )
print("## PYCORE_GEOMPATH ", PYCORE_GEOMPATH )
print("## PYCORE_OBJ_OUT  ", PYCORE_OBJ_OUT )
print("## PYCORE_BMP_OUT  ", PYCORE_BMP_OUT )
print("## M44_DISK_FILE   ", M44_DISK_FILE )
"""



##----------------------------------------------------
##   /usr/local/opt/python@3.10/bin/python3.10 ./imagecam.py  

#firstpass_bw(10, 1.5, 1.5, 1, 250, "images/in/art.jpg", "images/out", "output")

## (iteration , blur , contrast, bright, scaling(divs) , in, out )
#firstpass(10, 0, 1, 1, 250, "images/in/oil.png", "images/out", "output")

#firstpass(10, 1.5, 1.2, 1, 250, "images/in/er.jpg", "images/out", "output")

##----------------------------------------------------
##   /usr/local/opt/python@3.10/bin/python3.10 ./imagecam.py  
#secondpass("images/out/art.bmp", "images/out" , 8, False)

##----------------------------------------------------
#set the RGB values from last tool and run this 
#thirdpass( "images/out/commonbands.png", "images/out", "dxf" )
#thirdpass( "images/out/commonbands.png",  "images/out" , "geojson")

#thirdpass( "images/out/commonbands.png",  "images/out" , "dxf")
#thirdpass( "images/out/art.bmp",  "images/out" , "dxf", po_invert=True)



##***********************************************************##
##***********************************************************##
##***********************************************************##
##***********************************************************##
## DEFINE PYCORE COMMANDS (FROM PYGFX)




def loadkicad():
    kicad = pcbfile()

    #kicad.load_gcode('gcode/ngc/3D_Chips.ngc')

    #kicad.load_kicadpcb('gcode/kicad/sample1.kicad_pcb')
    #kicad.load_kicadpcb('gcode/kicad/zipper.kicad_pcb')
    #kicad.load_kicadpcb('gcode/kicad/simple.kicad_pcb')
    
    kicad.load_kicadpcb('gcode/kicad/cnc1.kicad_pcb')
    
    # kicad.show_geom()
    kicad.save_3d_obj(PYCORE_OBJ_OUT) 
 

def triangulate():
    obj = object3d()
    obj.load(PYCORE_OBJ_IN)
    obj.triangulate()
    obj.save(PYCORE_OBJ_OUT)


##------------------

def gen_normals():
    obj = object3d() 
    obj.load(PYCORE_OBJ_IN)

    vectors = [] 

    # iterate all faces, calculate a normal for it 
    # transfrom that normal vector to the 3D center of the face 
    for fid in range(len(obj.polygons)):
        nrml = obj.get_face_normal(fid) 
        cntr = obj.get_face_centroid(fid)
        vectors.append( (nrml, cntr) )

    obj2 = object3d()    
    obj2.vectorlist_to_obj(vectors)
    obj2.save(PYCORE_OBJ_OUT) 


##------------------

def object_rotate():
    obj = object3d()
    obj.load(PYCORE_OBJ_IN)
    #pts = [(2,2,2), (4,4,4), (8,8,8)]
    pts2 = obj.rotate_pts((45,45,45) )
    #print(pts2)
    obj.save(PYCORE_OBJ_OUT)


##------------------

def copy_sop():
    #be cautious of large number of polys. It gets slow real quick!

    obj = object3d() 
    obj.load(PYCORE_OBJ_IN)
    obj.copy_sop(span=(1,10), offset=(0,1,0), num=3, distance=.5)
    obj.save(PYCORE_OBJ_OUT)




##------------------


def modify_partial():
    obj = object3d()
    obj.load(PYCORE_OBJ_IN)
    
    geom = obj.sub_select_geom( span=(40,120) )


    newpts = obj.points

    obj2 = object3d() 
    obj2.insert_polygons(geom[0], newpts  )      
    obj2.rotate_pts((0,45,0))
    obj2.save(PYCORE_OBJ_OUT)    



##------------------
def circle_cube_pts():

    # PYCORE_OBJ_IN not needed 
    obj = object3d()
    obj.prim_circle(axis='z', pos=(0,0,0), spokes=20, dia=.5) 
    
    #obj.scale_pts((1,5,1))
        
    obj.triangulate(force=True)


    for f in range(1,len(obj.polygons)-1,1):
        if f%2==0:
            obj.extrude_face(f, .1)
    
    #obj.scale_pts((1,5,1))
    
    #obj.rotate_pts((35,35,0))

    #obj.extrude_face(len(obj.polygons)-2, -2)

    ## pts = obj.get_face_pts(3) 
    ## ct = 0
    ## for pt in pts:
    ##     tmp = object3d()
    ##     
    ##     #broken # tmp.prim_sphere(pos=pt, rot=(0,0,0), size=.05 )
    ##     #tmp.prim_locator_xyz(pos=pt, rot=(0,0,0), size=1)
    ##     #tmp.prim_cone('y', pos=pt, rot=(0,0,0), dia=.06, spokes=8)
    ##     tmp.prim_cube(size=.1, pos=pt, rot=(ct,ct,ct), pivot='world')
    ##     ct += 10
    ##     obj.insert(tmp)  

    #obj.triangulate(force=True)
    obj.save(PYCORE_OBJ_OUT)

##------------------
def sphericalcoords():
    """ test of spherical coordinates to a cartesian point 
        done in a nested loop to make a sphere
       
        PYCORE_OBJ_IN not needed 
    """

    obj = object3d()

    for theta in range(-180,180,20):
        print('## theta ', theta )
        for phi in range(-180,180,20):        
            sp = spherical(1.5, mu.dtr(theta), mu.dtr(phi) ) 
            pt=  sp.to_cartesian() 
            obj.prim_cube(pos=pt, size=.1, linecolor=(255,0,0), rot=(0,0,0), pivot='world')

    obj.save(PYCORE_OBJ_OUT) 


##------------------
def procedural_1():
    obj = object3d()


    obj.prim_circle(axis='y', pos=(0,0,0), spokes=10, dia=1.2) 
    obj.triangulate(force=True)

    for i in range(1,len(obj.polygons) ):   
        obj.extrude_face(i, 5)

    for theta in range(-180,180,20):
        print('## theta ', theta )
        for phi in range(-180,180,20):        
            sp = spherical(1.5, mu.dtr(theta), mu.dtr(phi) ) 
            pt = sp.to_cartesian() 
            obj.prim_cube(pos=pt, size=.1, linecolor=(255,0,0), rot=(0,0,0), pivot='world')

    obj.scale_pts( (1,.5,1) )
    obj.rotate_pts( (0, 45, 0) )
    
    #there is NO MOVE???
    #obj.move_pts( (0,-3,0) )

    for theta in range(-180,180,20):
        print('## theta ', theta )
        for phi in range(-180,180,20):        
            sp = spherical(1.5, mu.dtr(theta), mu.dtr(phi) ) 
            pt = sp.to_cartesian() 
            obj.prim_cube(pos=pt, size=.1, linecolor=(255,0,0), rot=(0,0,0), pivot='world')

    obj.save(PYCORE_OBJ_OUT) 



##------------------
def primitive(primtype):
    # PYCORE_OBJ_IN not needed 

    obj = object3d() 

    position = (0,0,0)
    rotation = (0,0,0)
    size = 1 
    axis = 'z'

    do_flush = True

    
    # obj.prim_line( axis=axis, pos=position, rot=rotation, size=size)
    # obj.save(PYCORE_OBJ_OUT)
    # if do_flush:
    #     obj._flush()

    if primtype == 'triangle':
        obj.prim_triangle( axis=axis, pos=position, rot=rotation, size=size)
        obj.save(PYCORE_OBJ_OUT)
        if do_flush:
            obj._flush()

    if primtype == 'quad':
        obj.prim_quad( axis=axis, pos=position, rot=rotation, size=size)
        obj.save(PYCORE_OBJ_OUT)
        if do_flush:
            obj._flush()

    if primtype == 'circle':
        obj.prim_circle( axis=axis, pos=position, dia=size) #rot=rotation
        obj.save(PYCORE_OBJ_OUT)
        if do_flush:
            obj._flush()

    if primtype == 'sphere':
        obj.prim_sphere(  pos=position, rot=rotation, size=size)
        obj.save(PYCORE_OBJ_OUT)
        if do_flush:
            obj._flush()

    if primtype == 'locator':
        obj.prim_locator(  pos=position, rot=rotation, size=size)
        obj.save(PYCORE_OBJ_OUT)
        if do_flush:
            obj._flush()

    if primtype == 'locatorxyz':
        obj.prim_locator_xyz(  pos=position, rot=rotation, size=size)
        obj.save(PYCORE_OBJ_OUT)
        if do_flush:
            obj._flush()

    if primtype == 'cone':
        obj.prim_cone( axis=axis, pos=position, dia=size) #rot=rotation
        obj.save(PYCORE_OBJ_OUT)
        if do_flush:
            obj._flush()


##------------------

def pt_transform():

    """ example of translate, rotate, scale of raw points 
        translate tools work with "ptgroups", or raw points
    """

    obj = object3d()
    obj.load(PYCORE_OBJ_IN)

    obj.points = obj.scale_pts( (1,1,1)  , pts=obj.points )   
    
    obj.points = obj.rotate_pts((0,0,45), pts=obj.points ) 
    
    #obj.points = obj.xform_pts( (0,2,0),  pts=obj.points ) 

    obj.save(PYCORE_OBJ_OUT)


##------------------
def face_extrude():
    """ brute force test of face extrude 
        extrudes all faces in a polygon object 
        also will display a cheapo test "progress bar"
        because it can be slow 
    """

    obj = object3d()
    obj.load(PYCORE_OBJ_IN)

    #end = len(obj.polygons)
    end = 20

    for i in range(1, 3 ):   
        obj.extrude_face(i, .5)
    
    #for i in range(1,100 ):   
    #    obj.extrude_face(i, .1)

    obj.save(PYCORE_OBJ_OUT)
    

##------------------
def m33_to_vectors(m33, transpose=False):
    outvecs = []

    if transpose:
        v1 = vec3(m33[0], m33[3], m33[6])
        v2 = vec3(m33[1], m33[4], m33[7])
        v3 = vec3(m33[2], m33[5], m33[8])  
    else:
        v1 = vec3(m33[0], m33[1], m33[2])
        v2 = vec3(m33[3], m33[4], m33[5])
        v3 = vec3(m33[6], m33[7], m33[8])

    outvecs = [v1,v2,v3]
    return outvecs 

def m44_to_vectors(m44, transpose=False):

    outvecs = []

    if transpose:
        v1 = vec3(m44[0], m44[4], m44[8])
        v2 = vec3(m44[1], m44[5], m44[9])
        v3 = vec3(m44[2], m44[6], m44[10])  
    else:
        v1 = vec3(m44[0], m44[1], m44[2])
        v2 = vec3(m44[4], m44[5], m44[6])
        v3 = vec3(m44[8], m44[9], m44[10])

    outvecs = [v1,v2,v3]
    return outvecs 


##------------------
def visualize_matrix_rotation():
    """ 
        turn a 3X3 matrix into 3 vectors and render them as lines 

        Use numpy to rotate a matrix around a vector 
        rotate a matrix around a vector with a scalar angle (degrees) 
        visualize the matrix as a 3D object 
    """
    if NUMPY_IS_LOADED:
        obj  = object3d()
        m33  = matrix33()
        
        def render_multi_matrices(lom):   
            outvecs = []
            for m in lom:      
                tmp = m33_to_vectors( m )
                outvecs.extend( tmp )
            return outvecs
 
        axis = vec3(3,3,0)
       
        m1 = m33.from_vec3( axis , 45 )
              
        rendervecs = render_multi_matrices( [m1,m33] )

        # viz the axis
        rendervecs.append(axis)
        obj.vectorlist_to_obj( rendervecs )
        obj.save(PYCORE_OBJ_OUT)  
    else:
        print("ERROR NUMPY IS DISABLED ")


##------------------
def visualize_perspective_matrix():
    obj = object3d()
    obj.load(PYCORE_OBJ_IN)    
    persp_m44 = matrix44() 

    
    #create a perspective matrix 
                                       # fov, aspect, znear, zfar):
    persp_m44 = persp_m44.buildPerspProjMat( 15, 1, 1, 2)
    obj.points = obj.apply_matrix_pts(obj.points, m44=persp_m44)
    
    #visualize the matrix in the object 
    rendervecs = m44_to_vectors( persp_m44 )
    obj.vectorlist_to_obj( rendervecs )    

    obj.save(PYCORE_OBJ_OUT) 


##------------------

#def visualize_matrices_multiplied():

##------------------
#def visualize_matrices_multiplied():


##------------------

#def lathe_profile(axis, pts, numturns ):


#def loft_two_profiles(pts1, pts2 ):

##------------------


def loadgcode():
    gkod = gcode()
    
    #gkod.load_gcode('ngc/3D_Chips.ngc')
    
    gkod.load_gcode('3d_obj/arcspiral.ngc')

    gkod.show_data()

    gkod.save_3d_object(PYCORE_OBJ_OUT) 




##------------------
def scratch_obj1():

    """ build a new polygon object in memory from points 
        
        This could be used for a million things like procedural modeling,
        custom file importers, exporters., etc  
        
        PYCORE_OBJ_IN not needed since we are generating the model           
    """ 

    geom  = [[],[]]
    geom2 = [[],[]]

    obj = object3d()
    #obj.load(PYCORE_OBJ_IN)

    #add new geom and auto increment the ids

    #3 sided polygons with color
    pts = [(1,1,1, 1,0,0),(0,1,1, 1,0,0),(-1,-1,1, 1,0,1),(2,-2,1, 0,1,0)]
    polys = [(1,2,3),  (4,3,1) ]

    geom = obj.insert_polygons(polys, pts, geom=geom) 
 
    #4 sided polygons
    pts = [(4,-4.3,-3),(5.5,-2.5,-2.1),(-2,2,-4),(4,-4.2,1)]
    polys = [(1,2,3,4) ]

    geom2 = obj.insert_polygons(polys, pts, geom=geom2) 

    # use insert to add geom to object 
    obj.insert(geom) 
    obj.insert(geom2) 

    obj.save( PYCORE_OBJ_OUT )

##------------------
def scratch_obj2():
    """ another demo of a model built from coordinates
        when you use asnew_shell , you can append to other models 
    """
    obj = object3d()
 
    #obj.load(PYCORE_OBJ_IN)

    #add new geom and auto increment the ids
    polys = [(1,2,3,4) ]
    pts = [(10,1,1),(0,11,1),(-1,-5,1),(2,-2,1)]
    obj.insert_polygons(polys, pts) 

    # #add new geom and auto increment the ids
    # pts = [(0,-3,-1),(2,-2,1),(3,-1,1)]
    # obj.insert_polygons([], pts)
    # #add polys without new points into same "shell"
    # obj.insert_polygons( [(1,2,3,4,5,6,7),(1,7,2)], None, ans=False)
    # #add new polygon in a new "shell" 
    # obj.insert_polygons( [(1,2,3,4)], [(3,3,3), (3,-4,5), (-4,-2.5,3.1), (6.2,-2.7,8)], ans=True)

    obj.save(PYCORE_OBJ_OUT)


##***********************************************************##
##***********************************************************##
##***********************************************************##
##***********************************************************##

def pyrender_ogl():
    ropr = simple_render()
    persp_m44 = matrix44()    
    persp_m44.load_file( M44_DISK_FILE )

    obj = object3d()
    obj.load(PYCORE_OBJ_IN)
    #obj.triangulate(force=True)  

    use_perpective = True 
    if use_perpective:
        #                          # fov, aspect, znear, zfar):    
        m44_2 = persp_m44.buildPerspProjMat( 100, 1, .1, 5)
        persp_m44 =  persp_m44 * m44_2
    
    obj.points = obj.apply_matrix_pts(obj.points, m44=persp_m44)

    # ropr.render_obj((100,0,255), 1, 1, 1, 1, 1000/abs(persp_m44.m[14]), object3d=obj)

    img_op = pixel_op()   
    img_op.load('textures/generated2.bmp') 
    #img_op.load('render_thing/images/uvmap.bmp') 

    #img_op.save("foobar.bmp")

    lightpos = (0, 2 ,-1)
    ropr = simple_render()

    ##----------
    #ropr.COLOR_MODE = 'flat'
    #ropr.COLOR_MODE = 'lighted'
    ropr.COLOR_MODE = 'lightedshaded'


    ropr.SHOW_VTXS        = False
    ropr.SHOW_FACE_CENTER = False
    ropr.SHOW_EDGES       = False    

    ##----------

    ropr.scanline(obj, 1000/abs(persp_m44.m[14]), lightpos=lightpos, texmap=img_op ) 
    ropr.save_image( PYCORE_BMP_OUT )




def bezier3d():
    obj = object3d()

    start = (1 ,  0, 0)
    ctrl1 = (.5,  0, 0)
    ctrl2 = ( 0, .5, 0)
    end   = (0 ,  1, 0)
    kurve = [start, ctrl1, ctrl2, end]

    start  = (0   ,  1   , 0)
    ctrl1  = (1.5 ,  1.0   , 0)
    ctrl2  = (1   ,  1.5 , 0)
    end    = (2   ,  2   , 0)
    kurve2 =[start, ctrl1, ctrl2, end]
    
    curves = [kurve, kurve2] 
    
    obj.draw_splines( 10, curves, drawctrls=True, drawhulls=True)

    obj.save(PYCORE_OBJ_OUT)




def lathe():
    """ needs to have the same num U and V to work """
    obj = object3d()

    # # simplest possible lathe example (must be square , num pts == num revolutions)
    # pts = [(.1,.1,0),(1,1,0),(2,2,0),(3,3,0)]
    # obj.lathe(pts, 4)


    # # using bezier curve function 
    # num = 23
    # start = (1 ,  0, 0)
    # ctrl1 = (.5,  0, 0)
    # ctrl2 = ( 0, .5, 0)
    # end   = (0 ,  1, 0)
    # curve = obj.cubic_bezier(num, start, ctrl1, ctrl2, end)
    # obj.lathe(curve, num)

    #------ 

    num = 10

    start  = (2 ,  0, 0)
    ctrl1  = (.5,  .5, 0)
    ctrl2  = (3, .5, 0)
    end    = (1 ,  1, 0)
    curve1 = [  start, ctrl1, ctrl2, end ]

    curvepts = obj.cubic_bezier(num, start, ctrl1, ctrl2, end )

    start  = (1   ,  1   , 0)
    ctrl1  = (2.5 ,  1.0   , 0)
    ctrl2  = (1   ,  1.5 , 0)
    end    = (2   ,  2   , 0)
    curve2 = [ start, ctrl1, ctrl2, end ]
  
    tmp      = obj.cubic_bezier(num, start, ctrl1, ctrl2, end )
   
    curvepts.extend(tmp)

    #obj.draw_splines( num, [curve1, curve2], drawctrls=True, drawhulls=True) 

    #obj.lathe2(curvepts, num*2)
    
    obj.lathe(curvepts, num*2)

    #for i in range(120):
    #    obj.extrude_face(i, -.4)

    #scal command kills COLOR - DOH !
    #obj.scale_pts( (1,2,1) )

    #------ 
    
    # num = 20
    # #cross_section = obj.calc_circle(pos=(1,0,0), rot=(0,180,0), dia=.5, axis='z', periodic=True, spokes=num-1)
    # cross_section = obj.calc_circle(pos=(1,0,0), dia=.5, axis='z', periodic=True, spokes=num-1)
    # cross_section.reverse()
    # obj.lathe(cross_section, num)


    obj.save(PYCORE_OBJ_OUT)
 






def test_pointgen():
    #OBJ.locate_pt_along3d( x1, y1, z1, x2, y2, z2, num):
    obj = object3d()
    grid = obj.test_data_grid(10,10,2)
    obj.print_grid( grid )

    #print( obj.get_grid_column(grid, 3) )
    #print( obj.get_grid_row(grid, 3) )
    
    #print( obj.indexer( ids=ids, span=None, unique=True, nth=None) )

    obj.show()

    ids = [1,2,3,4,5,6,7,8,9,10,11,12,13,14]
    # print( obj.indexer( ids, [1,10], False , 10) )
    print(obj.chunker(6, ids) )


#####################################################

"""
obj = object3d()
#obj.load('gnelscript/objects/cube.obj')
obj.load(PYCORE_OBJ_IN)

#obj.rotate_pts((45,45,45))
obj.scale_pts((.5,.5,.5))

obj.save(PYCORE_OBJ_OUT)
#obj.save("%s/%s"%(PYCORE_GEOMPATH, "cube.obj"))
""" 




## """ lookup and return the polygon indices and points for a single polygon 
##     reindex - if True  - renumber the new polygon indices startring at 1, 
##               if False - retain the oringial numbering 
##     geom - act on a geom obj passed in, or on self
## """


def extract():
    obj = object3d()
    obj.load(PYCORE_OBJ_IN)

    obj2 = object3d()

    ## you can use it just with a list of IDs 
    # geom = obj.get_face_geom( [22,23,1,2,3,4] )

    ## or you can use indexer for more power 
    pids = obj.indexer( ids=[30],span=[8,10])
    geom = obj.get_face_geom( pids, reindex=True )
    obj2.insert(geom)

    pids = obj.indexer( ids=[5,10])
    geom = obj.get_face_geom( pids, reindex=True )
    obj2.insert(geom)

    obj2.save(PYCORE_OBJ_OUT)
    


def kicad_test():
    """ experiment to parse a kicad pcb file and export it to gcode """


    tokens = kicadproj.split(os.sep) 
    projname = tokens[len(tokens)-1]
    pcbname = projname+'.kicad_pcb'

    kiparser = pcbfile()

    #kiparser.load_gcode('gcode/ngc/3D_Chips.ngc')
    kiparser.load_kicadpcb(kicadproj+os.sep+pcbname)

    # pts = kiparser.calc_circle(pos=(-2,5,0), dia=30, spokes=5)
    # kiparser.filled_polys.append(pts) 
    # pts = kiparser.calc_circle(pos=(1,1,0), dia=60, spokes=11)
    # kiparser.gr_polys.append(pts) 

    #kiparser.bufferinfo()
    #kiparser.showbuffers()
    #kiparser.show_geom()

    kiparser.export_ngc('cineballz.ngc')
    kiparser.save_3d_obj(PYCORE_OBJ_OUT)



def linuxcnctest():
    """ experiment to parse a gcode file """
    print("####### foo")
    gc_poly = gcode()
 
    cwd = os.getcwd()
    path = cwd+'/gnelscript/gcode/ngc/arcspiral.ngc'
    gc_poly.load_gcode(path)

    #gc_poly.load_gcode('gcode/ngc/arcspiral.ngc')
    #gc_poly.load_gcode('gcode/ngc/calibrate_5X1.ngc')
    #gc_poly.load_gcode('gcode/ngc/3D_Chips.ngc')

    gc_poly.show_data()



def linuxcnctest2():
    """ standalone gcode experiment """
    gcode = generate_gcode()
    gcode.milling_test_linear()
    
    #print( gcode.outfile )
    gcode.savengc("cineballs.ngc")



def gear_test():
    """ experiment to generate some gears """
    gearz = gear_generator()
    
    # build( shaftdia, dia, teeth_height, numteeth ):

    pts =  gearz.build(1, 2,.5, 12) 
    gearz.linegeom_fr_points(pts, color=(100,0,100), periodic=False )
    gearz.save(PYCORE_OBJ_OUT)

################################################

## parse commands coming in and run them
def runcommand():
    #SELECTION EXAMPLES 
    #slice_extract_and_makenew(PYCORE_OBJ_OUT)
    #extract_by_copy_hack(PYCORE_OBJ_OUT)
    #test_subsel_point_transform(PYCORE_OBJ_OUT)
    ## modify_a_subselect(PYCORE_OBJ_OUT) ##BROKEN 
    #select_polygons_spatially(  PYCORE_OBJ_OUT, (.25,.21,1), 5 )

    #VECTOR EXAMPLES 
    #rotate_around_vec(PYCORE_OBJ_OUT)
    #m33  = matrix33();render_m33_as_vec3s(m33, PYCORE_OBJ_OUT, transpose=False, vlist=None)
    #numpy_m33_fromvec(PYCORE_OBJ_OUT)
    #make_right_triangle(PYCORE_OBJ_OUT, 71.61)
    # unit_circle_viewer(PYCORE_OBJ_OUT) ##BROKEN 
    #visualize_cross_product( PYCORE_OBJ_OUT )
    #offset_between_2vecs( PYCORE_OBJ_OUT )
    #build_orthogonal_vector(PYCORE_OBJ_OUT)

    #WIP EXAMPLES 
    #place_obect_on_vector(PYCORE_OBJ_OUT)
    
    kicad_test()
    #scratch_obj2()
    #linuxcnctest()
    #linuxcnctest2()
    #gear_test()



    #extract()

    #face_extrude()

    #loadgcode()
    #loadkicad()
    #test_pointgen()

    #scratch_obj1()

    #bezier3d()
    #lathe()
    
    #visualize_matrix_rotation()
    #visualize_perspective_matrix()
    #gen_normals()
    
    #circle_cube_pts()
    #primitive('triangle')
    
    #procedural_1()
    

    #primitive('sphere')

    #pt_transform()
    #procedural_1()
    #modify_partial()
    #triangulate()
    #object_rotate()
    #copy_sop()
    pass







#####################################################



if __name__=="__main__":

    print("############ ", sys.argv )

    if sys.argv[2] == 'tcptest':
        tcpviz = tcpviz()
        #tcpviz.send_str('abcdefg')
        #tcpviz.close()

        #tcpviz.vz_locator() 
        tcpviz.vp_set_grid('1')

    if sys.argv[2] == 'tcpclose':
        tcpviz = tcpviz()
        tcpviz.close()

        #if sys.argv[1] =='close':
        #    tcp_send(0x00)
        #else:
        #    tcp_send(bytes(sys.argv[1], 'utf-8'))

    if sys.argv[2] == 'runcommand':
        runcommand()

    if sys.argv[2] == 'scanline':
        pyrender_ogl()    

    #if sys.argv[2] == 'normals':
    #    gen_normals()











