![alt text](https://github.com/keithlegg/gnelscript/blob/master/images/example/monkey_tex.png) 

# gnelscript


Procedural 3D Model, Render, and basic image composite in pure python. 

I am writing this to learn 3D math. It is for fun and to be used
as a platform to teach 3D math and graphics. 

Inspired by years of working with Maya 3D and Houdini. 
I love the "ops" idea from houdini, everything is an operator and
you chain them together.



    Organized by the following modules:
       -  math_ops         - math, vectors, matrices, the core logic of it all  
       -  point_ops        - points, polygons, and objects, the "brains" of geometry processing
       -  point_ops_2d     - 2d points, polygons
       -  obj2d            - data structure for 2D models.   Inherits all the other stuff  
       -  obj25d           - data structure for 2.5D models. Inherits all the other stuff  
       -  obj3d            - data structure for 3D models.   Inherits all the other stuff     
       -  raster_ops       - image manipulation with PIL, framebuffer  
       -  render           - brain dead simple 3D rendering, examples of using pygfx modules to build a 3D render  
       -  milling_ops      - process G-code files for CNC/CAM related stuff. Experimental. 
       -  kicad_ops        - parse the kicad pcb format for CNC/CAM related stuff. Experimental.
       -  gis_vector_ops   - GIS data import and export. Experimental.

    Example files are split up by type:
       -  objects_3d          - examples of 3D rendering and 2D vector processing 
       -  examples_render     - examples of 3D rendering 
       -  selection           - examples of extracting portions of polygon geometry  
       -  vector              - examples of 3D vector processing
    


    Stuff you can ignore for now :
       -  examples.dag         - directed acyclic graph tools - pygfx.dag is 14 years old and needs cleanup 
       -  examples.milling     - examples of G-code related tools. CNC and CAM in the future?
       -  examples.wip         - work in progess examples 
       -  examples.objects_2d  - nothing here worth looking at 
       -  pygfx.objects_2d        - 2d is a work in progress - it will behave like all the 3D code but in 2D 
       -  pygfx.unit_tests        - not done yet 
       -  pygfx.raytrace          - someones code I imported and am considering intergrating 
       -  pygfx.dag_ops           - directed acyclic graph, first stab at a scene graph (not working yet)
       -  pygfx.render_2d         - nothing here yet, future home of 2D rendering 



Raster tools structure:

    |-raster_op          - wrapper around PIL to work on whole images
        |                  obj.fb = framebuffer = PIL Image.Image object that is being operated on 
        |
        |---pixel_op     - operations at pixel level - this is the .fb framebuffer used almost everywhere



3D objects structure:

    can be used along side numpy or standalone
    inside the __init__ file there is a global NUMPY_IS_LOADED to enable/disable the components that use numpy 

    The object inheritance of a gnelscript 3D object is as follows:
    
    |-math_util - standalone math tools , the first prototype for math ops and all the others  

    |-math_ops           - just what it says. Vectors, Matricies, all the boring math stuff that everything is built on
       |-point_ops       - functions to work with coordinates directly
          |-polygon_ops     - deals with polygon indices, point indices 
              |             - geometry is stored as vertex arrays and face indexes with a flat data structure  
              |-object      - deals with groups of polygons, normals, polygon operators loading, saving, ect 
           
 

Designed to be a self contained python module that you import into other places.
It has some numpy functions that can be disabled with a global.



To enable numpy go to __init__ and set 

NUMPY_IS_LOADED=True


To use examples - just imort them like this 

    from gnelscript.examples.2d import *
    from gnelscript.examples.3d import *
    from gnelscript.examples.raster import *
    from gnelscript.examples.vector import *
    from gnelscript.examples.selection import *












