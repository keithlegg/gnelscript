
![alt text](https://github.com/keithlegg/gnelscript/blob/master/images/example/monkey_tex.png) 

# gnelscript


Render, 3D Model and composite images in pure python. 


I am writing this to learn 3D math. It is for fun and to be used
as a platform to teach 3D math and graphics. 

Inspired by years of working with Maya 3D and Houdini. 
I love the "ops" idea from houdini, everything is an operator and
you chain them together.



   Organized by the following modules:
   -  math_ops    - math, vectors, matrices, the core logic of it all  
   -  point_ops   - points, polygons, and objects, the "brains" of geometry processing
   -  obj3d       - data structure for 3D models. Inherits all the other stuff  
   -  raster_ops  - image manipulation with PIL, framebuffer  
   -  render      - brain dead simple 3D rendering, examples of using pygfx modules to build a 3D render  
   -  milling_ops - process G-code files for CNC/CAM related stuff. Experimental. 
   -  kicad_ops   - parse the kicad pcb format for CNC/CAM related stuff. Experimental.

   Example files are split up by type:
   -  examples            - examples of geometry processing
   -  examples_2d         - examples of 2D rendering and 2D vector processing 
   -  examples_render     - examples of 3D rendering 
   -  examples_selection  - examples of extracting portions of polygon geometry  
   -  examples_vector     - examples of 3D vector processing


   Stuff you can ignore for now :
   -  examples_milling  - examples of G-code related tools. CNC and CAM in the future?
   -  examples_wip      - work in progess examples 
   -  unit_tests        - not done yet 
   -  raytrace          - someones code I imported and am considering intergrating 
   -  dag_ops           - directed acyclic graph, first stab at a scene graph (not working yet)
   -  grid_ops          - nothing here yet, future home of mapping and grids 
   -  examples_2d       - nothing here worth looking at  
   -  render_2d         - nothing here yet, future home of 2D rendering 




Overall structure

    can be used along side numpy or standalone
    inside the __init__ file there is a global NUMPY_IS_LOADED to enable/disable the components that use numpy 

    The object inheritance of a gnelscript 3D object is as follows:
    
    |-math_ops           - just what it says. Vectors, Matricies, all the boring math stuff that everything is built on
       |-point_ops       - functions to work with coordinates directly
       |                 - geometry is stored as vertex arrays and face indexes with a flat structure to interface to outside data easily 
       |-polygon_ops     - deals with polygon indices, point indices 
           |-object3d    - deals with groups of polygons, loading, saving, ect 
           
 

    The elegance of the design is you can use operators directly, or use the objects without thinking about them.
    Every object has all operators already inside it, so things can be pluged together like legos. 
    This was very much inspired by the genius of Side Effects Houdini. 










