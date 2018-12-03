# pyrender2

"swiss army power tool" for 3D geometric data. 

Procedural modeling for polygons, 
Very basic 3D rendering in pure python.
Render and composite images with PIL. 


Lot of tools to process polygons like:

   - OBJ file load and save 
   - Sub select parts of a model 
   - Triangulation 
   - Polygon N sided face extrude 
   - visualize vectors and line geometry 





![alt text](https://github.com/keithlegg/pyrender2/blob/master/images/example/monkey.flat.png) 

I am writing this to learn 3D math. It is for fun and to be used
as a platform to teach 3D math and graphics. 

Inspired by years of working with Maya 3D and Houdini. 
I love the "ops" idea from houdini, everything is an operator and
you chain them together.







design goals:
   - simple. Only do what needs to be done and nothing more
   - flat code structure as possible. All components should be easy to find 
   - fun. Its a playground to lean vectors, matrices, math and more. 
   - only requires PIL and numpy 
   - numpy is loosely intergrated and can be disabled






   Organized by the following modules:
   -  math_ops    - math, vectors, matrices, the core logic of it all  
   -  point_ops   - points, polygons, and objects, the "brains" of geometry processing
   -  obj3d       - data structure for 3D models. Inherits all the other stuff  
   -  raster_ops  - image manipulation with PIL, framebuffer  
   -  render      - brain dead simple 3D rendering, examples of using pygfx modules to build a 3D render  

   Example files are split up by type:
   -  examples            - examples of geometry processing
   -  examples_2d         - examples of 2D rendering and 2D vector processing 
   -  examples_render     - examples of 3D rendering 
   -  examples_selection  - examples of extracting portions of polygon geometry  
   -  examples_vector     - examples of 3D vector processing

   Stuff you can ignore for now :
   -  examples_wip - work in progess examples 
   -  unit_tests   - not done yet 
   -  raytrace     - someones code I imported and am considering intergrating 
   -  dag_ops      - directed acyclic graph, first stab at a scene graph (not working yet)
   -  grid_ops     - nothing here yet, future home of mapping and grids 
   -  examples_2d  - nothing here worth looking at  
   -  render_2d    - nothing here yet, future home of 2D rendering 











