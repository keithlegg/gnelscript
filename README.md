# pyrender2

"swiss army power tool" for 3D geometric data. 

Procedural modeling for polygons, 
Very basic 3D rendering in pure python.
Render and composite images with PIL. 


Lot of tools to process polygons like:

OBJ file load and save 
Sub select parts of a model 
Triangulation 
Polygon N sided face extrude 
visualize vectors and line geometry 




I am writing this to learn 3D math. It is for fun and to be used
as a platform to teach 3D math and graphics. 

Inspired by years of working with Maya 3D and Houdini. 
I love the "ops" idea from houdini, everything is an operator and
you chain them together.



design goals:
   - simple. Only do what needs to be done and nothing more
   - flat code structure as possible. All components should be easy to find 
   - fun. Its a playground to lean vectors, matrices, math and more. 


only requires:
    - PIL 
    - numpy (numpy is loosely intergrated and could be pulled out)


   Organized by the following modules:
   -  unit_tests  - home for all the current examples of things you can do with the library
   -  math_ops    - math, vectors, matrices, the core logic of it all  
   -  point_ops   - points, polygons, and objects, the "brains" of geometry processing
   -  obj3d       - data structure for 3D models. Inherits all the other stuff  
   -  raster_ops  - image manipulation with PIL, framebuffer  
   -  render      - brain dead simple 3D rendering, examples of using pygfx modules to build a 3D render  

   Stuff you can ignore for now :
   -  raytrace    - someones code I imported and am considering intergrating 
   -  examples    - nothing here worth looking at 
   -  dag_ops     - directed acyclic graph, first stab at a scene graph (not working yet)
   -  grid_ops    - nothing here yet, future home of mapping and grids 
   -  examples_2d - nothing here worth looking at  
   -  render_2d   - nothing here yet, future home of 2D rendering 


The big changes from original pyrender:

  - cleaner, more organized code 
  - added 3X3 matrices 
  - matrix and vectors are their own types of object now
  - numpy integration, so it can interface to numpy arrays seamlessly 
  - better point operators 
  - better polygon operators 
  - working triangulation
  - improved object operators and primitive objects
  - object3D type 
  - "matrix" render - pass a 3X3, or 4X4 matrix to renderer directly 









