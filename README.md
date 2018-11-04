# pyrender2

Pure python Swiss Army Knife for 3D and 3D graphics. 

Procedural model, Render, composite images. 

I am writing this to learn 3D math and graphic better. It is for fun and a
tool to teach 3D math and graphics. It may even be useful for something.

Inspired by years of working with Maya 3D and Houdini. 
I love the "ops" idea from houdini, everything is an operator and
you chain them together.


only requires:
    - PIL 
    - numpy (numpy is loosely intergrated and could be pulled out)


Organized by the following modules:

   -  math_ops   - math, vectors, matrices, the core logic of it all  
   -  point_ops  - points, polygons, and objects 
   -  raster_ops - image manipulation with PIL, framebuffer  
   -  render     - brain dead simple 3D rendering, example using the other modules to render 3D 

   -  dag_ops    - directed acyclic graph, first stab at a scene graph (not working yet)
   -  grid_ops   - nothing here yet, future home of mapping and grids 
   -  render_2d  - nothing here yet, future home of 2D rendering 


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









