# pyrender2

Pure python graphics tools. 

Swiss Army Knife for 3D graphics. Render, model, composite. 

only requires PIL, and numpy 


A complete-ish rebuild of my oringial repository, pyrender. 
This code is a sandbox of many different tools I have worked with over the years.


Inspired by projects I have on over the last 10 years:

   virtual reality, game design, 3D animation, GIS and mapping,
   robotics, embdedded hardware, computer vision 


Organized by the following modules:

   -  math_ops   - math, vectors, matrices  
   -  point_ops  - points, polygons, and objects 
   -  raster_ops - image manipulation with PIL, framebuffer  
   -  render     - brain dead simple 3D rendering 

   -  dag_ops    - directed acyclic graph, first stab at a scene graph (not working yet)
   -  grid_ops   - nothing here yet, future home of mapping and grids 
   -  render_2d  - nothing here yet, future home of 2D rendering 


The big changes from original pyrender:

  - cleaner, more organized code 
  - added 3X3 matrices 
  - matrix and vectors are thier own types of object now
  - numpy integration, so it can interface to numpy arrays seamlessly 
  - point operators 
  - polygon operators 
  - improved object operators and primitive objects
  - object3D type 
  - "matrix" render - pass a 3X3, or 4X4 matrix to renderer directly 









