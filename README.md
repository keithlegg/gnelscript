# pyrender2

Pure python 3D tools. 

A complete-ish rebuild of my oringial repository, pyrender. 
This code is a sandbox of many different tools I have worked with over the years.

The big changes from oringial pyrender:

  - cleaner, more organized code 
  - added 3X3 matrices 
  - matrix and vectors are thier own types of object now
  - numpy integration, so it can interface to numpy arrays seamlessly 
  - point operators 
  - polygon operators 
  - improved object operators and primitive objects
  - object3D type 
  - "matrix" render - pass a 3X3, or 4X4 matrix to renderer directly 


Inspired by projects I have on over the last 10 years:

   virtual reality, game design, 3D animation, GIS and mapping
   robotics, embdedded hardware, computer vision 


Contains libraries for:

  -   2D vector 
  -   3D vector 
  -   3X3 matrix 
  -   4X4 matrix 



Organized by the following modules:

   -  dag_ops    - directed acyclic graph, first stab at a scene graph (not working yet)
   -  grid_ops   - nothing here yet
   -  math_ops   - matrix, vectors, and math stuff 
   -  raster_ops - image manipulation with PIL 





