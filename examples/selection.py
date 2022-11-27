


#from pygfx.render import *

from gnelscript.pygfx.math_ops import  *
from gnelscript.pygfx.point_ops import *
from gnelscript.pygfx.obj3d import  *


#####################################################
def extract_by_copy_hack(outfile):
    """ *slightly* higher level than raw geom 
        use the face lookup util to get the geom by face index 
    """

    obj = object3d()
    obj.load('3d_obj/sphere.obj')

    # duplicate all points with out thinking about it 
    all_pts = obj.points

    # extract two chunks of poly ids 
    polygr1 =  obj.get_pt_ids([0,1,2])  
    polygr2 =  obj.get_pt_ids([3,5])

    # make a new object and dump data into it
    obj2 = object3d() 
    obj2.points = all_pts  #move all points over - DEBUG add a clean func to remove unused

    for ply in polygr1:
        obj2.polygons.extend(ply)
    for ply in polygr2:
        obj2.polygons.extend(ply)        


    obj2.save(outfile)




#####################################################
def slice_extract_and_makenew(outfile):
    """ load two models, extract parts of them using the subselect tool 
        subselect grabs polys and points at the same time, with option to reindex
         
        reindex effectively makes a new object 

        weld them into a new model 
        fekkin awesome mate!  
    
        ##!! DEBUG DEPRECATION 
        PYCORE EXTRACT USES NEW get_face_geom - probably better 
           
           sub_select_geom should be rewritten to use that 
           maybe support mulitple groups like 
           [ [1-5],[9-12], etc ]

        ## or you can use indexer for more power 

        obj2.insert(geom)
        
        pids = obj.indexer( ids=[5,10])
        geom = obj.get_face_geom( pids, reindex=True )
        obj2.insert(geom)

    """
    
    obj = object3d() 
    obj.load('3d_obj/sphere2.obj')

    obj2 = object3d() 

    obj3 = object3d() 
    obj3.load('3d_obj/monkey.obj')

    ## weld two models together 
    #geom = obj.sub_select_geom( slice=(1,50), ids=[53,55,73], reindex=True)
    #geom2 = obj3.sub_select_geom( slice=(3,100) , ids=[101,105,148], reindex=True)    
    #obj2.insert_polygons(geom[0], geom[1]  ) 
    #obj2.insert_polygons(geom2[0], geom2[1]  )

    pids = obj.indexer( ids=[53,55,73],span=[1,50])
    geom = obj.get_face_geom( pids, reindex=True )


    pids = obj3.indexer( ids=[101,105,148],span=[3,100])
    geom2 = obj3.get_face_geom( pids, reindex=True )
    
    obj2.insert(geom)    
    obj2.insert(geom2)

    obj2.save( outfile )

#####################################################
def test_subsel_point_transform(outfile): 
    """ example of translate, rotate, scale of a point group 
        translate tools work with "ptgroups", or raw points 
    """
    obj = object3d()
    obj.load('3d_obj/monkey.obj')

    ptgrp = obj.get_pt_grp( span=(1,300) )
    xformed_pts = obj.scale_pts(2.5, ptgrp=ptgrp)   

    ptgrp = obj.get_pt_grp( span=(1,300) )
    xformed_pts = obj.rotate_pts((0,30,0),  ptgrp=ptgrp)

    ptgrp = obj.get_pt_grp( span=(1,100) )
    xformed_pts = obj.xform_pts( (0,2,0),  ptgrp=ptgrp) 

    obj.save(outfile)



#######################################################


def modify_a_subselect(outfile):
    """ UNFINSIHED ! """

    obj = object3d()
    #obj.load('3d_obj/sphere.obj')
    obj.load('3d_obj/cube.obj')

    geom = obj.sub_select_geom( span=[1,5]  , reindex=True )

    # rotate_pts( rot, pts=None, ptgrp=None):
    newpts = obj.rotate_pts(rot=(45,45,45), pts=geom[1])

    print(geom)
    geom2 = obj.sub_select_geom( span=[5,6]  , reindex=True )
    #newpts2 = obj.rotate_pts((-45,0,45), points=geom2[1])

    obj2 = object3d() 
    #obj2.insert_polygons(geom[0], newpts  )      
    obj2.insert_polygons(geom2[0], geom2[1], asnew_shell=False  ) 
    # obj2.insert_polygons(geom2[0], newpts2  , asnew_shell=False) 
    obj2.save(outfile)


#######################################################
def select_polygons_spatially( outfile, from_pt, dist ):
    """ 
        calculate the vector from face normal to a fixed point
        if the angle is within a threshold, select the face
        extract all those into a new model and export it  

        TO RUN: 
        
        select_polygons_spatially( (2, 1, -4) , 4.5)

    """
    obj = object3d() 
    obj.load('3d_obj/monkey.obj')

    #define the point to use 
    #from_pt = (2, 1, 4)

    # get the IDS of polygons near a point in space 
    fids = obj.select_by_location('polygons', from_pt, dist ) 

    # use sub select to grab the geometry from the face IDs 
    geom = obj.sub_select_geom( ids=fids , reindex=True)

    # dump that into a new file 
    obj2 = object3d() 

    # place a cube to mark the spot we are using as an anchor     
    obj2.prim_cube( size=.05, pos=from_pt, pivot='world' )  
    
    # save the geometry that was determined to be "near" the point 
    obj2.insert_polygons(geom[0], geom[1]  ) 
    obj2.save('3d_obj/selected.obj')
    #obj2.save(outfile)

    #pull the vectors out of a work buffer to visualize what the tool is doing     
    # #dump the vectors in the work buffer to a new object 
    obj3 = object3d() 
    obj3.vectorlist_to_obj( obj.vec_buffer )
    #obj3.save('3d_obj/dist_vectors.obj')
    obj3.save(outfile)


#select_polygons_spatially( (2, 1, -4) , 4.5)
