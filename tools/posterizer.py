

from gnelscript.tools.imagecam import * 
 




##----------------------------------------------------
##   /usr/local/opt/python@3.10/bin/python3.10 ./imagecam.py  


#firstpass_bw(10, 1.5, 1.5, 1, 250, "images/in/art.jpg", "images/out", "output")


## (iteration , blur , contrast, bright, scaling(divs) , in, out )
#firstpass(10, 0, 1, 1, 250, "images/in/oil.png", "images/out", "output")

#firstpass(3, 1.5, 1.2, 1, 250, "images/in/foo.jpg", "images/out", "output")



##----------------------------------------------------
##   /usr/local/opt/python@3.10/bin/python3.10 ./imagecam.py  
## secondpass("images/in/foo.jpg", "images/out" , 4, False)


##----------------------------------------------------
#set the RGB values from last tool and run this 
#thirdpass( "images/out/commonbands.png", "images/out", "dxf" )
#thirdpass( "images/out/commonbands.png",  "images/out" , "geojson")

#thirdpass( "images/out/foo.bmp",  "images/out" , "geojson", po_invert=False)

#thirdpass( "images/out/commonbands.png",  "images/out" , "geojson")




##----------------------------------------------------


def vector_magic():
    kiparser = generic_ngc()

    ##-- 

    #indexer(ids=None, span=None, unique=True, nth=None)
    ids = kiparser.indexer(span=[1,2])

    ##-- 

    kiparser.load_geojson('images/out/2.json', 0, getfids=None, getids=None)
    #kiparser.load_geojson('images/out/3.json', 0, getfids=None, getids=None)

    ##--
    bbox = kiparser.calc_bbox_pt(2, (5,5))
    pts = kiparser.cvt_2d_to_3d(kiparser.extents_fr_bbox(bbox, periodic=True))
    
    #debug - need to solve the clean_pts_str debacle?
    kiparser.gr_polys.append(pts)

    ##--
    bbox = kiparser.calc_bbox_pt(1.75, (-3,3))
    pts = kiparser.cvt_2d_to_3d(kiparser.extents_fr_bbox(bbox, periodic=True))
    
    #debug - need to solve the clean_pts_str debacle?    
    kiparser.gr_polys.append(pts)

    #print(pts)
    #kiparser.grply_inspect()
    #kiparser.cvt_grpoly_obj3d()
    #kiparser.save("3d_obj/foo.obj")

    #scale if you want to 
    gs = kiparser.global_scale
    xformed = []
    for ply in kiparser.gr_polys:
        xformed.append(kiparser.trs_points( ply, translate=(0,0,0), rotate=(0,0,0), scale=(gs,gs,gs) ))
    kiparser.gr_polys = xformed

    kiparser.calculate_paths()
    #kiparser.save_line_obj('3d_obj/foo.obj')
    kiparser.export_ngc("foo.ngc")


 



vector_magic()






# TODO 

""" 

 establish end orientation and scale of vectors related to CNC 

 build tools to scale, rotate, flip, etc if not already 

 re learn what tools I have already  





"""
