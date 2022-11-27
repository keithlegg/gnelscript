

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


    #indexer(ids=None, span=None, unique=True, nth=None):
    ids = kiparser.indexer(span=[1,2])
    kiparser.load_geojson('images/out/0.json', 0, getfids=None, getids=None)


    #kiparser.load_geojson('images/out/0.json', 0, getids=ids)

    #kiparser.load_geojson('images/out/0.json', 0)

    #kiparser.load_geojson('images/out/2.json', 0)
    #kiparser.load_geojson('images/out/3.json', 0)

    #pts = kiparser.cvt_2d_to_3d(kiparser.extents_fr_bbox([2,2,2,2], 3))
    #kiparser.gr_polys.append(pts) 

    #kiparser.gr_polys.append( [(-1, -1, 0), (-1, 5, 0), (5, 5, 0), (4, -4, 4), 
    #                           (5, 5, 5),   (6, 6, 6),  (7, 7, 7), (8, 8, 0)] ) 





    # kiparser.rh = .01     # retract height 
    # kiparser.ch = .1       # cut height 
    # kiparser.hp = (0,0,0) # home position 

    # calc_circle(self, pos=(0,0,0), rot=(0,0,0), dia=1, axis='z', periodic=True, spokes=23):

    # pts = kiparser.calc_circle(pos=(0,0,0), rot=(0,0,0), dia=1, axis='z', periodic=True, spokes=23)
    # kiparser.gr_polys.append(pts)


    # pts = kiparser.calc_circle(pos=(2,0,0), rot=(0,0,0), dia=.5, axis='z', periodic=True, spokes=8)
    # kiparser.grpts(pts)

    #pts = kiparser.locate_pt_along3d(0,0,0,5,5,5,8)
    #kiparser.gr_polys.extend([pts])



    print("gr_poly buffer has %s polygons "%len(kiparser.gr_polys) )
    print(len(kiparser.gr_polys))



    #if you want to export NGC/OBJ from geojson 

    #kiparser.calulate_paths( do_retracts=False)
    #kiparser.save_3d_obj('3d_obj/foo.obj', export_retracts=False)

    kiparser.calulate_paths()
    kiparser.save_3d_obj('3d_obj/foo.obj')

    kiparser.export_ngc("foo.ngc")

    #if you want to export OBJ without retracts from geojson
    #kiparser.calulate_paths(do_retract=False) 
    #kiparser.save_3d_obj('3d_obj/foo.obj')


 


vector_magic()






# TODO 

""" 

 establish end orientation and scale of vectors related to CNC 

 build tools to scale, rotate, flip, etc if not already 

 re learn what tools I have already  





"""
