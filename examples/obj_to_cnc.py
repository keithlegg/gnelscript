
from gnelscript.pygfx.vector_ops import *




##------------------------------------------------##
def scan_line_tool(filldensity, inobj, outfolder, outname):
    """ 
    DEBUG - experimental    
    """
    
    #size of circles showing a "hit"
    drawdia = .01 

    loader = object3d()
    loader.load(inobj)
    loader.triangulate()

    loader.rotate( 0,0, 30.5)

    g = loader.get_face_geom(0)
    
    print(g)

    vflo = vectorflow()
    polygon = [g[1][0], g[1][1], g[1][2]]

    #get extents of polygon
    bbox     = vflo.calc_2d_bbox(pts=polygon)
    res_x = bbox[2]-bbox[0]
    res_y = bbox[3]-bbox[1]

    hscanlines = int(res_y*filldensity)

    pts = vflo.scanline_ngon( 8, hscanlines, polygon )
    pop2 = pop2d()

    for row in pts:
        if row:
            if len(row)%2==0:
                #DEBUG = USE THE DISTANCE TO CALC THE NUM BETWEEN 
                p1 = vec3(row[0])
                p2 = vec3(row[1])
                scanlen = ((p1- p2).length)
                numchunks = int(scanlen*filldensity)

                newpts = pop2.locate_pt_along(row[0][0], row[0][1], row[1][0], row[1][1], numchunks)     
                for npt in newpts:
                    vflo.prim_circle(pos=(npt[0],npt[1],0), dia=drawdia, axis='z')
            else:
                print('ODD HITS %s'%len(row))
                for npt in row:
                    vflo.prim_circle(pos=(npt[0],npt[1],0), dia=drawdia, axis='z')

       
    #make origin lines
    #vflo.prim_line(axis='x', size=2, pos=(0,0,0))
    #vflo.prim_line(axis='y', size=2, pos=(0,0,0))

    vflo.save('%s.obj'%outname) 
    vflo.cvt_obj3d_grpoly()
    vflo.export_geojson_lines(outfolder, '%s.json'%outname)
    vflo.export_ngc(1, 0, .1, 2, '%s/%s.ngc'%(outfolder, outname), do3d=False)


##------------------------------------------------##
def vector_scanlines(infile, outfolder, outfname):
    #DEBUG WIP 
    
    import math 

    vflo = vectorflow() 
    vflo.load(infile )

    # prim_square(self,  pos=(0,0,0), rot=(0,0,0), size=1):
    # prim_triangle(self, pos=(0,0), rot=(0,0), size=1):

    #procedural object 
    for a in range(10): 
        px = math.sin(vflo.mu.dtr(a*30))
        py = math.cos(vflo.mu.dtr(a*30))
        vflo.prim_triangle(pos=(px,py,0), rot=(0,0,a*30), size=.2 )

    #scratch built object
    geom2 = [[],[]]
    polys = [(1,2,3,4) ]
    pts = vflo.calc_2d_bbox('z', aspts=True)
    vflo.insert_polygons(polys, pts ) 

    vflo.triangulate()


    vflo.cvt_obj3d_grpoly()
    vflo.export_geojson_lines(outfolder, '%s.json'%outfname         , periodic=True)
    vflo.export_sorted_extents(outfolder, '%s_extents.json'%outfname               )
    vflo.export_sorted_centroids(outfolder, '%s_centroid.json'%outfname            )

    vflo.export_ngc(1, 0, .1, 2, '%s/%s.ngc'%(outfolder, outfname), do3d=False)


##------------------------------------------------##
def make_vector_rastrum(infile, outfile, fid=0):
    """ 
        DEBUG - too slow to use - but interesting 

        take a 3d polygon(s) - 
            (for now keep it flat on Z axis,
             millin_ops scanlines will do a fully 3D triangle someday)

        get the 3d extents and centroid 
        generate a grid of dots to fire rays from 
        capture the intersections as lines  
    """

    vflo = vectorflow()
    vflo.load(infile)
    vflo.triangulate()
    
    vflo.get_face_geom(fid)

    ray = (vec3(0,0,-1), vec3(0, 0, 2))    #x axis ray 
    #ray = (vec3(0,-1,0), vec3(0, 2, 0))   #y axis ray 
    #ray = (vec3(-1, 0, 0), vec3(2, 0, 0)) #z axis ray 

    hits = vflo.ray_hit( ray[0],  ray[1])

    vflo.one_vec_to_obj(ray[1], pos=ray[0])
    vflo.save(outfile)


##------------------------------------------------##
def export_obj_asjson(infile, outfolder, outfname):
    import math 

    vflo = vectorflow() 
    vflo.load(infile )

    # prim_square(self,  pos=(0,0,0), rot=(0,0,0), size=1):
    # prim_triangle(self, pos=(0,0), rot=(0,0), size=1):

    #procedural object 
    for a in range(10): 
        px = math.sin(vflo.mu.dtr(a*30))
        py = math.cos(vflo.mu.dtr(a*30))
        vflo.prim_triangle(pos=(px,py,0), rot=(0,0,a*30), size=.2 )

    #scratch built object
    geom2 = [[],[]]
    polys = [(1,2,3,4) ]
    pts = vflo.calc_2d_bbox('z', aspts=True)
    vflo.insert_polygons(polys, pts ) 

    vflo.triangulate()


    vflo.cvt_obj3d_grpoly()
    vflo.export_geojson_lines(outfolder, '%s.json'%outfname         , periodic=True)
    vflo.export_sorted_extents(outfolder, '%s_extents.json'%outfname               )
    vflo.export_sorted_centroids(outfolder, '%s_centroid.json'%outfname            )

    vflo.export_ngc(1, 0, .1, 2, '%s/%s.ngc'%(outfolder, outfname), do3d=False)


##------------------------------------------------##
def prim_example( outfolder, outfname):
    vflo = vectorflow()
    
    #vflo.prim_circle('z', pos=(0,0,0), rot=(0,0,0), dia=1, spokes=9)

    vflo.prim_triangle('z', pos=(0,0,0), rot=(0,0,0), size=1, flip=False)

    vflo.rotate( 0,0,-37.5 )

    vflo.cvt_obj3d_grpoly()

    vflo.export_geojson_lines(outfolder, '%s.json'%outfname         , periodic=True)
    vflo.export_sorted_extents(outfolder, '%s_extents.json'%outfname               )
    vflo.export_sorted_centroids(outfolder, '%s_centroid.json'%outfname            )


##------------------------------------------------##
def spline_to_cnc(outfolder, outfname):
    vflo = vectorflow()

    start = (1 ,  0, 0)
    ctrl1 = (.5,  0, 0)
    ctrl2 = ( 0, .5, 0)
    end   = (0 ,  1, 0)
    kurve = [start, ctrl1, ctrl2, end]

    start  = (0   ,  1   , 0)
    ctrl1  = (1.5 ,  1.0   , 0)
    ctrl2  = (1   ,  1.5 , 0)
    end    = (2   ,  2   , 0)
    kurve2 =[start, ctrl1, ctrl2, end]

    curves = [kurve, kurve2] 

    #vflo.linegeom_fr_points( [start, ctrl1, ctrl2, end] )

    vflo.draw_splines( 10, curves, drawctrls=False, drawhulls=False, singlepoly=True)
    vflo.cvt_obj3d_grpoly()

    #vflo.export_geojson_lines(outfolder, '%s.json'%outfname         , periodic=True)
    vflo.export_sorted_extents(outfolder, '%s_extents.json'%outfname               )
    vflo.export_sorted_centroids(outfolder, '%s_centroid.json'%outfname            )

    vflo.export_ngc(1, 0, .1, 2, '%s/%s.ngc'%(outfolder, outfname), do3d=False)

    #vflo.save(outfile)   


##------------------------------------------------##

def scratch_obj_to_cnc(outfile):
    vflo = vectorflow()

    geom  = [[],[]]
    geom2 = [[],[]]

    #add new geom and auto increment the ids
    polys = [(1,2,3), (2,3,4) ]
    pts = [(1,1,1),(0,1,1),(-1,-1,1),(2,-2,1)]
    geom = vflo.insert_polygons(polys, pts, geom_obj=geom) 

    polys = [(1,2,3,4) ]
    pts = [(4,-4.3,-3),(1.5,-2.5,-2.1),(-2,2,-4),(4,-4.2,1)]
    geom2 = vflo.insert_polygons(polys, pts, geom_obj=geom2) 

    # use insert to add geom to object 
    vflo.insert(geom) 
    vflo.insert(geom2) 

    # see what we have done, or not done 
    vflo.show()
    
    #vflo.move(2,2,2)
    #Y axis in OBJ is Z axis on CNC - so - rotate Z  
    vflo.rotate(0,0,0)


    vflo.cvt_obj3d_grpoly()

    # do3d option is not done, but wow this shows potential 
    vflo.export_ngc(1, 0, .1, 2, outfile, do3d=False)

    #vflo.save("3d_obj/foo.obj")