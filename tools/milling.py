

from gnelscript.pygfx.obj3d import  *

from gnelscript.pygfx.milling_ops import *






##----------------------------------------##


def intersect_hemisphere( path , infile, heights):
    DEBUG_OUTLINES = False  

    for hgt in heights:

        cop = cam_op()
        
        #mush more to figure out here 
        #negative values for height dont work, you need to flip the normal 
        #DEBUG not all points get sorted
         
        normal = vec3(0,0,-1)
        origin = vec3(0,0,hgt)

        poly = cop.tm_meshplane_test(origin, normal, path , infile)
         
        ########## 
        pts   = [] 
        faces = []
        for i,pt in enumerate(poly): 
            pts.append(tuple(pt[0]))
            pts.append(tuple(pt[1]))
        for f in range(0,len(pts),2):
            faces.append([f+1,f+2])

        print('# hgt %s sorted %s pts   '%(hgt, len(pts)   ) )
        print('# hgt %s face has %s ids '%(hgt, len(faces) ) )
        
        if DEBUG_OUTLINES:
            obj = object3d()
            obj.points = pts 
            obj.polygons = faces    
            obj.save('output_%s.obj'%hgt) 

        #convert 3d points to a dict() 
        d = dict()
        for i,pt in enumerate(poly):
            d[i] = [tuple(pt[0]),tuple(pt[1])]

        # run the sort 
        pathids = cop.sort_linesegs( len(d), 0, d[0], d )

        # build a sorted list 
        sortedpts = []
        ply =[]
        for i,x in enumerate(pathids):
            sortedpts.append( d[x][0])
            ply.append(i+1)

        obj2 = object3d() 
        obj2.points = sortedpts 
        obj2.polygons = [ply]    
        obj2.save('sorted_%s.obj'%hgt)

##----------------------------------------##



"""
def merge_numeric_pairs( path , infile, heights):
    # this was the tool that I first solved the point sorting problem 
    
    DEBUG_OUTLINES = False  

    for hgt in heights:

        cop = cam_op()
        
        #mush more to figure out here 
        #negative values for height dont work, you need to flip the normal 
        #DEBUG not all points get sorted

        normal = vec3(0,0,-1)
        origin = vec3(0,0,hgt)

        poly = cop.tm_meshplane_test(origin, normal, path , infile)
         
        ########## 
        pts   = [] 
        faces = []
        for i,pt in enumerate(poly): 
            pts.append(tuple(pt[0]))
            pts.append(tuple(pt[1]))
        for f in range(0,len(pts),2):
            faces.append([f+1,f+2])

        print('# hgt %s sorted %s pts   '%(hgt, len(pts)   ) )
        print('# hgt %s face has %s ids '%(hgt, len(faces) ) )
        
        if DEBUG_OUTLINES:
            obj = object3d()
            obj.points = pts 
            obj.polygons = faces    
            obj.save('output_%s.obj'%hgt) 

        #convert 3d points to a dict() 
        d = dict()
        for i,pt in enumerate(poly):
            d[i] = [tuple(pt[0]),tuple(pt[1])]

        # run the sort 
        pathids = cop.sort_linesegs( len(d), 0, d[0], d )

        # build a sorted list 
        sortedpts = []
        ply =[]
        for i,x in enumerate(pathids):
            sortedpts.append( d[x][0])
            ply.append(i+1)

        obj2 = object3d() 
        obj2.points = sortedpts 
        obj2.polygons = [ply]    
        obj2.save('sorted_%s.obj'%hgt)
"""

##----------------------------------------##
def tm_multiplane(path, infile): 
    cop = cam_op()
    heights = [.1]
    normal = vec3(0,0,1)
    origin = vec3(0,0,0)

    planes = cop.tm_multiplane_test(heights, origin, normal, 5, path, infile, axis='y')
     
    vflo = vectorflow() 

    print(type(planes))

    #pts = cop.cvt_2d_to_3d(planes)
    #vflo.gr_polys = [cop.ltp(np.asarray(pts))]
    
    vflo.gr_polys = []    
    vflo._sort() 


    #vflo.export_geojson_polygon(path, 'heed.json') 
    vflo.export_geojson_lines(path, 'heed.json') 

    #print(planes[0][0][0])
    #print( type(planes[0]) ) 
    #print(planes.shape)



#tm_multiplane(GLOBAL_PROJ, '/3d_obj/head.obj') 





##----------------------------------------##
def tm_section(path, infile, numdivs): 
    """section ALMOST does what I want - but it requires a watertight mesh 

    """

    cop = cam_op()
    normal = vec3(0,1,1)
    origin = vec3(0,0,0)

    pts2d, m44 = cop.tm_section_test( origin, normal, 5, GLOBAL_PROJ, infile, axis='y')
    
    #poly_union = shp_mltply(pts2d)
    # print(type(poly_union))
    
    print('#tm_section - found %s polygons'%len(pts2d))


    ##--------------
    #method 1 - untested 
    # xe,ye = polydiff.exterior.xy
    # for LinearRing in polydiff.interiors:
    # xi,yi = LinearRing.xy

    ##--------------
    #method 2
    from shapely.geometry import mapping
    obj = mapping(pts2d[0])

    ##--------------
    #method 3 - untested
    #poly_union = shapely.geometry.MultiPolygon(pts2d)


    vflo = vectorflow()

    pts3d=  vflo.cvt_2d_to_3d(list(obj['coordinates'][0])) 

    vflo.gr_polys.append(pts3d) 
    vflo._sort()

    #vflo.export_geojson_polygon(path, 'trimesh.json') 
    vflo.export_geojson_lines (path, 'trimesh.json') 
    vflo.cvt_grpoly_obj3d(objtype='singlepoly')
    vflo.save('heed.obj')
  

    #return trimesh.path.segments.clean(out) 
    #ee = trimesh.path.path.Path(out)


#tm_section(GLOBAL_PROJ, '/3d_obj/monkey.obj', 3)


##----------------------------------------##
class gear_generator(object3d):
    def __init__(self):
        super().__init__()         
        #self.contact_angle
        self.tooth_spacing = .20 
        self.num_teeth = .5
        self.shaft_hole_dia = .4


    def build(self, shaftdia, dia, teeth_height, numteeth ):
        outer = self.calc_circle(pos=(0,0,0), rot=(0,0,0), dia=(dia+teeth_height), axis='z', periodic=True, spokes=36)
        inner = self.calc_circle(pos=(0,0,0), rot=(0,0,0), dia=dia               , axis='z', periodic=True, spokes=36)
        shaft = self.calc_circle(pos=(0,0,0), rot=(0,0,0), dia=shaftdia          , axis='z', periodic=True, spokes=36)

        outer.extend(inner)
        outer.extend(shaft)
        
        return outer 


##----------------------------------------##

class gear_mesh_generator(object3d):
    def __init__(self):
        self.gear1 = gear_generator()
        self.gear2 = gear_generator()




##----------------------------------------##
##----------------------------------------##
##----------------------------------------##

class placeholder(object):

    
    def tm_section_test(self, origin, normal, numdivs, path, infile, axis='y'):

        if type(origin)==vec3:
            origin = origin.aspt
                    
        if type(normal)==vec3:
            normal = normal.aspt

        mesh = trimesh.load_mesh('%s/%s'%(path,infile))

        slice2d = mesh.section(
            plane_normal=normal,
            plane_origin=origin,

        )
        
        if slice2d==None:
            print('### tm_section_test NO POLYGONS FOUND \n\n')
            return None 

        geom, m44 = slice2d.to_planar()
        pts = [poly for poly in geom.polygons_full]


        ##-------------------------------## 
        # def cross_section(mesh, plane_origin=[0,0,0], plane_normal=[1,0,0]):
        #     slice_ = mesh.section(plane_origin=plane_origin, 
        #                           plane_normal=plane_normal)
        #     # transformation matrix for to_planar 
        #      to_2D = trimesh.geometry.align_vectors(plane_normal, [0,0,-1])
        #     
        #     slice_2D, to_3D = slice_.to_planar(to_2D = to_2D)
        #     return slice_2D, to_3D
        
        #poly_union = shapely.geometry.MultiPolygon([poly for poly in slice_2D.polygons_full])

        return pts,m44

    ##-------------------------------##     
    def kdag_merge_points(self, pts):
        """ experiment with DG graphs in an attempt to merge points 
        """
        
        import networkx as nx
        G = nx.Graph()

    ##-------------------------------##     
    def nx_merge_points(self, pts):
        """ experiment with DG graphs in an attempt to merge points 
        """
        
        import networkx as nx
        G = nx.Graph()

    ##-------------------------------## 
    def boundary(self, mesh, close_paths=True):
        #https://stackoverflow.com/questions/76435070/how-do-i-use-python-trimesh-to-get-boundary-vertex-indices
        
        from collections import defaultdict 

        # Set of all edges and of boundary edges (those that appear only once).
        edge_set = set()
        boundary_edges = set()

        # Iterate over all edges, as tuples in the form (i, j) (sorted with i < j to remove ambiguities).
        # For each edge, three cases are possible:
        # 1. The edge has never been visited before. In this case, we can add it to the edge set and as a boundary
        #    candidate as well.
        # 2. The edge had already been visited once. We want to keep it into the set of all edges but remove it from the
        #    boundary set.
        # 3. The edge had already been visited at least twice. This is generally an indication that there is an issue with
        #    the mesh. More precisely, it is not a manifold, and boundaries are not closed-loops.
        for e in map(tuple, mesh.edges_sorted):
            if e not in edge_set:
                edge_set.add(e)
                boundary_edges.add(e)
            elif e in boundary_edges:
                boundary_edges.remove(e)
            else:
                raise RuntimeError(f"The mesh is not a manifold: edge {e} appears more than twice.")

        # Given all boundary vertices, we create a simple dictionary that tells who are their neighbours.
        neighbours = defaultdict(lambda: [])
        for v1, v2 in boundary_edges:
            neighbours[v1].append(v2)
            neighbours[v2].append(v1)

        # We now look for all boundary paths by "extracting" one loop at a time. After obtaining a path, we remove its edges
        # from the "boundary_edges" set. The algorithm terminates when all edges have been used.
        boundary_paths = []

        while len(boundary_edges) > 0:
            # Given the set of remaining boundary edges, get one of them and use it to start the current boundary path.
            # In the sequel, v_previous and v_current represent the edge that we are currently processing.
            v_previous, v_current = next(iter(boundary_edges))
            boundary_vertices = [v_previous]

            # Keep iterating until we close the current boundary curve (the "next" vertex is the same as the first one).
            while v_current != boundary_vertices[0]:
                # We grow the path by adding the vertex "v_current".
                boundary_vertices.append(v_current)

                # We now check which is the next vertex to visit.
                v1, v2 = neighbours[v_current]
                if v1 != v_previous:
                    v_current, v_previous = v1, v_current
                elif v2 != v_previous:
                    v_current, v_previous = v2, v_current
                else:
                    # This line should be un-reachable. I am keeping it only to detect bugs in case I made a mistake when
                    # designing the algorithm.
                    raise RuntimeError(f"Next vertices to visit ({v1=}, {v2=}) are both equal to {v_previous=}.")

            # Close the path (by repeating the first vertex) if needed.
            if close_paths:
                boundary_vertices.append(boundary_vertices[0])

            # "Convert" the vertices from indices to actual Cartesian coordinates.
            boundary_paths.append(mesh.vertices[boundary_vertices])

            # Remove all boundary edges that were added to the last path.
            boundary_edges = set(e for e in boundary_edges if e[0] not in boundary_vertices)

        # Return the list of boundary paths.
        return boundary_paths


    ##-------------------------------## 
    """
    # FROM https://lukeparry.uk/tag/slicing/

    # Generate a list of slices for each layer across the entire mesh
    zSlices = np.linspace(-5, 5, 200)
    k = 1000 # number of slices
    zMin, zMax = myMesh.bounds[:,2]
    zBox = np.linspace(zMin, zMax, k)


    tris = myMesh.triangles

    # Obtain the min and maximum Z values across the entire mesh
    zVals = tris[:, :, 2]
    triMin = np.min(zVals, axis=1)
    triMax = np.max(zVals, axis=1)


    # Below is a manual approach to sorting and collecting the triangles across each layer 

    if False:

        triSortIdx = np.argsort(triMinMax[1,:])
        triMinMaxSort = triMinMax[:,triSortIdx]

        startTime = time.time()

        sortTris = []

        iSects2 = []
        for i in range(len(zBox)):
            minInside = zBox[i].reshape(-1,1) > triMinMax[0, :].reshape(1, -1)
            maxInside = zBox[i].reshape(-1,1) < triMinMax[1, :].reshape(1, -1)
            iSects2.append(np.argwhere((minInside & maxInside).ravel()))

        print('endTime array base', time.time() - startTime)

    # An alteratnvie more compact way is to use binary search operator available within
    # numpy.searchsorted. This locates the bottom and top layer position 

    minIdx = np.searchsorted(zBox, triMin, side='left')
    maxIdx = np.searchsorted(zBox, triMax, side='left')

    # Attach the corresponding presorted triangles into the 
    iSects = [[] for i in range(len(zBox))]

    # The iterative part for assigning potential triangles for intersection on a triangle are performed here
    # Note: Process is very inefficient in native Python O(n*k) 
    for i in range(len(minIdx)):

        startLayer = minIdx[i]
        endLayer = maxIdx[i]
        for layer in iSects[startLayer:endLayer+1]:
            layer.append(i)


    #########################
    The pre-processing has been completed, now the final slincg may be complete. 
    Performing some micro-optimisations, further refactoring may be done to adapt the code previously present in Trimesh.

    The process of slicing or intersection, is simple. 
    Suprisingly, there are a no obvious references for this process. 
    Slicing a triangular mesh, relies firstly computing the potential intersection of each edge of a faceted mesh, 
    in this case it is generalised for an abritrary plane. Firstly, the vector between the mesh vertices 
     and the plane origin is calculated – red lines in the diagram. 

     The dot product is taken with the slicing plane normal nn. 
     The sign of the dot product indicates if the point lies above or below the plane – zero uniquely is the intersection.

    plane_origin = np.array([0,0,0])
    plane_normal = np.array([0,0,1])

    vertex_dots = np.dot(myMesh.vertices - plane_origin, plane_normal)

    """



    ##---------------------------------------------
    
 
    def ltp(self, linestrings):
        """
        ltp = linestrings_to_polygon

        Return valid polygon for unordered set of 3D linestrings in numpy array (n,2,3)
        
        """
        n_vertices, dim = linestrings.shape[0] * linestrings.shape[1], linestrings.shape[2]
        vertices = linestrings.reshape(n_vertices, dim)
        uniq, index, inv = np.unique(vertices.round(decimals=4), return_index=True, return_inverse=True, axis=0)
        assert uniq.shape[0] == vertices.shape[0] / 2
        linestrings_indices = inv.reshape(inv.shape[0]//2, 2)
        
        G = nx.from_edgelist([(a,b) for a,b in linestrings_indices])
        polygon = shapely.geometry.shp_ply(vertices[index][[a for a, b in nx.find_cycle(G)]])
        assert polygon.is_valid
        return polygon
    
  

    # if you do trimesh.load_path(sections) It will give you a Path3D object 
    # (or if you want to call the specific function, it's trimesh.path.io.misc.lines_to_path). 
    # If you transform the sections onto the plane either manually or through path.to_planar, 
    # It will then have path.polygons_full constructed for you.
    # As for the result of triangulate_polygon, that looks correct- the vertices are 2D because 
    # it's coming from a planar polygon, you just need to np.column_stack it with zeros to use it in a 3D mesh.
