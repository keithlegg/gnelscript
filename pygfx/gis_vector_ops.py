
from gnelscript import GEOJSON_IS_LOADED, SHAPELY_IS_LOADED

if GEOJSON_IS_LOADED:
    from geojson import dump 

    from geojson import Point as gjpt
    from geojson import Polygon as gjply
    from geojson import Feature as gjftr
    from geojson import LineString as gjln
    from geojson import FeatureCollection as gjfc
    from geojson import load as gjload

if SHAPELY_IS_LOADED:
    from shapely import buffer, BufferCapStyle, BufferJoinStyle
    from shapely import Point as shp_pt
    from shapely import Polygon as shp_ply
    from shapely import LineString as shp_ln

from gnelscript.pygfx.math_ops import  *
from gnelscript.pygfx.grid_ops import  *
from gnelscript.pygfx.raster_ops import  *

from gnelscript.pygfx.milling_ops import  *
#from gnelscript.pygfx.kicad_ops import  *


from gnelscript.pygfx.obj3d import  *
from gnelscript.pygfx.obj2d import  *


from gnelscript.pygfx.render import simple_render




class work_order(object3d):
    """container for cnc_op (defined in milling_ops.py)
    """

    def __init__(self):
        super().__init__() 

    def add(self):
        pass 



class vectorflow(object3d):
    """ 
        coordinate agnostic geometery processor 

        MAKE IT AS 3D AS POSSIBLE !!! - ONLY MAKE 2D WHEN REQUIRED  !!!

        partially supports kicad graphics, gdcode, json and obj files   

        It started out as an importer for kicad files (graphic polygons and ???) to GCODE tool.
        It turned into a tool turn geojson into GCODE  
        then the dag_ops tesselator was added for GCODE optimization and I got derailed making MC escher style vector renders  
        
        What it does: 
            it can load a whole geojeson file or selected indeces of a json file 
            it can append multiple geojson files together (beware the index changing though)

        ##-------------------------------##
        # IDEA 
        #1 - DAG - b tree 
        #2 - parent geom to nodes - add matrix scenegraph  
        #3-  scengraph for vector animation (attach a matrix to each node) 

    """

    def __init__(self):
        super().__init__()  
        
        self.outfile = []
        
        self.vf_gl_scale =  0.0393701 #NOT FULLY IMPLEMENTED - inch to mm 

        self.mu          = math_util()
        self.tesl        = tessellator()
        self.pop2d       = object2d() 
        self.pop3d       = object3d() 

        #self.kiparser      = pcbfile()
        #self.gcodearser    = gcode_object()

        # geometry buffers for JSON, NGC,sorting, processing, etc 
        self.gr_polys      = [] # list of list of points 
        self.gr_sort       = [] #               [[id, centroid, extents, len, points ]]  
        #self.gr_sort2      = [] # work buffer - [[id, centroid, extents, len, points ]]  

        self.ngc_buffer    = [] # list of list of points 
        self.jsonbuffer    = [] # list of list of points 

        # GEOM ARRAYS for export   
        self.ngc_to_obj    = [] # text buffer for obj 
        self.filled_polys  = [] # list of list of points 

        self.omit_ids      = [] #list of feature ids to NOT export (but leave in)

        self.rh = 1            # retract height 
        self.ch = .1           # cut height (top, start of cut)
        self.cdpi = .01        # cut depth per iteration on Z axis
        self.cmax = 1          # maximum cut depth on Z axis 

        self.hp = (0,0,self.rh) # home position 

        self.orig_minx = 0
        self.orig_miny = 0         
        self.orig_minz = 0
        self.orig_maxx = 0
        self.orig_maxy = 0
        self.orig_maxz = 0

        self.sort_minx = 0
        self.sort_miny = 0         
        self.sort_minz = 0
        self.sort_maxx = 0
        self.sort_maxy = 0
        self.sort_maxz = 0


    ##-------------------------------##
    def insert_gr_sort (self, points ):
        # add a new poly with all the fixins 
         
        #is it 2D or 3D?? 
        #self.checkdata(points)
        
        idx = len(self.gr_sort)
        
        cen = self.centroid(points)
        bbox = self.calc_3d_bbox(points)

        self.gr_sort.append([idx , cen, bbox, len(points), points])

        #               [[id, centroid, extents, len, points ]]  


    ##-------------------------------##
    ##-------------------------------##

    ##-------------------------------##
    """

    https://stackoverflow.com/questions/76435070/how-do-i-use-python-trimesh-to-get-boundary-vertex-indices

    def boundary(mesh, close_paths=True):
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



    def boundary_indices(self):
        
        # Get the indices of vertices that are on the boundaries of the mesh.
        # For each of n boundaries, return the boundary indices.
        # Returns:
        #     boundaries (ndarray): An (number of boundaries) x (variable) array of boundary face indices.
        
        # Proposed solution from StackExchange:
        # https://stackoverflow.com/questions/76435070/how-do-i-use-python-trimesh-to-get-boundary-vertex-indices/76907565#76907565
        connections = self.trimesh_obj.edges[trimesh.grouping.group_rows(
            self.trimesh_obj.edges_sorted, require_count=1)]

        # Start with the first vertex and then walk the graph of connections
        if (len(connections) == 0):
            # No boundaries!
            log.error("Mesh has no boundary edges!")
            raise ValueError("Mesh has no boundary edges!")

        # Tweak: Re-order the first entry, lowest first
        connections[0] = sorted(connections[0])

        # Use ChatGPT to reduce connections to lists of connected vertices
        adj_dict = {}
        for conn in connections:
            for vertex in conn:
                adj_dict.setdefault(vertex, []).extend(c for c in conn if c != vertex)

        groups = []
        visited = set()

        def dfs(vertex, group):
            group.append(vertex)
            visited.add(vertex)
            for conn_vertex in adj_dict[vertex]:
                if conn_vertex not in visited:
                    dfs(conn_vertex, group)

        for vertex in adj_dict:
            if vertex in visited:
                continue
            group = []
            dfs(vertex, group)
            groups.append(group)

        # Convert to numpy arrays
        boundaries = np.empty((len(groups)), dtype=object)
        for index, boundary in enumerate(groups):
            # Close the boundary by add the first element to the end
            boundary.append(boundary[0])
            new_array = np.asarray(boundary, dtype=int)
            boundaries[index] = new_array

        return boundaries

        """

    ##-------------------------------##        
    # def reindex_contiguous_segs (self, linesegs ):
    #     """ repair the messy trimesh plane geom 
    #         linesegs = list of list of line segments
    #                    only works if all segments are "watertight"  
    #     """
    #     fixed = [] 
    #     for seg in linesegs: 
    #         print(seg)


    ##-------------------------------##
    ##-------------------------------##
    #def scale_to_unit(self, inunits='mm', outunits='mm'):
    #    pass


    ##-------------------------------##
    def machine_size(self, units='inch', scale=None):

        mm_to_inch = 25.4

        # laser cutter 
        maxx = 1000
        maxy = 1200
        
        # bridgeport 
        #maxx = 1000
        #maxy = 1200

        if scale and scale!=0:
            maxx = maxx*scale
            maxy = maxy*scale              

        if units=='inch':
            return float(maxx/mm_to_inch),float(maxy/mm_to_inch) 

        if units=='mm':        
            return maxx, maxy

    @property
    def machine_center(self):
        maxx, maxy = self.machine_size() 
        return [ abs(maxx/2), abs(maxy/2) ]

    ##-------------------------------##        
    def export_machine_size(self, folder, name, type='laser_ngc'):
        """Trust the machine"""

        maxx, maxy = self.machine_size
        
        zheight = 0 

        tl = (0    , maxy , zheight)
        tr = (maxx , maxy , zheight)
        br = (maxx , 0    , zheight)
        bl = (0    , 0    , zheight)

        poly = [tl,tr,br,bl,tl]
        
        xx = vectorflow()
        xx.gr_polys.append(poly)
        xx._sort()
        if type=='laser_ngc':
            xx.export_ngc(1, 0, .1, 2, '%s/%s.ngc'%(folder,name), do_laser=True, do3d=False, do_retracts=False) 

    ##-------------------------------##
    def export_poly_rawpts(self, fid, filename):
        """export a single polygon for playing with  
        """
        
        outfile =[] 
        outfile.append(self.gr_sort[fid][4])           

        print('## EXPORTING file %s'%filename)
        with open(filename, 'w') as f:
            for l in outfile:
                f.write("%s\n" % l)

    ##-------------------------------##
    def polysize_info(self):
        biggest_x = 0 
        biggest_y = 0

        smallest_x = 100
        smallest_y = 100

        tmp = []
        for ply in self.gr_sort:
            bbox = ply[2]
            
            ply_sizex = abs(bbox[2]-bbox[0])
            ply_sizey = abs(bbox[3]-bbox[1])
            
            if ply_sizex > biggest_x:
                biggest_x = ply_sizex
            if ply_sizey > biggest_y:
                biggest_y = ply_sizey
            
            if ply_sizex < smallest_x:
                smallest_x = ply_sizex  
            if ply_sizey< smallest_y:
                smallest_y = ply_sizey  

        print('# done filtering polygons by size')
        print('# largest  :', [biggest_x ,biggest_y]  )        
        print('# smallest :', [smallest_x,smallest_y] )


    ##-------------------------------##        
    def filter_by_bbox(self, xsize, ysize, underover='bigger'):
        """ sifting screen for polygon bboxes. 
            filter out big or small based on threshold 
            
            ARGS:
                underover - 
                      bigger  : return only polygons bigger than the threshold
                      smaller : return only polygons smaller than the threshold

        """

        tmp = []
        for ply in self.gr_sort:
            bbox = ply[2]
            
            ply_sizex = abs(bbox[2]-bbox[0])
            ply_sizey = abs(bbox[3]-bbox[1])

            if underover=='bigger':
                if ply_sizex > xsize or ply_sizey>ysize:
                    tmp.append(ply)

            if underover=='smaller':
                if ply_sizex < xsize or ply_sizey<ysize:
                    tmp.append(ply)

        # debug hmmm - not sure on workflow here 
        #self.gr_sort = tmp
        return tmp       

    ##-------------------------------##
    def flush(self):
        #DEBUG - no worky, no testy 
        print('flushing geometry buffers')
        self.gr_polys      = []  
        self.gr_sort       = []   
        self.ngc_buffer    = []  
        self.jsonbuffer    = []  
        self.ngc_to_obj    = [] 
        self.filled_polys  = []  
        self.omit_ids      = [] 
        self.outfile       = []

    ##-------------------------------##
    def scanlines_to_segments(self, lines):
        """ DEBUG - NOT WORKING - need to skip over blank spaces (ODD pts)
            take output of scanline_ngon and draw the fills for each line 
              
        """
        zcoord = 0 

        for l in lines:
            #if length is EVEN - polygon has workable topolgy 
            if len(l)%2==0:
                #print(len(l))
                for i in range(0,len(l),2):
                    #print(i)
                    spt =  (l[i-1][0], l[i-1][1], zcoord)                         
                    ept =  (l[i][0]  , l[i][1]  , zcoord) 
                    self.pts_to_linesegment([spt,ept], periodic=False)

    ##-------------------------------##
    def scanline_ngon(self, numh, numv, drwply):
        """
            get extents of polygon 
            build up a grid of horizontal scanlines
            iterate each line segment testing for 2D intersection 
            return a list of those for each line
            if the number is EVEN assume its multiple lines 
            if the number is ODD is has poor topology 

        """
        #derive 2D BBOX from 3D function  
        bbox_pts = self.calc_2d_bbox( axis='z', pts=drwply, aspts=True)        
        bbox     = self.calc_2d_bbox( axis='z', pts=drwply)

        out_pts = [] 
        
        debug = False 

        #minx, miny , maxx, maxy
        res_x = bbox[2]-bbox[0]
        res_y = bbox[1]-bbox[3]
        
        xhalf = res_x/2 
        yhalf = res_y/2 
        ydivs = (res_y/numv)

        # render the bounding box (OBJ) 
        if debug:
            self.pts_to_linesegment(bbox_pts, periodic=False)
        
        # center needs to derive from calc_bbox
        bbox = self.calc_3d_bbox(pts=bbox_pts)

        center = (abs(bbox[0]+bbox[3])/2, abs(bbox[1]+bbox[4])/2)
        
        #to debug center calc  
        #print(center)
        #self.prim_circle(pos=(center[0], center[1], 0), dia=.2, axis='z')

        vecmath = vec2()     # use for math operations

        # build up line data for three sides of triangle
        s1 = ( drwply[0][0], drwply[0][1] )
        e1 = ( drwply[1][0], drwply[1][1] )
        #
        s2 = ( drwply[1][0], drwply[1][1])
        e2 = ( drwply[2][0], drwply[2][1])
        #
        s3 = ( drwply[2][0], drwply[2][1])
        e3 = ( drwply[0][0], drwply[0][1])


        #render the triangle  (OBJ)
        if debug:
            self.pts_to_linesegment(drwply, periodic=True)

        # define the scanline geometry, iterate each horizontal line of image
        for hscan in range(1,numv):

            thisline = []

            #ypos = center[1]+(hscan*ydivs)
            ypos = (center[1]-yhalf)+(hscan*ydivs)
            # build the scan vector 

            #DEBUG -need to extend beyond extents to catch all hits
            #buffer_dist = res_x
            s_hvec = (center[0]-xhalf, ypos)
            e_hvec = (center[0]+xhalf, ypos)
            
            #render the raster scan lines (OBJ) 
            if debug:
                self.pts_to_linesegment([s_hvec,e_hvec], periodic=False)

            #iterate all points, build line segments and brute force test intersetion with scanline 
            for i,pt in enumerate(drwply):
                if i>0:
                    #2d segment start and end 
                    s_v = drwply[i-1]
                    e_v = drwply[i]
                    
                    #run the test (slow but it works) 
                    hit = vecmath.intersect(s_hvec, e_hvec, (s_v[0],s_v[1]), (e_v[0],e_v[1]) )  
                    if hit:
                        thisline.append( hit )  
              

            #put each row in its own array  
            out_pts.append(thisline)

        return out_pts

    ##-------------------------------##
    def scanline_triangle(self, numh, numv, drwply):
        """
           build a series of lines across a triangle that represent a raster scan
           
           numh - int num width  
           numv - int num height 
           drwply   - tuple of tuple of 3 floats, xyz  

           returns [numh*numv]

        """
        
        #derive 2D BBOX from 3D function  
        bbox_pts = self.calc_2d_bbox( axis='z', pts=drwply, aspts=True)        
        bbox     = self.calc_2d_bbox( axis='z', pts=drwply)
                


        out_pts = [] 

        #minx, miny , maxx, maxy
        res_x = bbox[2]-bbox[0]
        res_y = bbox[1]-bbox[3]
        
        xhalf = res_x/2 
        yhalf = res_y/2 

        ydivs = (res_y/numv)

        # render the bounding box (OBJ) 
        #self.pts_to_linesegment(bbox_pts, periodic=False)
        
        # center needs to derive from calc_bbox
        bbox = self.calc_3d_bbox(pts=bbox_pts)

        center = (abs(bbox[0]+bbox[3])/2, abs(bbox[1]+bbox[4])/2)
        
        #to debug center calc  
        #print(center)
        #self.prim_circle(pos=(center[0], center[1], 0), dia=.2, axis='z')

        vecmath = vec2()     # use for math operations

        # build up line data for three sides of triangle
        s1 = ( drwply[0][0], drwply[0][1] )
        e1 = ( drwply[1][0], drwply[1][1] )
        #
        s2 = ( drwply[1][0], drwply[1][1])
        e2 = ( drwply[2][0], drwply[2][1])
        #
        s3 = ( drwply[2][0], drwply[2][1])
        e3 = ( drwply[0][0], drwply[0][1])


        #render the triangle  (OBJ)
        self.pts_to_linesegment(drwply, periodic=True)

        # define the scanline geometry, iterate each horizontal line of image
        for hscan in range(1,numv):

            thisline = []

            #ypos = center[1]+(hscan*ydivs)
            ypos = (center[1]-yhalf)+(hscan*ydivs)
            # build the scan vector 

            #DEBUG -need to extend beyond extents to catch all hits
            #buffer_dist = res_x
            s_hvec = (center[0]-xhalf, ypos)
            e_hvec = (center[0]+xhalf, ypos)
            
            #render the raster scan lines (OBJ) 
            self.pts_to_linesegment([s_hvec,e_hvec], periodic=False)

            # take the 3 edges of a triangle and determine if the horizontal scanline intersects any 
            i = vecmath.intersect(s_hvec, e_hvec, s1, e1) #left  side of triangle
            j = vecmath.intersect(s_hvec, e_hvec, s2, e2) #right side of triangle
            k = vecmath.intersect(s_hvec, e_hvec, s3, e3) #top   side of triangle


            #out_pts.append([s_hvec,e_hvec])
            #print(i,j,k)

            # debug tool - show the "hits" for the horizontal scanline
            if i: 
                #output.draw_fill_circle( i[0], i[1], 1, (255,0,0) ) 
                thisline.append( (i[0], i[1]) )
            if j: 
                #output.draw_fill_circle( j[0], j[1], 1, (0,255,0) ) 
                thisline.append( (j[0], j[1]) )
            if k: 
                #output.draw_fill_circle( k[0], k[1], 1, (0,0,255) )  
                thisline.append( (k[0], k[1]) )                

            #put each row in its own array  
            out_pts.append(thisline)

        return out_pts

    ##-------------------------------##
    def build_tesselation_sample(self):
        """ DAG + point generator = tesselator """
        print('# generating sample data from tesselator nodes')
        if len(self.tesl.nodes)==0:
            print('# error - no nodes in tesselator: exiting')
            return None 

        for cell in self.tesl.nodes:
            e1m = cell.getattrib('e1mid')[0]
            e2m = cell.getattrib('e2mid')[0]
            e3m = cell.getattrib('e3mid')[0]
            e4m = cell.getattrib('e4mid')[0]
            
            #BBOX center - not quite accurate, use the projected center below
            #cen = cell.getattrib('centroid') 
            tmp = vec2() 
            cen = tmp.intersect_2d_from_3D([e1m,e3m], [e2m,e4m])


            #mark the center with a shape (debug need to use multipolygon instead of linetype )
            #pts = self.pop2d.calc_circle_2d(cen[0],cen[1], .01, periodic=True, spokes=3)
            #cell.points.extend( pts )

            #if you connect the midpoint of edges, it subdivides 
            #cell.points.extend( [e1m,e3m] )
            #cell.points.extend( [e2m,e4m] )

            #lines radiating from center are basically a 2nd order divison 
            #cell.points.extend( [cen, e1m] )
            #cell.points.extend( [cen, e2m] )
            #cell.points.extend( [cen, e3m] )
            #cell.points.extend( [cen, e4m] )                        
            
            #this makes a quad diagonal inside the original 
            cell.points.extend( [e1m, e2m, e3m, e4m] ) 

    ##-------------------------------##
    def render_cells(self):
        """ bake the points stored in each node into the points to be dumped to a file 
        """

        for i,cell in enumerate(self.tesl.nodes):
            #cen = cell.getattrib('centroid')
            if cell.points:
                self.gr_polys.append(self.cvt_2d_to_3d(cell.points))
        
        self._sort()
 

    ##-------------------------------##
    def draw_dots_2d(self, dia, spokes, holes):
        """ render some points as circles
             inspired by hole cutter but just to see some coords in an image 
        """

        for h in holes:
            for pt in h:
                pts = self.pop3d.calc_circle(pos=(pt[0],pt[1],0), rot=(0,0,0), dia=dia, periodic=True, spokes=spokes)
                self.gr_polys.append(pts)
        self._sort()

    ##-------------------------------##
    def hole_cutter_2d(self, tooldia, holes):
        """ simple circle generator for chopping holes with gcode 
            
            ARGS: 
                holes = [ [dia, spokes, (x,y,z)/(x,y) ] ]

            NOTES:
                format of self.gr_sort = [[id, centroid, extents, len, points ]] 
        """

        for h in holes:
            print(h)
            if h:
                dia = h[0]
                spokes = h[1]
                pt = h[2]

                if len(pt)<3:
                    cen = [pt[0],pt[1],self.ch]
                if len(pt)==3:            
                    cen=pt   

                #pos=(0,0,0), rot=(0,0,0), dia=1, axis='z', periodic=True, spokes=23):
                pts = self.pop3d.calc_circle(pos=cen, rot=(0,0,0), dia=dia, periodic=True, spokes=spokes)
        
                self.gr_polys.append(pts)

        self._sort()

    ##-------------------------------##   
    def _set_cam_properties(self, rh, ch, cdpi, cmax):
        self.rh = rh           
        self.ch = ch            
        self.cdpi = cdpi       
        self.cmax = cmax          
    
    ##-------------------------------## 
    def _set_extents(self, bbox):
        """ set global extents for generating data 
            based on PIL coordinate which is [left, top, right, bottom] 
            [minx, miny, maxx, maxy]  
        """
        
        self.sort_minx = bbox[0]
        self.sort_miny = bbox[1]         
        self.sort_maxx = bbox[2]
        self.sort_maxy = bbox[3]

    ##-------------------------------##       
    def _omit_ids(self, ids=None, span=None, unique=True, nth=None):
        
        #DEBUG - not fully working - see indexer docs - nths with crash if you try it 
        pop = pop3d()
        ids = pop.indexer(ids, span, unique, nth)
        
        print(" ### omit ids: ", ids)

        self.omit_ids = ids

    ##-------------------------------##       
    def _make_periodic(self):
        """ DEBUG - SHOULD OPERATE ON GR_SORT 
            if data has polygon that are not closed - add the first point to the end to close 
        """
        for ply in self.gr_polys:
            first = ply[0]
            ply.append(first)

    ##-------------------------------##       
    def _sort(self):
        """ assemble data into [[id, centroid, extents, len, points ]] - put that in self.gr_sort   

            set extents of original data while running   
        """

        #print("indexing sort buffer ")
        self.gr_sort = []
        
        if len(self.gr_polys)==0:
            return None 

        #print(" gr_polys buffer has %s polys in it "%len(self.gr_polys) )
        for x,ply in enumerate(self.gr_polys):
            
            #print(type(ply))

            minx = 0
            miny = 0         
            maxx = 0
            maxy = 0


            try:
                #print('### len ', len(ply))
                for i,pt in enumerate(ply):
                    if i == 0:
                        minx=pt[0]
                        maxx=pt[0]
                        miny=pt[1]
                        maxy=pt[1]

                    if pt[0]<minx:
                        minx=pt[0]    
                    if pt[0]>maxx:
                        maxx=pt[0]  
                    if pt[1]<miny:
                        miny=pt[1]  
                    if pt[1]>maxy:
                        maxy=pt[1] 
            except: 
                print('_sort blew up !! ---------> ',pt)

            if minx<self.orig_minx:
                self.orig_minx=minx
            if miny<self.orig_miny:
                self.orig_miny=miny
            if maxx>self.orig_maxx:
                self.orig_maxx=maxx
            if maxy>self.orig_maxy:
                self.orig_maxy=maxy

            #DEBUG - consider 3D bbox and padding with spaces for future data
            self.gr_sort.append([x, self.centroid(ply) ,[minx,miny, maxx, maxy], len(ply), ply])
        

        # run initial extents calc on sorted polys         
        self.gl_extents()

        print("orig data extents %s %s %s %s "%(self.orig_minx, self.orig_maxx, self.orig_miny, self.orig_maxy) )

    ##-------------------------------------------## 
    ##-------------------------------------------## 
    def show_setup(self):
        print('##--------------------------##')
        print("retract height %s"%self.rh)
        print("cut height     %s"%self.ch)
        print("cut max        %s"%self.cmax)
        print("cut cdpi       %s"%self.cdpi)

        print('original extents  %s %s %s %s '%(self.orig_minx, self.orig_miny, self.orig_maxx, self.orig_maxy))
        print('sorted extents    %s %s %s %s '%(self.sort_minx, self.sort_miny, self.sort_maxx, self.sort_maxy)) 

    ##-------------------------------##
    def show_buffers(self, sid=None):
        """
        
            show evrything meaningful we can about our data 
            id is optional 
            - if passed use it to look up gr_poly and gr_sort, etc  

        """
        
        print('#################')
        if sid is None:
            #[[id, centroid, extents, len, points ]
            #for ply in self.gr_sort:
            #    print("%s %s %s "%(ply[0], ply[1], len(ply[2]) )) 

            print('size gr_polys     %s'%len(self.gr_polys))
            print('size gr_sort      %s'%len(self.gr_sort))

        else:
            grp = self.gr_polys[sid]
            grs = self.gr_sort[sid]
            print('gr_polys id %s size %s  '%( sid, len(grp) ) )
            print('gr_sort  id %s json id:%s extents %s len %s size: %s '%( sid, grs[0], grs[1], grs[2], grs[3] ) )

    ##---------------------- 
    def grply_inspect(self, index=None):
        """
            #return info about how many features, sub features, etc are in a file 
            
            args:
                index is positive int or iterable 

            todo:
                slice , index 
                get extents 
                get centroid 
                get as vectors 

        """
    
        print("##------------##")

        if index == None:
            print("gr_poly buffer has %s polygons "%len(self.gr_polys) )
            for i,p in enumerate(self.gr_polys):
                print("    feature %s has %s points"%(i, len(p)))
                print(p)

        else:
            pass 

    ##-------------------------------------------## 
    ##-------------------------------------------## 
    def gl_extents(self):
        """ global 2D/3D extents 

            we only work on gr_sort - gr_poly is a copy pf the orignial data

            if you run sort() - it automatically sets extents while it is sorting  
            if you want to re-run, use this 
        """

        #reset all values first to random point in model (first)
        #               [[id, centroid, extents, len, points ]] 
        randpt = self.gr_sort[0][4][0]
        self.sort_minx=randpt[0]
        self.sort_miny=randpt[1]
        self.sort_minz=randpt[2]
        self.sort_maxx=randpt[0]
        self.sort_maxy=randpt[1]
        self.sort_maxz=randpt[2]


        # #start with gr_sort - we will do OBJ below
        # [[id, centroid, extents, len, points ]]
        for row in self.gr_sort:
            ply = row[4] 

            minx = 0
            miny = 0         
            minz = 0  
            maxx = 0
            maxy = 0
            maxz = 0

            #print('### len ', len(ply))
            for i,pt in enumerate(ply):
                if i == 0:
                    minx=pt[0]
                    maxx=pt[0]
                    miny=pt[1]
                    maxy=pt[1]
                    
                    if len(pt)==3:
                        #some data will be 2D, some 3D
                        minz=pt[0]
                        maxz=pt[0]

                if pt[0]<minx:
                    minx=pt[0]    
                if pt[0]>maxx:
                    maxx=pt[0]  
                if pt[1]<miny:
                    miny=pt[1]  
                if pt[1]>maxy:
                    maxy=pt[1] 

                if len(pt)==3:
                    #some data will be 2D, some 3D
                    if pt[2]<minz:
                        minz=pt[2]  
                    if pt[2]>maxz:
                        maxz=pt[2] 

            if minx<self.sort_minx:
                self.sort_minx=minx
            if miny<self.sort_miny:
                self.sort_miny=miny
            if minz<self.sort_minz:
                self.sort_minz=minz

            if maxx>self.sort_maxx:
                self.sort_maxx=maxx
            if maxy>self.sort_maxy:
                self.sort_maxy=maxy
            if maxz>self.sort_maxz:
                self.sort_maxz=maxz

        if len(self.points):
            #do OBJ geometry next 
            bbox = self.calc_3d_bbox()
            #return [min_x, min_y, min_z, max_x, max_y, max_z ]
            for e in bbox:
                if bbox[0]<self.sort_minx:
                    self.sort_minx=bbox[0]
                if bbox[1]<self.sort_miny:
                    self.sort_miny=bbox[1]
                if bbox[2]<self.sort_minz:
                    self.sort_minz=bbox[2]

                if bbox[3]>self.sort_maxx:
                    self.sort_maxx=bbox[3]
                if bbox[4]>self.sort_maxy:
                    self.sort_maxy=bbox[4]
                if bbox[5]>self.sort_maxz:
                    self.sort_maxz=bbox[5]


        return [self.sort_minx, self.sort_miny, self.sort_minz, 
                self.sort_maxx, self.sort_maxy, self.sort_maxz ]

    ##-------------------------------------------## 
    def gl_centroid(self):
        """ global centroid - DEBUG - NOT TESTED
            we only work on gr_sort - gr_poly is a copy pf the original data
 
        """

        self.gl_extents()
        
        width  = abs(self.sort_maxx - self.sort_minx) 
        height = abs(self.sort_maxy - self.sort_miny)         
        
        return ( self.sort_maxx-(width/2), self.sort_maxy-(height/2) )
    
    ##-------------------------------------------## 
    def gl_width_height(self):
        self.gl_extents()
        
        width  = abs(self.sort_maxx - self.sort_minx) 
        height = abs(self.sort_maxy - self.sort_miny)         
        
        return (width, height )

    ##-------------------------------------------## 
    def get_extents_poly(self, zheight=0.0):

        self.gl_extents()
        
        width  = abs(self.sort_maxx - self.sort_minx) 
        height = abs(self.sort_maxy - self.sort_miny)   
        center = ( self.sort_maxx-(width/2), self.sort_maxy-(height/2) )
        
        print('## --------------- ')
        print('width  ', width  ) 
        print('height ', height )
        print('center ', center )
        
        tl = (self.sort_minx, self.sort_maxy, zheight)
        tr = (self.sort_maxx, self.sort_maxy, zheight)
        br = (self.sort_maxx, self.sort_miny, zheight)
        bl = (self.sort_minx, self.sort_miny, zheight)

        return [tl,tr,br,bl,tl]


    def gl_bbox(self):
        return self.cvt_3d_to_2d(self.get_extents_poly())

    ##-------------------------------------------## 
    def export_extents_ngc(self, folder, name, scale=None, type='laser_ngc'):
        """ 
        DEBUG - seems off - data is close but wrong 
        """
        poly = self.get_extents_poly()
        xx = vectorflow()
        xx.gr_polys.append(poly)
        xx._sort()
        if scale:
            xx.gl_scale(scale)

        if type=='laser_ngc':
            xx.export_ngc(1, 0, .1, 2, '%s/%s_extents.ngc'%(folder,name), do_laser=True, do3d=False, do_retracts=False)

    ##-------------------------------------------## 
    def gl_move_extents_corner(self, which='tl'):
        """
        shift all data to a corner   
        """
        shift = [0,0]

        bbox = self.gl_bbox()
        if which=='br':
            shift = bbox[0]
        if which=='bl':
            shift = bbox[1]
        if which=='tl':
            shift = bbox[2]
        if which=='tr':
            shift = bbox[3]

        self.gl_translate(shift[0], shift[1], 0)
    

         

    
    ##-------------------------------------------## 
    def gl_rotate(self, amt):
        """
            global rotate the entire dataset in 2d  
            this allows 3D moves but you probably want 2D - zero the Z axis for pos or just send it 2 coords
        """
        print("global rotate ", amt) 

        rotx = 0.0
        roty = 0.0
        rotz = 0.0
        
        if isinstance(amt,tuple) or isinstance(amt,list):
            rotx = amt[0]
            roty = amt[1]
            rotz = amt[2]
        else:
            #rotx = amt
            #roty = amt
            rotz = amt         

        pop = polygon_operator()

        for i,row in enumerate(self.gr_sort):
                # TRS all in one 
                #self.gr_sort[i][4] = pop.trs_points(self.gr_sort[i][4], rotate=(-rot[0], -rot[1], -rot[2] ))

                # ROUNDING WAS CAUSING GLITCHES IN THE NCVIEWER APP - MAY OR MAY NOT BE AN ISSUE 
                self.gr_sort[i][4] = pop.rotate_pts(rot=(-rotx, -roty, -rotz ), pts=self.gr_sort[i][4], doround=True)

        # dont forget to recalculate extents 
        self.gl_extents()

    ##-------------------------------------------## 
    def gl_scale(self, amt):        
        #derived from scale_pts
        # build a scale matrix 
        
        amtx = 1.0
        amty = 1.0
        amtz = 1.0
        
        if isinstance(amt,tuple) or isinstance(amt,list):
            amtx = amt[0]
            amty = amt[1]
            amtz = amt[2]
        else:
            amtx = amt
            amty = amt
            amtz = amt            

        sc_m33 = self.m33.identity
        sc_m33[0]  = amtx
        sc_m33[4]  = amty
        sc_m33[8]  = amtz    

        #for ply in self.gr_polys:
        #    print(ply[0])
        for ply in self.gr_sort:
            ply[4] = sc_m33.batch_mult_pts(ply[4])
        
        # dont forget to recalculate extents 
        self.gl_extents()
    
    ##-------------------------------------------## 
    def gl_move_center(self):
        """
            we only work on gr_sort - gr_poly is a copy pf the orignial data
        """
        print("global tansform to center ") 

        cen = self.gl_centroid()

        pop = polygon_operator()

        for i,row in enumerate(self.gr_sort):
            self.gr_sort[i][4] = pop.trs_points(self.gr_sort[i][4], (-cen[0], -cen[1], 0 ))

        # dont forget to recalculate extents 
        self.gl_extents()

    ##-------------------------------------------##  
    def gl_translate(self, xval, yval, zval):
        """
            we only work on gr_sort - gr_poly is a copy pf the orignial data
        """
        print("global tansform all pts ") 
        pop = polygon_operator()

        for i,row in enumerate(self.gr_sort):
            self.gr_sort[i][4] = pop.trs_points(self.gr_sort[i][4], (xval, yval, zval ))

        # dont forget to recalculate extents 
        self.gl_extents()

    ##-------------------------------------------## 
    def gl_move(self, pos):
        """
            DEBUG - NOT TESTED  
            global move the entire dataset in 2d (ignore Z axis) 
            this allows 3D moves but you probably want 2D - zero the Z axis for pos or just send it 2 coords
        """
        print("global move ", pos) 

        pop = polygon_operator()

        for i,row in enumerate(self.gr_sort):
            if len(pos)==2:
                self.gr_sort[i][4] = pop.trs_points(self.gr_sort[i][4], (-pos[0], -pos[1]) )
            if len(pos)==3:
                self.gr_sort[i][4] = pop.trs_points(self.gr_sort[i][4], (-pos[0], -pos[1], -pos[0] ))

        # dont forget to recalculate extents 
        self.gl_extents()

    ##-------------------------------##
    def _scrub(self, inp):
        """ clean up parsed characters from (kicad) file """

        out = inp.lower()
        out = out.strip()
        out = out.replace('(','')
        out = out.replace(')','')
        out = out.replace(' ','')
        out = out.replace('\"','')
        out = out.replace('\'','')
        return out

    ##-------------------------------##
    def _calculate_paths3d(self, do_laser=False, do_retracts=True, doround=True):
        """ 
            takes gr_sort buffer and makes an executable gcode file with retracts on Z for each polygon
           
            format of self.gr_sort = [[id, centroid, extents, len, points ]] 

        """
        pl = 6 #numeric rounding places 
        lastpt = (0,0,0)

        self.outfile.append('( exported with _calculate_paths3d )')
        self.outfile.append('(linear scale set to %s of internal coordinates)'%self.vf_gl_scale )
        self.outfile.append('  ')

        self.outfile.append('g20')                  #inches for unit 
        
        ##-----------------------------------------##

        #move to origin  
        self.outfile.append('g0 x%s y%s z%s f30'% ( self.hp[0], self.hp[1], self.rh) )   #rapid move to 0 
        self.ngc_to_obj.append( ( self.hp[0], self.hp[1], self.rh) ) 

        if do_retracts:
            self.outfile.append('g0z%s'% ( self.rh ) )  #retract in between cuts
   
        ##-----------------------------------------##
        self.outfile.append('  ')
        self.outfile.append('(exporting filled polygons )')

        ####
        self.outfile.append('  ')
        self.outfile.append('(exporting graphic polygons )')
        
        ##------------------------------

        # graphical polygons - build the gcode up with simple linear movements
        for row in self.gr_sort:
            
            # preprocess 
            export_ply = True 

            if row[0] in self.omit_ids:
                print("# omitting polygon ID %s"%row[0])
                export_ply = False


            #precache the rounded values so we dont confuse the logic of exporting below 
            if doround:
                strpoly = [] 
                for i,npt in enumerate(row[4]):
                    xc = f"{npt[0]:.18f}"
                    yc = f"{npt[1]:.18f}"
                    zc = f"{npt[2]:.18f}"
                    
                    if 'e' in xc or 'e' in yc or 'e' in zc:
                        print('ERROR EXPONENT FOUND ')
                        print(i, (xc,yc,zc))

                    strpoly.append((xc,yc,zc))

                gr_poly = strpoly 
            else:
                gr_poly = row[4]
            ##--

            if len(gr_poly) and export_ply:
                self.outfile.append('(exporting new polygon )')

                pt1 = gr_poly[0]

                ## first point with/without retracts 
                if do_retracts:
                    self.outfile.append( 'G0' )
                    #move to first point RH 
                    self.outfile.append('x%s y%s z%s'% (  pt1[0] , pt1[1], self.rh ) )               
                    self.ngc_to_obj.append( ( pt1[0], pt1[1], self.rh ) )             
                    self.outfile.append( 'G1' )
                else:
                    self.outfile.append( 'G0' )
                    self.outfile.append('x%s y%s z%s'% (  pt1[0] , pt1[1], pt1[2]) )               
                    self.ngc_to_obj.append( ( pt1[0], pt1[1], pt1[2] ) )             
                    self.outfile.append( 'G1' )

                ## iterate points in polygon 
                for i,pt in enumerate(gr_poly):
                    self.outfile.append( 'x%s y%s z%s'%( pt[0], pt[1],  pt[2] ) )
                    self.ngc_to_obj.append( (pt[0], pt[1],  pt[2]) )                   

                if do_retracts:
                    self.outfile.append( 'G0' )
                    gpt=gr_poly[0]

                    self.ngc_to_obj.append( (gpt[0], gpt[1], self.rh)  )
                    self.outfile.append( 'x%s y%s z%s'%( gpt[0], gpt[1], self.rh) )

                    # retract in between cuts
                    self.outfile.append('g0z%s'% ( self.rh ) )  


                self.outfile.append('  ')

        ##-----------------------------------------##        
        # self.outfile.append('(exporting segments )')

        ##-----------------------------------------##
        # rapid move at end 
        self.outfile.append('m2') #program end

    ##-------------------------------##
    def _calculate_paths2d(self, do_laser=False, do_retracts=True, laserpwm=300, doround=True):
        """ 
            SIMPLE SAUCE CAM PROGRAM 

            walks self.gr_poly and builds an OBJ and NGC file in self.outfile

            #milling commands 
            https://linuxcnc.org/docs/html/gcode/g-code.html
            Top 10 tasty GCODE commands:
                S - surface speed
                
                G20 (Use inch)
                G21 (Use mm)

                G90 (Set Absolute Coordinates)
                G0 -  Rapid Move
                G1 -  Linear Move
                
                M3, M4, M5  S ($)   Spindle Control
                M6 - 
                M9 - 
                M3 - 
                M2 - program end 

            #Laser commands  
                LASER ON/OFF could be M103/M105 or it could be M03/M05. depends on GRBL 
                M3 Constant Laser Power Mode
                M4 Dynamic Laser Power Mode
                S-MIN and S-MAX  Power Modulation

            #3d print commands

        """

        pl = 6 #numeric rounding places
        lastpt = (0,0,0)


        self.outfile.append('(exported with _calculate_paths2d )')
        self.outfile.append('(linear scale set to %s of internal coordinates)'%self.vf_gl_scale )
        self.outfile.append('  ')

        self.outfile.append('g20')                  #inches for unit 
        
        ##-----------------------------------------##

        #move to origin  
        self.outfile.append('g0 x%s y%s z%s f30'% ( self.hp[0], self.hp[1], self.rh) )   #rapid move to 0 
        self.ngc_to_obj.append( ( self.hp[0], self.hp[1], self.rh) ) 

        if do_retracts:
            self.outfile.append('g0z%s'% ( self.rh ) )  #retract in between cuts
   
        ##-----------------------------------------##
        self.outfile.append('  ')
        self.outfile.append('(exporting filled polygons )')

        # build the gcode up with simple linear movements 
        ##------------------------------

        """
        #DEBUG - FILL_POLYS ARE A LAYOVER FROM KICAD STUFF - UNTESTED 
        for fill_poly in self.filled_polys:
            if do_retracts:
                pt1 = fill_poly[0] 
                self.outfile.append('G0 x%s y%s z%s'% (  pt1[0], pt1[1], self.rh ) )  #first point at retract height   
                self.ngc_to_obj.append( ( self.hp[0], self.hp[1], self.rh) ) 
            
            for i,pt in enumerate(fill_poly):
                self.outfile.append( 'G1 x%s y%s z%s'%( pt[0], pt[1], self.ch ) )
                self.ngc_to_obj.append( ( pt[0], pt[1], self.ch ) )             
                lastpt =( pt[0], pt[1], self.ch )

            self.outfile.append( 'G0 x%s y%s z%s'%(lastpt[0], lastpt[1], self.ch) )
            self.ngc_to_obj.append( (lastpt[0], lastpt[1], self.ch)   ) 
            if do_retracts:
                self.outfile.append('g0z%s'% ( self.rh ) )  #retract in between cuts
                self.ngc_to_obj.append( (lastpt[0], lastpt[1], self.rh)   ) 
                self.outfile.append('  ')

        """
        ####

        self.outfile.append('  ')
        self.outfile.append('(exporting graphic polygons )')
        
        ##------------------------------
        # graphical polygons - build the gcode up with simple linear movements
        for row in self.gr_sort:
            
            # preprocess 
            export_ply = True 

            if row[0] in self.omit_ids:
                print("# omitting polygon ID %s"%row[0])
                export_ply = False
 
            #precache the rounded values so we dont confuse the logic of exporting below 
            if doround:
                strpoly = [] 
                for i,npt in enumerate(row[4]):
                    #xc = str(round(npt[0],pl))
                    #yc = str(round(npt[1],pl))
                    #zc = str(round(npt[2],pl))
                    xc = f"{npt[0]:.18f}"
                    yc = f"{npt[1]:.18f}"
                    zc = f"{npt[2]:.18f}"
                    
                    if 'e' in xc or 'e' in yc or 'e' in zc:
                        print('ERROR EXPONENT FOUND ')
                        print(i, (xc,yc,zc))

                    strpoly.append((xc,yc,zc))

                gr_poly = strpoly 
            else:
                gr_poly = row[4]

            ##--

            if len(gr_poly) and export_ply:
                self.outfile.append('(exporting new polygon )')

                pt1 = gr_poly[0]

                #### rapid move to first point (with ot without retract)
                if do_retracts and not do_laser:
                    self.outfile.append( 'G0' )    
                    self.outfile.append('x%s y%s z%s'% (  pt1[0] , pt1[1], self.rh ) )               
                    self.ngc_to_obj.append( ( pt1[0], pt1[1], self.rh ) )             
                    
                    self.outfile.append('x%s y%s z%s'% (  pt1[0] , pt1[1], self.ch ) )               
                    self.ngc_to_obj.append( ( pt1[0], pt1[1], self.rh ) ) 

                if not do_retracts or do_laser:
                    self.outfile.append( 'G0' )    
                    #move to first point RH 
                    self.outfile.append('x%s y%s z%s'% (  pt1[0] , pt1[1], self.ch ) )               
                    self.ngc_to_obj.append( ( pt1[0], pt1[1], self.ch ) )             

                
                ##draw the polygons   
                if do_laser:
                    #laser on
                    self.outfile.append( 'M3 S%s'%laserpwm )

                self.outfile.append( 'G1' )
                for i,pt in enumerate(gr_poly):
                    self.outfile.append( 'x%s y%s z%s'%( pt[0], pt[1], self.ch ) )
                    self.ngc_to_obj.append( (pt[0], pt[1], self.ch) )                   
                
 
                if do_laser:
                    self.outfile.append( 'M5 S0' )
                # move to last point at CH  
                #self.ngc_to_obj.append( (gr_poly[0][0], gr_poly[0][1], self.ch))   
                #self.outfile.append( 'x%s y%s z%s'%( (gr_poly[0][0], gr_poly[0][1], self.ch) ) )

                if do_retracts and not do_laser:
                    self.outfile.append( 'G0' )  
                    self.ngc_to_obj.append( (gr_poly[0][0], gr_poly[0][1], self.ch)  )
                    self.ngc_to_obj.append( (gr_poly[0][0], gr_poly[0][1], self.rh)  )
                    self.outfile.append( 'x%s y%s z%s'%( gr_poly[0][0], gr_poly[0][1], self.rh) )

                    #### retract in between cuts
                    self.outfile.append('g0z%s'% ( self.rh ) )  

                self.outfile.append('  ')

        ##-----------------------------------------##        
        # self.outfile.append('(exporting segments )')

        ##-----------------------------------------##
        # rapid move at end 
        self.outfile.append('m2') #program end

    ##-------------------------------##     
    ##-------------------------------##
    if GEOJSON_IS_LOADED:
        #def export_geojson_multiline(self, folder, name ):

        def export_geojson_lines(self, folder, name, periodic=False):
            # export all loaded polygons in gr_sort buffer as a single geojson polyline 
            features = []

            #[[id, centroid, extents, len, points ]]
            for i,s in enumerate(self.gr_sort):

                #make periodic - add the first point to the end 
                if periodic:
                    s[4].append(s[4][0])

                # [[id, centroid, extents, len, points ]]   
                features.append(gjftr(geometry=gjln(coordinates=s[4]), 
                                     properties={"id" : i 
                                                }
                                     ) 
                             )

            feature_collection = gjfc(features)
            filename= '%s/%s'%(folder,name)
            print('## EXPORTING file %s'%filename)

            with open(filename, 'w') as f:
                dump(feature_collection, f)

    ##-------------------------------##
        def export_geojson_polygon(self, folder, name):
            """
            DEBUG - WIP 
            DEBUG - this is mixing shapely objects with GEOJSON objects 
                the entire toolchain needs a rethink about geometry types - uhhhg 

                types are: Point, LineString, Polygon, MultiPoint, MultiLineString, MultiPolygon

            # export all loaded polygons in gr_sort buffer as a geojson polyline 
            """

            features = []

            # >>> Polygon([[(2.38, 57.322), (23.194, -20.28), (-120.43, 19.15), (2.38,   57.322)]])  
            # {"coordinates": [[[2.3..., 57.32...], [23.19..., -20.2...], [-120.4..., 19.1...]]], "type": "Polygon"}

            #[[id, centroid, extents, len, points ]]
            for i,s in enumerate(self.gr_sort):
                # [[id, centroid, extents, len, points ]]   
                #features.append(gjftr(geometry=shp_ply(self.cvt_3d_to_2d(s[4])), 
                features.append(gjftr(geometry=gjply(s[4]), 
                                     properties={"id" : i 
                                                }
                                     ) 
                             )

            feature_collection = gjfc(features)
            filename= '%s/%s'%(folder,name)
            print('## EXPORTING file %s'%filename)

            with open(filename, 'w') as f:
                dump(feature_collection, f)

    ##-------------------------------##
    def export_grid_gfx(self, name, folder , borders=True, centroids=True):
        """
            export cells as graphics 
        """
    
        self.gl_extents()

        if self.sort_minx==0 and self.sort_miny==0 and self.sort_maxx==0 and self.sort_maxy==0:
            raise ValueError('export_grid_gfx - extents are not set') 

        features = []
        for c in self.tesl.nodes:
            if centroids:
                #cell centroids  
                features.append(Feature(geometry=Point((c.coord_x+(c.width/2), c.coord_y+(c.height/2))), 
                                        properties={"id":c.name} 
                                       )
                            )

            if borders:
                #draw a square around the boundary of object 
                features.append(Feature(geometry=LineString(coordinates=c.boundary_pts), 
                              properties={"id" : 0 
                                         }
                              ) 
                      )

            ################
            # render multilines in points (nested arrays - broken up lines)            
            for ptgrp in c.points:
                features.append(Feature(geometry=LineString(coordinates=ptgrp), 
                                  properties={"id" : 0 
                                             }
                                  ) 
                          )

            #################

        feature_collection = FeatureCollection(features)
        with open('%s/%s_cells.json'%(folder, name), 'w') as f:
            dump(feature_collection, f)

    ##-------------------------------##
    def export_sorted_centroids(self, folder, name):
        """ export a centroid for each polygon as a geojson file """

        features = []

        #[[id, centroid, extents, len, points ]]
        for i,s in enumerate(self.gr_sort):
            features.append(Feature(geometry=Point((s[1][0],s[1][1])), properties={"id":i, "len" : len(s[4]) } ))

        feature_collection = FeatureCollection(features)
        with open('%s/%s'%(folder,name), 'w') as f:
            dump(feature_collection, f)

    ##-------------------------------##
    def cvt_grsort_todag(self):
        """ export a centroid for each polygon as a tesselation node (DAG) 
            
            There are more than one ways to calc a centroid: 
                you can always take the middle of a bbox, but it is not always the most accurate 
                if the polygon is square (4 sides) you can draw a line from the mid point of two opposing 
                edges, and intersect those forming a simple "projected" centroid 

                This function can do both, but keep in mind it only works with 4 sided 

        """

        features = []
        zdepth = 0

        #[[id, centroid, extents, len, points ]]
        for i,s in enumerate(self.gr_sort):
            #               [[id, centroid, extents, len, points ]] 
            xtn = s[2] 
            
            width = abs(xtn[2]-xtn[0])
            height = abs(xtn[3]-xtn[1])

            #fancy filtering based on polygons 
            #does it have 4 sides? 5== peridoic 4 pts 
            if(len(s[4])==5): 
                
                ## helper function - easier but limited 
                #self.tesl.new_cell_2d('ply_%s'%i, width,height, i,0,0, s[1][0],s[1][1],zdepth ) 
                
                ## construct the cell manually with more control 
                newc = cell('ply_%s'%i, width, height, i,0,0, s[1][0], s[1][1], s[1][2])
                
                #debug - data may not always be 3D - need to resolve this 
                newc.set_position([s[1][0], s[1][1], s[1][2]])
                
                newc.addattr('idx_x', i)
                #newc.addattr('idx_y', idx_y)
                #newc.addattr('idx_z', idx_z)

                newc.addattr('width', width)
                newc.addattr('height', height)

                ####
                
                # this is center of bbox - you can derive a more accurate center by bisecting the opposing edge midpoints of quad polygons
                newc.addattr('centroid', [s[1][0], s[1][1], s[1][2]])

                ## DEBUG midpoints of opposing edges can also derive an angle per cell  
                e1mid = self.locate_pt_along3d(s[4][0], s[4][1], 1)
                e2mid = self.locate_pt_along3d(s[4][1], s[4][2], 1)
                e3mid = self.locate_pt_along3d(s[4][2], s[4][3], 1)
                e4mid = self.locate_pt_along3d(s[4][3], s[4][0], 1)
                newc.addattr('e1mid', e1mid)
                newc.addattr('e2mid', e2mid)
                newc.addattr('e3mid', e3mid)
                newc.addattr('e4mid', e4mid)

                ####

                self.tesl.add(newc)


        ## loading/saving not implemented yet - need a wrapper for CELLS VS DAGNODES 
        #self.tesl.save_graph_file('%s/%s'%(folder,name))
        #self.tesl.save_tesselation('%s/%s'%(folder,name))

    ##-------------------------------##
    def export_global_extents(self, folder, name, ngc=False):
        """export only the global extents rectangle
           
           can export JSON or NGC 

           serves an an example of a simple gcode exporter 

           EXPORT_RASTER option adds raster lines across extents 

        """ 
         
        EXPORT_RASTER = False 
         
        #[[id, centroid, extents, len, points ]]
        features = []

        do_retracts = True 
        outfile = []

        filename= '%s/%s'%(folder,name)

        mx   = self.sort_minx
        my   = self.sort_miny
        maxx = self.sort_maxx
        maxy = self.sort_maxy
        height = maxy-my 
        width  = maxx-mx 

        numdivs = 100 
        step = height/numdivs

        # export gcode
        if ngc:
            lastpt = (0,0,0)
            
            outfile.append('(exported with export_global_extents )')
            outfile.append('(linear scale set to %s of internal coordinates)'%self.vf_gl_scale )
            outfile.append('  ')
            outfile.append('g20') #inches for unit 
            
            # move to origin  
            outfile.append('g0 x%s y%s z%s f30'% ( self.hp[0], self.hp[1], self.rh) )   #rapid move to 0 
            if do_retracts:
                outfile.append('g0z%s'% ( self.rh ) )  #retract in between cuts
            outfile.append('  ')
            outfile.append('(exporting filled polygons )')

            # no formatting (full precision)
            coords = self.extents_fr_bbox([mx,my,maxx,maxy], periodic=True)
            pt1 = coords[0]

            # first point at retract height 
            if do_retracts:
                outfile.append('x%s y%s z%s'% (  pt1[0] , pt1[1], self.rh ) )               
            
            # iterate points in polygon 
            outfile.append( 'G1' )
            for i,pt in enumerate(coords):
                outfile.append( 'x%s y%s z%s'%( pt[0], pt[1], self.ch ) )
            
            outfile.append( 'G0' )
            if do_retracts:
                outfile.append( 'x%s y%s z%s'%( coords[0][0], coords[0][1], self.rh) )
                outfile.append('g0z%s'% ( self.rh ) )  

            outfile.append('  ')

            print('## EXPORTING file %s'%filename)
            with open(filename, 'w') as f:
                for l in outfile:
                    f.write("%s\n" % l)

        #export json 
        else:  
            #sort extents  
            coords = self.extents_fr_bbox([mx,my,maxx,maxy], periodic=True)
            features.append(Feature(geometry=LineString(coordinates=coords), 
                                    properties={"id" : 0 
                                               }
                                    ) 
                            )
            
            # experiment to make raster lines across the data 
            if EXPORT_RASTER:
                for i in range(1,numdivs):
                    dist = i*step
                    coords = [[mx, my+dist],[maxx,my+dist]]
                    features.append(Feature(geometry=LineString(coordinates=coords), 
                            properties={"id" : i 
                                       }
                            ) 
                    )

            feature_collection = FeatureCollection(features)
            print('## EXPORTING file %s'%filename)
            with open(filename, 'w') as f:
                dump(feature_collection, f)

    ##-------------------------------##
    def export_sorted_extents(self, folder, name):
        """ export the extents of all loaded polygons as a geojson polyline 
        """

        features = []

        #[[id, centroid, extents, len, points ]]
        for i,s in enumerate(self.gr_sort):
            #    features.append(Feature(geometry=LineString(coordinates=self.extents_fr_bbox(s[2])), properties={"id":i, "len" : len(s[4]) } ))
            coords = self.extents_fr_bbox(s[2], periodic=True)

            features.append(Feature(geometry=LineString(coordinates=coords), 
                                 properties={"id" : 0 
                                            }
                                 ) 
                         )

        #sort extents  
        coords = self.extents_fr_bbox([self.sort_minx,self.sort_miny,self.sort_maxx,self.sort_maxy], periodic=True)
        features.append(Feature(geometry=LineString(coordinates=coords), 
                                properties={"id" : 0 
                                           }
                                ) 
                        )
                        
        feature_collection = FeatureCollection(features)
        filename= '%s/%s'%(folder,name)
        print('## EXPORTING file %s'%filename)

        with open(filename, 'w') as f:
            dump(feature_collection, f)

    ##-------------------------------##
    ##-------------------------------##
    def save_line_obj(self, name):
        """ dump a bunch of 3d points as a giant line to visualize in gnolmec 
            
            #do_retracts iterates the gr_polys buffer (not working yet)
            otherwise it iterates the ngc_to_obj buffer 

        """   
        
        #print("exporting %s polygons "%len(self.gr_polys))

        # make a series of connected lines from points 
        self.points = self.clean_pts_str(self.ngc_to_obj)

        for i,pt in enumerate(self.ngc_to_obj):
            if i>0 and i<len(self.ngc_to_obj):
                poly = (i,i+1)
                self.polygons.append( poly ) 

        #self.scale_pts(self.scale)
        self.save(name)

    ##-------------------------------##            
    def cvt_obj3d_grpoly(self, index=None):
        """ 
            convert 3d obj polys into gr_poly buffer to export to json, ect 
            we have to dump the z axis   
        """         
        print("converting %s polygons to grpoly format "%len(self.polygons))

        idx = 0  
    
        for ig, poly in enumerate(self.polygons):
            polygon = []
            for ip, pt in enumerate(poly):
                #polygons.append( (idx, idx+1) ) #2 point polygon indeces
                #print(self.points[pt])
                
                #subtract 1 beacuse OBJ data is not zero indexed 
                ptx = self.points[pt-1][0]
                pty = self.points[pt-1][1] 

                #allow for 2d and 3d data? 
                if len(self.points[pt-1])==3:               
                    ptz = self.points[pt-1][2] 
                    polygon.append((ptx, pty, ptz))
                else: 
                    polygon.append((ptx, pty, 0))

            # [[id, centroid, extents, len, points ]] 
            #self.gr_sort.append([ig, self.centroid((ptx,pty,0.0)) ,[minx,miny, maxx, maxy], len(ply), ply])
            self.gr_polys.append(polygon)
        
        self._sort()
 
    ##-------------------------------##
    def cvt_grpoly_3dobj(self, index=None):
        """DEBUG - sort of works but only with one polygon  """
        print("DEBUG cvt_grpoly_3dobj only works with 1 polygon")
        points = [] 
        pids = [] 

        for ply in self.gr_sort:
            idx=1
            #if index == ply[0]:
            for i,pt in enumerate(ply[4]):
                points.append(pt)
                if i < len(ply[4])-2:
                    pids.append(idx)               
                idx+=1
            
        self.points = points 
        self.polygons.append(pids) 
        
    def cvt_grpoly_obj3d(self, objtype='singlepoly'):
        """ 
            DEBUG - all these functions seem wonkey at best... 
            DEBUG - why not gr_sort? 
            load gr_polys into standard 3d data so we can run ops on it 
           
            ARGS:
                objtype = 'polyline' or 'singlepoly'

        """
        
        print("converting %s polygons "%len(self.gr_polys))

        idx = 0  
        
        for ig, grpoly in enumerate(self.gr_polys):
            polygons = []
            for ip, pt in enumerate(grpoly):

                idx = ip+len(self.points)
                
                if objtype == 'polyline': 
                    self.polygons.append([idx+1, idx+2])

                if objtype == 'singlepoly':
                    polygons.append( idx+1 )  

            if objtype == 'singlepoly':
                self.polygons.append(polygons)
            
            self.points.extend(self.clean_pts_str(grpoly))

        #print(self.rotation)
        #self.show() 

    ##-------------------------------##
    ##-------------------------------##
    ##-------------------------------------------##   
    def load_geojson(self, inputfile, getfids=None, zaxis=0):
        """ DEBUG - kind of works but we need to deal with multiple types of GEOM 
                    This was built with lines in mind originally 

            parse a geojson file - store points in arrays in GR_POLY buffer 
            
            Args:
            
                inputfile     - the file to load 
                getfids       - feature ids to load 
                                 can be single int or list of ints 
            
                zaxis         - value to set false zaxis to (2d to 3d) 

        """
        print('loading geojson file %s'%inputfile)

        plyidx = 1
        ptidx = 1 

        geojson_txt = None

        with open(inputfile) as f:
            gj = gjload(f)
        features = gj['features'] 
 
        fixedids = [] 

        ##---------------
        ##---------------
        #check that ids are in range, int type and in a list 
        if getfids:
            if type(getfids) == int:
                fixedids = [int(getfids)]
            else:
                for id in getfids:
                    if int(id)>=len(features):
                        print("id out of range %s"%id)
                        pass
                    else:
                        fixedids.append(int(id))

        ##---------------
        ##---------------
        # import all data ... 
        if getfids is None:
            for i,f in enumerate(features):
                # multipolygon type
                if f.geometry:
                    ##DEBUG ADDING THIS DATA TYPE 
                    if f.geometry.type == 'LineString':
                        tmp_poly = []
                        for coord in f.geometry.coordinates: 
                            #tmp_poly.append(coord)
                            tmp_poly.append( (coord[0], coord[1], zaxis) )

                        if tmp_poly:
                            #print(tmp_poly)
                            self.jsonbuffer.append(tmp_poly) 

                    if f.geometry.type == 'MultiPolygon':
                        for coord in f.geometry.coordinates:                    
                            for p,koord in enumerate(coord):
                                ptidx = 1
                                tmp_poly = [] 
                                for pt in koord:
                                    if type(pt[0])==float and type(pt[1])==float:
                                        tmp_poly.append( (pt[0], pt[1], zaxis) )
                                        ptidx += 1  
                                if tmp_poly:
                                    self.jsonbuffer.append(tmp_poly) 

                    # polygon type   
                    if f.geometry.type == 'Polygon':
                        for coord in f.geometry.coordinates:
                            ptidx = 1
                            tmp_poly = [] 
                            
                            print(coord)

                            for pt in coord:
                                if type(pt[0])==float and type(pt[1])==float:
                                    tmp_poly.append( (pt[0], pt[1], zaxis) )
                                    ptidx += 1  

                            #print("loaded %s points in polygon "%len(tmp_poly)) 
                            if tmp_poly:
                                self.jsonbuffer.append(tmp_poly) 
                    plyidx += 1 

        ##---------------        
        # or just cherry pick some features to import (below we also can pick the sub-features) 
        if getfids:
            for fid in fixedids:
                for i,f in enumerate(features):
                    #multipolygon - DEBUG - not fully tested yet 
                    if f.geometry.type == 'MultiPolygon':
                        if fid==i:                        
                            for coord in f.geometry.coordinates:                    
                                for p,koord in enumerate(coord):
                                    ptidx = 1
                                    tmp_poly = [] 
                                    for pt in koord:
                                        if type(pt[0])==float and type(pt[1])==float:
                                            tmp_poly.append( (pt[0], pt[1], zaxis) )
                                            ptidx += 1  
                                    if tmp_poly:
                                        self.jsonbuffer.append(tmp_poly) 
                    #polygon type
                    if f.geometry.type == 'Polygon':       
                        if fid==i:
                            for coord in f.geometry.coordinates:
                                ptidx = 1
                                tmp_poly = [] 
                                for pt in coord:
                                    if type(pt[0])==float and type(pt[1])==float:
                                        tmp_poly.append( (pt[0], pt[1], zaxis) )
                                        ptidx += 1  
                                #print("loaded %s points in polygon "%len(tmp_poly)) 
                                if tmp_poly:
                                    self.jsonbuffer.append(tmp_poly) 
                            plyidx += 1 

        ##---------------
        self.gr_polys = self.jsonbuffer 

        ##---------------
        #print('cloning data, sorting and calculating extents ')
        
        #make a copy of the data to work with - leaving the original in case we want to do something with it  
        self._sort()
        
        print("loaded %s polygons from %s "%(plyidx, inputfile)) 

    ##-------------------------------##
    def export_ngc(self, rh, ch, cdpi, cmax, filename, do3d=False, do_retracts=True, do_laser=False, laserpwm=400):
        """ rh         - retract height  
            ch         - cut height 
            cdpi       - cut depth per iteration 
            cmax       - cut depth max 
            filename 
            do3d=False
        """

        if ch==None:
            print("# exporting NGC file - cut height disabled (3d) ", filename)
        else:
            print("# exporting NGC file ", filename)
            self.ch = ch          # cut height (top, start of cut)
        
        self.rh = rh          # retract height 
        self.cdpi = cdpi      # cut depth per iteration on Z axis
        self.cmax = cmax      # maximum cut depth on Z axis 

        if do3d==True:
            print("## WARNING ## 3D export is not working yet")
            self._calculate_paths3d(do_laser=do_laser, do_retracts=do_retracts)
        else:
            self._calculate_paths2d(do_laser=do_laser, do_retracts=do_retracts,laserpwm=laserpwm)
        
        print("gr_sort buffer has %s polys in it. "%(len(self.gr_sort)))
        fobj = open( filename,"w") #encoding='utf-8'
        for line in self.outfile: 
            fobj.write(line+'\n')
        fobj.close()

    ##-------------------------------##
    ##-------------------------------##
    def make_grid(self, distance, folder, xcuts, ycuts, bbox=None):
        """ chop a square into smaller squares 
        """

        if bbox:
            self.tesl._set_extents(bbox) 
        else:
            self.tesl._set_extents([self.sort_minx, self.sort_miny, self.sort_maxx, self.sort_maxy]) 

        self.tesl.build_2d_cells(xcuts, ycuts, distance)
        

        #temporary export of geom I HAVE TO SEE THIS BEFORE I GO TO BED!
        features = []

        ##---
        # add attrs to derived DAG graph nodes 
        for c in self.tesl.nodes:
            cen = (c.coord_x+(c.width/2), c.coord_y+(c.height/2))
            c.addattr('centroid', cen )
            c.addattr('width' , c.width )
            c.addattr('height', c.height )

            if len(self.gr_sort):
                for sort in self.gr_sort:
                    # [[id, centroid, extents, len, points ]]
                    dist = self.mu.calc_line_length(cen[0], cen[1], sort[1][0], sort[1][1])
                    if float(dist) < float(distance):
                        
                        #DEBUG TODO - there is a problem of the same polygon getting added more than once 
                        #make sure we dont export multiple times 
                        # (check the number of points, then check thr actual points) 

                        tmp_pts=[cen, sort[1]]
                        features.append(Feature(geometry=LineString(coordinates=tmp_pts)) 
                                       )

        feature_collection = FeatureCollection(features)
        with open('%s/_distances.json'%(folder), 'w') as f:
            dump(feature_collection, f)

    ##-------------------------------##
    def mc_escher(self):
        """ DAG + point generator = tesselator 
            make a grid and draw points at each location 
        """
        for cell in self.tesl.nodes:
            cen = cell.getattrib('centroid')
            #print('# ', cell.name,' ', cen )
            
            pts = self.pop2d.calc_circle_2d(cen[0],cen[1], .75, periodic=True, spokes=3)
            cell.points.append( pts )

            pts2 = self.pop2d.batch_rotate_pts_2d( pts, cen, 180 )
            cell.points.append( pts2 )

            pts3 = self.pop2d.calc_circle_2d(cen[0],cen[1], 1.5, periodic=True, spokes=12)
            cell.points.append( pts3 )

    ##-------------------------------##
    def vec_render_obj(self, rx,ry,rz,  renderscale, folder, objfile):
        single_line = False

        framesize = self.machine_size 
        cen = self.machine_center

        ropr = simple_render() 
        obj = object3d()
        
        if type(objfile)==str:
            obj.load('%s/%s'%(folder,objfile))
        if type(objfile)==object3d:
            obj.polygons=objfile.polygons
            obj.points=objfile.points    


        #no rotation 
        #ropr.render_obj((100,0,255), 0, 0, 0, 1, renderscale, object3d=obj)
        
        #with rotation
        ropr.render_obj((100,0,255), rx, ry, rz, 1, renderscale, object3d=obj)

        #coords are in pixels - rather huge for a model 
        #ropr.vec_fr_buffer.scale_pts((.01,.01,.01))
    
        ropr.vec_fr_buffer.move_center()
      
        # used to test the vector render engine 
        #ropr.vec_fr_buffer.save('%s/test_%s.obj'%(folder,i) )

        out_pts = [] 

        # draw as a single line without breaks       
        if single_line:
            pts = [] 
            for pt in ropr.vec_fr_buffer.points:
                pts.append( (pt[0]+cen[0], pt[1]+cen[1] ) )  
            out_pts.append( pts )

        # draw as proper line segments         
        else:
            for ply in ropr.vec_fr_buffer.polygons:
                pts = []
                #in this setup it will only be two point polys - vector render only does that                    
                if len(ply)==2:
                    # pts are zero indexed hence the -1
                    pt1 = ropr.vec_fr_buffer.points[ply[0]-1] 
                    pt2 = ropr.vec_fr_buffer.points[ply[1]-1] 
                    pts.append( (pt1[0]+cen[0], pt1[1]+cen[1]) )
                    pts.append( (pt2[0]+cen[0], pt2[1]+cen[1]) )
                out_pts.append( pts )

        return out_pts 
        
    ##-------------------------------##
    def tess_vec_render(self, renderscale, path, objfile):
        """ make a grid and invoke render at each location 
            
            DEBUG - add mesh optimizer 
                - eliminate lines that overlap 
                - only render faces that face you (backface culling) 
                - eliminate lines too small/big 
        """

        # render as one line (with no breaks)
        single_line = False 
        
        ropr = simple_render() 
        obj = object3d()

        if type(objfile)==str:
            obj.load('%s/%s'%(path,objfile))
        if type(objfile)==object3d:
            obj.polygons=objfile.polygons
            obj.points=objfile.points    

        #obj.show()

        total = len(self.tesl.nodes)

        for i,cell in enumerate(self.tesl.nodes):
            cen = cell.getattrib('centroid')
            print('#rendering  %s@%s %s of %s '%(cell.name, cen, i,total) )

            #clear cache each time 
            ropr.vec_fr_buffer.points = [] 
            ropr.vec_fr_buffer.polygons = [] 

            ## color, rx, ry, rz, thick, scale 
            #ropr.render_obj((100,0,255), i*20, i*20, i*20,  1, renderscale, object3d=obj)
            
            #no rotation 
            ropr.render_obj((100,0,255), 0, 0, 0, 1, renderscale, object3d=obj)
            #with rotation
            #ropr.render_obj((100,0,255), i*30, -i*30, i*30, 1, renderscale, object3d=obj)

            #coords are in pixels - rather huge for a model 
            #ropr.vec_fr_buffer.scale_pts((.01,.01,.01))
        
            ropr.vec_fr_buffer.move_center()
          
            # used to test the vector render engine 
            #ropr.vec_fr_buffer.save('%s/test_%s.obj'%(path,i) )
    

            # draw as a single line without breaks       
            if single_line:
                pts = [] 
                for pt in ropr.vec_fr_buffer.points:
                    pts.append( (pt[0]+cen[0], pt[1]+cen[1] ) )  
                cell.points.append( pts )

            # draw as proper line segments         
            else:
                for ply in ropr.vec_fr_buffer.polygons:
                    pts = []
                    #in this setup it will only be two point polys - vector render only does that                    
                    if len(ply)==2:
                        # pts are zero indexed hence the -1
                        pt1 = ropr.vec_fr_buffer.points[ply[0]-1] 
                        pt2 = ropr.vec_fr_buffer.points[ply[1]-1] 
                        pts.append( (pt1[0]+cen[0], pt1[1]+cen[1]) )
                        pts.append( (pt2[0]+cen[0], pt2[1]+cen[1]) )
                    cell.points.append( pts )

    ##-------------------------------##
    def tess_objclone(self, pts=None, objfile=None):
        """ make a grid and insert points from an object into each tesselation cell  
            make sure to center the object for best results
            ( use object3D.move_center()  )

        """
        if objfile:
            if type(objfile)==str:
                self.pop2d.load(objfile)
            if type(objfile)==object3d: 
                self.pop2d.points = objfile.points        
                self.pop2d.polygons = objfile.polygons
        if pts:
            self.pop2d.points = pts 
            
        if not pts and not objfile:
            print("tess_objclone - no object or point data to work with ")
            return None 

        #self.pop2d.show()

        for cell in self.tesl.nodes:
            cen = cell.getattrib('centroid')
            #print('# ', cell.name,' ', cen )
        
            #DEBUG move goes into a broken loop if you try to use it - debug  
            #self.pop2d.move( cen[0],cen[1] )

            pts = []
            
            #append X Y coords to each cell (ignore Z)
            for pt in self.pop2d.points:
                pts.append( (pt[0]+cen[0], pt[1]+cen[1] ) )  

            cell.points.extend( pts )





##------------------------------##
##------------------------------##            
##------------------------------##




def arc_to_degree(NS, degrees, minutes, seconds, EW):
    """ arc minutes to decimal degree ( example n50d0'02"e ) """
    
    outdegrees = 0.0

    if NS =='n':
        outdegrees = degrees
        outdegrees = outdegrees + (minutes*.0166667) #1/60
        outdegrees = outdegrees + (seconds*.0166667*.0166667) #1/60
    if NS =='s':
        outdegrees = 180.0
        outdegrees = outdegrees + degrees
        outdegrees = outdegrees + (minutes*.0166667) #1/60
        outdegrees = outdegrees + (seconds*.0166667*.0166667) #1/60
    if EW =='w' and NS =='s':
        outdegrees = outdegrees * -1
    if EW =='e' and NS =='n':
        outdegrees = outdegrees * -1
  
    return outdegrees

