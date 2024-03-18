

from gnelscript.pygfx.obj3d import  *




def merge_numeric_pairs( path , infile, heights):

    for hgt in heights:

        cop = cam_op()
        normal = vec3(0,0,1)
        origin = vec3(0,0,hgt)
        poly = cop.tm_meshplane_test(origin, normal, 5, path , infile, axis='y')
         
        ########## 
        obj = object3d()
        pts   = [] 
        faces = []
        for i,pt in enumerate(poly): 
            pts.append(tuple(pt[0]))
            pts.append(tuple(pt[1]))
        for f in range(0,len(pts),2):
            faces.append([f+1,f+2])

        print('#hgt %s sorted %s pts   '%(hgt, len(pts)   ) )
        print('#hgt %s face has %s ids '%(hgt, len(faces) ) )

        obj.points = pts 
        obj.polygons = faces    
        obj.save('output_%s.obj'%hgt) 

        #convert 3d points to a dict() 
        d = dict()
        for i,pt in enumerate(poly):
            d[i] = [tuple(pt[0]),tuple(pt[1])]

        #run the sort 
        pathids = cop.create_path( len(d), 0, d[0], d )

        sortedpts = []

        #build a sorted list 
        ply =[]
        for i,x in enumerate(pathids):
            sortedpts.append( d[x][0])
            ply.append(i+1)

        obj2 = object3d() 
        obj2.points = sortedpts 
        obj2.polygons = [ply]    
        obj2.save('sorted_%s.obj'%hgt)




##------------------------------------##
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


#-------------------------------

class gear_mesh_generator(object3d):
    def __init__(self):
        self.gear1 = gear_generator()
        self.gear2 = gear_generator()