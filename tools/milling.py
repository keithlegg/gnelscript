

from gnelscript.pygfx.obj3d import  *




class gear_generator(object3d):
    def __init__(self):
        super().__init__()         
        #self.contact_angle
        self.tooth_spacing = 20 
        self.num_teeth = 5
        self.shaft_hole_dia = 4


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