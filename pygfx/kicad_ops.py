# ***** BEGIN GPL LICENSE BLOCK *****
#
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# ***** END GPL LICENSE BLOCK *****

import re 



class pcbfile():
    """ 
       A bit of history:
           the first tool I built 
           started with gcode export and kicad import 
           made this to load kicad and export gcode 
           cloned and renamed vectorflow() 
           cloned vectorflow into milling_op() 

    """

    def __init__(self):
        super().__init__()  

        #self.points        = []  NOPE! - conflicts with 3D OBJECT  
        #self.polygons      = []  NOPE! - conflicts with 3D OBJECT

        self.parse_depth = 0
        self.file_contents = []
        
        self.oneup_entity = None 
        self.cur_entity = None 
        self.cur_module = None 
        self.cur_module_pos = None 

        self.modules = []

        self.known_entities = ['module','gr_line','segment']

        self.kicad_units   = 'mils' #inch to mm -> scale 0.0393701
        self.object_units  = 'cm'

        self.outfile = []
        
        # GEOM ARRAYS for export   
        self.ngc_to_obj    = []

        self.line_segments = []  #list of list of points 
        self.ki_polygons   = []  #list of list of points 
        self.filled_polys  = []  #list of list of points 
        self.gr_polys      = []  #list of list of points 
        self.modules       = []  #list of list of [name, [(),()] ] 

        self.gr_polygon_buffer   = [] # buffer to dump into array of arrays (gr)
        self.fill_polygon_buffer = [] # buffer to dump into array of arrays (filled_polys)
        self.polygon_buffer      = [] # buffer to dump into array of arrays (modules)

        self.segment_buffer = [] # buffer to dump into array of arrays (line_segments)
        self.module_buffer  = [] # buffer to dump into array of arrays (modules)
 
        self.parsing_fillpolygon = False # stored state of parser object type
        self.parsing_grpolygon   = False # stored state of parser object type 
        self.parsing_polygon     = False # stored state of parser object type 

        self.module_depth    = 0     # stored state of parser object type
        #self.parsing_segment = False # stored state of parser object type

        self.rh = 1.0     # retract height 
        self.ch = .5      # cut height 
        self.hp = (0,0,0) # home position 

        self.scale =  -0.0393701 #mm to inch

    ##-------------------------------##
    def show_geom(self):
        #call the inherited polygon "show" method
        self.show()

    ##-------------------------------##
    def scrub(self, inp):
        """ clean up parsed characters from kicad file """

        out = inp.lower()
        out = out.strip()
        out = out.replace('(','')
        out = out.replace(')','')
        out = out.replace(' ','')
        out = out.replace('\"','')
        out = out.replace('\'','')
        return out

    ##-------------------------------##
    def fl_scrub(self, inp):
        out = self.scrub(inp)
        return float(out)

    ##-------------------------------##
    def bufferinfo(self):
        print('buffer sizes ')
        print('graphic polygons %s '%len(self.gr_polys     ) )
        print('filled  polygons %s '%len(self.filled_polys ) )
        print('        polygons %s '%len(self.ki_polygons     ) )                        

    ##-------------------------------##
    def showbuffers(self):
        print('graphic polygons %s '%self.gr_polys     )
        print('filled  polygons %s '%self.filled_polys )
        print('        polygons %s '%self.ki_polygons     ) 

    ##-------------------------------##
    def show_modules(self):
        for m in self.modules:
            print(m.name)
    
    ##-------------------------------##
    def save_3d_obj(self, name):
        """ dump a bunch of 3d poonts as a giant line to visualize in gnolmec"""   
        self.points = self.ngc_to_obj

        # make a series of connected lines from points  
        for i,pt in enumerate(self.ngc_to_obj):
            if i>0 and i<len(self.ngc_to_obj):
                poly = (i,i+1)
                self.polygons.append( poly ) 

        
        #self.scale_pts(self.scale)
        self.save(name)

    ##-------------------------------##
    def read_pcb(self,infile):
        #(gr_line (start 99.695 105.41) (end 119.38 149.86) (layer Dwgs.User) (width 0.15) (tstamp 6068B63B))
        #(gr_line (start 160.02 134.62) (end 99.695 105.41) (layer Dwgs.User) (width 0.15))
        #(gr_line (start 174.625 131.445) (end 160.02 134.62) (layer Dwgs.User) (width 0.15))
        #(gr_line (start 107.95 97.155) (end 174.625 131.445) (layer Dwgs.User) (width 0.15))
        #(gr_line (start 126.365 17.78) (end 107.95 97.155) (layer Dwgs.User) (width 0.15))
        #(gr_line (start 114.3 24.765) (end 126.365 17.78) (layer Dwgs.User) (width 0.15))
        #(gr_line (start 97.155 99.06) (end 114.3 24.765) (layer Dwgs.User) (width 0.15))
        #(gr_line (start 73.66 38.735) (end 97.155 99.06) (layer Dwgs.User) (width 0.15))
        #(gr_line (start 61.595 43.815) (end 73.66 38.735) (layer Dwgs.User) (width 0.15))
        #(gr_line (start 97.155 152.4) (end 61.595 43.815) (layer Dwgs.User) (width 0.15))
        #(gr_line (start 119.38 149.86) (end 97.155 152.4) (layer Dwgs.User) (width 0.15))
        #(via (at 139.7 76.835) (size 0.8) (drill 0.4) (layers F.Cu B.Cu) (net 0))
        #(via (at 153.035 88.9) (size 0.8) (drill 0.4) (layers F.Cu B.Cu) (net 0))
        #(via (at 149.225 99.695) (size 0.8) (drill 0.4) (layers F.Cu B.Cu) (net 0))
        pass

    ##-------------------------------##

    ##-------------------------------##
    def load_kicadpcb(self, zval, path, infname):
        """ DEBUG 

            This is just a simple test for loading a SINGLE polygon  
            needs tons of work 

            

            to build this properly we need to use recursion and sort all eleements in nested parens into blocks of text 


            This is not recursive, but clever enough to scan all the 
            entities in the file that have coordinates and store them in a list with a type identifier 

            all data needs to be 3d internally - add a zero value for Z where applicable 

        """ 

        var_module_name    = ''
        var_module_pos     = []
        var_module_lines   = []

        var_pad_name       = ''
        var_pad_size       = 0
        var_pad_xcoord     = 0
        var_pad_ycoord     = 0

        var_segment_start  = []
        var_segment_end    = []

        var_line_width     = 0
        var_line_start_xy  = []
        var_line_end_xy    = []



        f = open('%s/%s'%(path,infname), 'r')

        geometry = []

        for lc,line in enumerate(f):
            if '(' in line or ')' in line:

                newpoly =[] #[type, [coords] ]

                linetoked = line.split(' ') 
                cleanspaces = []
                for t in linetoked:
                    if t!='':
                        cleanspaces.append( self.scrub(t) )

                firstelem = self.scrub(cleanspaces[0]) 
                if firstelem!='':
                    newpoly.append(firstelem)

                for i,tok in enumerate(cleanspaces):
                    if tok !='':
                        if tok == 'xy':   
                            newpoly.append( (self.fl_scrub(cleanspaces[i+1])  , self.fl_scrub(cleanspaces[i+2]) ))  

                        if tok == 'start':   
                            newpoly.append( (self.fl_scrub(cleanspaces[i+1])  , self.fl_scrub(cleanspaces[i+2]) ))  
                        if tok == 'end':   
                            newpoly.append( (self.fl_scrub(cleanspaces[i+1])  , self.fl_scrub(cleanspaces[i+2]) ))  

                        if tok == 'center':   
                            newpoly.append( (self.fl_scrub(cleanspaces[i+1])  , self.fl_scrub(cleanspaces[i+2]) ))  

                if len(newpoly)==2:  
                    geometry.append(newpoly)

        ####################
        #DEBUG - only works with one poly right now, it will miss things that are longer than one line!
        
        print('######## num geometry elements read ', len(geometry ) )

        single_poly = [] 

        for p in geometry:
            if p[0]=='xy':
                pt = (p[1][0], p[1][1], zval )
                single_poly.append(pt)

        if single_poly:
            self.gr_polys.append(single_poly)



















