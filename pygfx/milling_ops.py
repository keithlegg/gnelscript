#!/usr/local/bin/python3



import math 
import re 

#from pygfx.point_ops import polygon_operator
#from pygfx.math_ops import vec3 

import os 

from gnelscript.pygfx.obj3d import object3d


#from pygfx.gcode.bridgeport import  parser_commands
from gnelscript.pygfx.gcode.linuxcnc import  parser_commands



"""
A  A axis of machine
B  B axis of machine
C  C axis of machine
D  Tool radius compensation number
F  Feed rate
G  General function (See table Modal Groups)
H  Tool length offset index
I  X offset for arcs and G87 canned cycles
J  Y offset for arcs and G87 canned cycles
K  Z offset for arcs and G87 canned cycles.
   Spindle-Motion Ratio for G33 synchronized movements.
L  generic parameter word for G10, M66 and others
M  Miscellaneous function (See table Modal Groups)
N  Line number
P  Dwell time in canned cycles and with G4.
   Key used with G10.
Q  Feed increment in G73, G83 canned cycles
R  Arc radius or canned cycle plane
S  Spindle speed
T  Tool selection
U  U axis of machine
V  V axis of machine
W  W axis of machine
X  X axis of machine
Y  Y axis of machine
Z  Z axis of machine

"""

# COMMENT = ';' # bridgeport uses this 

COMMENT = '('    # linuxcnc uses this 
PARAM   = '#<'   # parameter (variable) , followed with brackets 

#<xscale> = 1.0
#<yscale> = 1.0
#<zscale> = 1.0
#<fscale> = 10000.0
#<toolno> = 1
#<rpm>    = 1600

##------------------------------------------

## GCODE KNOWN DIALECTS 



#  X-Axis (table)  
#  Y-Axis (saddle)  
#  Z-Axis (quill)  


##------------------------------------------



##------------------------------------------

class generate_gcode(object):
    def __init__(self):
        pass
        

##------------------------------------------


class gcode(object3d):

    def __init__(self):
        super().__init__()  

        self.DEBUG_MODE = True
        # self.linear_units = 'in'
        # self.z_axis  ?? 

        #comments get dumped into this array (full or partial line)
        self.commented = [] # [index, string] 
        
        self.param_names  = [] # [name, value ]
        self.param_values = [] # [name, value ]

        #swappable dialects of gcode commands 
        self.dialect = parser_commands

        self.coord_words = ['N','G', 'X','Y','Z','U','V','W','I','J','K','R','P','F'] # ,'A','B','C','D']
        
        # simulated position of the cutting head
        self.POSX = 0
        self.POSY = 0
        self.POSZ = 0

        self.segments  = [] # [index, command, xyz_pos] 




    def show_data(self):
        print("## ## ## ## ## ## ## ## ## ##")

        for s in self.segments:
            print(s[1])

    def save_3d_object(self, filename):
        """ dump the XYZ positions into a polyline object """

        obj   = object3d() # container for 3D object 
        data = []

        line_data = []
        for pt in self.segments:
            line_data.append( pt[1] )

        obj.pts_to_linesegment(line_data)
        obj.save( filename )

    def contains_coord_words(self, string): 
        # check for any known coordinate words
        has_coords = False 
        for cw in self.coord_words:
            if cw in string:
                has_coords = True
        return  has_coords       

    def which_coord_words(self, string):
        words = []
        for cw2 in self.coord_words:
            if cw2 in string:
                words.append(cw2)
        return words

    def parse_params(self, string):
        """ discover and retain any parameters in the file """
        if PARAM in string:
            tmp = self.between_token_list(string, [PARAM, '>']) 

            if len(tmp)==2:
                if tmp[0] not in self.param_names:
                    if '=' in tmp[1]:
                        tmp2 = tmp[1].split('=')
                        self.param_names.append(tmp[0])
                        self.param_values.append(tmp2[1].replace(' ',''))
                else:
                    #IT IS ALREADY PARSED!
                    #deal with this here or in another function??
                    pass

    def between_token_list(self, string, tokens, return_match=False):
        """ give a list of tokens, return a list of betweens 
            
            string = 'y987a123b541c307d999'
            tokens = ['a','b','c','d']
            out = self.between_token_list( string, tokens)

            out will be:
                 ['123', '541', '307']

            unless return_match is True,
                then out will be:
                    [ ['a',123'], ['b',541'], ['c',307'] ]

        """

        #return re.split( tokens_str ,string)
        match = ''
        if len(tokens)==1:
            match = tokens[0] + '|'
        if len(tokens)>1:    
            for i,t in enumerate(tokens):
                if i<len(tokens)-1:
                    match = match + t + '|'
                if i==len(tokens)-1:
                    match = match + t 
        
        if return_match is False:
            tmp = re.split( match, string)
        else:
            match = '('+match+')'
            tmp = re.split( match, string)

        out = []
        for t in tmp:
            tmp2 = t.replace('\n','')
            tmp2 = tmp2.replace(';','')
            out.append(tmp2)

        if return_match is False:
            return out[1:]  
        else:
            out2 = [] 
            for i,tok in enumerate(out[1:]):
                if i%2==0:
                    out2.append( [tok, out[1:][i+1]] )
            return out2


    def save_gcode(self, filename):
        f = open( filename,"w", encoding='utf-8')
        
        #for seg in self.segments:
        #    #f.write(seg[])
        #    print( seg[0], seg[2])

    def load_gcode(self, filename):

        #if os.path.lexists(filename) == 0:
        #    self.scribe("%s DOES NOT EXIST !! "%filename )
        #    #raise
            
        if os.path.lexists(filename):
            f = open( filename,"r", encoding='utf-8')
            contents = f.readlines()

 
            #scan entire file for parametrs first 
            for lin in contents:
                self.parse_params(lin)

            print("#### params ", self.param_names )

            cleaned_contents = []

            for lin in contents:

                #if self.DEBUG_MODE:
                #print("#### LINE IS ", lin.replace("\n","") )

                # ignore commented out lines 
                if lin[0]==COMMENT or lin[0]=='\'':  
                    self.commented.append(lin)

                else:


                    # seperate the line index out and split from rest of command
                    lindex = self.between_token_list(lin, self.coord_words, return_match=True)
                    #n_idx = lindex[0]   
                    #comm = lin[len(str(n_idx)):] 
                    if lindex:
                        #print('############### ', lindex)

                        # execute the parameters(variables) if they exist 
                        for i,tok in enumerate(lindex):
                            for j,p in enumerate(self.param_names):
                                if p in tok[1]:
                                    # substitute the value with the name 
                                    tok[1] = (tok[1].replace('#<%s>'%self.param_names[j], self.param_values[j]) )
                                    try:
                                        # DONT LOAD ANY GCODE TEXTFILES WITH EXECUTABLE PYTHON IN THEM!
                                        # THIS IS A SECURITY HOLE FOR WANKERS TO EXPLOIT 
                                        
                                        # eval the damn thing and capture the result 
                                        #print( tok[0], eval(tok[1])[0])
                                        pass
                                        
                                    except:
                                        pass    

                        print( lindex )

                        ## #check for comments on the line but not at start   
                        ## if COMMENT in lin:
                        ##     tmp = lin.split(COMMENT)
                        ##     lin = tmp[0]
                        ##     if tmp[1] is not '\n':
                        ##         self.commented.append(tmp[1])
                        ## # check for known commands 
                        ## for key in self.dialect :
                        ##     if key in comm:
                        ##         print(" COM FOUND at INDEX %s ! "%n_idx , key, '---', self.dialect [key] )
                        ## # check for any known coordinate words
                        ## has_coords = self.contains_coord_words(comm) 
                        ## # if has coordinate words, determine which ones it has 
                        ## if has_coords is True:             
                        ##     line_contains = ( self.which_coord_words(comm) )
                        ##     output        = self.between_token_list(comm, line_contains)
                        ##     #print( line_contains, output )
                        ##     for i,token in enumerate(line_contains):
                        ##         if token=='X':
                        ##             self.POSX = float(output[i])
                        ##         if token=='Y':
                        ##             self.POSY = float(output[i])
                        ##         if token=='Z':
                        ##             self.POSZ = float(output[i])

                        #self.segments.append( [n_idx, [self.POSX, self.POSY, self.POSZ], comm ] )


"""
#!/usr/bin/env python
#----------------------------------------------------------------------------
# 29-Jun-2015 ShaneG
#
# Convert a collection of SVG paths into a Gcode file.
#----------------------------------------------------------------------------
from lxml import etree
from os.path import splitext, exists
from svg.path import Path, Line, parse_path
from optparse import OptionParser
from gcode import saveGCode

#--- Usage information
USAGE ='' 
# Usage:
#        %s [--cut depth] [--safe depth] filename
# 
# Where:
# 
#   --cut    depth    specify the target cut depth
#   --safe   depth    specify the safe distance for moves between cuts


# GCODE_PREFIX = ""
# G21 (Use mm)
# G90 (Set Absolute Coordinates)
# G0 Z%0.4f (Move clear of the board first)
# (end of prefix)
# 

GCODE_SUFFIX = ""
# (Operation complete)
# G0 X0 Y0 Z%0.4f
# M2
# %

def styleContains(style, values):
  #Determine if the given properties are set in a style definition.
  # First break up the style into properties
  props = dict()
  for entry in style.split(";"):
    k, v = [ e.strip() for e in entry.split(":") ]
    props[k] = v
  # Now check for a match
  for k in values.keys():
    if props.get(k.strip(), "").strip() <> values[k].strip():
      return False
  return True

def getScaling(doc):
  # Calculate the scaling from document points to mm
  wmm = float(doc.get("width")[:-2])
  hmm = float(doc.get("height")[:-2])
  vbox = list([ float(n) for n in doc.get("viewBox").split(" ") ])
  return wmm / vbox[2], hmm / vbox[3], vbox[3]

def getXY(point, sx, sy, h):
  return round(point.real * sx, 4), round((h - point.imag) * sy, 4)

def processPath(path, cut, safe, sx, sy, h, precision):
  # Process a single path, generating the required Gcode to implement it.
  gcode = list()
  tool = False
  lx, ly = None, None
  for segment in path:
    length = segment.length()
    if length == 0.0:
      continue
    gcode.append("(%s,%0.4fmm)" % (segment.__class__.__name__, length * sx))
    # Move the tool head if there is a break in the path
    x, y = getXY(segment.start, sx, sy, h)
    if ((x <> lx) or (y <> ly)):
      if tool:
        gcode.append("G0 Z%0.4f" % safe)
        tool = False
      gcode.append("G0 X%0.4f Y%0.4f" % (x, y))
    # Make sure the tool is down
    if not tool:
      gcode.append("G1 Z%0.4f" % cut)
      tool = True
    # Determine what sort of object we are processing
    consumed = False
    if segment.__class__ == Line:
      # Do the move to the end point
      x, y = getXY(segment.end, sx, sy, h)
      gcode.append("G1 X%0.4f Y%0.4f" % (x, y))
      lx, ly = x, y
      consumed = True
# TODO: Arcs can be done with G2/G3 commands
#    if segment.__class__ == Arc:
#      # Circles can be done with G2
#      if segment.radius.real == segment.radius.imag:
#        x, y = getXY(segment.end, sx, sy, h)
    # Fallback, interpolate as a sequence of straight lines
    if not consumed:
      # Assumes sx == sy so we get 1mm steps
      delta = precision / (sx * length)
      pos = delta
      while pos < 1.0:
        x, y = getXY(segment.point(pos), sx, sy, h)
        gcode.append("G1 X%0.4f Y%0.4f" % (x, y))
        pos = pos + delta
      x, y = getXY(segment.end, sx, sy, h)
      gcode.append("G1 X%0.4f Y%0.4f" % (x, y))
      lx, ly = x, y
      consumed = True
  # Bring the tool up again if needed
  if tool:
    gcode.append("G0 Z%0.4f" % safe)
    tool = False
  return gcode

#--- Main program
if __name__ == "__main__":
  # Set up program options
  parser = OptionParser()
  parser.add_option("-c", "--cut", action="store", type="float", dest="cut_depth")
  parser.add_option("-s", "--safe", action="store", type="float", dest="safe_depth")
  parser.add_option("-p", "--precision", action="store", type="float", dest="precision", default=1.0)
  options, args = parser.parse_args()
  # Check positional arguments
  if len(args) <> 1:
    print USAGE.strip() % argv[0]
    exit(1)
  # Make sure required arguments are present
  for req in ("cut_depth", "safe_depth"):
    if eval("options.%s" % req) is None:
      print "ERROR: Missing required argument '%s'" % req
      print USAGE.strip() % argv[0]
      exit(1)
  source = args[0]
  name, ext = splitext(source)
  if ext == "":
    source = name + ".svg"
  if not exists(source):
    print "ERROR: Could not find source file '%s'" % source
    exit(1)
  # Load the paths
  paths = list()
  svg = etree.parse(source)
  # Process any path elements with a 'd' attribute
  sx, sy, h = 1.0, 1.0, 0.0
  for path in svg.getiterator():
    if (path.tag == "svg") or path.tag.endswith("}svg"):
      sx, sy, h = getScaling(path)
    if (path.tag == "path") or path.tag.endswith("}path"):
      # Ignore 'red' construction lines
      s = path.get("style")
      if s is not None:
        if not styleContains(s, { 'stroke': "#000000" }):
          # Non-black strokes are construction lines
          continue
      d = path.get("d")
      if d is not None:
        paths.append(parse_path(d))
  print "Processing %d paths ..." % len(paths)
  results = list()
  for p in paths:
    results.extend(processPath(p, options.cut_depth, options.safe_depth, sx, sy, h, options.precision))
  # Write the file
  saveGCode(
    name + ".ngc",
    results,
    prefix = GCODE_PREFIX % options.safe_depth,
    suffix = GCODE_SUFFIX % options.safe_depth
    )


"""
