
#  http://linuxcnc.org/docs/html/gcode/overview.html

#-------------------------------------------
# A - A axis of machine
# B - B axis of machine
# C - C axis of machine
# D - Tool radius compensation number
# F - Feed rate
# G - General function (See table Modal Groups)
# H - Tool length offset index
# I - X offset for arcs and G87 canned cycles
# J - Y offset for arcs and G87 canned cycles
# K - Z offset for arcs and G87 canned cycles.
#     Spindle-Motion Ratio for G33 synchronized movements.
# L - generic parameter word for G10, M66 and others
# M - Miscellaneous function (See table Modal Groups)
# N - Line number
# P - Dwell time in canned cycles and with G4.
#     Key used with G10.
# Q - Feed increment in G73, G83 canned cycles
# R - Arc radius or canned cycle plane
# S - Spindle speed
# T - Tool selection
# U - U axis of machine
# V - V axis of machine
# W - W axis of machine
# X - X axis of machine
# Y - Y axis of machine
# Z - Z axis of machine

#-------------------------------------------

# 31-5000 - G code user parameters. These parameters are global in the G code file, and available for general use. Volatile.
# 5061-5069 - Coordinates of a G38 probe result (X, Y, Z, A, B, C, U, V & W). Coordinates are in the coordinate system in which the G38 took place. Volatile.
# 5070 - G38 probe result: 1 if success, 0 if probe failed to close. Used with G38.3 and G38.5. Volatile.
# 5161-5169 - "G28" Home for X, Y, Z, A, B, C, U, V & W. Persistent.
# 5181-5189 - "G30" Home for X, Y, Z, A, B, C, U, V & W. Persistent.
# 5211-5219 - "G92" offset for X, Y, Z, A, B, C, U, V & W. Persistent.
# 5210 - 1 if "G92" offset is currently applied, 0 otherwise. Persistent.
# 5211-5219 - G92 offset (X Y Z A B C U V W).
# 5220 - Coordinate System number 1 - 9 for G54 - G59.3. Persistent.
# 5221-5230 - Coordinate System 1, G54 for X, Y, Z, A, B, C, U, V, W & R. R denotes the XY rotation angle around the Z axis. Persistent.
# 5241-5250 - Coordinate System 2, G55 for X, Y, Z, A, B, C, U, V, W & R. Persistent.
# 5261-5270 - Coordinate System 3, G56 for X, Y, Z, A, B, C, U, V, W & R. Persistent.
# 5281-5290 - Coordinate System 4, G57 for X, Y, Z, A, B, C, U, V, W & R. Persistent.
# 5301-5310 - Coordinate System 5, G58 for X, Y, Z, A, B, C, U, V, W & R. Persistent.
# 5321-5330 - Coordinate System 6, G59 for X, Y, Z, A, B, C, U, V, W & R. Persistent.
# 5341-5350 - Coordinate System 7, G59.1 for X, Y, Z, A, B, C, U, V, W & R. Persistent.
# 5361-5370 - Coordinate System 8, G59.2 for X, Y, Z, A, B, C, U, V, W & R. Persistent.
# 5381-5390 - Coordinate System 9, G59.3 for X, Y, Z, A, B, C, U, V, W & R. Persistent.
# 5399 - Result of M66 - Check or wait for input. Volatile.
# 5400 - Tool Number. Volatile.
# 5401-5409 - Tool Offsets for X, Y, Z, A, B, C, U, V & W. Volatile.
# 5410 - Tool Diameter. Volatile.
# 5411 - Tool Front Angle. Volatile.
# 5412 - Tool Back Angle. Volatile.
# 5413 - Tool Orientation. Volatile.
# 5420-5428 - Current relative position in the active coordinate system including all offsets and in the current program units for X, Y, Z, A, B, C, U, V & W, volatile.
# 5599 - Flag for controlling the output of (DEBUG,) statements. 1=output, 0=no output; default=1. Volatile.
# 5600 - Toolchanger fault indicator. Used with the iocontrol-v2 component. 1: toolchanger faulted, 0: normal. Volatile.
# 5601 - Toolchanger fault code. Used with the iocontrol-v2 component. Reflects the value of the toolchanger-reason HAL pin if a fault occurred. Volatile.



# G17 use XY plane, 
# G20 inch mode, 
# G40 cancel diameter compensation, 
# G49 cancel length offset, 
# G54 use coordinate system 1, 
# G80 cancel canned cycles, 
# G90 absolute distance mode, 
# G94 feed/minute mode.


parser_commands =  {

    # 5220 - Coordinate System number 1 - 9 for G54 - G59.3. Persistent.
    # 5221-5230 - Coordinate System 1, G54 for X, Y, Z, A, B, C, U, V, W & R. R denotes the XY rotation angle around the Z axis. Persistent.
    # 5241-5250 - Coordinate System 2, G55 for X, Y, Z, A, B, C, U, V, W & R. Persistent.
    # 5261-5270 - Coordinate System 3, G56 for X, Y, Z, A, B, C, U, V, W & R. Persistent.
    # 5281-5290 - Coordinate System 4, G57 for X, Y, Z, A, B, C, U, V, W & R. Persistent.
    # 5301-5310 - Coordinate System 5, G58 for X, Y, Z, A, B, C, U, V, W & R. Persistent.
    # 5321-5330 - Coordinate System 6, G59 for X, Y, Z, A, B, C, U, V, W & R. Persistent.
    # 5341-5350 - Coordinate System 7, G59.1 for X, Y, Z, A, B, C, U, V, W & R. Persistent.
    # 5361-5370 - Coordinate System 8, G59.2 for X, Y, Z, A, B, C, U, V, W & R. Persistent.
    # 5381-5390 - Coordinate System 9, G59.3 for X, Y, Z, A, B, C, U, V, W & R. Persistent.

    "M8": "coolant_on", #?
    "M9": "coolant_off", #?
    "M5": "spindle_off", #?

    # G1 Code sets the Linear Interpolation mode. The format of the G1 command is:
    # G1X__Y__Z__F__;
    # Where X__Y__Z__ defines the endpoint of the move to be made. Simultaneous XYZ motion along a linear (straight
    # line) path occurs at the feedrate defined by the F__ word.    

    # "G0" : "rapid_traverse", #  G0/G1 - R__I__J__A__; Move to Pt
    "G00": "rapid_traverse",

    # "G1" : "linear_interpolation", 
    "G01": "linear_interpolation",

    # "G2": "circular_interpolation",        
    "G02": "circular_interpolation",

    "G4"  : "dwell",
    "G40" : "dwell",

    #"G12": "helical_1", # G12/G13 - A__Z__F__; Do Helix 
    #"G13": "helical_2", 

    #Circular interpolation is effective in the XY, YZ or ZX planes by preparatory codes as follows:
    "G17": "circle_XY_plane",
    "G18": "circle_ZX_plane",
    "G19": "circle_YZ_plane",

    #"G70": "inch_conversion",
    #"G71": "metric_conversion",
    "G20": "inch_conversion",
    #"G71": "metric_conversion",

    #"G80": "drill_cycle",

    #"G90": "abs_progr",  # data in a word represents the coordinate value of the point from part program zero.
    #"G91": "incr_progr", # data in a word is distance along designated axis from the existing position to the desired position

    #"G175": "inside_mill_cycle",   # G175 X__Y__Z__R__Z__Z__P__F__P__F__F__;
    #"G176": "outside_mill_circle", # G176 X__Y__Z__R__Z__Z__P__F__P__F__F__;

    # "X": "x_axis",
    # "Y": "y_axis",
    # "Z": "z_axis",  

    # X-Y Plane (G17)
    # G2X__Y__I__J__F__; CW
    # G3X__Y__I__J__F__; CCW
    
    # Z-X Plane (G18)
    # G2X__Z__I__K__F__; CW
    # G3X__Z__I__K__F__; CCW
    
    # Y-Z Plane (G 19)
    # G2Y__Z__J__K__F__; CW
    # G3Y__Z__J__K__F__; CCW

    #If G75 is active and G90 (absolute data input) is also active, I__, J__, K__ specifies the arc center coordinates with
    #respect to part program zero. If G91 (incremental data input) is active, I__, J__, K__ specifies the signed distance from
    #the arc startpoint to the center of the arc.
    
    #If G74 is active, I__, J__, K__ specifies the unsigned magnitude of the distance from the arc startpoint to the arc center
    #in G90 (Absolute) or G91 (Incremental) mode.

    # "I": "i_axis",
    # "J": "j_axis",
    # "K": "k_axis",  

    # "U": "u_axis",
    # "V": "v_axis",
    # "W": "w_axis", 

    #"M3": "end?", 
    "M30": "end",

    "E": "end",
}

