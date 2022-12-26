

# Builtin Modules:       bpy, bpy.data, bpy.ops, bpy.props, bpy.types, bpy.context, bpy.utils, bgl, blf, mathutils
# Convenience Imports:   from mathutils import *; from math import *
# Convenience Variables: C = bpy.context, D = bpy.data

# bpy.ops.mesh.primitive_cube_add(enter_editmode=False, align='WORLD', location=(0, 0, 1), scale=(1, 1, 1))
# bpy.ops.mesh.primitive_cylinder_add(enter_editmode=False, align='WORLD', location=(0, 0, 0), scale=(1, 1, 1))
# bpy.data.objects['Cylinder'].select_set(False); 

import sys

TOOLSPATH = '/Users/klegg/serv/gnolmec'

sys.path.append(TOOLSPATH)

BYCORE_OBJ_OUT= '%s/3d_obj/BYCORE.obj'%TOOLSPATH


from gnelscript.pygfx.math_ops import  *
from gnelscript.pygfx.point_ops import *
from gnelscript.pygfx.obj3d import  *

import bpy
import bmesh

import math
import mathutils

from bpy import context




""" 

open in terminal to get python output 

/Applications/Blender.app/Contents/MacOS/Blender


#-------
#notes from blender cam (not needed but relevant)

https://github.com/vilemduha/blendercam/blob/master/documentation/Blendercam%20Installation.md


If you are using a Blender with a bundled Python then packagers must be installed in the site-packages directory 
of the bundled python. For Blender 2.8 only Shapely is needed. 
To install it, open terminal, get to Blender directory and use PIP:

cd /Applications/Blender.app/Contents/Resources/3.3/python/bin  

./python3.10 -m ensurepip

./python3.10 -m pip install shapely
./python3.10 -m pip install pygeos


"""

##---------------------------------------


"""
    While blender stores the mesh data within object.data, this data is only valid in object mode, when you switch to edit mode a 
    bmesh copy of the mesh data is created, which replaces the object.data contents when you leave edit mode. 
    As you are using a duplicate mesh while editing any selection changes you make to object.data with python 
    will not effect the edit mesh and will be overwritten when exiting edit mode.


    Modes and Mesh Access
    When working with mesh data you may run into the problem where a script fails to run as expected in Edit-Mode. This is caused by Edit-Mode having its own data which is only written back to the mesh when exiting Edit-Mode.
    A common example is that exporters may access a mesh through obj.data (a bpy.types.Mesh) when the user is in Edit-Mode, where the mesh data is available but out of sync with the edit mesh.
    In this situation you can…
        Exit Edit-Mode before running the tool.

        Explicitly update the mesh by calling bmesh.types.BMesh.to_mesh.
        
        Modify the script to support working on the edit-mode data directly, see: bmesh.from_edit_mesh.
        
        Report the context as incorrect and only allow the script to run outside Edit-Mode.





    #--- 

    #get/set obejct mode 

    bpy.context.active_object.mode   # = 'OBJECT'

    bpy.ops.object.mode_set(mode='EDIT')

    bpy.context.active_object.mode   # = 'EDIT'


    #set EDIT/OBJECT mode  

    # bpy.ops.object.mode_set(mode = 'OBJECT')
    OR 
    # bpy.ops.object.editmode_toggle()

    #---

    # get polygons edited in EDIT mode 

    #objprior = context.active_object
    ## get data after modifiers ??
    #depsgraph = bpy.context.evaluated_depsgraph_get()
    #obj = objprior.evaluated_get(depsgraph)

    #---

    #iterate polygons selected in EDIT mode 

    # for ply in bpy.context.object.data.polygons:
    #     print(ply.select)



"""



##---------------------------------------



##---------------------------------------

"""
Setting bl_space_type = 'VIEW_3D' will place your panel within the 3dview. 
You then set bl_region_type = 'TOOLS' to have it show in the tools region 
while the less intuitive bl_region_type = 'UI' places it in the properties region.

If you want to place a panel in the properties editor, you would set bl_space_type = 'PROPERTIES' 
and bl_region_type = 'WINDOW', then you set bl_context to match the context you want it shown in  

- Object, Scene, World, Modifier...

Also remember that you need to register your panel class -
"""


class MyOperator(bpy.types.Operator):
    bl_idname = "scene.myoperator"
    bl_label = "export gnel"

    def execute(self, context):
        obj_to_gnel()
        return {"FINISHED"}

class View3DPanel:
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'TOOLS'
    
    #bl_region_type = 'UI'
    #bl_category = "Tool"

    @classmethod
    def poll(cls, context):
        return (context.object is not None)


class PanelOne(View3DPanel, bpy.types.Panel):
    bl_idname = "VIEW3D_PT_test_1"
    bl_label = "Tools"

    def draw(self, context):
        row = self.layout.row()
        row.operator("scene.myoperator")        
        self.layout.label(text="Small Class")

"""
class PanelTwo(View3DPanel, bpy.types.Panel):
    bl_idname = "VIEW3D_PT_test_2"
    bl_label = "toolz"
    def draw(self, context):
        self.layout.label(text="Also Small Class")
"""

bpy.utils.register_class(PanelOne)
#bpy.utils.register_class(PanelTwo)
bpy.utils.register_class(MyOperator)


##---------------------------------------

##---------------------------------------


#TO ADD TO GNELSCRIPT  

#cross section - start with a cube - then to a trapezoid cube, work up 
#intersection  
#projection  

#extents_3d
#extents_2d

#subdiv 

#spiral_2d( intersect_obj )

#spiral_fill( shape, intersect  ) 

#spiral_3d 







##---------------------------------------

def get_selection():
    """ get anything selected - object or component """
 
    context = bpy.context
    mode = context.active_object.mode 
    
    objects = []
    faces = []
    edges = []
    points = [] 
    
    if mode =='EDIT':
        print("# you are in edit mode ")
        ob = context.edit_object
        
        edit_mesh = ob.data
        # see: https://docs.blender.org/api/current/bmesh.types.html
        bm = bmesh.from_edit_mesh(edit_mesh)

        print("faces ", len(bm.faces))
        if len(bm.faces):
            for f in bm.faces:
                if f.select:
                    print(f.index)
                    print(f.normal)
                    print(f.loops)

        print("edges ", len(bm.edges))
        if len(bm.edges):
            for e in bm.edges: 
                if e.select:
                    print(e.index)
                    print(e.is_boundary)

        print("vertices ", len(bm.verts))
        if len(bm.verts):
            for v in bm.verts:            
                if v.select:
                    print(v.index)                    
                    print(v.co)

    if mode =='OBJECT':
        print("# you are in object mode ")
        obj = context.active_object

        faces = obj.data.polygons        
        print("faces ", len(faces))

        edges = obj.data.edges        
        print("edges ", len(edges))

        loops = obj.data.loops        
        print("loops ", len(loops))

        vertices = obj.data.vertices
        print("vertices ", len(vertices))


get_selection()

##---------------------------------------

def obj_to_gnel():
    """ example to export a polygon object selected in OBJECT mode """

    gn_polys = []
    gn_pts = []

    obj = context.active_object

    if obj:
        verts = obj.data.vertices
        for v in verts: 
            #blender vectors
            
            # matrix_parent_inverse
            # matrix_local
            # matrix_basis
            # -  Matrix access to location, rotation and scale (including deltas), 
            # -  before constraints and parenting are applied

            blvec = obj.matrix_world @ v.co
            #raw float data  
            gn_pts.append(blvec[:])

        #(reference a range of loops)
        faces = obj.data.polygons
        for f in faces:
            if f.select:
                #print(f.normal)
                #print(f.loop_start)
                #print(f.loop_total)
                gn_polys.append(f.vertices[:])

        #build a gnelscript object and save as OBJ 
        obj = object3d()
        geom = obj.insert_polygons(plyids=gn_polys, points=gn_pts, incrone=True) 

        obj.insert(geom) 
        obj.show()
        obj.save(BYCORE_OBJ_OUT)

   

##---------------------------------------

def obj_to_gnel2():
    """ example to iterate components in edit mode """

    gn_polys = []
    gn_pts = []

    context = bpy.context
    ob = context.edit_object
    me = ob.data

    bm = bmesh.from_edit_mesh(me)
    # list of selected faces
    selfaces = [f for f in bm.faces if f.select]
    if selfaces:
        print("%d faces selected" % len(selfaces))
    else:
        print("No Faces Selected")




##---------------------------------------

def sel_to_gnel():
    """ UNFINISHED 
        iterate mesh  
        mode and convert to gnelscript object 
    """

    ob = bpy.context.active_object
    me = ob.data

    selfaces =[]
    selverts =[]

    # check for edit mode
    editmode = False
    if ob.mode == 'EDIT':
        editmode =True
        # the following sets mode to object by default
        bpy.ops.object.mode_set()

    #(3 points in space)        
    for v in me.vertices:
        if v.select:
            #print(v.index)
            selverts.append(v)

    #(reference a range of loops)
    for f in me.polygons:
        if f.select:
            #print(f.index)
            selfaces.append(f)

    #(reference 2 vertices)
    for e in me.edges:
        print('edge ', e)  

    # (reference a single vertex and edge)
    for l in me.loops:
        print('loops ', e) 



##---------------------------------------
# 


def build():
    """ bmesh operators example  
        from : https://docs.blender.org/api/current/bmesh.ops.html
    """ 

    # Make a new BMesh
    bm = bmesh.new()

    # Add a circle XXX, should return all geometry created, not just verts.
    bmesh.ops.create_circle(
        bm,
        cap_ends=False,
        radius=0.2,
        segments=8)


    # Spin and deal with geometry on side 'a'
    edges_start_a = bm.edges[:]
    geom_start_a = bm.verts[:] + edges_start_a
    ret = bmesh.ops.spin(
        bm,
        geom=geom_start_a,
        angle=math.radians(180.0),
        steps=8,
        axis=(1.0, 0.0, 0.0),
        cent=(0.0, 1.0, 0.0))
    edges_end_a = [ele for ele in ret["geom_last"]
                   if isinstance(ele, bmesh.types.BMEdge)]
    del ret


    # Extrude and create geometry on side 'b'
    ret = bmesh.ops.extrude_edge_only(
        bm,
        edges=edges_start_a)
    geom_extrude_mid = ret["geom"]
    del ret


    # Collect the edges to spin XXX, 'extrude_edge_only' could return this.
    verts_extrude_b = [ele for ele in geom_extrude_mid
                       if isinstance(ele, bmesh.types.BMVert)]
    edges_extrude_b = [ele for ele in geom_extrude_mid
                       if isinstance(ele, bmesh.types.BMEdge) and ele.is_boundary]
    bmesh.ops.translate(
        bm,
        verts=verts_extrude_b,
        vec=(0.0, 0.0, 1.0))


    # Create the circle on side 'b'
    ret = bmesh.ops.spin(
        bm,
        geom=verts_extrude_b + edges_extrude_b,
        angle=-math.radians(180.0),
        steps=8,
        axis=(1.0, 0.0, 0.0),
        cent=(0.0, 1.0, 1.0))
    edges_end_b = [ele for ele in ret["geom_last"]
                   if isinstance(ele, bmesh.types.BMEdge)]
    del ret


    # Bridge the resulting edge loops of both spins 'a & b'
    bmesh.ops.bridge_loops(
        bm,
        edges=edges_end_a + edges_end_b)


    # Now we have made a links of the chain, make a copy and rotate it
    # (so this looks something like a chain)

    ret = bmesh.ops.duplicate(
        bm,
        geom=bm.verts[:] + bm.edges[:] + bm.faces[:])
    geom_dupe = ret["geom"]
    verts_dupe = [ele for ele in geom_dupe if isinstance(ele, bmesh.types.BMVert)]
    del ret

    # position the new link
    bmesh.ops.translate(
        bm,
        verts=verts_dupe,
        vec=(0.0, 0.0, 2.0))
    bmesh.ops.rotate(
        bm,
        verts=verts_dupe,
        cent=(0.0, 1.0, 0.0),
        matrix=mathutils.Matrix.Rotation(math.radians(90.0), 3, 'Z'))

    # Done with creating the mesh, simply link it into the scene so we can see it

    # Finish up, write the bmesh into a new mesh
    me = bpy.data.meshes.new("Mesh")
    bm.to_mesh(me)
    bm.free()


    # Add the mesh to the scene
    obj = bpy.data.objects.new("Object", me)
    bpy.context.collection.objects.link(obj)

    # Select and make active
    bpy.context.view_layer.objects.active = obj
    obj.select_set(True)




##---------------------------------------



##---------------------------------------


"""
#get point location and put a locator at thatr 

obj = context.active_object
v = obj.data.vertices[0]
co_final = obj.matrix_world @ v.co

# now we can view the location by applying it to an object
obj_empty = bpy.data.objects.new("Test", None)
context.collection.objects.link(obj_empty)
obj_empty.location = co_final
"""



#Another common operation is to get all vertex locations with the objects transform applied, in this case you can use list comprehension to get a list of transformed vertex locations...
#coords = [(obj.matrix_world @ v.co) for v in obj.data.vertices]

############################################

# this returns True if the first polygon is selected,
# or False in case it is not

#bpy.context.object.data.polygons[0].select

############################################

"""
import bpy
import numpy as np

from bpy import context

toggle = context.mode == 'EDIT_MESH'
if toggle:
    bpy.ops.object.editmode_toggle()

print("-" * 44)
dg = context.evaluated_depsgraph_get()
all_areas = []        
for ob in context.selected_objects:
    me = dg.objects[ob.name].data
    me.transform(ob.matrix_world) # globalize coords
    data = np.empty(len(me.polygons))
    me.polygons.foreach_get("select", data)
    selection = data.astype(bool)
    me.polygons.foreach_get("area", data)
    areas = data
    
    print(ob.name)
    print("Total ob area", np.sum(areas))
    print(np.count_nonzero(selection), "/", len(selection))
    all_areas.append(np.sum(areas[selection]))
    print("Selected area", all_areas[-1], "sq bu")
    #bpy.data.meshes.remove(me)
    print()

print("Total Area Selected")
print(sum(all_areas), "sq bu")
    
if toggle:
    bpy.ops.object.editmode_toggle() 

"""

############################################


#load object 
#bpy.ops.import_scene.obj(filepath='/Users/keithlegg/gnolmec/3d_obj/sphere.obj')

############################################




