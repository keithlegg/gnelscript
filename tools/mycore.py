import sys
sys.path.append('Users/keithlegg/gnolmec')





import maya.cmds as cmds

import maya.api.OpenMaya as om

#selected = cmds.ls(sl=True,long=True) or []
#for e in selected:
#   print(e)




import maya.cmds as cmds

pathOfFiles = "/Users/keithlegg/gnolmec/3d_obj"
fileType = "obj"

files = cmds.getFileList(folder=pathOfFiles, filespec='*.%s' % fileType)
if len(files) == 0:
    cmds.warning("No files found")
else:
    for f in files:
        cmds.file(pathOfFiles + f, i=True)



###############
# create object to manipulate
cmds.sphere( n='sphere1' )
# set rotation of sphere
cmds.xform( r=True, ro=(0, 90, 0) )
# change the rotate order but preserve the overall transformation
cmds.xform( p=True, roo='yzx' )

#######
# To find all non-manifold edges on a polygonal object called pPlane1
cmds.polyInfo( nme=True )
# Result: pPlane1.e[74] #
# To find all non-manifold vertices on a polygonal object called pPlane1
cmds.polyInfo( nmv=True )
# Result: pPlane1.vtx[38] pPlane1.vtx[49] #





#######################

#OLD?
from maya import OpenMaya as om

sel = om.MSelectionList()
om.MGlobal.getActiveSelectionList(sel)
comp_iter = om.MItSelectionList(sel, om.MFn.kMeshVertComponent)
component = om.MObject()
meshDagPath = om.MDagPath()


vlist = []
flist=[]
elist=[]

while not comp_iter.isDone():
	comp_iter.getDagPath(meshDagPath, component)
	if not component.isNull():
		vert_iter = om.MItMeshVertex(meshDagPath, component)
		while not vert_iter.isDone():
			vlist.append(vert_iter.index())
			vert_iter.next()
	comp_iter.next()

	comp_iter.getDagPath(meshDagPath, component)
	if not component.isNull():
		edge_iter = om.MItMeshEdge(meshDagPath, component)
		while not edge_iter.isDone():
			elist.append(edge_iter.index())
			edge_iter.next()
	comp_iter.next()
					
print(vlist)
print(elist)



#######################


from maya.OpenMaya import * #import in this space for better visualisation

mesh = "pSphere1" 
vertices = MIntArray()
indices = [4,2,1]
[vertices.append(x) for x in indices] #there might be a better solution for initializing a MIntArray with an tuple. 

comp = MFnSingleIndexedComponent()
compObj = comp.create(MFn.kMeshVertComponent)
comp.addElements(vertices)

# we query the MDagPath of the mesh
dag = MDagPath()
sel = MSelectionList()
sel.add(mesh)
sel.getDagPath(0, dag) 

# we then include a MSelection with the dag and its component(MObject)
sel2 = MSelectionList()
sel2.add(dag, compObj)
MGlobal.setActiveSelectionList(sel2)




#######################

def aimNode(node,target, myAimAxis=[1,0,0], myUpAxis=[0,1,0], myUpDir=[0,0,1]):
    """rotate a transform to aim a target point
    Args:
        node (str): transform node name
        target (list): world space point expresed as [float, float, float]
        myAimAxis (list, optional): the axis that point to the target. Defaults to [1,0,0].
        myUpAxis (list, optional): the axis up to stabilize the rotation. Defaults to [0,1,0].
        myUpDir (list, optional): [description]. Defaults to [0,0,1].
    """    
    selList = om.MSelectionList()
    selList.add(node)
    mdpThis = selList.getDagPath(0)
    mpTarget=om.MPoint(target)

    mmMatrix = mdpThis.inclusiveMatrix();
    mmMatrixInverse = mmMatrix.inverse();
    mpTargetInThisSpace=om.MPoint(mpTarget * mmMatrixInverse);
    mqToTarget=om.MQuaternion (om.MVector(myAimAxis).rotateTo(om.MVector(mpTargetInThisSpace)));
    # Apply that rotation to this node:
    mftThis=om.MFnTransform (mdpThis.transform());
    mftThis.rotateBy(mqToTarget, om.MSpace.kPreTransform);
    
    mmMatrix = mdpThis.inclusiveMatrix();
    mmMatrixInverse = mmMatrix.inverse();
    mvUp=om.MVector(myUpDir)
    mpWorldUpInThisSpace=om.MPoint (mvUp * mmMatrixInverse);
    # Get the quaternion that describes the rotation to make the local up vector point to the world up vector:
    mqToWorldUp=om.MQuaternion (om.MVector(om.MVector(myUpAxis)).rotateTo(om.MVector(mpWorldUpInThisSpace)));
    merToWorldUp=om.MEulerRotation (mqToWorldUp.asEulerRotation());
    merToWorldUp.y = 0.0;
    merToWorldUp.z = 0.0;
    # Apply that rotation to this node:
    mftThis.rotateBy(merToWorldUp, om.MSpace.kPreTransform);


#######################
def getLocalTranslation(node, pos, frame):
    """the the local values to a note to reach a world position on specific frame
    Args:
        node (str): node name
        pos (list): world space point expresed as [float, float, float]
        frame (int): frame to get the node position

    Returns:
        list: node local space point expresed as [float, float, float]
    """    
    # for some reason DGContex wont update properlly so
    cmds.currentTime(frame) # force frame to refresh matrix value
    parent = cmds.listRelatives(node, p=1, f=1)
    if not parent:
        return pos
    matrix = om.MMatrix(cmds.xform(parent, q=1, ws=1, m=1))
    point = om.MPoint(pos)
    point = om.MVector(point * matrix.inverse())
    return point

#######################
def getLocalRotation(node, rot, frame):
    """the the local values to a note to reach a world rotation on specific frame
    Args:
        node (str): node name
        rot (list): world space roation expresed as [float, float, float]
        frame (int): frame to get the node position

    Returns:
        list: node local space rotation expresed as degrees [float, float, float]
    """   
    # for some reason DGContex wont update properlly so
    cmds.currentTime(frame) # force frame to refresh matrix value
    parent = cmds.listRelatives(node, p=1, f=1)
    matrix = om.MMatrix(cmds.xform(parent, q=1, ws=1, m=1))
    trfMatrix = om.MTransformationMatrix()
    rotation = om.MVector([math.radians(a) for a in rot])
    euler = om.MEulerRotation(rotation)
    trfMatrix.rotateBy(euler, om.MSpace.kTransform)
    newMat = trfMatrix.asMatrix()
    newTrfMat = om.MTransformationMatrix(newMat *  matrix.inverse())
    return [math.degrees(a) for a in newTrfMat.rotation()]

#######################
def getAimVector(source, targetPosition):
    matrix = om.MMatrix(cmds.xform(source, q=1, ws=1, m=1))
    tPos = om.MVector(targetPosition)
    axes = [om.MVector(1,0,0), om.MVector(0,1,0), om.MVector(0,0,1), 
            om.MVector(-1,0,0), om.MVector(0,-1,0), om.MVector(0,0,-1)]
    distances = list()
    for each in axes:
        distances.append((om.MVector(each*matrix) - tPos).length())
    indx = distances.index(min(distances))  
    return axes[indx]  




def vector_to(source=None, target=None):
    """Calculate the distance between two nodes

    :param source: First node
    :param target: Second node
    :return: MVector (API2)
    """
    if source is None or target is None:
        # Default to selection
        selection = cmds.ls(sl=True, type='transform')
        if len(selection) != 2:
            raise RuntimeError('Select 2 transforms.')
        source, target = selection

    pos1 = cmds.xform(source, query=True, worldSpace=True, translation=True)
    pos2 = cmds.xform(target, query=True, worldSpace=True, translation=True)

    source = OpenMaya2.MPoint(pos1[0], pos1[1], pos1[2])
    target = OpenMaya2.MPoint(pos2[0], pos2[1], pos2[2])
    return target - source 







def create_arrow(jointName):
    curve = cmds.curve(
        name="%s_ForwardDirection" % jointName,
        degree=1,
        point=[
            (-1, 0, 0),
            (-1, 2, 0),
            (-2, 2, 0),
            (0, 4, 0),
            (2, 2, 0),
            (1, 2, 0),
            (1, 0, 0),
            (-1, 0, 0),
        ],
    )
    group = cmds.group()
    cmds.xform(objectSpace=True, pivots=(0, 0, 0))
    jointScale = cmds.jointDisplayScale(query=True)
    jointRadius = cmds.getAttr("%s.radius" % jointName)
    jointScale *= jointRadius
    cmds.xform(scale=(jointScale, jointScale, jointScale))

    return group 





    def sqGetPointLists(self, *args):
        cmds.select(self.baseCurve+".cv[*]")
        pointList = cmds.ls(selection=True, flatten=True)
        
        minX = 0
        maxX = 0
        sideA = 0
        sideB = 0
        for i in range(0, len(pointList)):
            pointPosX = cmds.xform(pointList[i], query=True, worldSpace=True, translation=True)[0]
            if pointPosX < minX:
                minX = pointPosX
                sideA = i
            elif pointPosX > maxX:
                maxX = pointPosX
                sideB = i
        if sideA > sideB:
            sideC = sideA
            sideA = sideB
            sideB = sideC
        
        pointListA = pointList[sideA:(sideB+1)]
        pointListB = pointList[sideB:]
        for j in range(0, (sideA+1)):
            pointListB.append(pointList[j])
        
        return pointListA, pointListB, sideA, sideB 








from maya.OpenMaya import * #import in this space for better visualisation

mesh = "pSphere1" 
vertices = MIntArray()
indices = [4,2,1]
[vertices.append(x) for x in indices] #there might be a better solution for initializing a MIntArray with an tuple. 

comp = MFnSingleIndexedComponent()
compObj = comp.create(MFn.kMeshVertComponent)
comp.addElements(vertices)

# we query the MDagPath of the mesh
dag = MDagPath()
sel = MSelectionList()
sel.add(mesh)
sel.getDagPath(0, dag) 

# we then include a MSelection with the dag and its component(MObject)
sel2 = MSelectionList()
sel2.add(dag, compObj)
MGlobal.setActiveSelectionList(sel2)








