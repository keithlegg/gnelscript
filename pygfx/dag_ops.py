
import codecs
import os
import re

from pathlib import Path

#from gnelscript.pygfx.point_ops import *


""" 

    THIS CODE IS IN NEED OF A REWRITE 

    TODO:

        make full on scene graph? 
        
        get rid of "obj" and "name" modes - just check input type 
        pass walk callback 
        cache data on nodes like:
            depth 
            ..?

        WHY THE HELL DO WE HAVE "PARENT_BUFFER_SCAN" 

        ADD ATTRS TO GRAPH ITSELF 

        ENFORCE RECURSION ERRORS - IE - prnt('chet','jones').prnt('jones','chet')
            -parent functon among others 

        clone branch 

        clone_branch as new tree

        graft_tree as branch

        delete branch  
        
        add_node_as_leaf 


"""

  


def new_stack (graphobj=None): 
    """ return a recursively nestable graph interface 
        in oder to "prime the pump" and launch a dagop stack
        we need to have an a dagop object with a data_tree passed in
        This function exists to bootstrap that process  
    """

    if graphobj is None:
        return dagop(data_graph('dd')) 

    if type(graphobj)==data_graph:
        return dagop(graphobj)

    if type(graphobj)==str:
        path = Path(graphobj)

        if path.is_file():
            DG = data_graph('dd')
            DG.load_graph_file(graphobj)
            return dagop(DG)
        
        else:
            DG  = data_graph(graphobj)
            return dagop(DG) 
      


       

##-------------------------------------##
class dagop(object):
    """ DAG operator- interface to the node and tree classes 
        self-referential, stackable and fun to play with 

        idea to explore ->  add a history component to this 
        so each graph returns a "stack" of operations that becomes 
        a modular instruction to reuse
    """
        
    def __init__(self, PKGRAPH):
        #merge with OTHER dagops? 
        #-perform branches and other fun things 
        self.DG = PKGRAPH

    def __add__( self, other ):
        #print type( other )
        raw_input()
        for newnod in other.__getNodes():
            self.DG.add( newnod )
        return self.DG
    
    def __get( self, nod ):
        self.DG.get( nod )
        return self.DG
        
    def __del( self, nod ):
        self.DG.delete( nod )
        return self.DG

    def __create( self, name, type ):
        node =  node_base( name )
        self.DG.add( node )
        return self.DG

    def __parent( self, cnode, pnode ):
        self.DG.parent(cnode,pnode)
        return self.DG

    def __addattr( self, nodename, attrname, attrval ):
        #self.DG.parent_name(cnode,pnode)
        return self.DG

    ##--------------------
    # non recursive functions 

    def __get_nodes( self ):
        return self.DG.getnodes()

    ##--------------------

    #DEBUG - NOT DONE - attach other tress as a branch 
    def attach_branch( self, insert_point, BRANCH ):
        #for node in BRANCH.walk():
            #add node to this branch
        #parent root node to inser_point 
        return self.DG

    #DEBUG - NOT DONE - save a branch as a tree
    def export_branch( self, insert_point, BRANCH ):
        #for node in BRANCH.walk():
            #add node to this branch with parenting, links, and attrs
        #parent root node to inser_point 
        return self.DG

    def delete_branch( self, branch_node ):
        #for node in BRANCH.walk():
            #add node to this branch
        #parent root node to inser_point 
        return self.DG

    ##--------------------
    #for debugging
    def pause(self):
       print('debug pause')
       raw_input()
       return dagop(self.DG)

    def load(self, graphfile):
       print('loading %s'%graphfile )
       self.DG.load_graph_file(graphfile)
       return dagop(self.DG)

    def save(self, graphfile):
       print('loading %s'%graphfile)
       self.DG.save_graph_file(graphfile)
       return dagop(self.DG)
      #######

    #DEBUG - NOT DONE -  
    def rotate(self, node, rotation):
        print('rotating node %s')%node
        #self.__createNode(name,type)
        return dagop(self.DG)

    #DEBUG - NOT DONE -
    def move(self, node, translation):
        print('moving node %s'%node)
        #self.__createNode(name,type)
        return dagop(self.DG)
     
    def adn(self, name, type):
        print('adding node %s type %s'%(name, type) )
        self.__create(name,type)
        return dagop(self.DG)

    def dln(self, name):
        print('deleteing node %s'%name)
        self.__del(name)
        return dagop(self.DG)

    def prnt(self, cname, pname):
        print ('parenting node %s to %s'%(cname, pname) )
        self.__parent(cname, pname)
        return dagop(self.DG)

    def reset(self):
       print('reseting ')
       self.DG.reset()
       return dagop(self.DG)

    def show( self, mode='full' ):
        self.DG.show( mode )
        return dagop(self.DG)

    ##################
    #these wont work with the stackable interface
    #because they dont return the right object 
    def get_nodes(self):
       return self.DG.getnodes()
    
    def gtn(self, nod):
       return self.DG.getnodes( nod )

##-------------------------------------##
class data_graph(object):
    def __init__( self, graphnamevar='dft_name' ):
        self.parse_object  = basic_parser() #to be continued (XML,etc) 
        self.graphname     = graphnamevar  #name of the graph
        self.graphtype     = 'dft_type' 
        self.nodes         = []  #the actual list of nodes in the graph
        self.textbuffer    = []  #context of file to write
        self.filetextname  = ''  #name for file to write
        self.textindent    = 3   #when writing a file out , an indentaion indicates attributes of a node
        ##########################
        #walkbuffer is also populated on file load
        self.walkbuffer    =[]   #A variable outside the memory of a method for storing a dag walk
        self.walkindent    =0
        self.walk_depth_buffer = [] #a list of the depths of all nodes
        self.PARENT_BUFFER_SCAN =[]
        self.last_created = '' #keep track of last node created 
         
    def __del__(self):
        del self.nodes
        del self.textbuffer

    def reset(self):
        self.nodes              = []  #the actual list of nodes in the graph
        self.textbuffer         = []  #context of file to write
        
    def setupwalk(self):
        self.walkbuffer    =[]   #for graph walking
        self.walkindent    =0    #for graph walking

    def getnodes(self):
      return self.nodes

    def setname(self, name):
      self.name=str(name)

    # def unlink(self,node1,node2):
    # def link_name(self,node1,node2):   
    # def copy_node(self,nodetocopy,nameofnew):
    # def is_root (self,node,parent):
    # def get_root (self,node,parent):
    # def is_child_of (self,node,parent):
    # def listchildren(self):
    # def unparent (self,node):
    # def get_downstream_datafolder(self,node):    

    def refresh(self):
      #A WAY TO CHECK THE "HEALTH" OF THE GRAPH
      # children may not be set for each node
      # becuase only parents are stored in the file format??
      # scan the nodes and fix any errors
      
      ##add in children attributes for each node ##
      for node in self.nodes:
         print('NODE IS '%node)

         fulnam = self.get_full_name( node )
         node.name = fulnam
         
         if node.hasparents:
            #print node.name +' HAS PARENTS! refrsh'
            parents = node.listparents()
            #print 'PARENTS ARE '
            #print parents
            
            #goto parent and check that current node is listed as child
            #if not add it in
            
         #if node.haschildren():
         #   print node.name +' HAS PARENTS! refrsh'
        
         ##check for nodes with same names          ##
         #check for empty spaces
         #check node health
         #check parenting
         #check attr health
         #check link health
         #duplicate names      


    def exists(self, nod):
        for node in self.nodes:
            if node is nod:
                return True  
            if type(nod) is node_base:
                if node.name==nod.name:
                    return True     
            if type(nod) is str:
                if node.name==nod:
                    return True 
            if type(nod) is str:#unicode:
                if node.name==nod:
                    return True 
        return False

    def get(self, nod):
        """ changed from find_node - 
        works with string or node - if it exists """
        if self.exists(nod):
            for node in self.nodes:
                #print("node %s type %s " % (node.name, node.nodetype) )
                if node is nod:
                    return node
                if type(nod) is node_base:
                    if node.name==nod.name:
                        return node     
                if type(nod) is str:
                    if node.name==nod:
                        return node 

    def find_any_attr(self, attrs, attrval):
        for node in self.nodes:
            for attr in attrs:
                if node.getattrib(attr)== attrval:
                   return node
        return None

    def find_attr(self, attrnam, attrval):
        for node in self.nodes:
            #print(node.getattrib(attrnam))
            if node.getattrib(attrnam)==attrval:
                return node
        return None
                                               
    def gettype(self, nod):
        """ type is not PYTHON type - type is a node attr """
        
        if type(nod) is node_base:
                return nod.gettype()
                
        if type(nod) is str:
            node = self.get(nod)
            if node !=None:
                return node.gettype()
        return None
   
    def settype(self, nod, ntype):
        """ type is not PYTHON type - type is a node attr """
        
        node = self.get( nod )
        if node !=None:
            return node.settype( ntype )
             
    def getnumber(self):
        return len(self.nodes)

    def getnodenames(self):
        output =[]
        for node in self.nodes:
            output.append(node.name)
        return output


    #PARENT (DIRECT OBJECT METHOD)
    def parent(self, node1, node2):
        node1 = self.get(node1)
        node2 = self.get(node2)            

        if node1 == None or node2 ==None:
            print('##### parent: nodes not found')
            return None
        
        if node1==node2:
            return None
        
        if not node1.hasparents:
            node1.hasparents = 1
            node1.parent.append('')
            node1.parent[0]= (node2.name)
            node2.haschildren = 1
            node2.children.append(node1.name)
 

        
    def copy_attrs(self, srcnod, destnod):
        node1 = self.get(srcnod)
        node2 = self.get(destnod)
        
        node2.ATTRS     = node1.ATTRS
        node2.ATTRVALS  = node1.ATTRVALS
        node2.ATTRTYPS  = node1.ATTRTYPS
 
    def node_iterate_getcurrent(self):
        return self.walkbuffer[self.walkindent]

    def node_iterate_getcurr_depth(self):
         #DEBUG THIS WONT RETURN A VALUE FOR THE LAST ONE
         #i think the list is not getting populated correctly
         if self.walkindent<=len(self.walk_depth_buffer)-1:
             #print self.walk_depth_buffer[self.walkindent]
             return self.walk_depth_buffer[self.walkindent]
             #this is a big hack , it is wrong!!
         if self.walkindent<=len(self.walk_depth_buffer):
             return 0
             #this is a big hack , it is wrong!!
             #fix the way list gets populated

    #debug this may need work
    def node_iterate_cur_parent(self,node):
        zeroindx= self.walkindent-1
        node = self.walkbuffer[zeroindx]
        return node.parent

    def node_iterate_next(self):
        out= self.walkbuffer[self.walkindent]
        self.walkindent=self.walkindent+1
        return out
       
    def walk_firstnode(self):
        firstnode = self.listnodes()
        temp = firstnode[0]
        self.walk(temp.name,'object')

    def get_node_attr( self, nod, attr ):
        node = self.get(nod)
        output=node.getattrib(attr)
        return output


    def get_nodes_with_attr(self, attr, mode='obj'):
      """ return nodes that contain an attr """
      output =[]

      for nodelist in self.nodes:
         allattrs= nodelist.listallattrs()
         if allattrs!=None:
           for attrname in allattrs:
              if attrname[0]==attr:
                 if (mode=='obj' or mode=='object'):  
                    output.append(nodelist)
                 if (mode=='name'):  
                    output.append(nodelist.name)
                                   
      return output

    def rename_node(self, oldname, newname):
       checkexsts = self.get(newname)
       if checkexsts!=None:
            #print 'Node ' + newname + ' already exists'
            return None
   
       nodobj =self.get(oldname)
       nodobj.name=newname
       parents   = nodobj.listparents()
       children = nodobj.listchildren()
       for parent in parents:
         pnode = self.get( parent)
         if pnode:
           pnode.children = ( [nodobj.name] )

       for child in children:
         cnode = self.get( child)
         if cnode:
           cnode.parent = ( [nodobj.name] )

    def set_attr(self, nodename, attr, value):
      node = self.get( nodename )
      if node ==None:
        #print 'ERROR NO NODE '+nodename
        return None

      node.setattrib(attr,value)

    def add_attr(self, nodename, attr, value):
      node = self.get(nodename)
      if node ==None:
        #print 'ERROR NO NODE '
        return None
        
      node.addattr(attr,value)

    
    def fixname(self, name):
      """ attempt to fix duplicate names """
      cnt = 1
      x = name
      if self.exists(name):
          while self.exists(x):
             x = name+str(cnt)     
             cnt+=1
      return x  
        
             
    def create(self, name):
      if self.exists(name): 
          newnode = node_base(self.fixname(name) )
          self.nodes.append(newnode)
                     
      if not self.exists(name):          
          newnode = node_base(name)
          self.nodes.append(newnode)
          self.last_created = newnode
          
      return newnode

    ################
    ## COULD MERGE THESE 

    def createnode_parent(self, name, parentnode):
      """ create a node and parent it to an existing node """
      pnod        = self.get(parentnode)
      checkexsts = self.get(name)
      if checkexsts!=None:
        #print 'node already exists : '+name
        return None

      if pnod ==None:
          #print 'ERROR NO NODE FOUND '
          return None
          
      if pnod !=None:
        newnode = node_base(name)
        self.nodes.append(newnode)
        self.parent( name, pnod.name)
        #parent_obj
        self.last_created = newnode
      return newnode

    def createnode_parent_attr(self, name, parentnode, ATTR, ATTRVAL):
      """ create a node with atrrs, parent to an existing node """
      PAR        = self.get(parentnode)
      checkexsts = self.get(name)
      if checkexsts!=None:
        print('node already exists error ')
        return None

      if PAR ==None:
          print( 'ERROR NO NODE FOUND ')
          return None

      if PAR !=None:
        newnode = node_base(name)
        newnode.addattr(ATTR,ATTRVAL)
        self.nodes.append(newnode)
        self.parent(name,PAR.name)
        self.last_created = newnode
      #return node so we can still work with it
      return newnode

    ################       

    def insert_node_parent(self, nodeobj, parname):
      PAR = self.get(parname)
      if PAR == None:
          print('ERROR NO NODE FOUND ')
          return None
      if PAR != None:
        self.nodes.append(nodeobj)
        self.parent(nodeobj.name,PAR.name)
      return None

    def delete(self, nodename):
        nod_obj = self.get(nodename)

        if nod_obj!= None:
            parents   = nod_obj.listparents()
            children = nod_obj.listchildren()
            for child in children:
                chobj = self.get(child)
                if chobj!= None:
                    chobj.parent = ( parents )
        temparry = []#
        for tempnode in self.nodes:
           if tempnode.name!=nodename:
              print('omitting %s'%tempnode.name)

              temparry.append(tempnode)
        self.nodes= temparry

    def write_file (self):
      fhandler = data_tree_file_io()
      fhandler.writefile_listlines(self.filetextname,self.textbuffer)
      print('#file saved: %s'%self.filetextname)

    def check_valid_node( self, query):
        for node in self.nodes:
            name = node.name
            if name == query:
                print('test test ', name)


    def is_parent_of (self, parentnode, node):
      self.PARENT_BUFFER_SCAN = [] 
      parentnode = self.get(parentnode)
      node = self.get(node)
         
      if self.exists(parentnode)==0:
        return False
      
      if node == parentnode or self.exists(node)==0:
         return False
         
      for item in  self.scan_parents(node):
         if item == parentnode:
            return True
      return False
      
   
    def list_parents_name (self, nodename):
       NOBJ = self.get(nodename)
       return self.list_parents(NOBJ)
       
    def list_parents ( self, node ):
      #setup walk first
      self.PARENT_BUFFER_SCAN = []
      self.scan_parents(node)
      return  self.PARENT_BUFFER_SCAN

    def scan_parents ( self, node ):
      """ same as walk_up but with different buffer """ 
      node = self.get(node)
      parents = node.listparents()
      if node != None:
          self.PARENT_BUFFER_SCAN.append(node)
      if parents:
        PARENTNODE =  self.get(parents[0])
        if PARENTNODE:
           self.scan_parents(PARENTNODE)
      return self.PARENT_BUFFER_SCAN[1:]
      
      
    def list_children (self, node):
        NOD=self.get(node)
        if NOD==None:
          print('ERROR datagraph.list_children no node %s'%node)
        if NOD!=None:
          return NOD.children
   
    def list_children_obj (self, node):
          CHILLINS = self.list_children(node)
          outout = []
          for child in CHILLINS:
             chillnodetmp = self.get(child)
             outout.append(chillnodetmp)
          return outout

    #get actual node object from name
    def find_rotate(self, nodename, xyzarray):
     NODEOBJ = self.get(nodename)
     if NODEOBJ !=None:
        #print 'debug rotating node '+nodename
        #print xyzarray
        NODEOBJ.setrotate(xyzarray)
        
    ###
    def find_move(self, nodename, xyzarray):
     NODEOBJ = self.get(nodename)
     if NODEOBJ !=None:
        print(xyzarray)
        NODEOBJ.xform_xyz(xyzarray)
      
    def save_graph_file (self, filename):
      self.textbuffer   =[]
      self.filetextname = str(filename)
      nodeindex = 0
      self.textbuffer.append('##GENERATED WITH kdag.data_graph ### \n')
      self.textbuffer.append('##GRAPH NAME IS  ' +self.graphname +'\n')

      for element in self.nodes:
         #if node has inputs display them on same line as name
         if element.hasinputs==1:
            inputsarray= element.listinputs()
            inputone = inputsarray[0]
            ##ERROR CHECK
            if len(element.listinputs()):
               self.textbuffer.append('\n#############################')
               self.textbuffer.append('#NODE ' + str(element.name) +' ( '+str(inputone)+' ) \n')
         #######... , other wise just display name
         if element.hasinputs==0:
            self.textbuffer.append('\n#############################')
            self.textbuffer.append('#NODE ' + str(element.name) +'\n')

         if element.hasoutputs==1:
            #print 'DEBUG HASOUTPUTS'
            outputs = element.listoutputs()
            for outp in outputs:
              #print '#OUTPUT '+ outp
              self.textbuffer.append((self.textindent*' ')+'#OUTPUT '+ outp)

         if element.hasxform==1:
            self.textbuffer.append((self.textindent*' ')+'#SET_POS ' + str(element.posx) +' '+ str(element.posy) +' '+ str(element.posz)  +'\n')
            self.textbuffer.append((self.textindent*' ')+'#SET_ROT ' + str(element.rotx) +' '+ str(element.roty) +' '+ str(element.rotz)  +'\n')

         if element.hasparents==1:
              self.textbuffer.append((self.textindent*' ')+'#SET_PARENT '+element.parent[0]+'\n')

            #self.textbuffer.append((self.textindent*' ')+'#SET_ROT ' + str(element.rotx) +' '+ str(element.roty) +' '+ str(element.rotz)  +'\n')
         if element.haschildren==1:
              childrennods = element.listchildren()
              for childscan in childrennods:
                   self.textbuffer.append((self.textindent*' ')+'#SET_CHILD '+childscan+'\n')

         #add type
         querytype= element.gettype()
         if querytype!='dft_base':
              self.textbuffer.append((self.textindent*' ')+'#SET_TYPE '+querytype+'\n')
              
         #ATTRS
         attrlist  = element.ATTRS     #NAMES OF ATTRS
         valuelist = element.ATTRVALS  #VALUES OF ATTRS
         attrindex =0
         for attr in attrlist:
             self.textbuffer.append((self.textindent*' ')+'#ATTR '+attr+' '+str(valuelist[attrindex])+' \n')
             attrindex = attrindex+1
         
      #save the file to disk 
      self.write_file()

    def load_graph_file  (self, filename):
      self.reset()
      fhandler = data_tree_file_io()
      if os.path.lexists(filename)==0:
         print('ERROR FILE %s DOES NOT EXIST'%filename)
         return None
      if os.path.lexists(filename):
         fhandler.readfilelines(filename)
         filelines = fhandler.filecontents_list
         #then do the load
         for line in filelines:
           words = line.split(' ')
           for word in words:
             if word != '':
               if word == '#GDATAFOLDER':
                 self.datafolder = words[1]
               if word == '#NODE':
                 #CREATE A NEW NODE OBJECT
                 newnode = node_base(words[1])
                 self.nodes.append(newnode)
                 self.walkbuffer.append(newnode)
                 #IF INPUT
                 if len(words)>2:
                   #print 'EXTRA INPUT '
                   #print words[3]
                   
                   #print 'DEBUG KDAG 659'
                   #print line
                   
                   inp_one = words[3]
                   #inp_one = words[5]
                   newnode.setinput (inp_one)
               if word == '#DATAFOLDER':
                newnode.datafolder = words[4]
               if word == '#COMMAND':
                commandtext = line.replace('#COMMAND','')
                newnode.COMMANDS.append(commandtext)
               if word == '#OUTPUT':
                newnode.hasoutputs = 1
                newnode.OUTPUTS.append( words[4] )
               if word == '#SET_PARENT':
                newnode.hasparents = 1
                newnode.parent.append( words[4] )
               if word == '#SET_CHILD':
                newnode.haschildren = 1
                childrentest = newnode.listchildren()
                checkint = 0
                for kid in childrentest:
                   if kid==words[4]:
                      checkint =checkint+1
                if checkint ==0:
                   newnode.children.append( words[4] )
                if checkint:
                   print('ERROR CLASH IN CHILDREN ')
               if word == '#SET_POS':
                newnode.hasxform = 1
                newnode.posx = float(words[4])
                newnode.posy = float(words[5])
                newnode.posz = float(words[6])
               if word == '#SET_ROT':
                newnode.hasxform = 1
                newnode.rotx = float(words[4])
                newnode.roty = float(words[5])
                newnode.rotz = float(words[6])
               if word == '#ATTR':
                newnode.ATTRS.append(    words[4])
                newnode.ATTRVALS.append( words[5])
               if word == '#SET_TYPE':
                  newnode.settype(    words[4])

    def walk_up_linked (self, rootnode, mode='obj'): #walk linked nodes instead of a hierarchy
        node = self.get(rootnode) #make sure you have the object
        if node == None:
            return 'ERROR '

        if mode =='name':
           #print '|'+(self.walkindent*'-')+'+'+ node.name
           self.walkbuffer.append( node.name )
           uplinked = node.listinputs()
           for uplinkednode in uplinked:
              self.walkindent=self.walkindent+1
              self.walk_up_linked(  self.get( uplinkednode ), mode )
              self.walkindent=self.walkindent-1

        if mode =='obj':
           #print '|'+(self.walkindent*'-')+'+'+ node.name
           self.walkbuffer.append( node )
           uplinked = node.listinputs()
           for uplinkednode in uplinked:
              self.walkindent=self.walkindent+1
              self.walk_up_linked(  self.get( uplinkednode ), mode )
              self.walkindent=self.walkindent-1
  
    def walklinked ( self, rootnode, mode='obj' ): #walk linked nodes instead of a hierarchy
         node = self.get(rootnode) #make sure you have the object
         if node == None:
            return 'ERROR '


         if mode =='name':
           print ('|'+(self.walkindent*'-')+'+'+ node.name)
           self.walkbuffer.append(node.name)
           childrentemp = node.listoutputs()
           for childd in childrentemp:
              self.walkindent=self.walkindent+1
              self.walklinked( self.get( childd ), mode)
              self.walkindent=self.walkindent-1

         if mode =='obj':
           print('|'+(self.walkindent*'-')+'+'+ node.name)
           self.walkbuffer.append(node)
           childrentemp = node.listoutputs()
           for childd in childrentemp:
              self.walkindent=self.walkindent+1
              self.walklinked( self.get( childd ), mode )
              self.walkindent=self.walkindent-1


         return self.walkbuffer

    def walk_path (self, mode='name'):
        out = ''
        if mode=='name':
            for x in self.walkbuffer:
                out = out + ('|'+ x )
            return out
        
        if mode=='filesys':
            for x in self.walkbuffer:
                out = out + ('/'+ x )
                print(x)
            return out
                           

    def trace_full_name(self, nod_obj, mode):
       """ trace the path up and convert to a name """
          
       self.setupwalk()
       if mode=='obj':
          self.walk_up( nod_obj, mode )
          self.walkbuffer.reverse()
          return self.walkbuffer
       if mode=='name':
          self.walk_up( nod_obj, mode )
          self.walkbuffer.reverse()   
          return self.walk_path( mode )
       if mode=='filesys':
          self.walk_up( nod_obj, 'name' )
          self.walkbuffer.reverse()
          return self.walk_path( mode ) 

    def walk(self, start_node, direction='down', mode='obj'): 
        """ interface to walk functions """
        
        self.setupwalk()
        if direction=='up':
            self.walk_up( start_node, mode)
        if direction=='down':
            self.walk_down( start_node, mode)
    
        return self.walkbuffer
         
    def walk_up(self, start_node, mode='obj'):    
        node = self.get( start_node )
        if node ==None:
          return None

        if mode =='name':
          #print '|'+(self.walkindent*'-')+'+'+ node.name
          self.walkbuffer.append( node.name )
          self.walk_depth_buffer.append( self.walkindent )
          parenttemp = node.parent #multiple []

          for parent in parenttemp:
              self.walkindent=self.walkindent+1
              self.walk_up(  self.get( parent ), mode )
              self.walkindent=self.walkindent-1

        if mode =='obj' or mode =='object':
          #print '|'+(self.walkindent*'-')+'+'+ node.name
          self.walkbuffer.append( node )
          self.walk_depth_buffer.append( self.walkindent )
          parenttemp = node.parent 
      
          for parent in parenttemp:
              self.walkindent=self.walkindent+1
              self.walk_up( self.get( parent ) , mode )
              self.walkindent=self.walkindent-1
 
   
    def walk_down(self, rootnode, mode='obj'):        #name string , or object mode
        node = self.get(rootnode) 
        
        #print('########## ', node.name)

        if node ==None:
          return None

        if mode =='obj' or mode =='object':
          print('|'+(self.walkindent*'-')+'+'+ node.name)
          
          self.walkbuffer.append(node)
          self.walk_depth_buffer.append( self.walkindent )
          childrentemp = node.children
      
          for childd in childrentemp:
              self.walkindent=self.walkindent+1
              self.walk_down( self.get( childd ), mode ) 
              self.walkindent=self.walkindent-1

        if mode =='name':
          #print '|'+(self.walkindent*'-')+'+'+ node.name
          self.walkbuffer.append( node.name )
          self.walk_depth_buffer.append( self.walkindent )
          childrentemp = node.children

          for childd in childrentemp:
              self.walkindent=self.walkindent+1
              self.walk_down( self.get( childd ), mode )
              self.walkindent=self.walkindent-1



    """
    def copy_branch( self, rootnode, newnodename, parent_to, spatial_offset, spatial_scale ):        #name string , or object mode
        #  optional parent (root) to node  
        mode = 'name'
        node = self.get(rootnode) #make sure you have the object
        if node ==None:
          return None

        if mode =='name':
          print '|'+(self.walkindent*'-')+'+'+ node.name
          self.walkbuffer.append(node.name)
          self.walk_depth_buffer.append(self.walkindent)
          childrentemp = node.children

          for childd in childrentemp:
              self.walkindent=self.walkindent+1
              self.walk(childd,mode)
              self.walkindent=self.walkindent-1

        if mode =='obj' or mode =='object':
          print '|'+(self.walkindent*'-')+'+'+ node.name
          self.walkbuffer.append(node)
          self.walk_depth_buffer.append(self.walkindent)
          childrentemp = node.children

          for childd in childrentemp:
              self.walkindent=self.walkindent+1
              self.walk(childd,mode)
              self.walkindent=self.walkindent-1
    """
    
    def link(self, node1, node2):
      #first check if both exist, then hook output of node1 to input node2
      if self.exists(node1)==0:
         print('ERROR NODE '+str(node1)+' DOES NOT EXIST ')
      if self.exists(node2)==0:
         print('ERROR NODE '+str(node2)+' DOES NOT EXIST ')
      ####
      for node in self.nodes:
        if node.name == node1.name:
            node.OUTPUTS.append(str(node2.name))
            node.hasoutputs =1
            
        if node.name == node2.name:
           node.INPUTS.append(str(node1.name))
           node.hasinputs =1
   
    def add(self, node):
        if self.exists(node):
            print('#kdag add: node %s already exists'% node)
            return False     
        if not self.exists(node):
           self.nodes.append(node)
           return True

    #SHOW ENTIRE CONTENTS OF A NODE 
    def dumpnode(self, nodename):
         nodtmp = self.get(nodename)
         if nodename ==None or nodtmp==None:
            return None
            
         print('#############################')
         print('#NODE NAME IS       '+ nodtmp.name)
         print('#NODE SHORTNAME IS  '+ nodtmp.shortname)
         print('#NODE TYPE IS : '+ nodtmp.nodetype)
         print('#HAS PARENTS '+ str( nodtmp.hasparents ) )
         if len(nodtmp.parent):
           parnts = ''
           for parent in nodtmp.parent:
                  parnts=parnts+parent
           print('#PARENT IS  %s'%parnts)       
         print('#HAS CHILDREN '+ str ( nodtmp.haschildren))
         if len(nodtmp.children):
           print('#CHILDREN ARE  ')
           for child in nodtmp.children:
                  print(child,' ')
         if len(nodtmp.INPUTS) :
           print('#NODE INPUTS ARE : ')
           print(nodtmp.INPUTS)
         if len(nodtmp.OUTPUTS):
           print('#NODE OUTPUTS ARE : ')
           print(nodtmp.OUTPUTS)
         if nodtmp.hasxform ==1:
           print('#TRANSLTION '+ ( str(nodtmp.posx)+' '+str(nodtmp.posy)+' '+str(nodtmp.posz) ) )
           print('#ROTATION '  + ( str(nodtmp.rotx)+' '+str(nodtmp.roty)+' '+str(nodtmp.rotz) ) )
         if nodtmp.listallattrs()!=None:
           print('#NODE ATTR DATA IS')
           print(nodtmp.listallattrs())
           
    def list_rootnodes(self, mode='obj'):
        out = []
        allnodes = self.listnodes()
        if allnodes==None:
            return None 
        for node in allnodes:
            parents = node.listparents()
            if (len(parents)==0 or parents==None ):
                if mode=='obj' or mode=='object':
                    out.append(node) 
                if mode=='name':
                    out.append(node.name) 
        return out
        
    def listnodes(self, mode='obj'): 
         out    = []
         for nodetmp in self.nodes:
             if mode == 'name':
                out.append(nodetmp.name)
             if mode == 'object' or mode =='obj':
                out.append(nodetmp)      
         return out
          
    def show_walk(self, mode='full'):
       if mode=='name':
          out = ''
          for each in self.walkbuffer:
             out= out+' '+each.name
          print(out)
          
       if mode=='full':
           print('#############################')
           print('#WALKBUFFER CONTAINS ' + str (len(self.walkbuffer)) + ' NODES ')       
           print(self.walkbuffer)
           
          
    #SHOW ENTIRE CONTENTS OF GRAPH (ALL NODES ,ATTRS ,ETC IN GRAPH )
    def show(self, mode='full'):
        if mode=='name':
            out = self.graphname + ': ' 
            for each in self.nodes:
                out= out+' '+each.name
            print(out)
              
        if mode=='full':
           print('#############################')
           print('#GRAPH NAME IS ' + self.graphname )
           #print('#GRAPH CONTAINS ' + str (len(self.nodes)) + ' NODES ')
           for each in self.nodes:
               print('\n')
               print('#############################')
               print('#NODE NAME IS       '+ each.name)
               #print('#NODE SHORTNAME IS  '+ each.shortname)
               #print('#NODE TYPE IS : '+ each.nodetype)
               #print( '#HAS PARENTS '+ str( each.hasparents ) )
           
               if len(each.parent):
                   parnts = ''
                   for parent in each.parent:
                          parnts=parnts+parent
                   print('#PARENT IS  %s'%parnts)

               #print ('#HAS CHILDREN '+ str ( each.haschildren))
               if len(each.children):
                   chlds = ''
                   for c in each.children:
                       chlds=chlds+' '+c
                   print('#CHILDREN ARE %s '%chlds)
               if len(each.INPUTS) :
                   print('#NODE INPUTS ARE : ')
                   print(each.INPUTS)
               if len(each.OUTPUTS):
                   print('#NODE OUTPUTS ARE : ')
                   print(each.OUTPUTS)
               # if each.hasxform ==1:
               #     print( '#TRANSLTION '+ ( str(each.posx)+' '+str(each.posy)+' '+str(each.posz) ) )
               #     print( '#ROTATION '  + ( str(each.rotx)+' '+str(each.roty)+' '+str(each.rotz) ) )
               # if each.listallattrs()!=None:
               #     print('#NODE ATTR DATA IS')
               #     print(each.listallattrs())

           print('#############################')


    def get_upstream_datafolder(self, node):
       upstream = node_base
       if node.hasinputs == 0:
           print('NO INPUTS')
       if node.hasinputs:
           inputnodnam = node.listinput()
           inputnod = self.get(inputnodnam)
           upstream  = inputnod
           return upstream.getdatafolder()
       return None #DEBUG

##-------------------------------------##
class node_base(object): 
    #def delete an attribute from a node
   
    def __init__( self, nam='dft_name' ):
        self.name         = nam         #FULL NAME OF NODE
        self.shortname    = nam         #CLIPPED PATH NAME  
        self.nodetype     = 'dft_base'  #WHAT IS IT?
        self.notes        = []          #TEXT ARRAY FOR NOTES
        self.parent       = []          #allow for multiple children and parents
        self.children     = []
        self.ATTRS        = []
        self.ATTRVALS     = []
        self.ATTRTYPS     = []         #this is optional , common idex to above attrs
        self.INPUTS       = []
        self.OUTPUTS      = []
        self.hasinputs    = 0
        self.hasoutputs   = 0          #IS IT CONNECTED TO ANY OTHER NODES?
        self.hasxform     = 0          #DOES IT HAVE A POSITION/ROTATION ?
        self.hasparents   = 0
        self.haschildren  = 0
        self.ORIGIN       = [0,0,0]
        self.posx         = 0
        self.posy         = 0
        self.posz         = 0
        self.rotx         = 0
        self.roty         = 0
        self.rotz         = 0
        self.scalex       = 0
        self.scaley       = 0
        self.scalez       = 0
    
    @property 
    def get_position(self):
        return (self.posx, self.posy, self.posz)

    @property 
    def get_rotation(self):
        return (self.rotx, self.roty, self.rotz)

    def set_position(self, xyz):
        self.hasxform = True
        self.posx = xyz[0]
        self.posy = xyz[1]
        self.posz = xyz[2]

    def set_rotation(self, xyz):
        self.hasxform = True      
        self.rotx = xyz[0]
        self.roty = xyz[1]
        self.rotz = xyz[2]

    def setname(self,name):
        self.name=str(name)

    def getname(self):
        return self.name

    def listparents(self):
        return self.parent

    def listchildren(self):
        return self.children

    def get_datafolder(self):
        return self.datafolder
   
    def settype(self,TYPE):
        self.nodetype = TYPE

    def gettype(self):
        return self.nodetype
   
    ## TODO ADD TYPES IN via self.ATTRTYPS[]
    def addattr (self, attrname, attrval):
        #DEBUG  TEMPORARILY DISABLED THIS CHECK   DEBUG
        #for at in self.ATTRS:
        #  if at == attrname:
        #     raise 'ERROR ALREADY EXISTS '
        self.ATTRS.append(    attrname )
        self.ATTRVALS.append(  attrval )

    #ADD FANCY SORTING LIKE BY ATTR VALUE <,>,=, etc
    #use compund datatypes ([attr,attrval]) ?
    def getdatafolder(self):
      return self.datafolder
   
    #print all attributes and thier values
    def listallattrs (self):
       if len(self.ATTRS)!=0:
         output = []
         #for x,y in map(None,self.ATTRS,self.ATTRVALS):
         #    output.append( [x,y] )
         for i,a in enumerate(self.ATTRS): 
             x = self.ATTRS[i]
             y = self.ATTRVALS[i]
             output.append( [x,y] )            
         return output
       else:
         return None

    def getallattrs (self):
        """ return the attributes names of a node """    
        output = []
        if len(self.ATTRS)!=0:
            for tempattr in self.ATTRS:
                output.append(tempattr)
        return output

    def getattrib (self, attrname):
      index=0    
      for attr in self.ATTRS :
         if attr == attrname:
           return self.ATTRVALS[index]
         index += 1  
      return None

    def setattrib (self, attrname, value):
      if self.hasattrib(attrname)==0:
         print('error setattr datagraph attr does not exist')
         
      if self.hasattrib(attrname):
        index = 0
        for attr in self.ATTRS :
           if attr == attrname:
             self.ATTRVALS[index]=value
           index = index+1

    def delete_attr_byname(self, attrtoremove):
      """  debug , delete an attribute from a node  """
      
      indexx = 0
      namebuffer = []
      valuebuffer = []
      #print 'removing attribute '+attrtoremove
      for attrname in self.ATTRS:
        if attrname !=attrtoremove :
           namebuffer.append( self.ATTRS[indexx] )
           valuebuffer.append( self.ATTRVALS[indexx] )
        indexx=indexx+1 #count remember

      self.ATTRS = namebuffer
      self.ATTRVALS = valuebuffer

    def hasattrib (self, attrname):
       """ check for attr exists ? """
       
       has = 0
       for attr in self.ATTRS:
          if attr ==attrname:
            has =1
            return has

    def hasAttribData (self):
       """ can be used as boolean or to count number of attrs """
       
       if self==None:
         return 0
       has = 0
       if (self.ATTRS)==0:
          return has
       if (self.ATTRS):
          has=has+1
          return has

    def delete_attrs(self):
       """ does this do anything usefull?? """
       
       self.ATTRS        = []
       self.ATTRVALS     = []
       self.ATTRTYPS     = []
       self.hasattrib = 0
        
    #returns [one,two]
    def listinput(self):
       if self.hasinputs==1:
          return self.INPUTS[0]
       else:
         return None

    def listinputs(self):
       out = []
       outone = ''
       outtwo = ''
       if self.hasinputs==1:
        for inputname in self.INPUTS:
          out.append(inputname)
        return out
       else:
         return None

    def listinputsindex(self, index):
       outvar = []
       count = 0

       for x in range( len( self.INPUTS) ):
         if len(self.INPUTS[x]) ==0 or len (self.INPUTS_TWO[x])==0:
           if x == index:
             print('ERROR INPUT VARS ARE EMPTY ')
             outvar.append( ['debug' ,'debug' ]) ##DEBUG NULL VALUE INSERTED
         else:
           if x == index:
             outvar.append( [str(self.INPUTS[x]) ,str(self.INPUTS_TWO[x]) ])
         return outvar

    def listoutput(self):
       return self.OUTPUTS[0]

    def listoutputs(self):
       out = []
       for nodeoutput in self.OUTPUTS:
         out.append(nodeoutput)
       return out

    def setinput(self, input1):
           self.hasinputs  = 1
           self.INPUTS =[]
           self.INPUTS.append(     str(input1) )

    def setinputs(self, input1, input2):
           self.hasinputs  = 1
           self.INPUTS =[]
           self.INPUTS.append(     str(input1) )
           self.INPUTS.append(     str(input2) )

    def setoutput(self, node):
       self.hasoutputs  = 1
       self.OUTPUTS =[]
       self.OUTPUTS.append( node )

    def show(self):
       print('NODE NAME IS '+ self.name)
       print('hasxform '+ str (self.hasxform ))
       ATTRS = self.listallattrs()
       if ATTRS !=None:
          for ATTR in ATTRS:
            print(ATTR)
       #
       #self.showTranslations()
       self.showParenting()
       
    def showTranslations(self):
        print('\n##############')
        print('NODE NAME IS '+ self.name)
        print('hasxform '+ str (self.hasxform ))
        print('##############')
        print('#TRANSLATION#')
        print(self.posx)
        print(self.posy)
        print(self.posz)
        print('##############')
        print('#ROTATION#')
        print(self.rotx)
        print(self.roty)
        print(self.rotz)
        print('##############\n\n')
        print('#ORIENTATION#')
        print('##############\n\n')
        print('#PARENTING#')
        print('#input/outputs')
        print(self.INPUTS)
        print(self.OUTPUTS)
        print('#has input/outputs')
        print(self.hasinputs)
        print(self.hasoutputs)
        print('#has xform')
        print(self.hasxform)
        print('#hasparents/haschildren')
        print(self.hasparents)
        print(self.haschildren)

    def showParenting(self):
        print('#hasparents/haschildren')
        print(self.hasparents)
        print(self.haschildren)

    def getxform(self):
       out = []
       out.append([self.posx, self.posy ,self.posz ])
       temp = str(out[0])
       temp.replace(' ','')
       return temp
       
    def getrotate(self):
       out = []
       out.append([self.rotx, self.roty ,self.rotz ])
       return out

    def getorient(self):
       out = []
       out.append([self.orientx, self.orienty ,self.orientz ])
       return out

    def scale_entire_graph(self, amount):
        allnodes = self.listnodes('object')
        for node in allnodes:
            node.scalexform(amount) 
                   
    def scalexform(self, amount):
        self.posx =  self.posx * amount
        self.posy =  self.posy * amount
        self.posz =  self.posz * amount
   
    def setxform(self, fbt):
       out = []
       if fbt==None:
         return

       buffer = [fbt[0],fbt[1] ,fbt[2]]
       out.append( buffer)
       self.xform(out)
       self.hasxform =1
       
    def setrotate(self, fbt):
       if fbt==None:
         return
       out = []

       out.append([fbt[0],fbt[1] ,fbt[2]])
       self.rotate(fbt)
       self.hasrotate =1

    def xform(self, fbt):
       """  translations are optional, expects an array of arrays - [ [xyz],[xyz] ] """
       
       self.hasxform = 1
       pos = fbt[0]
       self.posx = pos[0]
       self.posy = pos[1]
       self.posz = pos[2]
     
    def xform_xyz(self, xyz):
       self.hasxform = 1
       pos = xyz
       self.posx = xyz[0]
       self.posy = xyz[1]
       self.posz = xyz[2]

    def rotate(self, fbt):
       self.hasxform = 1 
       self.rotx = fbt[0]
       self.roty = fbt[1]
       self.rotz = fbt[2]

##-------------------------------------##
class data_tree_file_io (object):

    def __init__(self):
        self.filecontents_list =[]
        self.filecontents      =''

    def writefile_listlines (self, path, listvar):
        file_object = codecs.open ( path, "w", encoding='utf-8')
        #file_object.writelines (listvar)
        for line in listvar:
            #file_object.write(line+'\n')
            
            #DEBUG 
            #file_object.write(unicode(line)+'\n')    
            #print(line)        
            file_object.write(str(line)+'\n')
        file_object.close()
 
    def readfilelines (self, path):
        if os.path.lexists(path) == 0:
            print(path, " DOES NOT EXIST !! " )
        if os.path.lexists(path):
            if len(path) == 0:
                print("filename not defined")
            else:
                f = codecs.open( path,"r", encoding='utf-8')
                self.filecontents = f.readlines()
                for x in self.filecontents :
                    #lines = x.split(" ")
                    nonewline = x.split('\n')
                    self.filecontents_list.append(nonewline[0])

##-------------------------------------##
class basic_parser( object ):
    
    def __init__(self):
        self.pnodes               = [] 
        self.walkbuffer           = []
        self.tag_data             = []
        self.attr_known_tags      = []
        self.attr_unknown_tags    = [] 

    def all_between(self, line, left, right):
        regex_obj = re.compile('(\\' +left + '.+?\\'+right+')')    
        if regex_obj.findall(line):
            return regex_obj.findall(line)
            
    def all_between_curlys(self, line):
        regex_obj = re.compile("(\{.+?\})")
        if regex_obj.findall(line):
            return regex_obj.findall(line)
             
    def all_between_parens(self, line):
        regex_obj = re.compile("(\(.+?\))")
        if regex_obj.findall(line):
            return regex_obj.findall(line)
        
    def crop_leftmost(self, line, slpitchar):
        regex_obj = re.compile('^[^'+slpitchar+']*')
        if regex_obj.findall(line):
            return regex_obj.findall(line)

    def crop_rightmost(self, line, slpitchar):
       regex_obj = re.compile('[^'+slpitchar+']*$')
       if regex_obj.findall(line):
         return regex_obj.findall(line) 

    def invcrop_rightmost(self, line, slpitchar):
       regex_obj = re.compile('['+slpitchar+'].+')
       if regex_obj.findall(line):
         return regex_obj.findall(line)

    # string $rghtExp = (".+("+$char+")");
    def invcrop_leftmost(self, line, slpitchar):
       regex_obj = re.compile('.+['+slpitchar+']')
       if regex_obj.findall(line):
         return regex_obj.findall(line)

    def stringify(self, listobj):
      stringdata = ''
      for item in listobj:
         if item !='' and item !=' ':
           stringdata=(stringdata+' '+item)
      return stringdata

    def getnext(self, splitup, trigger):
       """ splits up a line by whitespace and retrieves next element """
       
       #splitup = line.split(' ') #or if you fancy a string , riddle me this
       count =0 
       out = ''
       kleanup = splitup.strip()
       ### 
       size = len(splitup)
       for item in kleanup:
          if item == trigger:
            if size>count+1:
             out= kleanup[count+1] 
             return out #return first occurance
          count =count+1
       #return out.lstrip(' ')
       return None
 
    def testparse (self, filename):
      #self.reset()
      fhandler = data_tree_file_io()
      if os.path.lexists(filename)==0:
         raise 'ERROR FILE DOES NOT EXIST '
      if os.path.lexists(filename):
         fhandler.readfilelines(filename)
         filelines = fhandler.filecontents_list

         #then do the load
         for line in filelines:
           words = line.split(' ')
           ###
           restofline = ''
           #go through each word on ecah line
           for word in words:
               #look for words  (tags) that begin with a #
               if word[:1]=='#':
                 #but not ones that begin with two ##'s or single #'s
                 if word[1:2]!='#' and word[1:2]!='':
                   tag =  word

                   ################################ 
                   if tag=='#NODE':
                     nextword = self.getnext(words,'#NODE')
                     newnode = node_base(nextword) 
                     self.pnodes.append(newnode)
                     #self.data_graph.nodes.append(newnode) 
                     #SET INPUTS#  newnode.setinput (inp_one)
                   ################################ 
                   if tag=='#SET_TYPE':
                     nextword = self.getnext(words,'#SET_TYPE')
                     newnode.settype(    nextword)

                   ################################ 
                   if tag=='#OUTPUT':
                     newnode.hasoutputs = 1
                     nextword = self.getnext(words,'#OUTPUT')
                     newnode.OUTPUTS.append( nextword )  
                   ################################ 
                   if tag=='#SET_PARENT':
                     newnode.hasparents = 1
                     nextword = self.getnext(words,'#SET_PARENT')
                     newnode.parent.append( nextword ) 
                   ################################ 
                   if tag=='#SET_CHILD':
                     newnode.haschildren = 1
                     nextword = self.getnext(words,'#SET_CHILD')
                     #add check for exsisting before adding DEBUG
                     childrentest = newnode.listchildren()
                     checkint = 0
                     #DEBUG WHAT DOES THIS DO, CHECK DUPES ?? July 26, 2011
                     for kid in childrentest:
                        if kid==nextword:
                          checkint =checkint+1
                     if checkint ==0:
                        newnode.children.append( nextword )
                     if checkint:
                        print('ERROR CLASH IN CHILDREN  ')

                   ################################
                   if tag=='#SET_POS':
                     newnode.hasxform = 1
                     #
                     posx = self.getnext(words,'#SET_POS')
                     newnode.posx=float(posx)
                     #
                     posy = self.getnext(words,posx)
                     newnode.posx=float(posy)
                     #
                     posz = self.getnext(words,posy)
                     newnode.posx=float(posz)

                   ################################
                   if tag=='#SET_ROT':
                     newnode.hasxform = 1
                     #
                     rotx = self.getnext(words,'#SET_ROT')
                     newnode.rotx=float(rotx)
                     #
                     roty = self.getnext(words,rotx)
                     newnode.roty=float(roty)
                     #
                     rotz = self.getnext(words,roty)
                     newnode.rotz=float(rotz)

                   ################################
                   if tag=='#ATTR':
                     nextword = self.getnext(words,'#ATTR')
                     tempstring = self.stringify(words)
                     attrval = self.crop_rightmost(tempstring,'ATTR')
                     ######################
                     newnode.ATTRS.append(    nextword)
                     tempp = attrval[0]
                     temp2=tempp.replace(nextword,'')
                     temp2=temp2.replace(' ','')
                     newnode.ATTRVALS.append( temp2 )







