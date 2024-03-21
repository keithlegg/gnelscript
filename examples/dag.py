

from gnelscript.pygfx.dag_ops import *

from gnelscript import GLOBAL_PROJ






##------------------------------------##
def dag_sort():
    #failed test of spatial indexing with ye olde DAG tool 
    
    dag = data_graph() 

    nn = node_base('nod2')
    nn.addattr('sp', 3) 
    nn.addattr('ep', 4)
    nn.addattr('d', 0)  
    dag.add(nn)

    nn = node_base('nod4')
    nn.addattr('sp', 5) 
    nn.addattr('ep', 4)
    nn.addattr('d', 0)  
    dag.add(nn)

    nn = node_base('nod1')
    nn.addattr('sp', 3) 
    nn.addattr('ep', 2) 
    nn.addattr('d', 0) 
    dag.add(nn)

    nn = node_base('nod3')
    nn.addattr('sp', 1) 
    nn.addattr('ep', 2)
    nn.addattr('d', 0) 
    dag.add(nn)

    #dag.show()

    def flip(node): 
        s = node.getattrib('sp')
        e = node.getattrib('ep')
        node.setattrib('ep', s)
        node.setattrib('sp', e)

    def snod(nod1,nod2):
        print('%s s %s e %s '%(nod1.name, nod1.getattrib('sp'), nod1.getattrib('ep') ) )
        print('%s s %s e %s '%(nod2.name, nod2.getattrib('sp'), nod2.getattrib('ep') ) )    
        print('-')

    #for x in dag.nodes:
    for f in dag.nodes: 
        for nod in dag.nodes: 
            #snod(f,nod)
            
            #if f.getattrib('d')==0:
            if f.getattrib('sp') == nod.getattrib('ep'):
                print('##  match s>e %s s  %s e'%(f.name,nod.name) ) 
                nod.setattrib('d',1)
                dag.parent(nod, f) 

            #if nod.getattrib('d')==0:
            if f.getattrib('ep') == nod.getattrib('sp'):
                print('##  match e<s %s e  %s s'%(f.name,nod.name) ) 
                f.setattrib('d',1)
                dag.parent(f, nod)

            if f is not nod and nod.getattrib('d')==0: 
                #if no match - flip the order of input/outputs  
                #if nod.getattrib('d')==0:
                if f.getattrib('sp') == nod.getattrib('sp') or f.getattrib('ep') == nod.getattrib('ep'):
                    print('##  flip %s e  %s s'%(f.name,nod.name) ) 
                    flip(nod)

    ##WONDER WALK - BAMBOO STEAMER
    start  = dag.find_any_attr(['sp','ep'],1)

    print('-->  start %s'%start.name)

    nodes = dag.walk(start)
    for n in nodes:
        print(n.name, n.getattrib('sp'), n.getattrib('ep'), n.getattrib('d') )


        
def stackops():
    dg = new_stack()

    # stack commands from the interface
    dg.adn('frank', 'beer').adn('frankx', 'beer').adn('franky', 'beer').prnt('frankx','frank').prnt('franky','frank')

    #OR run commands directly on the graph  
    nodes = dg.DG.list_children_obj('frank')
    print(nodes)

    dg.DG.save_graph_file(GLOBAL_PROJ)



def example_stack():
    dg = new_stack(GLOBAL_PROJ)
    print(dg)


def example_stack2():
    dag = data_graph()
    dg = new_stack(GLOBAL_PROJ)

    for i,nod in enumerate(dg.DG.walk('frank')):
        print(' '*i,nod.name)
        print(' '*i,'|')


def example_stack3():
    dg = new_stack()
    print(dg)
        





#print( dg.DG.getnodes() )
#print( dg.DG.trace_full_name(dg.DG.get('frankx') ))
# stack commands from the interface

#dg = new_stack()
#dg.adn('bob', 'beer').adn('jones', 'beer').adn('chet', 'beer').prnt('chet','jones').prnt('jones','bob')
#print(dg.DG.walk('bob'))

#OR run commands directly on the graph  
#$nodes =dg.DG.save_graph_file('%s/foo.klsys'%GLOBAL_PROJ)

#kl = new_stack()





