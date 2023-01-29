

from gnelscript.pygfx.dag_ops import *

from gnelscript import GLOBAL_PROJ





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





