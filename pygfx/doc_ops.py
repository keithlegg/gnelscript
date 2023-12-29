
""" simple tool to self document """


#import inspect 
#func.__code__.co_argcount
#func.__code__.co_varnames
#func2.__defaults__


OMIT_PY_OBJS = ['__class__','__delattr__','__dict__','__dir__','__doc__','__eq__','__format__',
                '__ge__','__getattribute__','__gt__','__hash__','__init__','__init_subclass__',
                '__le__','__lt__','__module__','__ne__','__new__','__reduce__','__reduce_ex__',
                '__repr__','__setattr__','__sizeof__','__str__','__subclasshook__','__weakref__']




def gn_dir(module):
    print('##-----------------------##')
    print('gnelscript module: %s'%module  ) 
    indent = ''
    for d in dir(module):
        indent = '    '
        if d not in OMIT_PY_OBJS:
            print('%s%s'%(indent,d))



def gn_dir_type(module):
    for a in module.__dict__:
        print(a)

