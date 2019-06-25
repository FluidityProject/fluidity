#TODO: use the rewrite to swap depending whether mpi is there are not..

import numpy
import sys

import H5hut as H5hut_rewrite
from H5hut import *


__h5hut_table__ = {}

#
# Rewrite functions for reading and writing datasets and attributes for
# the types
# - Int32
# - Int64
# - Float32
# - Float64
#
# These functions have either three or four arguments:
# - functions to read or write a dataset have three arguments
# - functions to read an attribute have three arguments
# - functions to write an attribute have four attributes (except the functions
#   to write a string attribute, which are not handled here)
#
# The third arguments is an array in all cases.
#
funcx = """
def {0}(*args, **kwargs): 
    import numpy
    if not len(args) == 3:
        print 'wrong number of arguments'
        return -2

    if type (args[2]) is numpy.ndarray:
        dtype = args[2].dtype
        size = args[2].size
    elif isinstance (args[2], (numpy.str)):
        dtype = numpy.str
        size = len (args[2])
    elif isinstance (args[2], (list, tuple)):
        dtype = numpy.dtype (type(args[2][0]))
        size = len (args[2])
    else:
        print 'third argument must be an array'
        return -2

    if 'Attrib' in '{0}' and 'Write' in '{0}' and not dtype == type(''):
        return __h5hut_table__['{0}'][dtype] (args[0], args[1], args[2], size, **kwargs)
    else:
        return __h5hut_table__['{0}'][dtype](*args, **kwargs)
"""

def __update_types__():
    import numpy
    import sys
    import inspect

    __h5hut_api__ = dir(H5hut_rewrite)

    __h5hut_types__ = [ "Int32", "Int64", "Float32", "Float64", "String" ]
    __numpy_types__ = [ numpy.dtype(numpy.int32), 
                        numpy.dtype(numpy.int64), 
                        numpy.dtype(numpy.float32), 
                        numpy.dtype(numpy.float64),
                        numpy.str]

    #
    # loop over all H5hut C-functions and above types.
    # Replace function if type is in name.
    #
    for __i__ in __h5hut_api__:
        for __j__ in __h5hut_types__:
            if __j__ in __i__:
                key = __i__[:-len(__j__)]

                numpy_key = __numpy_types__[__h5hut_types__.index(__j__)]
                func = H5hut_rewrite.__dict__[__i__]

                if key not in __h5hut_table__:
                    __h5hut_table__[key] = {}

                __h5hut_table__[key][numpy_key] = func

                del H5hut_rewrite.__dict__[__i__]
                del sys.modules["H5hut"].__dict__[__i__]

    for __new_func__ in __h5hut_table__.keys():
        types = __h5hut_table__[__new_func__]
        #args = inspect.getargspec(types[0])
        #print __new_func__
        result = funcx.format(__new_func__) 
        exec(result)

        #H5hut_rewrite.__dict__[__new_func__] = globals()[__new_func__]
        sys.modules["H5hut"].__dict__[__new_func__] = locals()[__new_func__]

__update_types__()
