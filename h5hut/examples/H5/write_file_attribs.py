#!/usr/bin/env python

from H5hut import *
import numpy as np

FNAME =        "example_file_attribs.h5"

ATTR_STRING =  "FileAttrString"
ATTR_INT32 =   "FileAttrInt32"
ATTR_INT64 =   "FileAttrInt64"
ATTR_FLOAT32 = "FileAttrFloat32"
ATTR_FLOAT64 = "FileAttrFloat64"

string_value = "This is a string attribute attached to the file."
int32_value = np.array ([0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144], dtype='int32')
int64_value = np.array ([42, 43, 44, 45], dtype='int64')
float32_value = np.array ([2.71828, ], dtype='float32')
float64_value = np.array ([3.14159265358979323846264338327950288419716939937510,],
                            dtype='float64')

f = H5OpenFile (FNAME, H5_O_WRONLY, H5_PROP_DEFAULT)

H5WriteFileAttrib (f, ATTR_STRING,  string_value)
H5WriteFileAttrib (f, ATTR_INT32,   int32_value)
H5WriteFileAttrib (f, ATTR_INT64,   int64_value)
H5WriteFileAttrib (f, ATTR_FLOAT32, float32_value)
H5WriteFileAttrib (f, ATTR_FLOAT64, float64_value)

H5CloseFile (f)
