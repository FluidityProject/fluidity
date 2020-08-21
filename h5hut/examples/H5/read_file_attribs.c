/*
  Copyright (c) 2006-2016, The Regents of the University of California,
  through Lawrence Berkeley National Laboratory (subject to receipt of any
  required approvals from the U.S. Dept. of Energy) and the Paul Scherrer
  Institut (Switzerland).  All rights reserved.

  License: see file COPYING in top level of source distribution.
*/

#include "H5hut.h"
#include "examples.h"

#include <stdlib.h>

#define FNAME           "example_file_attribs.h5"
#define VERBOSITY       H5_VERBOSE_ERROR
#define DEBUG_MSK	0

#define ATTR_STRING     "FileAttrString"
#define ATTR_INT32      "FileAttrInt32"
#define ATTR_INT64      "FileAttrInt64"
#define ATTR_FLOAT32    "FileAttrFloat32"
#define ATTR_FLOAT64    "FileAttrFloat64"

int
main (
	int argc,
	char** argv
	) {
	MPI_Init (&argc, &argv);
        H5AbortOnError ();
        H5SetVerbosityLevel (VERBOSITY);
	H5SetDebugMask (DEBUG_MSK);

        // if file properties is set to default, MPI_COMM_WORLD will be used
        h5_file_t f = H5OpenFile (FNAME, H5_O_RDONLY, H5_PROP_DEFAULT);
        h5_size_t len;

        H5GetFileAttribInfoByName (f, ATTR_STRING, NULL, &len);
        char* attr_string = (char*)malloc (len+1);
        H5ReadFileAttribString (f, ATTR_STRING, attr_string);
        printf ("%s: %s\n", ATTR_STRING, attr_string);
        free (attr_string);

        H5GetFileAttribInfoByName (f, ATTR_INT32, NULL, &len);
        int32_t* attr_int32 = (int32_t*)malloc (sizeof(*attr_int32)*len);
        H5ReadFileAttribInt32 (f, ATTR_INT32, attr_int32);
        printf ("%s:", ATTR_INT32);
        for (int i = 0; i < len; i++) {
                printf (" %d", attr_int32[i]);
        }
        printf ("\n");
        free (attr_int32);

        H5GetFileAttribInfoByName (f, ATTR_INT64, NULL, &len);
        int64_t* attr_int64 = (int64_t*)malloc (sizeof(*attr_int64)*len);
        H5ReadFileAttribInt64 (f, ATTR_INT64, attr_int64);
        printf ("%s:", ATTR_INT64);
        for (int i = 0; i < len; i++) {
                printf (" %lld", (long long int)attr_int64[i]);
        }
        printf ("\n");
        free (attr_int64);

        H5GetFileAttribInfoByName (f, ATTR_FLOAT32, NULL, &len);
        h5_float32_t* attr_float32 = (h5_float32_t*)malloc (sizeof(*attr_float32)*len);
        H5ReadFileAttribFloat32 (f, ATTR_FLOAT32, attr_float32);
        printf ("%s:", ATTR_FLOAT32);
        for (int i = 0; i < len; i++) {
                printf (" %f", attr_float32[i]);
        }
        printf ("\n");
        free (attr_float32);

        H5GetFileAttribInfoByName (f, ATTR_FLOAT64, NULL, &len);
        h5_float64_t* attr_float64 = (h5_float64_t*)malloc (sizeof(*attr_float64)*len);
        H5ReadFileAttribFloat64 (f, ATTR_FLOAT64, attr_float64);
        printf ("%s:", ATTR_FLOAT64);
        for (int i = 0; i < len; i++) {
                printf (" %f", attr_float64[i]);
        }
        printf ("\n");
        free (attr_float64);

        // cleanup
	H5CloseFile (f);
	MPI_Finalize ();
	return 0;
}
