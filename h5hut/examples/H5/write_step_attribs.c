/*
  Copyright (c) 2006-2015, The Regents of the University of California,
  through Lawrence Berkeley National Laboratory (subject to receipt of any
  required approvals from the U.S. Dept. of Energy) and the Paul Scherrer
  Institut (Switzerland).  All rights reserved.

  License: see file COPYING in top level of source distribution.
*/

#include "H5hut.h"
#include "examples.h"

#define FNAME           "example_step_attribs.h5"

#define ATTR_STRING     "StepAttrString"
#define ATTR_INT32      "StepAttrInt32"
#define ATTR_INT64      "StepAttrInt64"
#define ATTR_FLOAT32    "StepAttrFloat32"
#define ATTR_FLOAT64    "StepAttrFloat64"

#define asize(array)    (sizeof(array)/sizeof(array[0]))

int
main (
	int argc,
	char** argv
	) {
        char* string_value = "This is a string attribute bound to this step.";
        int32_t int32_value[] = {0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144};
        int64_t int64_value[] = {42, 43, 44, 45};
        h5_float32_t float32_value[] = {2.71828};
        h5_float64_t float64_value[] = {3.14159265358979323846264338327950288419716939937510};

	MPI_Init (&argc, &argv);

        H5AbortOnError ();

        // if file properties is set to default, MPI_COMM_WORLD will be used
        h5_file_t f = H5OpenFile (FNAME, H5_O_WRONLY, H5_PROP_DEFAULT);
        H5SetStep (f, 1);
        H5WriteStepAttribString  (f, ATTR_STRING,  string_value);
        H5WriteStepAttribInt32   (f, ATTR_INT32,   int32_value,   asize(int32_value));
        H5WriteStepAttribInt64   (f, ATTR_INT64,   int64_value,   asize(int64_value));
        H5WriteStepAttribFloat32 (f, ATTR_FLOAT32, float32_value, asize(float32_value));
        H5WriteStepAttribFloat64 (f, ATTR_FLOAT64, float64_value, asize(float32_value));

	H5CloseFile (f);
	MPI_Finalize ();
	return 0;
}
