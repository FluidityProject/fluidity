/*
  Copyright (c) 2006-2015, The Regents of the University of California,
  through Lawrence Berkeley National Laboratory (subject to receipt of any
  required approvals from the U.S. Dept. of Energy) and the Paul Scherrer
  Institut (Switzerland).  All rights reserved.

  License: see file COPYING in top level of source distribution.
*/

#include "H5hut.h"
#include "examples.h"

// name of output file
const char* fname = "example_field.h5";

// H5hut verbosity level
const h5_int64_t h5_verbosity = H5_VERBOSE_DEFAULT;

int
main (
        int argc,
        char* argv[]
        ){

        // initialize MPI & H5hut
        MPI_Init (&argc, &argv);
        H5AbortOnError ();
        H5SetVerbosityLevel (h5_verbosity);
	//H5SetDebugMask (-1);

        // open file and create first step
        h5_file_t file = H5OpenFile (fname, H5_O_RDWR, H5_PROP_DEFAULT);
        H5SetStep (file, 0); 

	if (!H5BlockHasField (file, "data")) {
		printf ("Doesn't have field data with name 'data' in step#0\n");
		goto done;
	}

	h5_int32_t attrib[1] = { 42 };
	H5BlockWriteFieldAttribInt32 (
		file,
		"data",
		"The answer",
		attrib,
		sizeof (attrib) / sizeof (*attrib));
	h5_float64_t origin[3] = { 0.0, 0.0, 1.0 };
	H5Block3dSetFieldOrigin (
		file,
		"data",
		origin[0], origin[1], origin[2]);
	h5_float64_t spacing[3] = { 1.0, 2.0, 3.0 };
	H5Block3dSetFieldSpacing (
		file,
		"data",
		spacing[0], spacing[1], spacing[2]);
done:
        // done
        H5CloseFile(file);
	MPI_Finalize ();
        return 0;
}
