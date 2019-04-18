/*
  Copyright (c) 2006-2017, The Regents of the University of California,
  through Lawrence Berkeley National Laboratory (subject to receipt of any
  required approvals from the U.S. Dept. of Energy) and the Paul Scherrer
  Institut (Switzerland).  All rights reserved.

  License: see file COPYING in top level of source distribution.
*/

#include "H5hut.h"
#include "examples.h"

#include <inttypes.h>

// name of output file
const char* fname = "example_write_field.h5";

// H5hut verbosity level
const h5_int64_t h5_verbosity = H5_VERBOSE_DEFAULT;

int
main (
        int argc,
        char* argv[]
        ){

        // initialize MPI & H5hut
        MPI_Init (&argc, &argv);
        MPI_Comm comm = MPI_COMM_WORLD;
        int comm_size = 1;
        MPI_Comm_size (comm, &comm_size);
	int comm_rank = 0;
        MPI_Comm_rank (comm, &comm_rank);
        H5AbortOnError ();
        H5SetVerbosityLevel (h5_verbosity);
	//H5SetDebugMask (-1);

        // open file and create first step
        h5_file_t file = H5OpenFile (fname, H5_O_RDONLY, H5_PROP_DEFAULT);
        H5SetStep (file, 0); 

	if (!H5BlockHasFieldData (file)) {
		goto done;
	}
	printf ("Has field data in step#0\n");

	if (!H5BlockHasField (file, "data")) {
		goto done;
	}
	printf ("Has field data with name 'data' in step#0\n");

	h5_size_t field_rank;
	h5_size_t field_dims[3];
	h5_size_t elem_rank;
	h5_int64_t type;
	H5BlockGetFieldInfoByName (
		file,
		"data",
		&field_rank,
		field_dims,
		&elem_rank,
		&type);
	char* stype = "unknown";
	if (type == H5_INT64_T) {
		stype = "H5_INT64_T";
	} else if (type == H5_INT32_T) {
		stype = "H5_INT32_T";
	} else if (type == H5_FLOAT64_T) {
		stype = "H5_FLOAT64_T";
	} else if (type == H5_FLOAT32_T) {
		stype = "H5_FLOAT32_T";
	} else if (type == H5_STRING_T) {
		stype = "H5_STRING_T";
	}
	printf ("rank of field:       %" PRId64 "\n", field_rank);
	printf ("dims of field:      [%" PRId64 ", %" PRId64 ", %" PRId64"]\n",
		field_dims[0], field_dims[1], field_dims[2]);
	printf ("rank of field data: %" PRId64 "\n", elem_rank);
	printf ("type of field data: '%s'\n", stype);
done:
        // done
        H5CloseFile(file);
	MPI_Finalize ();
        return 0;
}
