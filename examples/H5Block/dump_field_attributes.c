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
const char* fname = "example_field.h5";

// H5hut verbosity level
const h5_int64_t h5_verbosity = H5_VERBOSE_DEFAULT;

static inline void
dump_int64_attrib (
	h5_file_t file,
	const char* const field_name,
	const char* const attrib_name,
	h5_size_t attrib_nelems
	) {
	h5_int64_t attrib_data[attrib_nelems];
	H5BlockReadFieldAttribInt64 (
		file,
		field_name,
		attrib_name,
		attrib_data);
	printf ("Attribute: '%s'\n", attrib_name);
	printf ("    Type: H5_INT64_T\n");
	printf ("    Data: %" PRId64, attrib_data[0]);
	for (size_t i = 1; i < attrib_nelems; i++) {
		printf (", %" PRId64, attrib_data[i]);
	}
	printf ("\n");
}

static inline void
dump_int32_attrib (
	h5_file_t file,
	const char* const field_name,
	const char* const attrib_name,
	h5_size_t attrib_nelems
	) {
	h5_int32_t attrib_data[attrib_nelems];
	H5BlockReadFieldAttribInt32 (
		file,
		field_name,
		attrib_name,
		attrib_data);
	printf ("Attribute: '%s'\n", attrib_name);
	printf ("    Type: H5_INT32_T\n");
	printf ("    Data: %ld", (long)attrib_data[0]);
	for (size_t i = 1; i < attrib_nelems; i++) {
		printf (", %ld", (long)attrib_data[i]);
	}
	printf ("\n");
}

static inline void
dump_float64_attrib (
	h5_file_t file,
	const char* const field_name,
	const char* const attrib_name,
	h5_size_t attrib_nelems
	) {
	h5_float64_t attrib_data[attrib_nelems];
	H5BlockReadFieldAttribFloat64 (
		file,
		field_name,
		attrib_name,
		attrib_data);
	printf ("Attribute: '%s'\n", attrib_name);
	printf ("    Type: H5_FLOAT64_T\n");
	printf ("    Data: %2f", attrib_data[0]);
	for (size_t i = 1; i < attrib_nelems; i++) {
		printf (", %2f", attrib_data[i]);
	}
	printf ("\n");
}

static inline void
dump_float32_attrib (
	h5_file_t file,
	const char* const field_name,
	const char* const attrib_name,
	h5_size_t attrib_nelems
	) {
	h5_float32_t attrib_data[attrib_nelems];
	H5BlockReadFieldAttribFloat32 (
		file,
		field_name,
		attrib_name,
		attrib_data);
	printf ("Attribute: '%s'\n", attrib_name);
	printf ("    Type: H5_FLOAT32_T\n");
	printf ("    Data: %2f", attrib_data[0]);
	for (size_t i = 1; i < attrib_nelems; i++) {
		printf (", %2f", attrib_data[i]);
	}
	printf ("\n");
}

static inline void
dump_string_attrib (
	h5_file_t file,
	const char* const field_name,
	const char* const attrib_name,
	h5_size_t attrib_nelems
	) {
	char attrib_data[attrib_nelems];
	H5BlockReadFieldAttribString (
		file,
		field_name,
		attrib_name,
		attrib_data);
	printf ("Attribute: '%s'\n", attrib_name);
	printf ("    Type: H5_STRING_T\n");
	printf ("    Data: %s", attrib_data);
}

void
dump_attrib (
	h5_file_t file,
	const char* const field_name,
	const char* const attrib_name,
	h5_int64_t attrib_type,
	h5_size_t attrib_nelems
	) {
	if (attrib_type == H5_INT64_T) {
		dump_int64_attrib (file, field_name, attrib_name, attrib_nelems);
	} else if (attrib_type == H5_INT32_T) {
		dump_int32_attrib (file, field_name, attrib_name, attrib_nelems);
	} else if (attrib_type == H5_FLOAT64_T) {
		dump_float64_attrib (file, field_name, attrib_name, attrib_nelems);
	} else if (attrib_type == H5_FLOAT32_T) {
		dump_float32_attrib (file, field_name, attrib_name, attrib_nelems);
	} else if (attrib_type == H5_STRING_T) {
		dump_string_attrib (file, field_name, attrib_name, attrib_nelems);
	}
}

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
        h5_file_t file = H5OpenFile (fname, H5_O_RDWR, H5_PROP_DEFAULT);
        H5SetStep (file, 0); 

	// test wheter field exists
	const char field_name[] = "data";
	if (!H5BlockHasField (file, field_name)) {
		printf ("Doesn't have field data with name 'data' in step#0\n");
		goto done;
	}

	// get number of attributes attached to field
	h5_ssize_t n_attribs = H5BlockGetNumFieldAttribs (
		file,
		field_name);
	printf ("Field has %" PRId64 " attributes attached.\n",
		n_attribs);

	// dump all attached attributes
	for (h5_size_t i = 0; i < n_attribs; i++) {
		char attrib_name[128];
		h5_size_t sizeof_attrib_name = sizeof (attrib_name);
		h5_int64_t attrib_type;
		h5_size_t attrib_nelems;
		H5BlockGetFieldAttribInfo (
				file,
				field_name,
				i,
				attrib_name, sizeof_attrib_name,
				&attrib_type,
				&attrib_nelems);

		dump_attrib (file,
			     field_name,
			     attrib_name, attrib_type, attrib_nelems);
	}
done:
        // done
        H5CloseFile(file);
	MPI_Finalize ();
        return 0;
}
