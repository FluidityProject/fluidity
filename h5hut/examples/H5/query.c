/*
  Copyright (c) 2006-2015, The Regents of the University of California,
  through Lawrence Berkeley National Laboratory (subject to receipt of any
  required approvals from the U.S. Dept. of Energy) and the Paul Scherrer
  Institut (Switzerland).  All rights reserved.

  License: see file COPYING in top level of source distribution.
*/

#include "H5hut.h"

#include "examples.h"

#define FNAME1          "example_file_attribs.h5"
#define FNAME2          "example_step_attribs.h5"

/*
  Due to the way types are defined in H5hut, we cannot use "switch() {}"
*/
const char*
type2string (
        h5_int64_t type
        ) {
        if (type == H5_FLOAT64_T) 
                return "H5_FLOAT64_T";
        if (type == H5_FLOAT32_T)
                return "H5_FLOAT32_T";
        if (type == H5_INT64_T)
                return "H5_INT64_T";
        if (type == H5_INT32_T)
                return "H5_INT32_T";
        if (type == H5_STRING_T)
                return "H5_STRING_T";
        return "unknown type";
}

static inline void
print_header (
        h5_int64_t n
        ) {
        if (n > 0) {
                printf ("\t%-6s %-30s %-15s %-10s\n", "idx", "name", "type", "dim");
        }
}

static inline void
print_query_result (
        h5_int64_t i,
        const char* const name,
        h5_int64_t type,
        h5_int64_t dim
        ) {
        printf ("\t%-6lld %-30s %-15s %-10lld\n", (long long)i, name, type2string(type), (long long)dim);
}

void
query_file_attribs (
        h5_int64_t f
        ) {
        char name[H5_MAX_NAME_LEN];
        h5_int64_t type;
        h5_size_t dim;

        // query # of file attributes
        h5_int64_t n = H5GetNumFileAttribs (f);
        printf ("\tNumber of file attributes: %lld\n", (long long)n);

        // output name and type of all file attribute
        print_header (n);
        for (h5_int64_t i = 0; i < n; i++) {
                H5GetFileAttribInfo (f, i, name, sizeof(name), &type, &dim);
                print_query_result (i, name, type, dim);
        }
        printf ("\n");
}

void
query_step_attribs (
        h5_int64_t f,
        h5_int64_t step
        ) {
        char name[H5_MAX_NAME_LEN];
        h5_int64_t type;
        h5_size_t dim;

        H5SetStep (f, step);

        // query # of step attributes
        h5_int64_t n = H5GetNumStepAttribs (f);
        printf ("\tNumber of step attributes: %lld\n", (long long)n);

        // output name and type of all step attribute
        print_header (n);
        for (h5_int64_t i = 0; i < n; i++) {
                H5GetStepAttribInfo (f, i, name, sizeof(name), &type, &dim); 
                print_query_result (i, name, type, dim);
        }
}
        
void
query_file (
        const char* const fname
        ) {
        printf ("\nFile: %s\n", fname);
        // if file properties is set to default, MPI_COMM_WORLD will be used
        h5_file_t f = H5OpenFile (fname, H5_O_RDONLY, H5_PROP_DEFAULT);

        // query and output file attribs
        query_file_attribs (f);

        // query # of steps, if > 0: go to first step, query and output step attribs
        h5_int64_t n = H5GetNumSteps (f);
        printf ("\tNumber of steps: %lld\n", (long long)n);

        if (n > 0) {
                // go to first step 
                h5_int64_t i = -1;
                while (!H5HasStep (f, ++i));

                query_step_attribs (f, i);
        }
	H5CloseFile (f);
}

int
main (
	int argc,
	char** argv
	) {
	MPI_Init (&argc, &argv);
        H5AbortOnError ();

        query_file (FNAME1);
        query_file (FNAME2);

	MPI_Finalize ();
	return 0;
}
