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

// number of particles we are going to write per core
const ssize_t dim_x = 8;
const ssize_t dim_y = 8;
const ssize_t dim_z = 32;


static inline ssize_t idx (
	ssize_t i, ssize_t i_dim,
	ssize_t j, ssize_t j_dim,
	ssize_t k
	) {
	return (i + j*i_dim + k*i_dim*j_dim);
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

	// slice field in Z direction
	ssize_t dim = dim_z;
	ssize_t n_slices = dim / comm_size;
	ssize_t remaining_slices = dim % comm_size;
	ssize_t start;
	if (comm_rank < remaining_slices) {
		n_slices++;
		start = comm_rank * n_slices;
	} else {
		start = comm_rank*n_slices + remaining_slices;
	}
	ssize_t end = start + n_slices - 1;

	// setting our view
	ssize_t i_start = 0;
	ssize_t i_end = dim_x - 1;
	ssize_t j_start = 0;
	ssize_t j_end = dim_y - 1;
	ssize_t k_start = start;
	ssize_t k_end = end;
	
	// create fake data
	ssize_t i_dim = i_end - i_start + 1;
	ssize_t j_dim = j_end - j_start + 1;
        h5_int64_t data[(i_end-i_start+1) * (j_end-j_start+1) * (k_end-k_start+1)];
	for (int k = k_start; k <= k_end; k++) {
		for (int j = j_start; j <= j_end; j++) {
			for (int i = i_start; i <= i_end; i++) {
				ssize_t _idx = idx (i-i_start, i_dim, j-j_start, j_dim, k-k_start);
				data[_idx] = (h5_int64_t)idx (i, dim_x, j, dim_y, k);
			}
		}
	}

        // open file and create first step
        h5_file_t file = H5OpenFile (fname, H5_O_WRONLY, H5_PROP_DEFAULT);
        H5SetStep (file, 0); 

	// set view on data for this core
	H5Block3dSetView (file, i_start, i_end, j_start, j_end, k_start, k_end);
	
        // write data
        H5Block3dWriteScalarFieldInt64 (file, "data", data);

        // done
        H5CloseFile(file);
	MPI_Finalize ();
        return 0;
}
