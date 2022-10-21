/*
  Copyright (c) 2006-2015, The Regents of the University of California,
  through Lawrence Berkeley National Laboratory (subject to receipt of any
  required approvals from the U.S. Dept. of Energy) and the Paul Scherrer
  Institut (Switzerland).  All rights reserved.

  License: see file COPYING in top level of source distribution.
*/

#include <stdlib.h>

#include "H5hut.h"
#include "examples.h"

// name of output file
const char* fname = "example_setview.h5";

// H5hut verbosity level
const h5_int64_t h5_verbosity = H5_VERBOSE_DEFAULT;
const h5_int64_t h5_debug_mask = 0;

// we are going to write multiple consecutive blocks
//const h5_int64_t num_blocks = 32;
//const h5_int64_t num_particles_per_block = 1048576*8;

int
main (
	int argc,
	char* argv[]
        ){
	if (argc < 3) {
		fprintf (stderr, "Usage: %s <number_of_blocks> <sizeof_block>\n", argv[0]);
		exit (1);
	}

	char* endptr = NULL;
	long long n = strtoll (argv[1], &endptr, 10);
	if (*endptr != 0) {
		fprintf (stderr, "first argument (number of blocks) is not a unsigned integer!\n");
		exit (1);
	}
	if (n < 1) {
		fprintf (stderr, "first argument (number of block) must be >= 1!\n");
		exit (1);
	}
	if (n == LLONG_MAX) {
		fprintf (stderr, "first argument (number of block) to large!\n");
		exit (1);
	}
	h5_int64_t num_blocks = (h5_int64_t)n;
	n = strtoll (argv[2], &endptr, 10);
	if (*endptr != 0) {
		fprintf (stderr, "second argument (sizeof blocks) is not a unsigned integer!\n");
		exit (1);
	}
	if (n < 1024) {
		fprintf (stderr, "second argument (sizeof block) must be >= 1024!\n");
		exit (1);
	}
	if (n == LLONG_MAX) {
		fprintf (stderr, "second argument (sizeof block) to large!\n");
		exit (1);
	}
	h5_int64_t num_particles_per_block = (h5_int64_t)n;

        // initialize MPI & H5hut
        MPI_Init (&argc, &argv);
        MPI_Comm comm = MPI_COMM_WORLD;
        int comm_rank = 0;
        MPI_Comm_rank (comm, &comm_rank);
        H5AbortOnError ();
        H5SetVerbosityLevel (h5_verbosity);
        H5SetDebugMask (h5_debug_mask);

	h5_prop_t prop = H5CreateFileProp ();
	H5SetPropFileAlign (prop, 1048576*8);
	H5SetPropFileMPIOIndependent (prop, &comm);
	//H5SetPropFileMPIOCollective (prop, &comm);

        // open file and create first step
        h5_file_t file = H5OpenFile (fname, H5_O_WRONLY, prop);
	//H5PartSetChunkSize (file, 1048576*1);
        H5SetStep (file, 0);

	/*
	  If we want to write consecutive blocks, the 'view' can be defined
	  with H5PartSetview(). Otherwise we have to define the total number
	  of particles with H5PartSetNumParticles().
	 */
        const h5_int64_t offset = comm_rank * num_blocks * num_particles_per_block;
	H5PartSetView (
		file,
		offset,
		offset + num_blocks*num_particles_per_block -1);

        // write multiple consecutive blocks
        for (int i = 0; i < num_blocks; i++) {
		// create fake data
		//h5_int32_t data[num_particles_per_block];
		h5_int64_t *data;
		data = calloc (num_particles_per_block, sizeof(*data));
		for (int j = 0; j < num_particles_per_block; j++) {
			data[j] = j + i*num_particles_per_block + offset;
		}
		
                // set the "view" to select a subset of the dataset
                H5PartSetView (
                        file,
                        offset + i*num_particles_per_block,
                        offset + (i+1)*num_particles_per_block - 1);
                // write data
                H5PartWriteDataInt64 (file, "data", data);
        }

        // done
        H5CloseFile(file);
	MPI_Finalize();
        return 0;
}
