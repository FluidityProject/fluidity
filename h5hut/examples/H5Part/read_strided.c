/*
  Copyright (c) 2006-2015, The Regents of the University of California,
  through Lawrence Berkeley National Laboratory (subject to receipt of any
  required approvals from the U.S. Dept. of Energy) and the Paul Scherrer
  Institut (Switzerland).  All rights reserved.

  License: see file COPYING in top level of source distribution.
*/

#include "H5hut.h"
#include "examples.h"

#include <stdlib.h>

// name of input file
const char* fname = "example_strided.h5";

// H5hut verbosity level
const h5_int64_t h5_verbosity = H5_VERBOSE_DEFAULT;
//const h5_int64_t h5_verbosity = H5_DEBUG_ALL;

int
main (
        int argc, char* argv[]
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

        // open file and go to first step
        h5_file_t file = H5OpenFile (fname, H5_O_RDONLY, H5_PROP_DEFAULT);
        H5SetStep (file, 0);

        // compute number of particles this process has to read
        h5_ssize_t num_particles_total = H5PartGetNumParticles (file);
        h5_ssize_t num_particles = num_particles_total / comm_size;
        if (comm_rank+1 == comm_size)
                num_particles += num_particles_total % comm_size;

	printf ("[proc %d]: particles in view: %lld\n", comm_rank, (long long)num_particles);
	printf ("[proc %d]: total number of particles: %lld\n",
		comm_rank, (long long unsigned)num_particles_total);
	
        // set number of particles and memory stride
        H5PartSetNumParticlesStrided (file, num_particles, 6);

        // read data
        h5_float64_t* data = calloc (6*num_particles, sizeof (*data));
        H5PartReadDataFloat64 (file, "x",  data+0);
        H5PartReadDataFloat64 (file, "y",  data+1);
        H5PartReadDataFloat64 (file, "z",  data+2);
        H5PartReadDataFloat64 (file, "px", data+3);
        H5PartReadDataFloat64 (file, "py", data+4);
        H5PartReadDataFloat64 (file, "pz", data+5);

	// print dataset "x"
        for (int i = 0; i < num_particles*6; i+=6) {
                printf ("[proc %d]: local index = %d, value = %6.3f\n",
                        comm_rank, i, data[i]);
        }

	// cleanup
	free (data);
        H5CloseFile (file);
	MPI_Finalize ();
        return 0;
}
