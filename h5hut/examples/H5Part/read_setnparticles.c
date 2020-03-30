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
const char* fname = "example_setnparticles.h5";

// H5hut verbosity level
const h5_int64_t h5_verbosity = H5_VERBOSE_DEFAULT;

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

	// set number of particles
        H5PartSetNumParticles (file, num_particles);

        // read and print data
        h5_int32_t* data = calloc (num_particles, sizeof (*data));
        H5PartReadDataInt32 (file, "data", data);
        for (int i = 0; i < num_particles; i++) {
                printf ("[proc %d]: local index = %d, value = %d\n",
                        comm_rank, i, data[i]);
        }

        // cleanup
	free (data);
        H5CloseFile (file);
	MPI_Finalize ();
        return 0;
}
