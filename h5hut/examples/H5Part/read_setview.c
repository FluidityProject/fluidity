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
const char* fname = "example_setview.h5";

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

        // compute and set a "canonical" view:
	// all cores get almost the same number of particles
        h5_int64_t num_particles_total = H5PartGetNumParticles (file);
        h5_int64_t num_particles = num_particles_total / comm_size;
        h5_int64_t remainder = num_particles_total % comm_size;
        h5_int64_t start = comm_rank * num_particles;

        // adjust number of local particles
        if (comm_rank < remainder)
                num_particles++;

        // adjust start
        if (comm_rank < remainder) 
                start += comm_rank;
        else
                start += remainder;
        
        // Note:
	// setting end = start - 1 forces the selection of zero particles!
        h5_int64_t end = start + num_particles - 1;
        
        printf ("[proc %d]: set view to [%lld..%lld]\n", comm_rank, (long long)start, (long long)end);
        H5PartSetView (file, start, end);

	// read and print data
        h5_int32_t* data = calloc (num_particles, sizeof (*data));
        H5PartReadDataInt32 (file, "data", data);
        for (int i = 0; i < num_particles; i++) {
                printf ("[proc %d]: global index = %lld; local index = %d, value = %d\n",
                        comm_rank, (long long)(start+i), i, data[i]);
        }

	// cleanup
	free (data);
        H5CloseFile (file);
	MPI_Finalize ();
        return 0;
}
