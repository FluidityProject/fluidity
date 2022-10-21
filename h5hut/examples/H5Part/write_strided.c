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
const char* fname = "example_strided.h5";

// H5hut verbosity level
const h5_int64_t h5_verbosity = H5_VERBOSE_DEFAULT;

// number of particles we are going to write per core
const h5_int64_t num_particles = 99;

int
main (
        int argc,
	char* argv[]
        ){
	
        // initialize MPI & H5hut
        MPI_Init (&argc, &argv);
        MPI_Comm comm = MPI_COMM_WORLD;
        int comm_rank = 0;
        MPI_Comm_rank (comm, &comm_rank);
        H5AbortOnError ();
        H5SetVerbosityLevel (h5_verbosity);

        // open file and create first step
        h5_file_t file = H5OpenFile (fname, H5_O_WRONLY, H5_PROP_DEFAULT);
        H5SetStep (file, 0);

        // create fake data
        h5_float64_t data[6*num_particles];
        h5_int64_t id[num_particles];
        for (int i = 0; i < num_particles; i++) {
                data [6*i + 0] = 0.0 + i + num_particles * comm_rank;
                data [6*i + 1] = 0.1 + i + num_particles * comm_rank;
                data [6*i + 2] = 0.2 + i + num_particles * comm_rank;
                data [6*i + 3] = 0.3 + i + num_particles * comm_rank;
                data [6*i + 4] = 0.4 + i + num_particles * comm_rank;
                data [6*i + 5] = 0.5 + i + num_particles * comm_rank;
                id [i] = i + num_particles * comm_rank;
        }

        // define number of items this processor will write and set the
        // in-memory striding
        H5PartSetNumParticlesStrided (file, num_particles, 6);

        // write strided data
        H5PartWriteDataFloat64 (file, "x",  data+0);
        H5PartWriteDataFloat64 (file, "y",  data+1);
        H5PartWriteDataFloat64 (file, "z",  data+2);
        H5PartWriteDataFloat64 (file, "px", data+3);
        H5PartWriteDataFloat64 (file, "py", data+4);
        H5PartWriteDataFloat64 (file, "pz", data+5);

        // disable striding to write the ID's
        H5PartSetNumParticles (file, num_particles);
        H5PartWriteDataInt64 (file, "id", id);

        // cleanup
	H5CloseFile (file);
	MPI_Finalize ();
	return 0;
}
