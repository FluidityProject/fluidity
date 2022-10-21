/*
  Copyright (c) 2006-2015, The Regents of the University of California,
  through Lawrence Berkeley National Laboratory (subject to receipt of any
  required approvals from the U.S. Dept. of Energy) and the Paul Scherrer
  Institut (Switzerland).  All rights reserved.

  License: see file COPYING in top level of source distribution.
*/

/*
  Note:
  Running this example on more than one core is possible but the result
  might not be what you expect. Please read the HDF5 documentation about 
  the VFD core driver.
*/
#include "H5hut.h"
#include "examples.h"

#include <stdlib.h>

// name of input file
const char* fname = "example_core_vfd.h5";

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
        int comm_rank = 0;
        MPI_Comm_rank (comm, &comm_rank);
        H5AbortOnError ();
        H5SetVerbosityLevel (h5_verbosity);

        // open file and create first step
        h5_prop_t prop = H5CreateFileProp ();
        H5SetPropFileCoreVFD (prop, 0);
        h5_file_t file = H5OpenFile (fname, H5_O_RDONLY, prop);
        H5SetStep (file, 0);

	// with core cfd we read the hole file on all cores!
        h5_int64_t num_particles = H5PartGetNumParticles (file);
        printf ("[proc %d]: particles in view: %lld\n", comm_rank, (long long)num_particles);

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
        return H5_SUCCESS;
}

