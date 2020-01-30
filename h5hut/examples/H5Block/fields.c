/*
  Copyright (c) 2006-2017, The Regents of the University of California,
  through Lawrence Berkeley National Laboratory (subject to receipt of any
  required approvals from the U.S. Dept. of Energy) and the Paul Scherrer
  Institut (Switzerland).  All rights reserved.

  License: see file COPYING in top level of source distribution.
*/

#include <H5hut.h>
#include "examples.h"

#define XSIZE		8
#define YSIZE		8
#define ZSIZE		8
#define DATASIZE	XSIZE*YSIZE*ZSIZE

#define VERBOSITY       H5_VERBOSE_DEFAULT
#define FNAME           "example_fields.h5"

int
main (
        int argc,
        char** argv
        ) {
        h5_int64_t verbosity = VERBOSITY;

        h5_float64_t ex[DATASIZE];
        h5_float64_t ey[DATASIZE];
        h5_float64_t ez[DATASIZE];
        h5_float64_t q[DATASIZE];

        // initialize MPI & H5hut
        int comm_rank = 0;
        int comm_size = 1;
        MPI_Init (&argc, &argv);
        MPI_Comm comm = MPI_COMM_WORLD;
        MPI_Comm_rank (comm, &comm_rank);
        MPI_Comm_size (comm, &comm_size);

        H5AbortOnError ();
        H5SetVerbosityLevel (verbosity);

        // open file and go to step#0
        h5_file_t file = H5OpenFile (FNAME, H5_O_RDONLY, H5_PROP_DEFAULT);
        H5SetStep (file, 0);

        H5Block3dSetView (file,
                           comm_rank*XSIZE, (comm_rank+1)*XSIZE - 1,
                           0, YSIZE - 1,
                           0, ZSIZE - 1);
        H5Block3dWriteScalarFieldFloat64(file, "Q", q);
        H5Block3dWriteVector3dFieldFloat64(file, "E", ex, ey, ez);
        H5CloseFile(file);

        MPI_Finalize();
        return H5_SUCCESS;
}

