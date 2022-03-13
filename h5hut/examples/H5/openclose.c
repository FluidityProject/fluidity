/*
  Copyright (c) 2006-2018, The Regents of the University of California,
  through Lawrence Berkeley National Laboratory (subject to receipt of any
  required approvals from the U.S. Dept. of Energy) and the Paul Scherrer
  Institut (Switzerland).  All rights reserved.

  License: see file COPYING in top level of source distribution.
*/

#include "H5hut.h"

#include "examples.h"

#if !defined(VERBOSITY)
#define VERBOSITY       H5_VERBOSE_ERROR
#endif

/*
  Use H5_DEBUG_ALL to get full debug output
*/
#if !defined(DEBUGMASK)
//#define DEBUGMASK	H5_DEBUG_NONE
#define DEBUGMASK	H5_DEBUG_ALL
#endif

int
main (
	int argc,
	char** argv
	) {
	MPI_Comm comm = MPI_COMM_WORLD;

	MPI_Init (&argc, &argv);
	H5SetVerbosityLevel (VERBOSITY);
	H5SetDebugMask (DEBUGMASK);
        h5_prop_t prop = H5CreateFileProp ();
        H5SetPropFileMPIOCollective (prop, &comm);
        h5_file_t f = H5OpenFile ("testfile.h5", H5_O_WRONLY, prop);
        H5CloseProp (prop);
	H5CloseFile (f);
	H5Finalize ();
	MPI_Finalize ();
	
	return 0;
}
