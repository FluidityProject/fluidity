/*
  Copyright (c) 2006-2015, The Regents of the University of California,
  through Lawrence Berkeley National Laboratory (subject to receipt of any
  required approvals from the U.S. Dept. of Energy) and the Paul Scherrer
  Institut (Switzerland).  All rights reserved.

  License: see file COPYING in top level of source distribution.
*/

#include <H5hut.h>
#include <fcntl.h>
#include <stdlib.h>
#include <sys/stat.h>

#define FNAME           "attach_file.h5"
#define ATTACHMENT      "attach_file"
#define VERBOSITY       H5_VERBOSE_ERROR
#define DEBUG_MSK	H5_DEBUG_ALL

#include "examples.h"

int
main (
        int argc,
        char* argv[]
        ) {
	MPI_Init (&argc, &argv);

	H5SetErrorHandler (H5AbortErrorhandler);
        H5SetVerbosityLevel (VERBOSITY);
	H5SetDebugMask (DEBUG_MSK);

	h5_file_t f = H5OpenFile (FNAME, H5_O_WRONLY, H5_PROP_DEFAULT);
	H5AddAttachment (f, ATTACHMENT);
	H5CloseFile (f);
	f = H5OpenFile (FNAME, H5_O_RDWR, H5_PROP_DEFAULT);
	h5_ssize_t num_attachments = H5GetNumAttachments (f);
	printf ("Number of attachments: %lld\n", (long long int)num_attachments);
	int i;
	char fname[FILENAME_MAX];
	h5_size_t fsize;
	for (i=0; i < num_attachments; i++) {
		H5GetAttachmentInfoByIdx (f, i, fname, sizeof(fname), &fsize);
		printf (
			"Attachment %d: Name: %s, Size: %llu\n",
			i, fname, (long long unsigned)fsize);
		H5GetAttachment (f, fname);
		H5DeleteAttachment (f, fname);
	}
	H5CloseFile (f);
	return 0;
}
