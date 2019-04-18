/*
  Copyright (c) 2006-2017, The Regents of the University of California,
  through Lawrence Berkeley National Laboratory (subject to receipt of any
  required approvals from the U.S. Dept. of Energy) and the Paul Scherrer
  Institut (Switzerland).  All rights reserved.

  License: see file COPYING in top level of source distribution.
*/

#include "H5hut.h"
#include "examples.h"

#include <stdlib.h>

#define NPROCS  8
#define DEFAULT_VERBOSITY       H5_VERBOSE_DEFAULT

struct H5BlockPartition {
	h5_int64_t	i_start;
	h5_int64_t	i_end;
	h5_int64_t	j_start;
	h5_int64_t	j_end;
	h5_int64_t	k_start;
	h5_int64_t	k_end;
};

struct H5BlockPartition Layout1[1] = {
	{ 0, 63, 0, 63, 0, 511 }
};

struct H5BlockPartition Layout8[8] = {
		{  0,63,  0,63,   0, 63},
		{  0,63,  0,63,  64,127},
		{  0,63,  0,63, 128,191},
		{  0,63,  0,63, 192,255},
		{  0,63,  0,63, 256,319},
		{  0,63,  0,63, 320,383},
		{  0,63,  0,63, 384,447},
		{  0,63,  0,63, 448,511}
};

struct H5BlockPartition Layout8G[8] = {
		{  0,63,  0,63,   0, 64},
		{  0,63,  0,63,  63,128},
		{  0,63,  0,63, 127,192},
		{  0,63,  0,63, 191,256},
		{  0,63,  0,63, 255,320},
		{  0,63,  0,63, 319,384},
		{  0,63,  0,63, 383,448},
		{  0,63,  0,63, 447,511}
};

struct H5BlockPartition Layout16[16] = {
		{  0,63,  0,31,   0, 63},
		{  0,63, 32,63,   0, 63},
		{  0,63,  0,31,  64,127},
		{  0,63, 32,63,  64,127},
		{  0,63,  0,31, 128,191},
		{  0,63, 32,63, 128,191},
		{  0,63,  0,31, 192,255},
		{  0,63, 32,63, 192,255},
		{  0,63,  0,31, 256,319},
		{  0,63, 32,63, 256,319},
		{  0,63,  0,31, 320,383},
		{  0,63, 32,63, 320,383},
		{  0,63,  0,31, 384,447},
		{  0,63, 32,63, 384,447},
		{  0,63,  0,31, 448,511},
		{  0,63, 32,63, 448,511}
};

struct H5BlockPartition Layout16G[16] = {
		{  0,63,  0,32,   0, 64},
		{  0,63, 31,63,   0, 64},
		{  0,63,  0,32,  63,128},
		{  0,63, 31,63,  63,128},
		{  0,63,  0,32, 127,192},
		{  0,63, 31,63, 127,192},
		{  0,63,  0,32, 191,256},
		{  0,63, 31,63, 191,256},
		{  0,63,  0,32, 255,320},
		{  0,63, 31,63, 255,320},
		{  0,63,  0,32, 319,384},
		{  0,63, 31,63, 319,384},
		{  0,63,  0,32, 383,448},
		{  0,63, 31,63, 383,448},
		{  0,63,  0,32, 447,511},
		{  0,63, 31,63, 447,511}
};


struct H5BlockPartition Layout32[32] = {
		{  0,31,  0,31,   0, 63},
		{  0,31, 32,63,   0, 63},
		{ 32,63,  0,31,   0, 63},
		{ 32,63, 32,63,   0, 63},
		{  0,31,  0,31,  64,127},
		{  0,31, 32,63,  64,127},
		{ 32,63,  0,31,  64,127},
		{ 32,63, 32,63,  64,127},
		{  0,31,  0,31, 128,191},
		{  0,31, 32,63, 128,191},
		{ 32,63,  0,31, 128,191},
		{ 32,63, 32,63, 128,191},
		{  0,31,  0,31, 192,255},
		{  0,31, 32,63, 192,255},
		{ 32,63,  0,31, 192,255},
		{ 32,63, 32,63, 192,255},
		{  0,31,  0,31, 256,319},
		{  0,31, 32,63, 256,319},
		{ 32,63,  0,31, 256,319},
		{ 32,63, 32,63, 256,319},
		{  0,31,  0,31, 320,383},
		{  0,31, 32,63, 320,383},
		{ 32,63,  0,31, 320,383},
		{ 32,63, 32,63, 320,383},
		{  0,31,  0,31, 384,447},
		{  0,31, 32,63, 384,447},
		{ 32,63,  0,31, 384,447},
		{ 32,63, 32,63, 384,447},
		{  0,31,  0,31, 448,511},
		{  0,31, 32,63, 448,511},
		{ 32,63,  0,31, 448,511},
		{ 32,63, 32,63, 448,511}
};

struct H5BlockPartition Layout32G[32] = {
		{  0,32,  0,32,   0, 64},
		{  0,32, 31,63,   0, 64},
		{ 31,63,  0,32,   0, 64},
		{ 31,63, 31,63,   0, 64},
		{  0,32,  0,32,  63,128},
		{  0,32, 31,63,  63,128},
		{ 31,63,  0,32,  63,128},
		{ 31,63, 31,63,  63,128},
		{  0,32,  0,32, 127,192},
		{  0,32, 31,63, 127,192},
		{ 31,63,  0,32, 127,192},
		{ 31,63, 31,63, 127,192},
		{  0,32,  0,32, 191,256},
		{  0,32, 31,63, 191,256},
		{ 31,63,  0,32, 191,256},
		{ 31,63, 31,63, 191,256},
		{  0,32,  0,32, 255,320},
		{  0,32, 31,63, 255,320},
		{ 31,63,  0,32, 255,320},
		{ 31,63, 31,63, 255,320},
		{  0,32,  0,32, 319,384},
		{  0,32, 31,63, 319,384},
		{ 31,63,  0,32, 319,384},
		{ 31,63, 31,63, 319,384},
		{  0,31,  0,31, 383,448},
		{  0,31, 31,63, 383,448},
		{ 31,63,  0,31, 383,448},
		{ 31,63, 31,63, 383,448},
		{  0,32,  0,32, 447,511},
		{  0,32, 31,63, 447,511},
		{ 31,63,  0,32, 447,511},
		{ 31,63, 31,63, 447,511}
};

#define _calc_index( i, i_dims, j, j_dims, k, k_dims ) \
		(i + j*i_dims + k*i_dims*j_dims)

static h5_int64_t
_write_data (
	h5_file_t f,
	int myproc,
        struct H5BlockPartition* view
	) {

	h5_int64_t i, j, k, idx;
	h5_int64_t herr;
	h5_float64_t *data;
	h5_int64_t i_dims = view->i_end - view->i_start + 1;
	h5_int64_t j_dims = view->j_end - view->j_start + 1;
	h5_int64_t k_dims = view->k_end - view->k_start + 1;

	printf ( "Writing scalar field data to step #%lld\n", (long long)H5GetStep (f));

	data = (h5_float64_t*)malloc ( i_dims * j_dims * k_dims * sizeof (h5_float64_t) );
	for ( i = 0; i < i_dims; i++ ) {
		for ( j = 0; j < j_dims; j++ ) {
			for ( k = 0; k < k_dims; k++ ) {
				idx = _calc_index (
					i, i_dims,
					j, j_dims,
					k, k_dims );
				*(data + idx) = k
					+ 1000*j
					+ 100000*i
					+ 10000000*myproc;
			}
		}
	}

	herr = H5Block3dSetView (
		f,
		view->i_start, view->i_end,
		view->j_start, view->j_end,
		view->k_start, view->k_end );
	if ( herr < 0 ) return herr;

	herr = H5Block3dWriteScalarFieldFloat64 ( f, "TestField", data );
	if ( herr < 0 ) return herr;

	free ( data );

	return 1;
}

static h5_int64_t
_write_attributes (
	h5_file_t f,
	const int myproc
	) {
        printf ("Writing attributes to field '%s' in step #%lld\n",
		"TestField", (long long)H5GetStep (f));
	h5_int64_t herr = H5BlockWriteFieldAttribString (
		f,
		"TestField",
		"TestString",
		"42" );
	if ( herr < 0 ) return -1;

	h5_int32_t i4_val[1] = { 42 };
	herr = H5BlockWriteFieldAttribInt32 (
		f,
		"TestField",
		"TestInt32",
		i4_val, 1 );
	if ( herr < 0 ) return -1;

	h5_int64_t i8_val[1] = { 42 };
	herr = H5BlockWriteFieldAttribInt64 (
		f,
		"TestField",
		"TestInt64",
		i8_val, 1 );
	if ( herr < 0 ) return -1;

	h5_float32_t r4_val[1] = { 42.0 };
	herr = H5BlockWriteFieldAttribFloat32 (
		f,
		"TestField",
		"TestFloat32",
		r4_val, 1 );
	if ( herr < 0 ) return -1;

	h5_float64_t r8_val[1] = { 42.0 };
	herr = H5BlockWriteFieldAttribFloat64 (
		f,
		"TestField",
		"TestFloat64",
		r8_val, 1 );
	if ( herr < 0 ) return -1;

	herr = H5Block3dSetFieldOrigin ( f, "TestField", 1.0, 2.0, 3.0 );
	if ( herr < 0 ) return -1;

	herr = H5Block3dSetFieldSpacing ( f, "TestField", 2.0, 3.0, 4.0 );
	if ( herr < 0 ) return -1;

	return H5_SUCCESS;
}

static h5_int64_t
_write_file (
	const char *fname,
	const int myproc,
	MPI_Comm comm,
	struct H5BlockPartition *layout
	) {
		printf ("PROC[%d]: Open file \"%s\" for writing ...\n",
		myproc, fname );
        h5_file_t file = H5OpenFile (fname, H5_O_WRONLY, H5_PROP_DEFAULT);
        H5SetStep (file, 0);
	
	_write_data (file, myproc, layout);
	_write_attributes (file, myproc);
	
	H5CloseFile (file);
	
	return H5_SUCCESS;
}

static h5_int64_t
_read_data (
	h5_file_t f,
	int myproc,
	struct H5BlockPartition *layout
	) {

	h5_int64_t i, j, k, idx;
	h5_float64_t *data;
	h5_int64_t i_dims = layout->i_end - layout->i_start + 1;
	h5_int64_t j_dims = layout->j_end - layout->j_start + 1;
	h5_int64_t k_dims = layout->k_end - layout->k_start + 1;

	printf ("Reading Step #%lld\n", (long long)H5GetStep (f));

	data = (h5_float64_t*)malloc (i_dims * j_dims * k_dims * sizeof (*data));

	H5Block3dSetView (
		f,
		layout->i_start, layout->i_end,
		layout->j_start, layout->j_end,
		layout->k_start, layout->k_end );
	H5Block3dReadScalarFieldFloat64 ( f, "TestField", data );

	for (i = 0; i < i_dims; i++) {
		for (j = 0; j < j_dims; j++) {
			for (k = 0; k < k_dims; k++) {
				idx = _calc_index (
					i, i_dims,
					j, j_dims,
					k, k_dims );

				h5_float64_t value = k
					+ 1000 * j
					+ 100000 * i
					+ 10000000 * myproc;
				if (*(data + idx) != value) {
					printf (
						"PROC[%d]: "
						"value missmatch for (%lld,%lld,%lld); is: %f;"
						" should be: %f\n",
						myproc,
						(long long)i, (long long)j, (long long)k,
						*( data + idx ), value );
					return -1;
				}
			}
		}
	}

	free (data);

	return H5_SUCCESS;
}

static h5_int64_t
_read_attributes (
	h5_file_t f,
	const int myproc,
	MPI_Comm comm
	) {
	h5_int64_t timestep = 0;
	
	H5SetStep (f, timestep);

	char sval[16];
	H5BlockReadFieldAttribString (
		f,
		"TestField",
		"TestString",
		sval );
	if (strcmp (sval, "42") != 0) {
		printf ("Error reading string attribute: "
			"Value is \"%s\" and should be \"42\"\n", sval);
	}

	h5_int32_t i4_val[1];
	H5BlockReadFieldAttribInt32 (
		f,
		"TestField",
		"TestInt32",
		i4_val );
	if (i4_val[0] != 42) {
		printf ("Error reading int32 attribute: "
			"Value is %lld and should be 42\n",
			(long long) i4_val[0]);
	}

	h5_int64_t i8_val[1];
	H5BlockReadFieldAttribInt64 (
		f,
		"TestField",
		"TestInt64",
		i8_val );
	if (i8_val[0] != 42) {
		printf ("Error reading int64 attribute: "
			"Value is %lld and should be 42\n",
			(long long) i8_val[0]);
	}

	h5_float32_t r4_val[1];
	H5BlockReadFieldAttribFloat32 (
		f,
		"TestField",
		"TestFloat32",
		r4_val );
	if (r4_val[0] != 42.0) {
		printf ("Error reading float64 attribute: "
			"Value is %f and should be 42.0\n",
			r4_val[0]);
	}

	h5_float64_t r8_val[1];
	H5BlockReadFieldAttribFloat64 (
		f,
		"TestField",
		"TestFloat64",
		r8_val );
	if (r8_val[0] != 42.0) {
		printf ("Error reading float64 attribute: "
			"Value is %f and should be 42.0\n",
			r8_val[0]);
	}

	h5_float64_t x_origin;
	h5_float64_t y_origin;
	h5_float64_t z_origin;
	h5_float64_t x_spacing;
	h5_float64_t y_spacing;
	h5_float64_t z_spacing;

	H5Block3dGetFieldOrigin (
		f, "TestField",
		&x_origin,
		&y_origin,
		&z_origin );

	if (x_origin != 1.0 || y_origin != 2.0 || z_origin != 3.0) {
		printf (
			"Error reading field origin: Read values (%f,%f,%f)\n",
			x_origin, y_origin, z_origin);
	}
	H5Block3dGetFieldSpacing (
		f, "TestField",
		&x_spacing,
		&y_spacing,
		&z_spacing );
	if (x_spacing != 2.0 || y_spacing != 3.0 || z_spacing != 4.0) {
		printf (
			"Error reading field spacing: Read values (%f,%f,%f)\n",
			x_spacing, y_spacing, z_spacing);
	}

	return H5_SUCCESS;
}

static h5_int64_t
_read_file (
	const char *fname,
	const int myproc,
	MPI_Comm comm,
	struct H5BlockPartition *layout
	) {
	printf ("PROC[%d]: Open file \"%s\" for reading ...\n",
		myproc, fname );

	h5_file_t file = H5OpenFile (fname, H5_O_RDONLY, H5_PROP_DEFAULT);
        H5SetStep (file, 0);

	_read_data (file, myproc, layout);
        _read_attributes (file, myproc, comm);

	H5CloseFile (file);
	
	return H5_SUCCESS;
}

int
main (
	int argc,
	char **argv
	) {
        h5_int64_t verbosity = DEFAULT_VERBOSITY;
	char *fname;
	int opt_with_ghosts = 0;
	int opt_read = 0;
	int opt_write = 0;
	struct H5BlockPartition *layout;

        if (argc == 1) {
                fprintf ( stderr,
                          "Usage: %s -w|-r [-g]\n",
                          argv[0] );
                return 1;
        }
	while (--argc) {
		if (strcmp (argv[argc], "-r") == 0)
			opt_read = 1;
		else if (strcmp (argv[argc], "-w") == 0)
			opt_write = 1;
		else if (strcmp (argv[argc], "-g") == 0)
			opt_with_ghosts = 1;
		else {
			fprintf (stderr,
				 "Illegal option %s\n\n"
				 "Usage: %s -w -r -g\n",
				 argv[argc], argv[0]);
			return 1;
		}
	}

        // initialize MPI & H5hut
        int comm_rank = 0;
        int comm_size = 1;
        MPI_Init (&argc, &argv);
        MPI_Comm comm = MPI_COMM_WORLD;
        MPI_Comm_rank (comm, &comm_rank);
        MPI_Comm_size (comm, &comm_size);

        H5AbortOnError ();
        H5SetVerbosityLevel (verbosity);

	switch (comm_size) {
	case 1:
		fname  = "blockfile1.h5";
		layout = &Layout1[comm_rank];
		break;
	case 8:
		if (opt_with_ghosts) {
			fname  = "blockfile8G.h5";
			layout = &Layout8G[comm_rank];
		} else {
			fname  = "blockfile8.h5";
			layout = &Layout8[comm_rank];
		}
		break;
	case 16:
		if (opt_with_ghosts) {
			fname  = "blockfile16G.h5";
			layout = &Layout16G[comm_rank];
		} else {
			fname  = "blockfile16.h5";
			layout = &Layout16[comm_rank];
		}
		break;
	case 32:
		if (opt_with_ghosts) {
			fname  = "blockfile32G.h5";
			layout = &Layout32G[comm_rank];
		} else {
			fname  = "blockfile32.h5";
			layout = &Layout32[comm_rank];
		}
		break;
	default:
		printf ( "Run this test on %d, %d, %d or %d processor(s)!\n",
			 1, 8, 16, 32);
		return 1;
	}

	if (opt_write) {
		_write_file (fname, comm_rank, comm, layout);
	} else if (opt_read) {
		_read_file (fname, comm_rank, comm, layout);
	}

	MPI_Finalize();
	return 0;
}
