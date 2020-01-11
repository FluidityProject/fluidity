/*
  Copyright (c) 2006-2015, The Regents of the University of California,
  through Lawrence Berkeley National Laboratory (subject to receipt of any
  required approvals from the U.S. Dept. of Energy) and the Paul Scherrer
  Institut (Switzerland).  All rights reserved.

  License: see file COPYING in top level of source distribution.
*/

#include "h5_private.h"
#include "h5core/h5_log.h"

#include "h5core/h5u_io.h"

/*==================Writing data ============*/
#define h5pt_writedata_r8 FC_MANGLING (					\
                h5pt_writedata_r8,                                      \
                H5PT_WRITEDATA_R8 )
h5_int64_t
h5pt_writedata_r8 (
	const h5_int64_t* const fh,
	const char* const name,
	const h5_float64_t* const data,
	const int l_name
        ) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_int64_t,
                      "fh=%p, name='%s', data=%p, l_name=%d",
                      (h5_file_p)f, name, data, l_name);
	char *name2 = h5_strdupfor2c ( name, l_name );
	h5_int64_t herr = h5u_write_dataset (
		f, name2, (void*)data, H5_FLOAT64_T );
	free ( name2 );
	H5_API_RETURN(herr);
}

#define h5pt_writedata_r4 FC_MANGLING (					\
                h5pt_writedata_r4,                                      \
                H5PT_WRITEDATA_R4 )
h5_int64_t
h5pt_writedata_r4 (
	const h5_int64_t* const fh,
	const char* const name,
	const h5_float32_t* const data,
	const int l_name
        ) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_int64_t,
                      "fh=%p, name='%s', data=%p, l_name=%d",
                      (h5_file_p)f, name, data, l_name);
	char *name2 = h5_strdupfor2c ( name, l_name );
	h5_int64_t herr = h5u_write_dataset (
		f, name2, (void*)data, H5_FLOAT32_T );
	free ( name2 );
	H5_API_RETURN(herr);
}

#define h5pt_writedata_i8 FC_MANGLING (					\
                h5pt_writedata_i8,                                      \
                H5PT_WRITEDATA_I8 )
h5_int64_t
h5pt_writedata_i8 (
	const h5_int64_t* const fh,
	const char* const name,
	const h5_int64_t* const data,
	const int l_name
        ) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_int64_t,
                      "fh=%p, name='%s', data=%p, l_name=%d",
                      (h5_file_p)f, name, data, l_name);
	char *name2 = h5_strdupfor2c ( name, l_name );
	h5_int64_t herr = h5u_write_dataset (
		f, name2, (void*)data, H5_INT64_T );
	free ( name2 );
	H5_API_RETURN(herr);
}

#define h5pt_writedata_i4 FC_MANGLING (					\
                h5pt_writedata_i4,                                      \
                H5PT_WRITEDATA_I4 )
h5_int64_t
h5pt_writedata_i4 (
	const h5_int64_t* const fh,
	const char* const name,
	const h5_int32_t* const data,
	const int l_name
        ) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_int64_t,
                      "fh=%p, name='%s', data=%p, l_name=%d",
                      (h5_file_p)f, name, data, l_name);
	char *name2 = h5_strdupfor2c ( name, l_name );
	h5_int64_t herr = h5u_write_dataset (
		f, name2, (void*)data, H5_INT32_T );
	free ( name2 );
	H5_API_RETURN(herr);
}


/*==================Reading data ============*/
#define h5pt_readdata_r8 FC_MANGLING (					\
                h5pt_readdata_r8,                                       \
                H5PT_READDATA_R8 )
h5_int64_t
h5pt_readdata_r8 (
	const h5_int64_t* const fh,
	const char* const name,
	h5_float64_t* const data,
	const int l_name
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_int64_t,
                      "fh=%p, name='%s', data=%p, l_name=%d",
                      (h5_file_p)f, name, data, l_name);
	char *name2 = h5_strdupfor2c ( name, l_name );
	h5_int64_t herr = h5u_read_dataset (
		f, name2, data, H5_FLOAT64_T );
	free ( name2 );
	H5_API_RETURN(herr);
}

#define h5pt_readdata_r4 FC_MANGLING (					\
                h5pt_readdata_r4,                                       \
                H5PT_READDATA_R4 )
h5_int64_t
h5pt_readdata_r4 (
	const h5_int64_t* const fh,
	const char* const name,
	h5_float32_t* const data,
	const int l_name
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_int64_t, 
                      "fh=%p, name='%s', data=%p, l_name=%d",
                      (h5_file_p)f, name, data, l_name);
	char *name2 = h5_strdupfor2c ( name, l_name );
	h5_int64_t herr = h5u_read_dataset (
		f, name2, data, H5_FLOAT32_T );
	free ( name2 );
	H5_API_RETURN(herr);
}

#define h5pt_readdata_i8 FC_MANGLING (					\
                h5pt_readdata_i8,                                       \
                H5PT_READDATA_I8 )
h5_int64_t
h5pt_readdata_i8 (
	const h5_int64_t* const fh,
	const char* const name,
	h5_int64_t* const data,
	const int l_name
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_int64_t,
                      "fh=%p, name='%s', data=%p, l_name=%d",
                      (h5_file_p)f, name, data, l_name);
	char *name2 = h5_strdupfor2c ( name, l_name );
	h5_int64_t herr = h5u_read_dataset (
		f, name2, data, H5_INT64_T );

	free ( name2 );
	H5_API_RETURN(herr);
}

#define h5pt_readdata_i4 FC_MANGLING (					\
                h5pt_readdata_i4,                                       \
                H5PT_READDATA_I4 )
h5_int64_t
h5pt_readdata_i4 (
	const h5_int64_t* const fh,
	const char* const name,
	h5_int32_t* const data,
	const int l_name
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_int64_t,
                      "hf=%p, name='%s', data=%p, l_name=%d",
                      (h5_file_p)f, name, data, l_name);
	char *name2 = h5_strdupfor2c ( name, l_name );
	h5_int64_t herr = h5u_read_dataset (
		f, name2, data, H5_INT32_T );
	free ( name2 );
	H5_API_RETURN(herr);
}
