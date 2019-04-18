/*
  Copyright (c) 2006-2015, The Regents of the University of California,
  through Lawrence Berkeley National Laboratory (subject to receipt of any
  required approvals from the U.S. Dept. of Energy) and the Paul Scherrer
  Institut (Switzerland).  All rights reserved.

  License: see file COPYING in top level of source distribution.
*/

#include "FC.h"
#include "h5_private.h"
#include "h5core/h5_log.h"
#include "h5core/h5b_io.h"

#define h5bl_3d_write_scalar_field_r8 FC_GLOBAL (	\
                h5bl_3d_write_scalar_field_r8,		\
                H5BL_3D_WRITE_SCALAR_FIELD_R8 )
h5_err_t
h5bl_3d_write_scalar_field_r8 (
	const h5_int64_t* const fh,
	const char* const name,
	const h5_float64_t* const buffer,
	const int l_name
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_err_t,
		      "fh=%p, name='%*s', buffer=%p",
		      (h5_file_p)f, l_name, name, buffer);
	char *name2 =  h5_strdupfor2c ( name, l_name );
	h5_err_t herr = h5b_write_scalar_data (
		f, name2, (void*)buffer, H5_FLOAT64_T );
	free ( name2 );
	H5_API_RETURN(herr);
}

#define h5bl_3d_read_scalar_field_r8 FC_GLOBAL (	\
                h5bl_3d_read_scalar_field_r8,		\
                H5BL_3D_READ_SCALAR_FIELD_R8 )
h5_err_t
h5bl_3d_read_scalar_field_r8 (
	const h5_int64_t* const fh,
	const char* const name,
	h5_float64_t* const buffer,
	const int l_name
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_err_t,
		      "fh=%p, name='%s', buffer=%p, l_name=%d",
		      (h5_file_p)f, name, buffer, l_name);
	char *name2 =  h5_strdupfor2c ( name,  l_name );
	h5_err_t herr = h5b_read_scalar_data (
		f, name2, buffer, H5_FLOAT64_T );
	free ( name2 );
	H5_API_RETURN(herr);
}

#define h5bl_3d_write_vector3d_field_r8 FC_GLOBAL (	\
                h5bl_3d_write_vector3d_field_r8,        \
                H5BL_3D_WRITE_VECTOR3D_FIELD_R8 )
h5_err_t
h5bl_3d_write_vector3d_field_r8 (
	const h5_int64_t* const fh,
	const char* const name,
	const h5_float64_t* const x_buf,
	const h5_float64_t* const y_buf,
	const h5_float64_t* const z_buf,
	const int l_name
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_err_t,
		      "fh=%p, name='%s', x_buf=%p, y_buf=%p, z_buf=%p, l_name=%d",
		      (h5_file_p)f, name, x_buf, y_buf, z_buf, l_name);
	char *name2 =  h5_strdupfor2c ( name,  l_name );
	h5_err_t herr = h5b_write_vector3d_data (
		f, name2,
		(void*)x_buf, (void*)y_buf, (void*)z_buf, H5_FLOAT64_T );
	free ( name2 );
	H5_API_RETURN(herr);
}

#define h5bl_3d_read_vector3d_field_r8 FC_GLOBAL (       \
                h5bl_3d_read_vector3d_field_r8,		   \
                H5BL_3D_READ_VECTOR3D_FIELD_R8 )
h5_err_t
h5bl_3d_read_vector3d_field_r8 (
	const h5_int64_t* const fh,
	const char* const name,
	h5_float64_t* const x_buf,
	h5_float64_t* const y_buf,
	h5_float64_t* const z_buf,
	const int l_name
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_err_t,
		      "fh=%p, name='%s', x_buf=%p, y_buf=%p, z_buf=%p, l_name=%d",
		      (h5_file_p)f, name, x_buf, y_buf, z_buf, l_name);
	char *name2 =  h5_strdupfor2c ( name,  l_name );
	h5_err_t herr = h5b_read_vector3d_data (
		f, name2,
		(void*)x_buf, (void*)y_buf, (void*)z_buf, H5_FLOAT64_T );
	free ( name2 );
	H5_API_RETURN(herr);
}

#define h5bl_3d_write_scalar_field_r4 FC_GLOBAL (        \
                h5bl_3d_write_scalar_field_r4,		   \
                H5BL_3D_WRITE_SCALAR_FIELD_R4 )
h5_err_t
h5bl_3d_write_scalar_field_r4 (
	const h5_int64_t*const fh,
	const char* const name,
	const h5_float32_t* const buffer,
	const int l_name
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_err_t,
		      "fh=%p, name='%s', buffer=%p, l_name=%d",
		      (h5_file_p)f, name, buffer, l_name);
	char *name2 =  h5_strdupfor2c ( name, l_name );
	h5_err_t herr = h5b_write_scalar_data (
		f, name2, (void*)buffer, H5_FLOAT32_T );
	free ( name2 );
	H5_API_RETURN(herr);
}

#define h5bl_3d_read_scalar_field_r4 FC_GLOBAL (	\
                h5bl_3d_read_scalar_field_r4,		\
                H5BL_3D_READ_SCALAR_FIELD_R4 )
h5_err_t
h5bl_3d_read_scalar_field_r4 (
	const h5_int64_t* const fh,
	const char* const name,
	h5_float32_t* const buffer,
	const int l_name
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_err_t,
		      "fh=%p, name='%s', buffer=%p, l_name=%d",
		      (h5_file_p)f, name, buffer, l_name);
	char *name2 =  h5_strdupfor2c ( name,  l_name );
	h5_err_t herr = h5b_read_scalar_data (
		f, name2, buffer, H5_FLOAT32_T );
	free ( name2 );
	H5_API_RETURN(herr);
}

#define h5bl_3d_write_vector3d_field_r4 FC_GLOBAL (	\
                h5bl_3d_write_vector3d_field_r4,        \
                H5BL_3D_WRITE_VECTOR3D_FIELD_R4 )
h5_err_t
h5bl_3d_write_vector3d_field_r4 (
	const h5_int64_t* const fh,
	const char* const name,
	const h5_float32_t* const x_buf,
	const h5_float32_t* const y_buf,
	const h5_float32_t* const z_buf,
	const int l_name
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_err_t,
		      "fh=%p, name='%s', x_buf=%p, y_buf=%p, z_buf=%p, l_name=%d",
		      (h5_file_p)f, name, x_buf, y_buf, z_buf, l_name);
	char *name2 =  h5_strdupfor2c ( name,  l_name );
	h5_err_t herr = h5b_write_vector3d_data (
		f, name2,
		(void*)x_buf, (void*)y_buf, (void*)z_buf, H5_FLOAT32_T );
	free ( name2 );
	H5_API_RETURN(herr);
}

#define h5bl_3d_read_vector3d_field_r4 FC_GLOBAL (	\
                h5bl_3d_read_vector3d_field_r4,         \
                H5BL_3D_READ_VECTOR3D_FIELD_R4 )
h5_err_t
h5bl_3d_read_vector3d_field_r4 (
	const h5_int64_t* const fh,
	const char* const name,
	h5_float32_t* const x_buf,
	h5_float32_t* const y_buf,
	h5_float32_t* const z_buf,
	const int l_name
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_err_t,
		      "fh=%p, name='%s', x_buf=%p, y_buf=%p, z_buf=%p, l_name=%d",
		      (h5_file_p)f, name, x_buf, y_buf, z_buf, l_name);
	char *name2 =  h5_strdupfor2c ( name,  l_name );
	h5_err_t herr = h5b_read_vector3d_data (
		f, name2,
		(void*)x_buf, (void*)y_buf, (void*)z_buf, H5_FLOAT32_T );
	free ( name2 );
	H5_API_RETURN (herr);
}

#define h5bl_3d_write_scalar_field_i8 FC_GLOBAL (	\
                h5bl_3d_write_scalar_field_i8,          \
                H5BL_3D_WRITE_SCALAR_FIELD_I8 )
h5_err_t
h5bl_3d_write_scalar_field_i8 (
	const h5_int64_t* const fh,
	const char* const name,
	const h5_int64_t* const buffer,
	const int l_name
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_err_t,
		      "fh=%p, name='%s', buffer=%p, l_name=%d",
		      (h5_file_p)f, name, buffer, l_name);
	char *name2 =  h5_strdupfor2c ( name, l_name );
	h5_err_t herr = h5b_write_scalar_data (
		f, name2, (void*)buffer, H5_INT64_T );
	free ( name2 );
	H5_API_RETURN (herr);
}

#define h5bl_3d_read_scalar_field_i8 FC_GLOBAL (	\
                h5bl_3d_read_scalar_field_i8,		\
                H5BL_3D_READ_SCALAR_FIELD_I8 )
h5_err_t
h5bl_3d_read_scalar_field_i8 (
	const h5_int64_t* const fh,
	const char* const name,
	h5_int64_t* const buffer,
	const int l_name
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_err_t,
		      "fh=%p, name='%s', buffer=%p, l_name=%d",
		      (h5_file_p)f, name, buffer, l_name);
	char *name2 =  h5_strdupfor2c ( name,  l_name );
	h5_err_t herr = h5b_read_scalar_data (
		f, name2, buffer, H5_INT64_T );
	free (name2);
	H5_API_RETURN (herr);
}

#define h5bl_3d_write_vector3d_field_i8 FC_GLOBAL (	\
                h5bl_3d_write_vector3d_field_i8,        \
                H5BL_3D_WRITE_VECTOR3D_FIELD_I8 )
h5_err_t
h5bl_3d_write_vector3d_field_i8 (
	const h5_int64_t* const fh,
	const char* const name,
	const h5_int64_t* const x_buf,
	const h5_int64_t* const y_buf,
	const h5_int64_t* const z_buf,
	const int l_name
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_err_t,
		      "fh=%p, name='%s', x_buf=%p, y_buf=%p, z_buf=%p, l_name=%d",
		      (h5_file_p)f, name, x_buf, y_buf, z_buf, l_name);
	char *name2 =  h5_strdupfor2c ( name,  l_name );
	h5_err_t herr = h5b_write_vector3d_data (
		f, name2,
		(void*)x_buf, (void*)y_buf, (void*)z_buf, H5_INT64_T );
	free ( name2 );
	H5_API_RETURN(herr);
}

#define h5bl_3d_read_vector3d_field_i8 FC_GLOBAL (	\
                h5bl_3d_read_vector3d_field_i8,         \
                H5BL_3D_READ_VECTOR3D_FIELD_I8 )
h5_err_t
h5bl_3d_read_vector3d_field_i8 (
	const h5_int64_t *const fh,
	const char* const name,
	h5_int64_t* const x_buf,
	h5_int64_t* const y_buf,
	h5_int64_t* const z_buf,
	const int l_name
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_err_t,
		      "fh=%p, name='%s', x_buf=%p, y_buf=%p, z_buf=%p, l_name=%d",
		      (h5_file_p)f, name, x_buf, y_buf, z_buf, l_name);
	char *name2 =  h5_strdupfor2c ( name,  l_name );
	h5_err_t herr = h5b_read_vector3d_data (
		f, name2,
		(void*)x_buf, (void*)y_buf, (void*)z_buf, H5_INT64_T );
	free ( name2 );
	H5_API_RETURN(herr);
}

#define h5bl_3d_write_scalar_field_i4 FC_GLOBAL (	\
                h5bl_3d_write_scalar_field_i4,          \
                H5BL_3D_WRITE_SCALAR_FIELD_I4 )
h5_err_t
h5bl_3d_write_scalar_field_i4 (
	const h5_int64_t* const fh,
	const char* const name,
	const h5_int32_t* const buffer,
	const int l_name
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_err_t,
		      "fh=%p, name='%s', buffer=%p, l_name=%d",
		      (h5_file_p)f, name, buffer, l_name);
	char *name2 =  h5_strdupfor2c ( name, l_name );
	h5_err_t herr = h5b_write_scalar_data (
		f, name2, (void*)buffer, H5_INT32_T );
	free ( name2 );
	H5_API_RETURN(herr);
}

#define h5bl_3d_read_scalar_field_i4 FC_GLOBAL (	\
                h5bl_3d_read_scalar_field_i4,		\
                H5BL_3D_READ_SCALAR_FIELD_I4 )
h5_err_t
h5bl_3d_read_scalar_field_i4 (
	const h5_int64_t* const fh,
	const char* const name,
	h5_int32_t* const buffer,
	const int l_name
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_err_t,
		      "fh=%p, name='%s', buffer=%p, l_name=%d",
		      (h5_file_p)f, name, buffer, l_name);
	char *name2 =  h5_strdupfor2c ( name,  l_name );
	h5_err_t herr = h5b_read_scalar_data (
		f, name2, buffer, H5_INT32_T );
	free ( name2 );
	H5_API_RETURN(herr);
}

#define h5bl_3d_write_vector3d_field_i4 FC_GLOBAL (	\
                h5bl_3d_write_vector3d_field_i4,        \
                H5BL_3D_WRITE_VECTOR3D_FIELD_I4 )
h5_err_t
h5bl_3d_write_vector3d_field_i4 (
	const h5_int64_t* const fh,
	const char* const name,
	const h5_int32_t* const x_buf,
	const h5_int32_t* const y_buf,
	const h5_int32_t* const z_buf,
	const int l_name
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_err_t,
		      "fh=%p, name='%s', x_buf=%p, y_buf=%p, z_buf=%p, l_name=%d",
		      (h5_file_p)f, name, x_buf, y_buf, z_buf, l_name);
	char *name2 =  h5_strdupfor2c ( name,  l_name );
	h5_err_t herr = h5b_write_vector3d_data (
		f, name2,
		(void*)x_buf, (void*)y_buf, (void*)z_buf, H5_INT32_T );
	free ( name2 );
	H5_API_RETURN(herr);
}

#define h5bl_3d_read_vector3d_field_i4 FC_GLOBAL (	\
                h5bl_3d_read_vector3d_field_i4,         \
                H5BL_3D_READ_VECTOR3D_FIELD_I4 )
h5_err_t
h5bl_3d_read_vector3d_field_i4 (
	const h5_int64_t* const fh,
	const char* const name,
	h5_int32_t* const x_buf,
	h5_int32_t* const y_buf,
	h5_int32_t* const z_buf,
	const int l_name
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_err_t,
		      "fh=%p, name='%s', x_buf=%p, y_buf=%p, z_buf=%p, l_name=%d",
		      (h5_file_p)f, name, x_buf, y_buf, z_buf, l_name);
	char *name2 =  h5_strdupfor2c ( name,  l_name );
	h5_err_t herr = h5b_read_vector3d_data (
		f, name2,
		(void*)x_buf, (void*)y_buf, (void*)z_buf, H5_INT32_T );
	free ( name2 );
	H5_API_RETURN(herr);
}
