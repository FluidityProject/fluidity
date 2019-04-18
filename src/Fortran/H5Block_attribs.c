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
#include "h5core/h5b_attribs.h"

/*
   __ _ _   _  ___ _ __ _   _ 
  / _` | | | |/ _ \ '__| | | |
 | (_| | |_| |  __/ |  | |_| |
  \__, |\__,_|\___|_|   \__, |
     |_|                |___/
*/

#define h5bl_getnfieldattribs FC_GLOBAL (				\
                h5bl_getnfieldattribs,                                  \
                H5BL_GETNFIELDATTRIBS)
h5_int64_t
h5bl_getnfieldattribs (
	const h5_int64_t* const fh,
	const char* const name,
	const int l_name
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_int64_t,
		      "fh=%p, name='%*s'",
                      (h5_file_p)f, l_name, name);
	char* name2 = h5_strdupfor2c ( name, l_name );
	h5_int64_t herr = h5b_get_num_field_attribs (f, name2);
	free (name2);
	H5_API_RETURN (herr);
}

#define h5bl_getfieldattribinfo FC_GLOBAL (				\
                h5bl_getfieldattribinfo,                                \
                h5bl_getfieldattribinfo)
h5_int64_t
h5bl_getfieldattribinfo (
	const h5_int64_t* const fh,
	const char* const field_name,
	const h5_int64_t* const attrib_idx,
	char* const attrib_name,
	h5_int64_t* const attrib_nelem,
	const int l_field_name,
	const int l_attrib_name
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_err_t,
		      "fh=%p, field_name='%*s', attrib_idx=%lld, "
                      "attrib_name=%p, attrib_nelem=%p",
		      (h5_file_p)f,
		      l_field_name, field_name,
		      (long long)*attrib_idx, 
		      attrib_name, attrib_nelem);

	char *field_name2 = h5_strdupfor2c ( field_name, l_field_name );
	h5_int64_t attrib_type;
	h5_int64_t herr = h5b_get_field_attrib_info_by_idx (
                f,
                field_name2, *attrib_idx - 1,
		attrib_name, l_attrib_name,
		&attrib_type,
		(h5_size_t*)attrib_nelem );

	h5_strc2for ( attrib_name, l_attrib_name );

	free (field_name2);
	H5_API_RETURN (herr);
}

/*
  _    __    
 (_)  / /__  
 | | / / _ \ 
 | |/ / (_) |
 |_/_/ \___/ 
*/
static inline h5_int64_t
write_field_attrib (
        const h5_file_t fh,
	const char* field_name,
	const int l_field_name,
	const char* attrib_name,
	const int l_attrib_name,
        const h5_int64_t attrib_type,
	const void* attrib_value,
	const hsize_t attrib_nelems
        ) {
	char *field_name2 = h5_strdupfor2c (field_name, l_field_name);
	char *attrib_name2 = h5_strdupfor2c (attrib_name, l_attrib_name);
	h5_int64_t h5err = h5b_write_field_attrib (
                fh, field_name2,
                attrib_name2, attrib_type,
                attrib_value, attrib_nelems);
	free (field_name2);
	free (attrib_name2);
        return h5err;
}

static inline h5_int64_t
read_field_attrib (
        const h5_file_t fh,
	const char* field_name,
	const int l_field_name,
	const char* attrib_name,
	const int l_attrib_name,
        const hid_t attrib_type,
	void* attrib_value
	) {
	char *field_name2 = h5_strdupfor2c (field_name, l_field_name);
	char *attrib_name2 = h5_strdupfor2c (attrib_name, l_attrib_name);
	h5_int64_t h5err = h5b_read_field_attrib (
                fh, field_name2,
                attrib_name2, attrib_type, attrib_value);
	free (field_name2);
	free (attrib_name2);
	return h5err;
}

/*
      _        _             
  ___| |_ _ __(_)_ __   __ _ 
 / __| __| '__| | '_ \ / _` |
 \__ \ |_| |  | | | | | (_| |
 |___/\__|_|  |_|_| |_|\__, |
                       |___/ 
*/
#define h5bl_writefieldattrib_string FC_GLOBAL (			\
                h5bl_writefieldattrib_string,                           \
                H5BL_WRITEFIELDATTRIB_STRING)
h5_int64_t
h5bl_writefieldattrib_string (
	const h5_int64_t* const fh,
	const char* const field_name,
	const char* const attrib_name,
	const char* const attrib_value,
	const int l_field_name,
	const int l_attrib_name,
	const int l_attrib_value
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_int64_t,
		      "fh=%p, field_name='%.*s', "
                      "attrib_name='%.*s' attrib_value='%.*s'",
		      (h5_file_p)f,
                      l_field_name, field_name,
                      l_attrib_name, attrib_name,
                      l_attrib_value, attrib_value);
	char* attrib_value2 = h5_strdupfor2c (attrib_value, l_attrib_value);
	h5_int64_t h5err = write_field_attrib (
                f,
                field_name, l_field_name,
                attrib_name, l_attrib_name,
                H5_STRING_T,
                attrib_value2, strlen(attrib_value2)+1 );
	free (attrib_value2);
	H5_API_RETURN (h5err);
}

#define h5bl_readfieldattrib_string FC_GLOBAL (			\
                h5bl_readfieldattrib_string,                            \
                H5BL_READFIELDATTRIB_STRING)
h5_err_t
h5bl_readfieldattrib_string (
	const h5_int64_t* const fh,
	const char* const field_name,
	const char* const attrib_name,
	char* const attrib_value,
	const int l_field_name,
	const int l_attrib_name,
	const int l_attrib_value
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_int64_t,
		      "f=%p, field_name='%.*s', attrib_name='%.*s' attrib_value='%p'",
		      (h5_file_p)f,
                      l_field_name, field_name,
                      l_attrib_name, attrib_name,
                      attrib_value);
        h5_int64_t h5err = read_field_attrib (
                f,
		field_name, l_field_name,
                attrib_name, l_attrib_name,
		H5_STRING_T, attrib_value);

	h5_strc2for (attrib_value, l_attrib_value);
	H5_API_RETURN (h5err);
}

/*
                 _ 
  _ __ ___  __ _| |
 | '__/ _ \/ _` | |
 | | |  __/ (_| | |
 |_|  \___|\__,_|_|
*/
#define h5bl_writefieldattrib_r8 FC_GLOBAL (			    \
                h5bl_writefieldattrib_r8,                           \
                H5BL_WRITEFIELDATTRIB_R8)
h5_int64_t
h5bl_writefieldattrib_r8 (
	const h5_int64_t* const fh,
	const char* const field_name,
	const char* const attrib_name,
	const h5_float64_t* const attrib_value,
	const h5_int64_t* const attrib_nelems,
	const int l_field_name,
	const int l_attrib_name
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_int64_t,
		      "fh=%p, field_name='%.*s', attrib_name='%.*s', "
                      "attrib_value=%p, attrib_nelems=%lld",
		      (h5_file_p)f,
                      l_field_name, field_name,
                      l_attrib_name, attrib_name,
                      attrib_value, (long long)*attrib_nelems);
	H5_API_RETURN (write_field_attrib (
                               f,
                               field_name, l_field_name,
                               attrib_name, l_attrib_name,
                               H5_FLOAT64_T,
                               attrib_value, *attrib_nelems));
}

#define h5bl_readfieldattrib_r8 FC_GLOBAL (			    \
                h5bl_readfieldattrib_r8,                            \
                H5BL_READFIELDATTRIB_R8)
h5_err_t
h5bl_readfieldattrib_r8 (
	const h5_int64_t* const fh,
	const char* const field_name,
	const char* const  attrib_name,
	h5_float64_t* const attrib_value,
	const int l_field_name,
	const int l_attrib_name
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_int64_t,
		      "fh=%p, field_name='%.*s', attrib_name='%.*s', "
                      "attrib_value=%p",
		      (h5_file_p)f,
                      l_field_name, field_name,
                      l_attrib_name, attrib_name,
                      attrib_value);
        H5_API_RETURN (read_field_attrib (
                               f,
                               field_name, l_field_name,
                               attrib_name, l_attrib_name,
                               H5_FLOAT64_T,
                               attrib_value));
}

#define h5bl_writefieldattrib_r4 FC_GLOBAL (			    \
                h5bl_writefieldattrib_r4,                           \
                H5BL_WRITEFIELDATTRIB_R4)
h5_int64_t
h5bl_writefieldattrib_r4 (
	const h5_int64_t* const fh,
	const char* const field_name,
	const char* const attrib_name,
	const h5_float32_t* const attrib_value,
	const h5_int64_t* const attrib_nelems,
	const int l_field_name,
	const int l_attrib_name
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_int64_t,
		      "fh=%p, field_name='%.*s', attrib_name='%.*s', "
                      "attrib_value=%p, attrib_nelems=%lld",
		      (h5_file_p)f,
                      l_field_name, field_name,
                      l_attrib_name, attrib_name,
                      attrib_value, (long long)*attrib_nelems);
	H5_API_RETURN (write_field_attrib (
                               f,
                               field_name, l_field_name,
                               attrib_name, l_attrib_name,
                               H5_FLOAT32_T,
                               attrib_value, *attrib_nelems));
}

#define h5bl_readfieldattrib_r4 FC_GLOBAL (			    \
                h5bl_readfieldattrib_r4,                            \
                H5BL_READFIELDATTRIB_R4)
h5_err_t
h5bl_readfieldattrib_r4 (
	const h5_int64_t* const fh,
	const char* const field_name,
	const char* const attrib_name,
	h5_float32_t* const attrib_value,
	const int l_field_name,
	const int l_attrib_name
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_int64_t,
		      "fh=%p, field_name='%.*s', attrib_name='%.*s', "
                      "attrib_value=%p",
		      (h5_file_p)f,
                      l_field_name, field_name,
                      l_attrib_name, attrib_name,
                      attrib_value);
        H5_API_RETURN (read_field_attrib (
                               f,
                               field_name, l_field_name,
                               attrib_name, l_attrib_name,
                               H5_FLOAT32_T,
                               attrib_value));
}

/*
  _       _                       
 (_)_ __ | |_ ___  __ _  ___ _ __ 
 | | '_ \| __/ _ \/ _` |/ _ \ '__|
 | | | | | ||  __/ (_| |  __/ |   
 |_|_| |_|\__\___|\__, |\___|_|   
                  |___/
*/
#define h5bl_writefieldattrib_i8 FC_GLOBAL (			    \
                h5bl_writefieldattrib_i8,                           \
                H5BL_WRITEFIELDATTRIB_I8)
h5_int64_t
h5bl_writefieldattrib_i8 (
	const h5_int64_t* const fh,
	const char* const field_name,
	const char* const attrib_name,
	const h5_int64_t* const attrib_value,
	const h5_int64_t* const attrib_nelems,
	const int l_field_name,
	const int l_attrib_name
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_int64_t,
		      "fh=%p, field_name='%.*s', attrib_name='%.*s', "
                      "attrib_value=%p, attrib_nelems=%lld",
		      (h5_file_p)f,
                      l_field_name, field_name,
                      l_attrib_name, attrib_name,
                      attrib_value, (long long)*attrib_nelems);
	H5_API_RETURN (write_field_attrib (
                               f,
                               field_name, l_field_name,
                               attrib_name, l_attrib_name,
                               H5_INT64_T,
                               attrib_value, *attrib_nelems));
}

#define h5bl_readfieldattrib_i8 FC_GLOBAL (			    \
                h5bl_readfieldattrib_i8,                            \
                H5BL_READFIELDATTRIB_I8)
h5_err_t
h5bl_readfieldattrib_i8 (
	const h5_int64_t* const fh,
	const char* const field_name,
	const char* const attrib_name,
	h5_int64_t* const attrib_value,
	const int l_field_name,
	const int l_attrib_name,
	const int l_attrib_value
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_int64_t,
		      "fh=%p, field_name='%.*s', attrib_name='%.*s', "
                      "attrib_value=%p",
		      (h5_file_p)f,
                      l_field_name, field_name,
                      l_attrib_name, attrib_name,
                      attrib_value);
        H5_API_RETURN (read_field_attrib (
                               f,
                               field_name, l_field_name,
                               attrib_name, l_attrib_name,
                               H5_INT64_T,
                               attrib_value));
}

#define h5bl_writefieldattrib_i4 FC_GLOBAL (			    \
                h5bl_writefieldattrib_i4,                           \
                H5BL_WRITEFIELDATTRIB_I4 )
h5_int64_t
h5bl_writefieldattrib_i4 (
	const h5_int64_t* const fh,
	const char* const field_name,
	const char* const attrib_name,
	const h5_int32_t* const attrib_value,
	const h5_int64_t* const attrib_nelems,
	const int l_field_name,
	const int l_attrib_name
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_int64_t,
		      "fh=%p, field_name='%.*s', attrib_name='%.*s', "
                      "attrib_value=%p, attrib_nelems=%lld",
		      (h5_file_p)f,
                      l_field_name, field_name,
                      l_attrib_name, attrib_name,
                      attrib_value, (long long)*attrib_nelems);
	H5_API_RETURN (write_field_attrib (
                               f,
                               field_name, l_field_name,
                               attrib_name, l_attrib_name,
                               H5_INT32_T,
                               attrib_value, *attrib_nelems));
}

#define h5bl_readfieldattrib_i4 FC_GLOBAL (			    \
                h5bl_readfieldattrib_i4,                            \
                H5BL_READFIELDATTRIB_I4)
h5_err_t
h5bl_readfieldattrib_i4 (
	const h5_int64_t* const fh,
	const char* const field_name,
	const char* const attrib_name,
	h5_int32_t* const attrib_value,
	const int l_field_name,
	const int l_attrib_name
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_int64_t,
		      "fh=%p, field_name='%.*s', attrib_name='%.*s', "
                      "attrib_value=%p",
		      (h5_file_p)f,
                      l_field_name, field_name,
                      l_attrib_name, attrib_name,
                      attrib_value);
        H5_API_RETURN (read_field_attrib (
                               f,
                               field_name, l_field_name,
                               attrib_name, l_attrib_name,
                               H5_INT32_T,
                               attrib_value));
}

#define h5bl_get_fieldorigin FC_GLOBAL (	\
		h5bl_get_fieldorigin,		\
		H5BL_GET_FIELDORIGIN)

h5_int64_t
h5bl_get_fieldorigin (
	const h5_int64_t* const fh,
	const char* const field_name,
	h5_float64_t* const x_origin,
	h5_float64_t* const y_origin,
	h5_float64_t* const z_origin,
	const int l_field_name
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_int64_t,
		      "fh=%p, field_name='%.*s'"
                      ", x_origin=%p"
                      ", y_origin=%p"
                      ", z_origin=%p",
		      (h5_file_p)f,
                      l_field_name, field_name,
		      x_origin, y_origin, z_origin);
	h5_float64_t origin[3];
	TRY (read_field_attrib (
		     f,
		     field_name, l_field_name,
		     H5BLOCK_FIELD_ORIGIN_NAME, sizeof (H5BLOCK_FIELD_ORIGIN_NAME),
		     H5_FLOAT64_T,
		     origin));
	*x_origin = origin[0];
	*y_origin = origin[1];
	*z_origin = origin[2];

	H5_API_RETURN (H5_SUCCESS);
}

#define h5bl_set_fieldorigin FC_GLOBAL (	\
		h5bl_set_fieldorigin,		\
		H5BL_SET_FIELDORIGIN)

h5_int64_t
h5bl_set_fieldorigin (
	const h5_int64_t* const fh,
	const char* const field_name,
	h5_float64_t* const x_origin,
	h5_float64_t* const y_origin,
	h5_float64_t* const z_origin,
	const int l_field_name
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_int64_t,
		      "fh=%p, field_name='%.*s'"
                      ", x_origin=%g"
                      ", y_origin=%g"
                      ", z_origin=%g",
		      (h5_file_p)f,
                      l_field_name, field_name,
		      *x_origin, *y_origin, *z_origin);
	h5_float64_t origin[3] = { *x_origin, *y_origin, *z_origin };
	TRY (write_field_attrib (
		     f,
		     field_name, l_field_name,
		     H5BLOCK_FIELD_ORIGIN_NAME, sizeof (H5BLOCK_FIELD_ORIGIN_NAME),
		     H5_FLOAT64_T,
		     origin, 3));

	H5_API_RETURN (H5_SUCCESS);
}

#define h5bl_get_fieldspacing FC_GLOBAL (		\
		h5bl_get_fieldspacing,			\
		H5BL_GET_FIELDSPACING)

h5_int64_t
h5bl_get_fieldspacing (
	const h5_int64_t* const fh,
	const char* const field_name,
	h5_float64_t* const x_spacing,
	h5_float64_t* const y_spacing,
	h5_float64_t* const z_spacing,
	const int l_field_name
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_int64_t,
		      "fh=%p, field_name='%.*s'"
                      ", x_spacing=%p"
                      ", y_spacing=%p"
                      ", z_spacing=%p",
		      (h5_file_p)f,
                      l_field_name, field_name,
		      x_spacing, y_spacing, z_spacing);
	h5_float64_t spacing[3];
	TRY (read_field_attrib (
		     f,
		     field_name, l_field_name,
		     H5BLOCK_FIELD_ORIGIN_NAME, sizeof (H5BLOCK_FIELD_ORIGIN_NAME),
		     H5_FLOAT64_T,
		     spacing));
	*x_spacing = spacing[0];
	*y_spacing = spacing[1];
	*z_spacing = spacing[2];

	H5_API_RETURN (H5_SUCCESS);
}

#define h5bl_set_fieldspacing FC_GLOBAL (		\
		h5bl_set_fieldspacing,			\
		H5BL_SET_FIELDSPACING)

h5_int64_t
h5bl_set_fieldspacing (
	const h5_int64_t* const fh,
	const char* const field_name,
	h5_float64_t* const x_spacing,
	h5_float64_t* const y_spacing,
	h5_float64_t* const z_spacing,
	const int l_field_name
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_int64_t,
		      "fh=%p, field_name='%.*s'"
                      ", x_spacing=%g"
                      ", y_spacing=%g"
                      ", z_spacing=%g",
		      (h5_file_p)f,
                      l_field_name, field_name,
		      *x_spacing, *y_spacing, *z_spacing);
	h5_float64_t spacing[3] = { *x_spacing, *y_spacing, *z_spacing };
	TRY (read_field_attrib (
		     f,
		     field_name, l_field_name,
		     H5BLOCK_FIELD_ORIGIN_NAME, sizeof (H5BLOCK_FIELD_SPACING_NAME),
		     H5_FLOAT64_T,
		     spacing));

	H5_API_RETURN (H5_SUCCESS);
}


static inline h5_int64_t
set_field_coords (
	const h5_file_t f,
	int rank,
	const char* field_name,
	const int l_field_name,
        const char* attrib_name,
	const int l_attrib_name,
	const h5_float64_t* coords,
	const h5_int64_t n_coords
	) {
	char *field_name2 = h5_strdupfor2c (field_name, l_field_name);
	char *attrib_name2 = h5_strdupfor2c (attrib_name, l_attrib_name);

	h5_int64_t h5err = h5b_set_3d_field_coords (
		f, rank,
		field_name2, attrib_name2,
		coords, n_coords);

	free (field_name2);
	free (attrib_name2);

	return (h5err);
}

static inline h5_int64_t
get_field_coords (
	const h5_file_t f,
	int rank,
	const char* field_name,
	const int l_field_name,
        const char* attrib_name,
	const int l_attrib_name,
	h5_float64_t* const coords,
	const h5_int64_t n_coords
	) {
	char *field_name2 = h5_strdupfor2c (field_name, l_field_name);
	char *attrib_name2 = h5_strdupfor2c (attrib_name, l_attrib_name);

	h5_int64_t h5err = h5b_get_3d_field_coords (
		f, rank,
		field_name2, attrib_name2,
		coords, n_coords);

	free (field_name2);
	free (attrib_name2);

	return (h5err);
}

#define h5bl_set_fieldxcoords FC_GLOBAL (	 \
		h5bl_set_fieldxcoords,		 \
		H5BL_SET_FIELDXCOORDS)
h5_int64_t
h5bl_set_fieldxcoords (
	const h5_int64_t* const fh,
	const char* const field_name,
	h5_float64_t* const coords,
	const h5_int64_t* n_coords,
	const int l_field_name
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_err_t,
		      "fh=%p, field_name='%.*s'"
                      "coords=%p, n_coords=%llu",
                      (h5_file_p)f,
                      l_field_name, field_name,
                      coords, (long long unsigned)n_coords);
        H5_API_RETURN (set_field_coords (
                               f, 0,
			       field_name, l_field_name,
			       H5BLOCK_FIELD_XCOORD_NAME, sizeof(H5BLOCK_FIELD_XCOORD_NAME),
                               coords, *n_coords));
}

#define h5bl_get_fieldxcoords FC_GLOBAL (	 \
		h5bl_get_fieldxcoords,		 \
		H5BL_GET_FIELDXCOORDS)
h5_int64_t
h5bl_get_fieldxcoords (
	const h5_int64_t* const fh,
	const char* const field_name,
	h5_float64_t* const coords,
	const h5_int64_t* n_coords,
	const int l_field_name
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_err_t,
		      "fh=%p"
		      ", field_name='%.*s'"
                      ", coords=%p, n_coords=%llu",
                      (h5_file_p)f,
                      l_field_name, field_name,
                      coords, (long long unsigned)n_coords);
        H5_API_RETURN (get_field_coords (
                               f, 0,
			       field_name, l_field_name,
			       H5BLOCK_FIELD_XCOORD_NAME, sizeof (H5BLOCK_FIELD_XCOORD_NAME),
                               coords, *n_coords));
}

#define h5bl_set_fieldycoords FC_GLOBAL (	 \
		h5bl_set_fieldycoords,		 \
		H5BL_SET_FIELDYCOORDS)
h5_int64_t
h5bl_set_fieldycoords (
	const h5_int64_t* const fh,
	const char* const field_name,
	h5_float64_t* const coords,
	const h5_int64_t* n_coords,
	const int l_field_name
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_err_t,
		      "fh=%p"
		      ", field_name='%.*s'"
                      ", coords=%p, n_coords=%llu",
                      (h5_file_p)f,
                      l_field_name, field_name,
                      coords, (long long unsigned)n_coords);
        H5_API_RETURN (set_field_coords (
                               f, 1,
			       field_name, l_field_name,
			       H5BLOCK_FIELD_YCOORD_NAME, sizeof (H5BLOCK_FIELD_YCOORD_NAME),
                               coords, *n_coords));
}

#define h5bl_get_fieldycoords FC_GLOBAL (	 \
		h5bl_get_fieldycoords,		 \
		H5BL_GET_FIELDyCOORDS)
h5_int64_t
h5bl_get_fieldycoords (
	const h5_int64_t* const fh,
	const char* const field_name,
	h5_float64_t* const coords,
	const h5_int64_t* n_coords,
	const int l_field_name
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_err_t,
		      "fh=%p"
		      ", field_name='%.*s'"
                      ", coords=%p, n_coords=%llu",
                      (h5_file_p)f,
                      l_field_name, field_name,
                      coords, (long long unsigned)n_coords);
        H5_API_RETURN (get_field_coords (
                               f, 1,
			       field_name, l_field_name,
			       H5BLOCK_FIELD_YCOORD_NAME, sizeof (H5BLOCK_FIELD_YCOORD_NAME),
                               coords, *n_coords));
}


#define h5bl_set_fieldzcoords FC_GLOBAL (	 \
		h5bl_set_fieldzcoords,		 \
		H5BL_SET_FIELDZCOORDS)
h5_int64_t
h5bl_set_fieldzcoords (
	const h5_int64_t* const fh,
	const char* const field_name,
	h5_float64_t* const coords,
	const h5_int64_t* n_coords,
	const int l_field_name
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_err_t,
		      "fh=%p"
		      ", field_name='%.*s'"
                      ", coords=%p, n_coords=%llu",
                      (h5_file_p)f,
                      l_field_name, field_name,
                      coords, (long long unsigned)n_coords);
        H5_API_RETURN (set_field_coords (
                               f, 2,
			       field_name, l_field_name, 
			       H5BLOCK_FIELD_ZCOORD_NAME, sizeof (H5BLOCK_FIELD_ZCOORD_NAME),
                               coords, *n_coords));
}

#define h5bl_get_fieldzcoords FC_GLOBAL (	 \
		h5bl_get_fieldzcoords,		 \
		H5BL_GET_FIELDZCOORDS)
h5_int64_t
h5bl_get_fieldzcoords (
	const h5_int64_t* const fh,
	const char* const field_name,
	h5_float64_t* const coords,
	const h5_int64_t* n_coords,
	const int l_field_name
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_err_t,
		      "fh=%p"
		      ", field_name='%.*s'"
                      "coords=%p, n_coords=%llu",
                      (h5_file_p)f,
                      l_field_name, field_name,
                      coords, (long long unsigned)n_coords);
        H5_API_RETURN (get_field_coords (
                               f, 2,
			       field_name, l_field_name,
			       H5BLOCK_FIELD_ZCOORD_NAME, sizeof (H5BLOCK_FIELD_ZCOORD_NAME),
                               coords, *n_coords));
}
