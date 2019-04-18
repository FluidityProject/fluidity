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
#include "h5core/h5_file_attribs.h"
#include "h5core/h5_step_attribs.h"

/*
   __ _ _              _   _        _ _           _            
  / _(_) | ___    __ _| |_| |_ _ __(_) |__  _   _| |_ ___  ___ 
 | |_| | |/ _ \  / _` | __| __| '__| | '_ \| | | | __/ _ \/ __|
 |  _| | |  __/ | (_| | |_| |_| |  | | |_) | |_| | ||  __/\__ \
 |_| |_|_|\___|  \__,_|\__|\__|_|  |_|_.__/ \__,_|\__\___||___/

   __ _ _   _  ___ _ __ _   _ 
  / _` | | | |/ _ \ '__| | | |
 | (_| | |_| |  __/ |  | |_| |
  \__, |\__,_|\___|_|   \__, |
     |_|                |___/
*/

#define h5_getnfileattribs FC_GLOBAL(		\
                h5_getnfileattribs,		\
                H5_GETNFILEATTRIBS)
h5_int64_t
h5_getnfileattribs (
	const h5_int64_t* const fh
	) {
	h5_file_t f = h5_filehandlefor2c(fh);
	H5_API_ENTER (h5_int64_t, "fh=%p", (h5_file_p)f);
        H5_API_RETURN (h5_get_num_file_attribs (f));
}

#define h5_getfileattribinfo FC_GLOBAL(	 \
                h5_getfileattribinfo,            \
                H5_GETFILEATTRIBINFO)
h5_int64_t
h5_getfileattribinfo (
	const h5_int64_t* const fh,
	const h5_int64_t* attrib_idx,
	char* attrib_name,
        h5_int64_t* attrib_type,
	h5_int64_t* attrib_nelem,
	const int l_attrib_name
        ) {
	h5_file_t f = h5_filehandlefor2c(fh);
	H5_API_ENTER (h5_int64_t,
		      "fh=%p, "
		      "attrib_idx=%lld, "
		      "attrib_name=%p, "
		      "attrib_type=%p, "
		      "attrib_nelem=%p",
		      (h5_file_p)fh,
		      (long long)*attrib_idx,
		      attrib_name, attrib_type, attrib_nelem);
        h5_int64_t h5err = h5_get_file_attrib_info_by_idx (
                f,
                *attrib_idx - 1,
                attrib_name, l_attrib_name,
                attrib_type,
                (h5_size_t*)attrib_nelem);
	h5_strc2for (attrib_name, l_attrib_name);
	H5_API_RETURN (h5err);
}

#define h5_getfileattribinfo_by_name FC_GLOBAL(	 \
                h5_getfileattribinfo_by_name,            \
                H5_GETFILEATTRIBINFO_BY_NAME)
h5_int64_t
h5_getfileattribinfo_by_name (
	const h5_int64_t* const fh,
	const char* const _name,
        h5_int64_t* const _type,
	h5_int64_t* const _nelem,
	const int l_name
        ) {
	h5_file_t f = h5_filehandlefor2c(fh);
	H5_API_ENTER (h5_int64_t,
		      "fh=%p, "
		      "name=%.*s, "
		      "type=%p, "
		      "nelem=%p",
		      (h5_file_p)fh,
		      l_name, _name, _type, _nelem);
	char* name = h5_strdupfor2c (_name, l_name);
        h5_int64_t h5err = h5_get_file_attrib_info_by_name (
                f,
                name,
                _type,
                (h5_size_t*)_nelem);
	H5_API_RETURN (h5err);
}

/*
   __ _ _              _   _        _ _           _            
  / _(_) | ___    __ _| |_| |_ _ __(_) |__  _   _| |_ ___  ___ 
 | |_| | |/ _ \  / _` | __| __| '__| | '_ \| | | | __/ _ \/ __|
 |  _| | |  __/ | (_| | |_| |_| |  | | |_) | |_| | ||  __/\__ \
 |_| |_|_|\___|  \__,_|\__|\__|_|  |_|_.__/ \__,_|\__\___||___/

  _    __    
 (_)  / /__  
 | | / / _ \ 
 | |/ / (_) |
 |_/_/ \___/ 
*/
static inline h5_int64_t
write_file_attrib (
        const h5_file_t f,
	const char* name,
	const int l_name,
        const hid_t type,
	const void* buffer,
	const hsize_t l_buffer
        ) {
	char *name2 = h5_strdupfor2c (name, l_name);
	h5_int64_t herr = h5_write_file_attrib (f, name2, type, buffer, l_buffer );
	free (name2);
        return herr;
}

static inline h5_int64_t
read_file_attrib (
	const h5_file_t f,
	const char* name,
	const int l_name,
        const hid_t type,
	void* const buffer
	) {
	char* name2 = h5_strdupfor2c ( name, l_name );
	h5_int64_t herr = h5_read_file_attrib (f, name2, type, buffer);
	free (name2);
	return herr;
}

#define h5_writefileattrib_string FC_GLOBAL (				\
                h5_writefileattrib_string,                              \
                H5_WRITEFILEATTRIB_STRING)
h5_int64_t
h5_writefileattrib_string (
	h5_int64_t *const fh,
	const char *name,
	const char *buffer,
	const int l_name,
	const int l_buffer
	) {
	h5_file_t f = h5_filehandlefor2c(fh);
	H5_API_ENTER (h5_int64_t,
		      "fh=%p, name='%.*s', buffer='%.*s'",
		      (h5_file_p)f, l_name, name, l_buffer, buffer);
	char *buffer2 = h5_strdupfor2c (buffer, l_buffer);
	h5_int64_t herr = write_file_attrib (
		f, name, l_name, H5_STRING_T, buffer2, strlen(buffer2)+1 );
	free (buffer2);
	H5_API_RETURN (herr);
}

#define h5_readfileattrib_string FC_GLOBAL (				\
                h5_readfileattrib_string,                               \
                H5_READFILEATTRIB_STRING)
h5_int64_t
h5_readfileattrib_string (
	h5_int64_t *const fh,
	const char *name,
	char *buffer,
	const int l_name,
	const int l_buffer
	) {
	h5_file_t f = h5_filehandlefor2c(fh);
	H5_API_ENTER (h5_int64_t,
		      "fh=%p, name='%.*s', buffer='%.*s'",
		      (h5_file_p)f, l_name, name, l_buffer, buffer);
	h5_int64_t herr = read_file_attrib (f, name, l_name, H5_STRING_T, buffer);
	h5_strc2for (buffer, l_buffer);
	H5_API_RETURN (herr);
}

#define h5_writefileattrib_r8 FC_GLOBAL (	\
                h5_writefileattrib_r8,          \
                H5_WRITEFILEATTRIB_R8)
h5_int64_t
h5_writefileattrib_r8 (
	const h5_int64_t *const fh,
	const char *name,
	const h5_int64_t *buffer,
	const h5_int64_t *nelem,
	const int l_name
	) {
	h5_file_t f = h5_filehandlefor2c(fh);
	H5_API_ENTER (h5_int64_t,
		      "fh=%p, name='%.*s', buffer=%p, nelem=%lld",
		      (h5_file_p)f, l_name, name, buffer, (long long)*nelem);
        H5_API_RETURN (write_file_attrib(
                               f,
                               name, l_name,
                               H5_FLOAT64_T,
                               buffer, (hsize_t)*nelem));
}

#define h5_readfileattrib_r8 FC_GLOBAL (	\
                h5_readfileattrib_r8,		\
                H5_READFILEATTRIB_R8 )
h5_int64_t
h5_readfileattrib_r8 (
	h5_int64_t *const fh,
	const char *name,
	h5_int64_t *buffer,
	const int l_name
	) {
	h5_file_t f = h5_filehandlefor2c(fh);
	H5_API_ENTER (h5_int64_t,
		      "(h5_file_p)fh=%p, name='%.*s', buffer=%p",
		      (h5_file_p)f, l_name, name, buffer);
        H5_API_RETURN (read_file_attrib(
		f,
                name, l_name, 
                H5_FLOAT64_T,
                (void*)buffer));
}

#define h5_writefileattrib_r4 FC_GLOBAL (	\
                h5_writefileattrib_r4,		\
                H5_WRITEFILEATTRIB_R4 )
h5_int64_t
h5_writefileattrib_r4 (
	const h5_int64_t *const fh,
	const char *name,
	const h5_int32_t *buffer,
	const h5_int64_t *nelem,
	const int l_name
	) {
	h5_file_t f = h5_filehandlefor2c(fh);
	H5_API_ENTER (h5_int64_t,
		      "fh=%p, name='%.*s', buffer=%p, nelem=%lld",
		      (h5_file_p)f, l_name, name, buffer, (long long)*nelem);
        H5_API_RETURN (write_file_attrib(
                               f,
                               name, l_name,
                               H5_FLOAT32_T,
                               buffer, (hsize_t)*nelem));
}

#define h5_readfileattrib_r4 FC_GLOBAL (	\
                h5_readfileattrib_r4,		\
                H5_READFILEATTRIB_R4 )
h5_int64_t
h5_readfileattrib_r4 (
	const h5_int64_t *const fh,
	const char *name,
	h5_int32_t *buffer,
	const int l_name
	) {
	h5_file_t f = h5_filehandlefor2c(fh);
	H5_API_ENTER (h5_int64_t,
		      "f=%p, name='%.*s', buffer=%p",
		      (h5_file_p)f, l_name, name, buffer);
        H5_API_RETURN (read_file_attrib(
		f,
                name, l_name, 
                H5_FLOAT32_T,
                buffer));
}

#define h5_writefileattrib_i8 FC_GLOBAL (	\
                h5_writefileattrib_i8,          \
                H5_WRITEFILEATTRIB_I8)
h5_int64_t
h5_writefileattrib_i8 (
	const h5_int64_t *const fh,
	const char *name,
	const h5_int64_t *buffer,
	const h5_int64_t *nelem,
	const int l_name
	) {
	h5_file_t f = h5_filehandlefor2c(fh);
	H5_API_ENTER (h5_int64_t,
		      "fh=%p, name='%.*s', buffer=%p, nelem=%lld",
		      (h5_file_p)f, l_name, name, buffer, (long long)*nelem);
        H5_API_RETURN (write_file_attrib(
                               f,
                               name, l_name,
                               H5_INT64_T,
                               buffer, (hsize_t)*nelem));
}

#define h5_readfileattrib_i8 FC_GLOBAL (	\
                h5_readfileattrib_i8,		\
                H5_READFILEATTRIB_I8 )
h5_int64_t
h5_readfileattrib_i8 (
	const h5_int64_t *const fh,
	const char *name,
	h5_int64_t *buffer,
	const int l_name
	) {
	h5_file_t f = h5_filehandlefor2c(fh);
	H5_API_ENTER (h5_int64_t,
		      "fh=%p, name='%.*s', buffer=%p",
		      (h5_file_p)f, l_name, name, buffer);
        H5_API_RETURN (read_file_attrib(
		f,
                name, l_name, 
                H5_INT64_T,
                buffer));
}

#define h5_writefileattrib_i4 FC_GLOBAL (	\
                h5_writefileattrib_i4,          \
                H5_WRITEFILEATTRIB_I4 )
h5_int64_t
h5_writefileattrib_i4 (
	h5_int64_t *const fh,
	const char *name,
	const h5_int32_t *buffer,
	const h5_int64_t *nelem,
	const int l_name
	) {
	h5_file_t f = h5_filehandlefor2c(fh);
	H5_API_ENTER (h5_int64_t,
		      "fh=%p, name='%.*s', buffer[12]=%d, %x, %x, nelem=%lld",
		      (h5_file_p)f, l_name, name, buffer[12], buffer[13], buffer[14], (long long)*nelem);
        H5_API_RETURN (write_file_attrib(
                               f,
                               name, l_name,
                               H5_INT32_T,
                               buffer, (hsize_t)*nelem));
}

#define h5_readfileattrib_i4 FC_GLOBAL (	\
                h5_readfileattrib_i4,           \
                H5_READFILEATTRIB_I4 )
h5_int64_t
h5_readfileattrib_i4 (
	h5_int64_t *const fh,
	const char *name,
	h5_int32_t *buffer,
	const int l_name
	) {
	h5_file_t f = h5_filehandlefor2c(fh);
	H5_API_ENTER (h5_int64_t,
		      "fh=%p, name='%.*s', buffer=%p",
		      (h5_file_p)f, l_name, name, buffer);
        H5_API_RETURN (read_file_attrib(
		f,
                name, l_name, 
                H5_INT32_T,
                buffer));
}

/*
      _                     _   _        _ _           _            
  ___| |_ ___ _ __     __ _| |_| |_ _ __(_) |__  _   _| |_ ___  ___ 
 / __| __/ _ \ '_ \   / _` | __| __| '__| | '_ \| | | | __/ _ \/ __|
 \__ \ ||  __/ |_) | | (_| | |_| |_| |  | | |_) | |_| | ||  __/\__ \
 |___/\__\___| .__/   \__,_|\__|\__|_|  |_|_.__/ \__,_|\__\___||___/
             |_|                                                    
 
   __ _ _   _  ___ _ __ _   _ 
  / _` | | | |/ _ \ '__| | | |
 | (_| | |_| |  __/ |  | |_| |
  \__, |\__,_|\___|_|   \__, |
     |_|                |___/
*/
#define h5_getnstepattribs FC_GLOBAL( \
                h5_getnstepattribs,  \
                H5_GETNSTEPATTRIBS)
h5_int64_t
h5_getnstepattribs (
	const h5_int64_t* fh
	) {
	h5_file_t f = h5_filehandlefor2c(fh);
	H5_API_ENTER (h5_int64_t,
                      "fh=%p",
                      (h5_file_p)f);
        H5_API_RETURN (h5_get_num_iteration_attribs (f));
}

#define h5_getstepattribinfo FC_GLOBAL(	 \
		h5_getstepattribinfo,		 \
		H5_GETSTEPATTRIBINFO)
h5_int64_t
h5_getstepattribinfo (
	const h5_int64_t* fh,
	const h5_int64_t* attrib_idx,
	char* attrib_name,
        h5_int64_t* attrib_type,
	h5_int64_t* attrib_nelem,
	const int l_attrib_name
        ) {
	h5_file_t f = h5_filehandlefor2c(fh);
	H5_API_ENTER (h5_int64_t,
		      "f=%p, "
		      "attrib_idx=%lld, "
		      "attrib_name=%p, "
		      "attrib_type=%p, "
		      "attrib_nelem=%p",
		      (h5_file_p)f,
		      (long long)*attrib_idx,
		      attrib_name, attrib_type, attrib_nelem);
        h5_int64_t h5err = h5_get_iteration_attrib_info_by_idx (
                f,
                *attrib_idx - 1,
                attrib_name, l_attrib_name,
                attrib_type,
                (h5_size_t*)attrib_nelem);
	h5_strc2for (attrib_name, l_attrib_name);
	H5_API_RETURN (h5err);
}

#define h5_getstepattribinfo_by_name FC_GLOBAL(	 \
                h5_getstepattribinfo_by_name,            \
                H5_GETSTEPATTRIBINFO_BY_NAME)
h5_int64_t
h5_getstepattribinfo_by_name (
	const h5_int64_t* const fh,
	const char* const _name,
        h5_int64_t* const _type,
	h5_int64_t* const _nelem,
	const int l_name
        ) {
	h5_file_t f = h5_filehandlefor2c(fh);
	H5_API_ENTER (h5_int64_t,
		      "fh=%p, "
		      "name=%.*s, "
		      "type=%p, "
		      "nelem=%p",
		      (h5_file_p)fh,
		      l_name, _name, _type, _nelem);
	char* name = h5_strdupfor2c (_name, l_name);
        h5_int64_t h5err = h5_get_iteration_attrib_info_by_name (
                f,
                name,
                _type,
                (h5_size_t*)_nelem);
	H5_API_RETURN (h5err);
}

/*
      _                     _   _        _ _           _            
  ___| |_ ___ _ __     __ _| |_| |_ _ __(_) |__  _   _| |_ ___  ___ 
 / __| __/ _ \ '_ \   / _` | __| __| '__| | '_ \| | | | __/ _ \/ __|
 \__ \ ||  __/ |_) | | (_| | |_| |_| |  | | |_) | |_| | ||  __/\__ \
 |___/\__\___| .__/   \__,_|\__|\__|_|  |_|_.__/ \__,_|\__\___||___/
             |_|                                                    

  _    __    
 (_)  / /__  
 | | / / _ \ 
 | |/ / (_) |
 |_/_/ \___/ 
*/

static inline h5_int64_t
write_iteration_attrib (
        const h5_file_t fh,
	const char* name,
	const int l_name,
        const hid_t type,
	const void* buffer,
	const hsize_t l_buffer
        ) {
	char *name2 = h5_strdupfor2c (name, l_name);
	h5_int64_t herr = h5_write_iteration_attrib (
		fh, name2, type, buffer, l_buffer );
	free (name2);
        return herr;
}

static inline h5_int64_t
read_iteration_attrib (
	const h5_file_t fh,
	const char* name,
	const int l_name,
        const hid_t type,
	void* const buffer
	) {
	char* name2 = h5_strdupfor2c ( name, l_name );
	h5_int64_t herr = h5_read_iteration_attrib (fh, name2, type, buffer);
	free (name2);
	return herr;
}

#define h5_writestepattrib_string FC_GLOBAL (				\
                h5_writestepattrib_string,                              \
                H5_WRITESTEPATTRIB_STRING)
h5_int64_t
h5_writestepattrib_string (
	const h5_int64_t *const fh,
	const char *name,
	const char *buffer,
	const int l_name,
	const int l_buffer
	) {
	h5_file_t f = h5_filehandlefor2c(fh);
	H5_API_ENTER (h5_int64_t,
		      "f=%p, name='%.*s', buffer='%.*s'",
		      (h5_file_p)f, l_name, name, l_buffer, buffer);
	char *buffer2 = h5_strdupfor2c (buffer, l_buffer);
	h5_int64_t herr = write_iteration_attrib (
		f, name, l_name, H5_STRING_T, buffer2, strlen(buffer2)+1 );
	free (buffer2);
	H5_API_RETURN (herr);
}

#define h5_readstepattrib_string FC_GLOBAL (				\
                h5_readstepattrib_string,                               \
                H5_READSTEPATTRIB_STRING)
h5_int64_t
h5_readstepattrib_string (
	const h5_int64_t *const fh,
	const char *name,
	char *buffer,
	const int l_name,
	const int l_buffer
	) {
	h5_file_t f = h5_filehandlefor2c(fh);
	H5_API_ENTER (h5_int64_t,
		      "f=%p, name='%.*s', buffer='%.*s'",
		      (h5_file_p)f, l_name, name, l_buffer, buffer);
	h5_int64_t herr = read_iteration_attrib (
		f, name, l_name, H5_STRING_T, buffer);
	h5_strc2for (buffer, l_buffer);
	H5_API_RETURN (herr);
}

#define h5_writestepattrib_r8 FC_GLOBAL (	\
                h5_writestepattrib_r8,          \
                H5_WRITESTEPATTRIB_R8)
h5_int64_t
h5_writestepattrib_r8 (
	const h5_int64_t *const fh,
	const char *name,
	const h5_float64_t *buffer,
	const h5_int64_t *nelem,
	const int l_name
	) {
	h5_file_t f = h5_filehandlefor2c(fh);
	H5_API_ENTER (h5_int64_t,
		      "f=%p, name='%.*s', buffer=%p, nelem=%lld",
		      (h5_file_p)f, l_name, name, buffer, (long long)*nelem);
        H5_API_RETURN (write_iteration_attrib(
                               f,
                               name, l_name,
                               H5_FLOAT64_T,
                               buffer, (hsize_t)*nelem));
}

#define h5_readstepattrib_r8 FC_GLOBAL (	\
                h5_readstepattrib_r8,		\
                H5_READSTEPATTRIB_R8 )
h5_int64_t
h5_readstepattrib_r8 (
	const h5_int64_t *const fh,
	const char *name,
	h5_float64_t *buffer,
	const int l_name
	) {
	h5_file_t f = h5_filehandlefor2c(fh);
	H5_API_ENTER (h5_int64_t,
		      "f=%p, name='%.*s', buffer=%p",
		      (h5_file_p)f, l_name, name, buffer);
        H5_API_RETURN (read_iteration_attrib(
		f,
                name, l_name, 
                H5_FLOAT64_T,
                (void*)buffer));
}

#define h5_writestepattrib_r4 FC_GLOBAL (	\
                h5_writestepattrib_r4,		\
                H5_WRITESTEPATTRIB_R4 )
h5_int64_t
h5_writestepattrib_r4 (
	h5_int64_t *const fh,
	const char *name,
	const h5_float64_t *buffer,
	const h5_int64_t *nelem,
	const int l_name
	) {
	h5_file_t f = h5_filehandlefor2c(fh);
	H5_API_ENTER (h5_int64_t,
		      "f=%p, name='%.*s', buffer=%p, nelem=%lld",
		      (h5_file_p)f, l_name, name, buffer, (long long)*nelem);
        H5_API_RETURN (write_iteration_attrib(
                               f,
                               name, l_name,
                               H5_FLOAT32_T,
                               buffer, (hsize_t)*nelem));
}

#define h5_readstepattrib_r4 FC_GLOBAL (	\
                h5_readstepattrib_r4,		\
                H5_READSTEPATTRIB_R4 )
h5_int64_t
h5_readstepattrib_r4 (
	const h5_int64_t *const fh,
	const char *name,
	h5_float64_t *buffer,
	const int l_name
	) {
	h5_file_t f = h5_filehandlefor2c(fh);
	H5_API_ENTER (h5_int64_t,
		      "f=%p, name='%.*s', buffer=%p",
		      (h5_file_p)f, l_name, name, buffer);
        H5_API_RETURN (read_iteration_attrib(
		f,
                name, l_name, 
                H5_FLOAT32_T,
                buffer));
}

#define h5_writestepattrib_i8 FC_GLOBAL (	\
                h5_writestepattrib_i8,          \
                H5_WRITESTEPATTRIB_I8)
h5_int64_t
h5_writestepattrib_i8 (
	h5_int64_t *const fh,
	const char *name,
	const h5_float64_t *buffer,
	const h5_int64_t *nelem,
	const int l_name
	) {
	h5_file_t f = h5_filehandlefor2c(fh);
	H5_API_ENTER (h5_int64_t,
		      "f=%p, name='%.*s', buffer=%p, nelem=%lld",
		      (h5_file_p)f, l_name, name, buffer, (long long)*nelem);
        H5_API_RETURN (write_iteration_attrib(
                               f,
                               name, l_name,
                               H5_INT64_T,
                               buffer, (hsize_t)*nelem));
}

#define h5_readstepattrib_i8 FC_GLOBAL (	\
                h5_readstepattrib_i8,		\
                H5_READSTEPATTRIB_I8 )
h5_int64_t
h5_readstepattrib_i8 (
	const h5_int64_t *const fh,
	const char *name,
	h5_float64_t *buffer,
	const int l_name
	) {
	h5_file_t f = h5_filehandlefor2c(fh);
	H5_API_ENTER (h5_int64_t,
		      "f=%p, name='%.*s', buffer=%p",
		      (h5_file_p)f, l_name, name, buffer);
        H5_API_RETURN (read_iteration_attrib(
		f,
                name, l_name, 
                H5_INT64_T,
                buffer));
}

#define h5_writestepattrib_i4 FC_GLOBAL (	\
                h5_writestepattrib_i4,          \
                H5_WRITESTEPATTRIB_I4 )
h5_int64_t
h5_writestepattrib_i4 (
	const h5_int64_t *const fh,
	const char *name,
	const h5_float64_t *buffer,
	const h5_int64_t *nelem,
	const int l_name
	) {
	h5_file_t f = h5_filehandlefor2c(fh);
	H5_API_ENTER (h5_int64_t,
		      "f=%p, name='%.*s', buffer=%p, nelem=%lld",
		      (h5_file_p)f, l_name, name, buffer, (long long)*nelem);
        H5_API_RETURN (write_iteration_attrib(
                               f,
                               name, l_name,
                               H5_INT32_T,
                               buffer, (hsize_t)*nelem));
}

#define h5_readstepattrib_i4 FC_GLOBAL (	\
                h5_readstepattrib_i4,           \
                H5_READSTEPATTRIB_I4 )
h5_int64_t
h5_readstepattrib_i4 (
	const h5_int64_t *const fh,
	const char *name,
	h5_float64_t *buffer,
	const int l_name
	) {
	h5_file_t f = h5_filehandlefor2c(fh);
	H5_API_ENTER (h5_int64_t,
		      "f=%p, name='%.*s', buffer=%p",
		      (h5_file_p)f, l_name, name, buffer);
        H5_API_RETURN (read_iteration_attrib(
		f,
                name, l_name, 
                H5_INT32_T,
                buffer));
}
