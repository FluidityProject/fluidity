/*
  Copyright (c) 2006-2015, The Regents of the University of California,
  through Lawrence Berkeley National Laboratory (subject to receipt of any
  required approvals from the U.S. Dept. of Energy) and the Paul Scherrer
  Institut (Switzerland).  All rights reserved.

  License: see file COPYING in top level of source distribution.
*/

#include "h5_private.h"
#include "h5core/h5_log.h"
#include "h5core/h5b_model.h"
#include "h5core/h5b_io.h"

#define h5bl_hasfielddata FC_MANGLING (				  \
                h5bl_hasfielddata,				  \
                H5BL_HASFIELDDATA )
h5_int64_t
h5bl_hasfielddata (
	const h5_int64_t* const fh
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_int64_t,
                      "fh=%p",
                      (h5_file_p)f);
	H5_API_RETURN (h5b_has_field_data ( f ));
}


#define h5bl_3d_hasview FC_MANGLING (				\
                h5bl_hasview,                                   \
                H5BL_HASVIEW )
h5_int64_t
h5bl_3d_hasview (
	const h5_int64_t* const fh
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_int64_t,
                      "fh=%p",
                      (h5_file_p)f);
	H5_API_RETURN (h5b_3d_has_view ( f ));
}


#define h5bl_3d_setview FC_MANGLING (					\
                h5bl_3d_setview,                                        \
                H5BL_3D_SETVIEW )
h5_int64_t
h5bl_3d_setview (
	const h5_int64_t* const fh,
	const h5_int64_t* const i_start,	/*!< start index of i */
	const h5_int64_t* const i_end,  	/*!< end index of i */
	const h5_int64_t* const j_start,	/*!< start index of j */
	const h5_int64_t* const j_end,          /*!< end index of j */
	const h5_int64_t* const k_start,	/*!< start index of k */
	const h5_int64_t* const k_end		/*!< end index of k */
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_int64_t,
		      "fh=%p, "
		      "i_start=%lld, i_end=%lld, "
		      "j_start=%lld, j_end=%lld, "
		       "k_start=%lld, k_end=%lld",
		      (h5_file_p)f,
		      (long long)*i_start, (long long)*i_end,
		      (long long)*j_start, (long long)*j_end,
		      (long long)*k_start, (long long)*k_end);
	H5_API_RETURN(h5b_3d_set_view (
		f,
		*i_start-1, *i_end-1,
		*j_start-1, *j_end-1,
		*k_start-1, *k_end-1,
                0));
}

#define h5bl_3d_getview FC_MANGLING (					\
                h5bl_3d_getview,                                        \
                H5BL_3D_GETVIEW )
h5_int64_t
h5bl_3d_getview (
	const h5_int64_t* const fh,
	h5_int64_t* const i_start,
	h5_int64_t* const i_end,
	h5_int64_t* const j_start,
	h5_int64_t* const j_end,
	h5_int64_t* const k_start,
	h5_int64_t* const k_end
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_int64_t,
		      "fh=%p, "
		      "i_start=%p, i_end=%p, "
		      "j_start=%p, j_end=%p, "
		      "k_start=%p, k_end=%p",
		      (h5_file_p)f,
		      i_start, i_end,
		      j_start, j_end,
		      k_start, k_end);
	TRY (h5b_3d_get_view (
                     f,
                     (h5_size_t*)i_start, (h5_size_t*)i_end,
                     (h5_size_t*)j_start, (h5_size_t*)j_end,
                     (h5_size_t*)k_start, (h5_size_t*)k_end ));
        *i_start += 1;
        *i_end +=   1;
        *j_start += 1;
        *j_end +=   1;
        *k_start += 1;
        *k_end +=   1;
        H5_API_RETURN (H5_SUCCESS);
}



#define h5bl_3d_getreducedview FC_MANGLING (				\
                h5bl_3d_getreducedview,                                 \
                H5BL_3D_GETREDUCEDVIEW )
h5_int64_t
h5bl_3d_getreducedview (
	const h5_int64_t* const fh,
	h5_int64_t* const i_start, 
	h5_int64_t* const i_end,
	h5_int64_t* const j_start,
	h5_int64_t* const j_end,
	h5_int64_t* const k_start,
	h5_int64_t* const k_end
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_int64_t,
		      "fh=%p, "
		      "i_start=%p, i_end=%p, "
		      "j_start=%p, j_end=%p, "
		      "k_start=%p, k_end=%p",
		      (h5_file_p)f,
		      i_start, i_end,
		      j_start, j_end,
		      k_start, k_end);
	TRY (h5b_3d_get_reduced_view (
                     f,
                     (h5_size_t*)i_start, (h5_size_t*)i_end,
                     (h5_size_t*)j_start, (h5_size_t*)j_end,
                     (h5_size_t*)k_start, (h5_size_t*)k_end));
        *i_start += 1;
        *i_end +=   1;
        *j_start += 1;
        *j_end +=   1;
        *k_start += 1;
        *k_end +=   1;
        H5_API_RETURN (H5_SUCCESS);
}

#define h5bl_3d_setchunk FC_MANGLING (				\
                h5bl_3d_setchunk,                               \
                H5BL_3D_SETCHUNK )
h5_int64_t
h5bl_3d_setchunk (
	const h5_int64_t* const fh,
	const h5_int64_t* const i,
	const h5_int64_t* const j,
	const h5_int64_t* const k
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_int64_t,
		      "fh=%p, i=%lld, j=%lld, k=%lld",
		      (h5_file_p)f, (long long)*i, (long long)*j, (long long)*k);
	H5_API_RETURN(h5b_3d_set_chunk (f, *i, *j, *k));
}

#define h5bl_3d_getchunk FC_MANGLING (				\
                h5bl_3d_getchunk,                               \
                H5BL_3D_GETCHUNK )
h5_int64_t
h5bl_3d_getchunk (
	const h5_int64_t* const fh,
	const char* const field_name,
	h5_int64_t* const i,
	h5_int64_t* const j,
	h5_int64_t* const k,
	const int l_field_name
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_int64_t,
		      "fh=%p, i=%p, j=%p, k=%p",
		      (h5_file_p)f, i, j, k);
	char *field_name2 = h5_strdupfor2c (field_name, l_field_name);
	h5_int64_t h5err = h5b_3d_get_chunk (f, field_name2, (h5_size_t*)i, (h5_size_t*)j, (h5_size_t*)k);
	free (field_name2);
	H5_API_RETURN (h5err);
}


#if defined(H5_HAVE_PARALLEL)
#define h5bl_3d_setgrid FC_MANGLING (		\
		h5bl_3d_setgrid,		\
		h5bl_3d_setgrid)
h5_int64_t
h5bl_3d_setgrid (
	const h5_int64_t* const fh,
	const h5_int64_t* const i,
	const h5_int64_t* const j,
	const h5_int64_t* const k
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_int64_t,
		      "fh=%p, i=%lld, j=%lld, k=%lld",
		      (h5_file_p)f, (long long)*i, (long long)*j, (long long)*k);
	H5_API_RETURN(h5b_3d_set_grid (f, *i, *j, *k));
}

#define h5bl_3d_getgrid FC_MANGLING (	        \
		h5bl_3d_getgrid,		\
		H5BL_3D_GETGRID)
h5_int64_t
h5bl_3d_getgrid (
	const h5_int64_t* const fh,
	const h5_int64_t* const proc,
	h5_int64_t* const i,
	h5_int64_t* const j,
	h5_int64_t* const k
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_int64_t,
		      "fh=%p, proc=%lld, i=%p, j=%p, k=%p",
		      (h5_file_p)f, (long long)proc, i, j, k);
	H5_API_RETURN(h5b_3d_get_grid_coords (f, (int)*proc, i, j, k));
}

#define h5bl_3d_setdims FC_MANGLING (		\
		h5bl_3d_setdims,                \
		H5BL_3D_SETDIMS)
h5_int64_t
h5bl_3d_setdims (
	const h5_int64_t* const fh,
	const h5_int64_t* const i,
	const h5_int64_t* const j,
	const h5_int64_t* const k
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_int64_t,
		      "fh=%p, i=%lld, j=%lld, k=%lld",
		      (h5_file_p)f, (long long)*i, (long long)*j, (long long)*k);
	H5_API_RETURN(h5b_3d_set_dims (f, *i, *j, *k));
}

#endif

#define h5bl_3d_sethalo FC_MANGLING (		\
         	h5bl_3d_sethalo,                \
		H5BL_3D_SETHALO)
h5_int64_t
h5bl_3d_sethalo (
	const h5_int64_t* const fh,
	const h5_int64_t* const i,
	const h5_int64_t* const j,
	const h5_int64_t* const k
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_int64_t,
		      "fh=%p, i=%lld, j=%lld, k=%lld",
		      (h5_file_p)f, (long long)*i, (long long)*j, (long long)*k);
	H5_API_RETURN(h5b_3d_set_halo (f, *i, *j, *k));
}

#define h5bl_getnumfields FC_MANGLING (					\
                h5bl_getnumfields,                                      \
                H5BL_GETNUMFIELDS )
h5_int64_t
h5bl_getnumfields (
	const h5_int64_t* const fh
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_int64_t,
                      "fh=%p",
                      (h5_file_p)f);
	H5_API_RETURN (h5b_get_num_fields (f));
}

#define h5bl_getfieldinfo FC_MANGLING (					\
                h5bl_getfieldinfo,                                      \
                H5BL_GETFIELDINFO )
h5_int64_t
h5bl_getfieldinfo (
	const h5_int64_t* const fh,
	const h5_int64_t *idx,
	char* const name,
	h5_size_t* const field_rank,
	h5_size_t* const field_dims,
	h5_size_t* const elem_rank,
	h5_int64_t* const type,
	const int l_name
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_int64_t,
		      "fh=%p, idx=%lld, "
		      "name=%*s,"
		      "field_rank=%p, field_dims=%p, elem_rank=%p, type=%p",
		      (h5_file_p)f, (long long)*idx, l_name, name,
		      field_rank, field_dims, elem_rank, type);
	h5_int64_t herr = h5b_get_field_info (
		f, *idx - 1, name, (h5_size_t)l_name,
		field_rank, field_dims, elem_rank, type );
	h5_strc2for ( name, l_name );
	H5_API_RETURN(herr);
}

#define h5bl_getfieldinfobyname FC_MANGLING (	\
		h5bl_getfieldinfobyname,        \
		H5BL_GETFIELDINFOBYNAME)
h5_int64_t
h5bl_getfieldinfobyname (
	const h5_int64_t* const fh,
	const char* const field_name,
	h5_size_t* const field_rank,
	h5_size_t* const field_dims,
	h5_size_t* const elem_rank,
	h5_int64_t* const type,
	const int l_field_name
	) {
	h5_file_t f = h5_filehandlefor2c (fh);
	H5_API_ENTER (h5_int64_t,
		      "fh=%p"
		      ", field_name=%*s,"
		      ", field_rank=%p, field_dims=%p, elem_rank=%p, type=%p",
		      (h5_file_p)f, l_field_name, field_name,
		      field_rank, field_dims, elem_rank, type);
	char *field_name2 = h5_strdupfor2c (field_name, l_field_name);
	h5_int64_t herr = h5b_get_field_info_by_name (
		f, field_name2,
		field_rank, field_dims, elem_rank, type );
	free (field_name2);
	H5_API_RETURN(herr);
}
