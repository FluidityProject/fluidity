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
#include "h5core/h5_file.h"

#include <hdf5.h>

#define h5_createprop_file FC_GLOBAL(    \
                h5_createprop_file,	   \
                H5_CREATEPROP_FILE)
h5_int64_t
h5_createprop_file (
        void
        ) {
        H5_API_ENTER (h5_int64_t, "%s", "");
        H5_API_RETURN ((h5_int64_t)h5_create_prop (H5_PROP_FILE));
}

#if defined(H5_HAVE_PARALLEL)
#define h5_setprop_file_mpio FC_GLOBAL( \
                h5_setprop_file_mpio,	  \
                H5_SETPROP_FILE_MPIO)
h5_int64_t
h5_setprop_file_mpio (
        h5_int64_t* _prop,
	MPI_Fint* _comm
        ) {
        H5_API_ENTER (h5_int64_t,
                      "prop=%lld, comm=%lld",
                      (long long int)*_prop, (long long int)*_comm);
        h5_prop_t prop = (h5_prop_t)*_prop;
        MPI_Comm comm = MPI_Comm_f2c (*_comm);
        H5_API_RETURN ((h5_int64_t)h5_set_prop_file_mpio_collective (prop, &comm));
}

#define h5_setprop_file_mpio_collective FC_GLOBAL( \
                h5_setprop_file_mpio_collective,     \
                H5_SETPROP_FILE_MPIO_COLLECTIVE)
h5_int64_t
h5_setprop_file_mpio_collective (
        h5_int64_t* _prop,
	MPI_Fint* _comm
        ) {
        H5_API_ENTER (h5_int64_t,
                      "prop=%lld, comm=%lld",
                      (long long int)*_prop, (long long int)*_comm);
        h5_prop_t prop = (h5_prop_t)*_prop;
        MPI_Comm comm = MPI_Comm_f2c (*_comm);
        H5_API_RETURN ((h5_int64_t)h5_set_prop_file_mpio_collective (prop, &comm));
}

#define h5_setprop_file_mpio_independent FC_GLOBAL( \
                h5_setprop_file_mpio_independent,     \
                H5_SETPROP_FILE_MPIO_INDEPENDENT)
h5_int64_t
h5_setprop_file_mpio_independent (
        h5_int64_t* _prop,
	MPI_Fint* _comm
        ) {
        H5_API_ENTER (h5_int64_t,
                      "prop=%lld, comm=%lld",
                      (long long int)*_prop, (long long int)*_comm);
        h5_prop_t prop = (h5_prop_t)*_prop;
        MPI_Comm comm = MPI_Comm_f2c (*_comm);
        H5_API_RETURN ((h5_int64_t)h5_set_prop_file_mpio_independent (prop, &comm));
}

#if H5_VERSION_LE(1,8,12)
#define h5_setprop_file_mpio_posix FC_GLOBAL( \
                h5_setprop_file_mpio_posix,  \
                H5_SETPROP_FILE_MPIO_POSIX)
h5_int64_t
h5_setprop_file_mpio_posix (
        h5_int64_t* _prop,
	MPI_Fint* _comm
        ) {
        H5_API_ENTER (h5_int64_t,
                      "prop=%lld, comm=%lld",
                      (long long int)*_prop, (long long int)*_comm);
        h5_prop_t prop = (h5_prop_t)*_prop;
        MPI_Comm comm = MPI_Comm_f2c (*_comm);
        H5_API_RETURN ((h5_int64_t)h5_set_prop_file_mpio_posix (prop, &comm));
}
#endif

#endif

#define h5_setprop_file_corevfd FC_GLOBAL( \
                h5_setprop_file_corevfd,     \
                H5_SETPROP_FILE_COREVFD)
h5_int64_t
h5_setprop_file_corevfd (
        h5_int64_t* _prop,
	h5_int64_t* increment
        ) {
        H5_API_ENTER (h5_int64_t,
                      "prop=%lld, increment=%lld",
                      (long long int)*_prop, (long long int)*increment);
        h5_prop_t prop = (h5_prop_t)*_prop;
        H5_API_RETURN ((h5_int64_t)h5_set_prop_file_core_vfd (prop, *increment));
}

#define h5_setprop_file_align FC_GLOBAL (  \
                h5_setprop_file_align,	     \
                H5_SETPROP_FILE_ALIGN)
h5_int64_t
h5_setprop_file_align (
        h5_int64_t* _prop,
        h5_int64_t* align
        ) {
        H5_API_ENTER (h5_err_t,
                      "prop=%lld, align=%lld",
                      (long long int)*_prop, (long long int)*align);
        h5_prop_t prop = (h5_prop_t)*_prop;
        H5_API_RETURN (h5_set_prop_file_align (prop, *align));
}

#define h5_setprop_file_throttle FC_GLOBAL (		  \
                h5_setprop_file_throttle,                 \
                H5_SETPROP_FILE_THROTTLE)

h5_int64_t
h5_setprop_file_throttle (
        h5_int64_t* _prop,
        h5_int64_t* throttle
        ) {
        H5_API_ENTER (
                h5_err_t,
                "prop=%lld, throttle=%lld",
                (long long int)*_prop, (long long int)*throttle);
        h5_prop_t prop = (h5_prop_t)*_prop;
        H5_API_RETURN (h5_set_prop_file_throttle (prop, *throttle));
}

#define h5_closeprop FC_GLOBAL (		\
                h5_closeprop,                   \
                H5_CLOSEPROP)
h5_int64_t
h5_closeprop (
        h5_int64_t* _prop
        ) {
        H5_API_ENTER (h5_err_t,
                      "prop=%lld",
                      (long long int)*_prop);
        h5_prop_t prop = (h5_prop_t)*_prop;
        H5_API_RETURN (h5_close_prop (prop));
}

#define h5_openfile FC_GLOBAL( \
                h5_openfile,  \
                H5_OPENFILE)
h5_int64_t
h5_openfile (
	const char* _fname,
        h5_int64_t* _mode,
        h5_int64_t* _props,
	const int _len_fname
        ) {
        H5_API_ENTER (h5_int64_t,
                      "fname = %*s, mode=%lld, props=%lld",
                      _len_fname, _fname, (long long int)*_mode, (long long int)*_props);
        char* fname = h5_strdupfor2c (_fname, _len_fname);
        h5_int64_t mode = *_mode;
        h5_prop_t props = (h5_prop_t)*_props;
        h5_file_t f = h5_open_file2 (fname, mode, props);
        free (fname);
        H5_API_RETURN ((h5_int64_t)f);
}

#define h5_closefile FC_GLOBAL(                  \
                h5_closefile,			   \
                H5_CLOSEFILE)
h5_int64_t
h5_closefile (
	const h5_int64_t *f
	) {
	h5_file_t fh = h5_filehandlefor2c(f);
	H5_API_ENTER (h5_int64_t, "f=%p", (h5_file_p)fh);
	H5_API_RETURN (h5_close_file (fh));
}

#define h5_checkfile FC_GLOBAL(                  \
                h5_checkfile,			   \
                H5_CHECKFILE)
h5_int64_t
h5_checkfile (
	const h5_int64_t *f
	) {
	h5_file_t fh = h5_filehandlefor2c(f);
	H5_API_ENTER (h5_int64_t, "f=%p", (h5_file_p)fh);
	H5_API_RETURN (h5_check_filehandle (fh));
}

#define h5_flushfile FC_GLOBAL(                      \
                h5_flushfile,			       \
                H5_FLUSHFILE)
h5_int64_t
h5_flushfile (
	const h5_int64_t* f
	) {
	h5_file_t fh = h5_filehandlefor2c(f);
	H5_API_ENTER (h5_int64_t, "f=%p", (h5_file_p)fh);
	H5_API_RETURN (h5_flush_file (fh));
}

#define h5_flushstep FC_GLOBAL(		    \
                h5_flushstep,                       \
                H5_FLUSHSTEP)
h5_int64_t
h5_flushstep (
	const h5_int64_t* f
	) {
	h5_file_t fh = h5_filehandlefor2c(f);
	H5_API_ENTER (h5_int64_t, "f=%p", (h5_file_p)fh);
	H5_API_RETURN (h5_flush_iteration (fh));
}

#define h5_finalize FC_GLOBAL(		   \
                h5_finalize,                       \
                H5_FINALIZE)
h5_int64_t
h5_finalize (
	void
	) {
	H5_API_ENTER (h5_int64_t, "%s", "");
	H5_API_RETURN (h5_close_h5hut ());
}



/* debug output */
#define h5_set_verbosity_level FC_GLOBAL(	\
                h5_set_verbosity_level,         \
                H5_SET_VERBOSITY_LEVEL)
h5_int64_t
h5_set_verbosity_level (
	const h5_int64_t *level
	) {

	H5_API_ENTER (h5_int64_t, "level=%lld", (long long)*level);
	H5_API_RETURN(h5_set_loglevel (*level));
}

#define h5_abort_on_error FC_GLOBAL( \
                h5_abort_on_error,     \
                H5_ABORT_ON_ERROR)
h5_int64_t
h5_abort_on_error (
        void
        ) {
	H5_API_ENTER (h5_int64_t, "%s", "");
        h5_set_loglevel (1);
        H5_API_RETURN (h5_set_errorhandler (h5_abort_errorhandler));
}

#define h5_get_error_number FC_GLOBAL(	\
                h5_get_error_number,		\
                H5_GET_ERROR_NUMBER)
h5_int64_t
h5_get_error_number (
        void
        ) {
	H5_API_ENTER (h5_int64_t, "%s", "");
        H5_API_RETURN (h5_get_errno ());
}
