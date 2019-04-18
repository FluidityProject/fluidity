%module H5hutpy
%{
#define SWIG_FILE_WITH_INIT
#if defined(PARALLEL_IO)
#include <mpi.h>
#endif
#include <stdint.h>
#include "h5core/h5_types.h"
//#include "H5.h"
#include "H5hut.h"
%}

%import <stdint.i>

%include numpy.i
%include cstring.i
%include cpointer.i

%apply unsigned long int { h5_prop_t };
%apply unsigned long int { h5_file_t };

///////////////////////////////////////////////////////////////////////////////
%typemap(out) h5_err_t H5HasFileAttrib {
    if($1)
        $result = (PyObject *)Py_True;
    else
        $result = (PyObject *)Py_False;
}


%cstring_bounded_output(char* const attrib_name, 256);
extern h5_err_t H5GetFileAttribInfo (
	const h5_file_t,
	const h5_size_t,
	char* const attrib_name,
	const h5_size_t l_attrib_name=256,
	h5_int64_t *OUTPUT,
	h5_size_t *OUTPUT);

extern h5_err_t H5GetFileAttribInfoByName (
	const h5_file_t,
	const char* const name,
	h5_int64_t *OUTPUT, h5_size_t *OUTPUT);

extern h5_err_t H5GetFileAttribName (
	const h5_file_t,
	const h5_id_t,
	char* const attrib_name,
	const h5_size_t l_attrib_name=256);

extern h5_err_t H5ReadFileAttribString (
	const h5_file_t f,
	const char* const name,
	char* const buffer);

extern h5_err_t H5GetStepAttribInfo (
	const h5_file_t,
	const h5_size_t,
	char* const attrib_name,
	const h5_size_t l_attrib_name=256,
	h5_int64_t *OUTPUT,
	h5_size_t *OUTPUT);

extern h5_err_t H5GetStepAttribInfoByName (
	const h5_file_t,
	const char* const name,
	h5_int64_t *OUTPUT, h5_size_t *OUTPUT);

extern h5_err_t H5GetStepAttribName (
	const h5_file_t,
	const h5_id_t,
	char* const attrib_name,
	const h5_size_t l_attrib_name=256);

%clear char* const attrib_name;

///////////////////////////////////////////////////////////////////////////////
%cstring_bounded_output(char* const dataset_name, 256);
extern h5_err_t H5PartGetDatasetInfo (
	const h5_file_t,
	const h5_id_t,
	char* const dataset_name, const h5_size_t l_dataset_name=256,
	h5_int64_t *OUTPUT, h5_size_t *OUTPUT);

extern h5_err_t H5PartGetDatasetInfoByName (
	const h5_file_t,
	const char* const name,
	h5_int64_t *OUTPUT, h5_size_t *OUTPUT);

extern h5_err_t H5PartGetDatasetName (
	const h5_file_t,
	const h5_id_t,
	char* const dataset_name,
	const h5_size_t l_dataset_name=256);
%clear char* const dataset_name;

%typemap(out) h5_err_t H5PartHasDataset {
    if($1)
        $result = (PyObject *)Py_True;
    else
        $result = (PyObject *)Py_False;
}

%typemap(out) h5_err_t H5PartHasView {
    if($1)
        $result = (PyObject *)Py_True;
    else
        $result = (PyObject *)Py_False;
}

%apply (unsigned long long* IN_ARRAY1) { h5_size_t* };
%apply (unsigned int* IN_ARRAY1) { h5_uint32_t* }; //uint32_t
%apply (unsigned long long* IN_ARRAY1) { h5_uint64_t* }; //uint64_t
%apply (int* IN_ARRAY1) { h5_int32_t* }; //int32_t
%apply (long long* IN_ARRAY1) { h5_int64_t* }; //int64_t
%apply (float* IN_ARRAY1) { h5_float32_t* };
%apply (double* IN_ARRAY1) { h5_float64_t* };


%init %{
import_array();
%}

#ifdef (PARALLEL_IO)
%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);
%typemap(in) MPI_Comm* {
    MPI_Comm *ptr = (MPI_Comm *)0;
    int res = SWIG_AsPtr_MPI_Comm($input, &ptr);
    if (!SWIG_IsOK(res) || !ptr) {
      SWIG_exception_fail(SWIG_ArgError((ptr ? res : SWIG_TypeError)),
			  "in method '" "$symname" "', argument " "$argnum"" of type '" "MPI_Comm""'");
    }
    $1 = ptr;
    if (SWIG_IsNewObj(res)) free((char*)ptr);
}
#endif

%ignore h5_report_errorhandler;
%ignore h5_abort_errorhandler;
%ignore h5priv_vprintf;
%ignore h5_verror;
%ignore H5ReportErrorhandler;
%ignore H5AbortErrorhandler;

%include "h5core/h5_types.h"

%include "H5_attachments.h"

%rename(H5OpenFile) H5OpenFile2 (const char* const, const h5_int64_t, const h5_prop_t);
%include "H5_file.h"

%include "H5_model.h"
%include "H5_file_attribs.h"
%include "H5_step_attribs.h"
%include "H5_log.h"

%include "H5Block_attribs.h"
%include "H5Block_io.h"
%include "H5Block_model.h"

%include "H5Part_io.h"
%include "H5Part_model.h"

%clear  h5_size_t*;
%clear  h5_int64_t*;
