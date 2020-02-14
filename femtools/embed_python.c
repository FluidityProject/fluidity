/*  Copyright (C) 2006 Imperial College London and others.

Please see the AUTHORS file in the main source directory for a full list
of copyright holders.

Prof. C Pain
Applied Modelling and Computation Group
Department of Earth Science and Engineering
Imperial College London

amcgsoftware@imperial.ac.uk

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation,
version 2.1 of the License.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
USA
*/



#include "confdefs.h"
#include "string.h"

#ifdef HAVE_PYTHON
#include "Python.h"
#if PY_MAJOR_VERSION >= 3
#define PyInt_FromLong PyLong_FromLong
#define PyInt_AsLong PyLong_AsLong
#define PyString_Size PyUnicode_GET_SIZE
#define PyString_AsString PyUnicode_AsUTF8
#endif
#endif
#ifdef HAVE_NUMPY
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"
#endif

int eval_user_func(char *function, int function_len, PyObject *pLocals, PyObject **pFunc) {
  PyObject *pMain, *pGlobals, *pCode;
  char *function_c;

  // the function string passed down from Fortran needs terminating,
  // so make a copy and fiddle with it (remember to free it)
  function_c = (char *)malloc(function_len+1);
  memcpy(function_c, function, function_len);
  function_c[function_len] = 0;

  // Get a reference to the main module and global dictionary
  pMain = PyImport_AddModule("__main__");
  pGlobals = PyModule_GetDict(pMain);

  // Execute the user's code and clean up allocated string.
  pCode = PyRun_String(function_c, Py_file_input, pGlobals, pLocals);
  free(function_c);
  Py_DECREF(pCode);

  // Extract the function from the code.
  *pFunc = PyDict_GetItemString(pLocals, "val");
  if (*pFunc == NULL) {
      printf("Couldn't find a 'val' function in your Python code.\n");
      return 1;
  }

  // Check for errors in executing user code.
  if (PyErr_Occurred()) {
    PyErr_Print();
    return 1;
  }

  return 0;
}

static inline void set_pos_tuple(int i, int dim, PyObject *pPos, double *x, double *y, double *z) {
  PyObject *px;
  px = PyFloat_FromDouble(x[i]);
  PyTuple_SetItem(pPos, 0, px);
  if (dim > 1) {
    px = PyFloat_FromDouble(y[i]);
    PyTuple_SetItem(pPos, 1, px);
  }
  if (dim > 2) {
    px = PyFloat_FromDouble(z[i]);
    PyTuple_SetItem(pPos, 2, px);
  }
}

void set_field_from_python(char *function, int function_len, int dim, int nodes, int *natt,
			   double *x, double *y, double *z, double t, double *dt,
			   int *stat, int *result_dim, void **result,
			   int (*set_result)(int, int *, void *, void **, PyObject *),
			   void *result_user_data)
{
#ifndef HAVE_PYTHON
  int i;
  strncpy(function, "No Python support!\n", (size_t)function_len);
  for (i = 0; i < function_len; i++)
  {
    if (function[i] == '\0')
      function[i] = ' ';
  }

  *stat=1;
  return;
#else

  PyObject *pFunc, *pLocals, *pResult, *pArgs, *pKwArgs, *pPos, *pT;
  int i, res;

  // Global and local namespace dictionaries for our code.
  pLocals = PyDict_New();

  // load the user's function as a Python object -- borrows locals
  if (eval_user_func(function, function_len, pLocals, &pFunc)) {
    *stat = 1;
    return;
  }

  // Create Python objects for function arguments
  pT = PyFloat_FromDouble(t);
  pPos = PyTuple_New(dim);
  pArgs = PyTuple_New(2);
  // note that PyTuple_SetItem steals refs, so we don't have to decrement
  // the refcounts of pT or pPos
  PyTuple_SetItem(pArgs, 0, pPos);
  PyTuple_SetItem(pArgs, 1, pT);
  pKwArgs = PyDict_New();

  // particle routines need dt too
  if (dt != NULL) {
    PyObject *pdT = PyFloat_FromDouble(*dt);
    PyDict_SetItemString(pKwArgs, "dt", pdT);
    Py_DECREF(pdT);
  }
  // array routines need the array length
  if (natt != NULL) {
    PyObject *pnAtt = PyLong_FromLong(*natt);
    PyDict_SetItemString(pKwArgs, "n", pnAtt);
    Py_DECREF(pnAtt);
  }

  // Check for a Python error in the function call
  if (PyErr_Occurred()) {
    PyErr_Print();
    *stat = 1;
    return;
  }

  // populate position tuple
  for (i = 0; i < nodes; i++) {
    set_pos_tuple(i, dim, pPos, x, y, z);
    pResult = PyObject_Call(pFunc, pArgs, pKwArgs);

    // Check for a Python error in the function call
    if (PyErr_Occurred()) {
      PyErr_Print();
      *stat = 1;
      return;
    }

    res = set_result(i, result_dim, result_user_data, result, pResult);
    if (res) {
      *stat = res;
      return;
    }
    // Check for a Python error in result.
    if (PyErr_Occurred()) {
      PyErr_Print();
      *stat = 1;
      return;
    }

    Py_DECREF(pResult);
  }

  // clean up
  Py_DECREF(pArgs);
  Py_DECREF(pKwArgs);
  Py_DECREF(pLocals);
  PyGC_Collect();
  *stat = 0;
  #endif
}

int set_scalar_result_double(int i, int *dim, void *data, void **_result, PyObject *pResult)
{
  double *result = *(double**)_result;
  result[i] = PyFloat_AsDouble(pResult);
  return 0;
}

int set_scalar_result_double_array(int i, int *dim, void *data, void **_result, PyObject *pResult)
{
  int natt = *((int *)data);
  // cast from pointer to 1d array, to 2d array (from fortran -- result(natt, npart))
  double (*result)[natt] = (double (*)[natt])(*(double**)_result);

  if (PyObject_Length(pResult) != natt) {
    fprintf(stderr, "Error: length of object returned from python (%d) does not match"
	    " the expected number of attributes (%d)\n",
	    (int)PyObject_Length(pResult), natt);
    return 1;
  }

  // copy to result array for this particle
  for (int j = 0; j < natt; j++) {
    PyObject *px = PySequence_GetItem(pResult, j);
    result[i][j] = PyFloat_AsDouble(px);
    Py_DECREF(px);
  }

  return 0;
}

#define set_scalar_field_from_python F77_FUNC(set_scalar_field_from_python, SET_SCALAR_FIELD_FROM_PYTHON)
void set_scalar_field_from_python(char *function, int function_len, int dim,
                                  int nodes, double *x, double *y, double *z, double t,
				  double *result, int *stat)
{
  set_field_from_python(function, function_len, dim, nodes, NULL, x, y, z, t, NULL, stat,
			NULL, (void**)&result, set_scalar_result_double, NULL);
}

void set_scalar_particles_from_python(char *function, int function_len, int dim, int ndete,
				      double *x, double *y, double *z, double t, double dt,
				      double *result, int *stat)
{
  set_field_from_python(function, function_len, dim, ndete, NULL, x, y, z, t, &dt, stat,
			NULL, (void**)&result, set_scalar_result_double, NULL);
}

void set_scalar_particles_from_python_array(char *function, int function_len, int dim, int npart,
					    int natt, double *x, double *y, double *z, double t, double dt,
					    double *result, int *stat)
{
  set_field_from_python(function, function_len, dim, npart, &natt, x, y, z, t, &dt, stat,
			NULL, (void**)&result, set_scalar_result_double_array, &natt);
}

int set_scalar_result_integer(int i, int *dim, void *data, void **_result, PyObject *pResult)
{
  int *result = *(int**)_result;
  result[i] = PyLong_AsLong(pResult);
  return 0;
}

#define set_integer_array_from_python F77_FUNC(set_integer_array_from_python, SET_INTEGER_ARRAY_FROM_PYTHON)
void set_integer_array_from_python(char* function, int function_len, int dim,
                                   int nodes, double *x, double *y, double *z, double t,
                                   int* result, int* stat)
{
  set_field_from_python(function, function_len, dim, nodes, NULL, x, y, z, t, NULL, stat,
			NULL, (void**)&result, set_scalar_result_integer, NULL);
}

// populate an array of output arrays per dimension, e.g. double *result[] = {result_x, ...}
int set_vector_result_double(int i, int *dim, void *data, void **result, PyObject *pResult)
{
  if (PyObject_Length(pResult) != *dim) {
      fprintf(stderr, "Error: length of object returned from python (%d) does not match the allocated dimension of the vector field (%d).\n",
              (int)PyObject_Length(pResult), *dim);
      return 1;
  }

  for (int d = 0; d < *dim; d++) {
    PyObject *px = PySequence_GetItem(pResult, d);
    ((double**)result)[d][i] = PyFloat_AsDouble(px);
    Py_DECREF(px);
  }

  return 0;
}

// similar to set_vector_result_double, but expects a contiguous result array -- in fortran: result(dim, ndete)
int set_vector_result_double_contiguous(int i, int *dim, void *data, void **_result, PyObject *pResult)
{
  if (PyObject_Length(pResult) != *dim) {
      fprintf(stderr, "Error: length of object returned from python (%d) does not match the allocated dimension of the vector field (%d).\n",
              (int)PyObject_Length(pResult), *dim);
      return 1;
  }

  // cast from pointer to 1d array, to 2d array for easier indexing
  double (*result)[*dim] = (double (*)[*dim])(*(double**)_result);

  PyObject *px;
  for (int d = 0; d < *dim; d++) {
    px = PySequence_GetItem(pResult, d);
    result[i][d] = PyFloat_AsDouble(px);
    Py_DECREF(px);
  }

  return 0;
}

// expects a contiguous result array for an attribute array -- result(dim*natt, ndete)
// pResult should be a numpy array
int set_vector_result_double_array(int i, int *dim, void *data, void **_result, PyObject *pResult)
{
  PyArrayObject *pArray = (PyArrayObject *)PyArray_ContiguousFromObject(pResult, NPY_DOUBLE, 2, 2);
  int natt = *((int *)data);
  double (*result)[*dim*natt] = (double (*)[*dim*natt])(*(double**)_result);

  if (PyArray_DIMS(pArray)[0] != natt || PyArray_DIMS(pArray)[1] != *dim) {
    fprintf(stderr, "Error: function did not return (%dx%d) array, got (%ldx%ld)\n",
      natt, *dim, PyArray_DIMS(pArray)[0], PyArray_DIMS(pArray)[1]);
    return 1;
  }

  for (int j = 0; j < natt; j++) {
    for (int d = 0; d < *dim; d++) {
      result[i][j*(*dim) + d] = *((double*)PyArray_GETPTR2(pArray, j, d));
    }
  }

  return 0;
}

#define set_vector_field_from_python F77_FUNC(set_vector_field_from_python, SET_VECTOR_FIELD_FROM_PYTHON)
void set_vector_field_from_python(char *function, int function_len, int dim,
                                  int nodes, double *x, double *y, double *z, double t,
                                  int result_dim, double *result_x, double *result_y, double *result_z,
                                  int *stat)
{
  double *results[] = {result_x, result_y, result_z};

  set_field_from_python(function, function_len, dim, nodes, NULL, x, y, z, t, NULL, stat,
			&result_dim, (void**)results, set_vector_result_double, NULL);
}

//#define set_vector_particles_from_python F77_FUNC(set_vector_particles_from_python, SET_VECTOR_PARTICLES_FROM_PYTHON)
void set_vector_particles_from_python(char *function, int function_len, int dim,
				      int ndete, double *x, double *y, double *z, double t, double dt,
				      double *result, int *stat)
{
  set_field_from_python(function, function_len, dim, ndete, NULL, x, y, z, t, &dt, stat,
			&dim, (void**)&result, set_vector_result_double_contiguous, NULL);
}

void set_vector_particles_from_python_array(char *function, int function_len, int dim,
				      int ndete, int natt, double *x, double *y, double *z, double t, double dt,
				      double *result, int *stat)
{
  set_field_from_python(function, function_len, dim, ndete, &natt, x, y, z, t, &dt, stat,
			&dim, (void**)&result, set_vector_result_double_array, &natt);
}

int set_tensor_result_double(int i, int *dim, void *data, void **result, PyObject *pResult)
{
  // get a 2d array from result
  PyArrayObject *pArray = (PyArrayObject *)PyArray_ContiguousFromObject(pResult, NPY_DOUBLE, 2, 2);
  if (PyErr_Occurred()) {
    PyErr_Print();
    return 1;
  }
  if (PyArray_DIMS(pArray)[0] != dim[0] || PyArray_DIMS(pArray)[1] != dim[1]) {
    fprintf(stderr, "Error: dimensions of array returned from python ([%d, %d]) do not match allocated dimensions of the tensor_field ([%d, %d])).\n",
	    (int) PyArray_DIMS(pArray)[0], (int) PyArray_DIMS(pArray)[1], dim[0], dim[1]);
    return 1;
  }

  for (int ii = 0; ii < dim[0]; ii++) {
    for (int jj = 0; jj < dim[1]; jj++) {
      (*(double**)result)[i*(dim[0]*dim[1]) + jj*dim[0] + ii] = *((double*)PyArray_GETPTR2(pArray, ii, jj));
    }
  }
  Py_DECREF(pArray);

  return 0;
}

int set_tensor_result_double_array(int i, int *dim, void *data, void **result, PyObject *pResult)
{
  PyArrayObject *pArray = (PyArrayObject *)PyArray_ContiguousFromObject(pResult, NPY_DOUBLE, 3, 3);
  if (PyErr_Occurred()) {
    PyErr_Print();
    return 1;
  }
  int natt = *((int *)data);

  if (PyArray_DIMS(pArray)[0] != natt || PyArray_DIMS(pArray)[1] != dim[0] || PyArray_DIMS(pArray)[2] != dim[1]) {
    fprintf(stderr, "Error: function did not return (%dx%dx%d) array, got (%ldx%ldx%ld)\n",
      natt, dim[0], dim[1], PyArray_DIMS(pArray)[0], PyArray_DIMS(pArray)[1], PyArray_DIMS(pArray)[2]);
    return 1;
  }

  for (int j = 0; j < natt; j++) {
    for (int ii = 0; ii < dim[0]; ii++) {
      for (int jj = 0; jj < dim[1]; jj++) {
        (*(double**)result)[i*(natt*dim[0]*dim[1]) + j * (dim[0]*dim[1]) + jj*dim[0] + ii] =
          *((double*)PyArray_GETPTR3(pArray, j, ii, jj));
      }
    }
  }

  Py_DECREF(pArray);

  return 0;
}

void set_tensor_particles_from_python(char *function, int function_len, int dim,
                                      int npart, double *x, double *y, double *z,
                                      double t, double dt, double *result, int *stat)
{
  int result_dims[] = {dim, dim};
  set_field_from_python(function, function_len, dim, npart, NULL, x, y, z, t, &dt, stat,
      result_dims, (void**)&result, set_tensor_result_double, NULL);
}

void set_tensor_particles_from_python_array(char *function, int function_len, int dim,
                                      int npart, int natt, double *x, double *y, double *z,
                                      double t, double dt, double *result, int *stat)
{
  int result_dims[] = {dim, dim};
  set_field_from_python(function, function_len, dim, npart, &natt, x, y, z, t, &dt, stat,
      result_dims, (void**)&result, set_tensor_result_double_array, &natt);
}

#define set_tensor_field_from_python F77_FUNC(set_tensor_field_from_python, SET_TENSOR_FIELD_FROM_PYTHON)
void set_tensor_field_from_python(char *function, int function_len, int dim,
                                  int nodes, double *x, double *y, double *z, double t,
				  int *result_dim, double *result, int *stat)
{
#ifndef HAVE_NUMPY
  int i;
  strncpy(function, "No Numpy support!\n", (size_t)function_len);
  for (i=0; i < function_len; i++)
  {
    if (function[i] == '\0')
      function[i] = ' ';
  }
  *stat = 1;
  return;
#else
  import_array();

  set_field_from_python(function, function_len, dim, nodes, NULL, x, y, z, t, NULL, stat,
			result_dim, (void**)&result, set_tensor_result_double, NULL);

#endif
}

static inline void set_dict_from_fields(int i, int dim, int *nfields, int name_len, char *_field_names, double *_field_vals,
			  PyObject *pNames)
{
  int all_fields = nfields[0] + dim*nfields[1] + dim*dim*nfields[2];
  double (*field_vals)[all_fields] = (double (*)[all_fields])_field_vals;
  char (*field_names)[name_len] = (char (*)[name_len])_field_names;

  // for creation of numpy arrays of the right size
  npy_intp dims[] = {dim, dim};

  // set values for fields dictionary
  for (int j = 0; j < nfields[0]; j++) {
    // scalar fields
    PyObject *pField = PyFloat_FromDouble(field_vals[i][j]);
    PyDict_SetItemString(pNames, field_names[j], pField);
    Py_DECREF(pField);
  }

  for (int j = 0; j < nfields[1]; j++) {
    // vector fields
    PyObject *pVec = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    int offset = nfields[0] + j*dim;

    for (int k = 0; k < dim; k++) {
      ((double*)PyArray_DATA((PyArrayObject*)pVec))[k] = field_vals[i][k + offset];
    }

    PyDict_SetItemString(pNames, field_names[j + nfields[0]], pVec);
    Py_DECREF(pVec);
  }

  for (int j = 0; j < nfields[2]; j++) {
    // tensor fields
    PyObject *pTens = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    int offset = nfields[0] + nfields[1]*dim + j*dim*dim;

    for (int k = 0; k < dim; k++) {
      for (int l = 0; l < dim; l++) {
	((double*)PyArray_DATA((PyArrayObject*)pTens))[k + l*dim] = field_vals[i][k*dim + l + offset];
      }
    }

    PyDict_SetItemString(pNames, field_names[j + nfields[0] + nfields[1]], pTens);
    Py_DECREF(pTens);
  }
}

static void set_dict_from_old_atts(int i, int dim, int *natts, int name_len, char *_att_names,
				   int *att_dims, double *_att_vals, PyObject *pNames)
{
  int all_atts = natts[0] + dim*natts[1] + dim*dim*natts[2];
  double (*att_vals)[all_atts] = (double (*)[all_atts])_att_vals;
  char (*att_names)[name_len] = (char (*)[name_len])_att_names;

  npy_intp dims[] = {0, dim, dim};

  int ai = 0; // attribute index

  for (int j = 0; j < natts[0]; j++) {
    // scalar old attributes
    PyObject *pAtt;
    int ad = att_dims[j];

    if (ad == 0) {
      // non-array
      pAtt = PyFloat_FromDouble(att_vals[i][ai++]);
    } else {
      dims[0] = ad;
      pAtt = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
      for (int d = 0; d < ad; d++) {
	*((double*)PyArray_GETPTR1((PyArrayObject*)pAtt, d)) = att_vals[i][ai + d];
      }
      ai += ad;
    }
    PyDict_SetItemString(pNames, att_names[j], pAtt);
    Py_DECREF(pAtt);
  }
  for (int j = 0; j < natts[1]; j++) {
    // vector old attributes
    PyObject *pAtt;
    int ad = att_dims[j + natts[0]];

    if (ad == 0) {
      // non-array
      pAtt = PyArray_SimpleNew(1, dims+1, NPY_DOUBLE);
      for (int k = 0; k < dim; k++) {
	*((double*)PyArray_GETPTR1((PyArrayObject*)pAtt, k)) = att_vals[i][ai + k];
      }
      ai += dim;
    } else {
      dims[0] = ad;
      pAtt = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
      for (int d = 0; d < ad; d++) {
	for (int k = 0; k < dim; k++) {
	  *((double*)PyArray_GETPTR2((PyArrayObject*)pAtt, d, k)) = att_vals[i][ai + k];
	}
	ai += dim;
      }
    }
    PyDict_SetItemString(pNames, att_names[j + natts[0]], pAtt);
    Py_DECREF(pAtt);
  }
  for (int j = 0; j < natts[2]; j++) {
    // tensor old attributes
    PyObject *pAtt;
    int ad = att_dims[j + natts[0] + natts[1]];

    if (ad == 0) {
      // non-array
      pAtt = PyArray_SimpleNew(2, dims+1, NPY_DOUBLE);
      for (int k = 0; k < dim; k++) {
	for (int l = 0; l < dim; l++) {
	  *((double*)PyArray_GETPTR2((PyArrayObject*)pAtt, l, k)) = att_vals[i][ai + l];
	}
	ai += dim;
      }
    } else {
      dims[0] = ad;
      pAtt = PyArray_SimpleNew(3, dims, NPY_DOUBLE);
      for (int d = 0; d < ad; d++) {
	for (int k = 0; k < dim; k++) {
	  for (int l = 0; l < dim; l++) {
	    *((double*)PyArray_GETPTR3((PyArrayObject*)pAtt, d, l, k)) = att_vals[i][ai + l];
	  }
	  ai += dim;
	}
      }
    }
    PyDict_SetItemString(pNames, att_names[j + natts[0] + natts[1]], pAtt);
    Py_DECREF(pAtt);
  }
}

void set_field_from_python_fields(char *function, int function_len, int dim, int nodes, int *natt,
				  double *x, double *y, double *z, double t, double dt,
				  int name_len, int *nfields, char *field_names, double *field_vals,
				  int *old_nfields, char *old_field_names, double *old_field_vals,
				  int *old_natts, char *old_att_names, int *old_att_dims, double *old_atts,
				  int *stat, int *result_dim, void **result,
				  int (*set_result)(int, int *, void *, void **, PyObject*),
				  void *result_user_data)
{
#ifndef HAVE_PYTHON
  int i;
  strncpy(function, "No Python support!\n", (size_t)function_len);
  for (i = 0; i < function_len; i++)
  {
    if (function[i] == '\0')
      function[i] = ' ';
  }

  *stat=1;
  return;
#else
  PyObject *pLocals, *pFunc, *pNames, *pT, *pdT, *pPos, *pArgs, *pKwArgs, *pResult;

  import_array();

  // load the user's function as a Python object -- borrows locals
  pLocals = PyDict_New();
  if (eval_user_func(function, function_len, pLocals, &pFunc)) {
    *stat = 1;
    return;
  }

  // create objects to hold arguments
  pNames = PyDict_New();
  pT = PyFloat_FromDouble(t);
  pdT = PyFloat_FromDouble(dt);
  pPos = PyTuple_New(dim);

  // build function argument tuple
  pArgs = PyTuple_New(4);
  PyTuple_SetItem(pArgs, 0, pPos);
  PyTuple_SetItem(pArgs, 1, pT);
  PyTuple_SetItem(pArgs, 2, pdT);
  PyTuple_SetItem(pArgs, 3, pNames);

  pKwArgs = PyDict_New();
  if (natt != NULL) {
    PyObject *pnAtt = PyLong_FromLong(*natt);
    PyDict_SetItemString(pKwArgs, "n", pnAtt);
    Py_DECREF(pnAtt);
  }

  for (int i = 0; i < nodes; i++) {
    // set position
    set_pos_tuple(i, dim, pPos, x, y, z);
    set_dict_from_fields(i, dim, nfields, name_len, field_names, field_vals, pNames);
    set_dict_from_fields(i, dim, old_nfields, name_len, old_field_names, old_field_vals, pNames);
    set_dict_from_old_atts(i, dim, old_natts, name_len, old_att_names, old_att_dims, old_atts, pNames);

    if (PyErr_Occurred()) {
      PyErr_Print();
      *stat = 1;
      return;
    }

    pResult = PyObject_Call(pFunc, pArgs, pKwArgs);
    if (PyErr_Occurred()) {
      PyErr_Print();
      *stat = 1;
      return;
    }

    int res = set_result(i, result_dim, result_user_data, result, pResult);
    if (res) {
      *stat = res;
      return;
    }

    Py_DECREF(pResult);
  }

  Py_DECREF(pArgs);
  Py_DECREF(pKwArgs);
  Py_DECREF(pLocals);
  PyGC_Collect();
  *stat = 0;
#endif
}

void set_scalar_particles_from_python_fields(char *function, int function_len, int dim, int nodes,
					     double *x, double *y, double *z, double t, double dt,
					     int name_len, int *nfields, char *field_names, double *field_vals,
					     int *old_nfields, char *old_field_names, double *old_field_vals,
					     int *old_natts, char *old_att_names, int *old_att_dims, double *old_atts,
					     double *result, int *stat)
{
  set_field_from_python_fields(function, function_len, dim, nodes, NULL,
			       x, y, z, t, dt, name_len,
			       nfields, field_names, field_vals,
			       old_nfields, old_field_names, old_field_vals,
			       old_natts, old_att_names, old_att_dims, old_atts,
			       stat, NULL, (void**)&result, set_scalar_result_double, NULL);
}

void set_scalar_particles_from_python_fields_array(char *function, int function_len, int dim, int nodes, int natt,
						   double *x, double *y, double *z, double t, double dt,
						   int name_len, int *nfields, char *field_names, double *field_vals,
						   int *old_nfields, char *old_field_names, double *old_field_vals,
						   int *old_natts, char *old_att_names, int *old_att_dims, double *old_atts,
						   double *result, int *stat)
{
  set_field_from_python_fields(function, function_len, dim, nodes, &natt,
			       x, y, z, t, dt, name_len,
			       nfields, field_names, field_vals,
			       old_nfields, old_field_names, old_field_vals,
			       old_natts, old_att_names, old_att_dims, old_atts,
			       stat, NULL, (void**)&result, set_scalar_result_double_array, &natt);
}

void set_vector_particles_from_python_fields(char *function, int function_len, int dim, int nodes,
					     double *x, double *y, double *z, double t, double dt,
					     int name_len, int *nfields, char *field_names, double *field_vals,
					     int *old_nfields, char *old_field_names, double *old_field_vals,
					     int *old_natts, char *old_att_names, int *old_att_dims, double *old_atts,
					     double *result, int *stat)
{
  set_field_from_python_fields(function, function_len, dim, nodes, NULL,
			       x, y, z, t, dt, name_len,
			       nfields, field_names, field_vals,
			       old_nfields, old_field_names, old_field_vals,
			       old_natts, old_att_names, old_att_dims, old_atts,
			       stat, &dim, (void**)&result,
			       set_vector_result_double_contiguous, NULL);
}

void set_vector_particles_from_python_fields_array(char *function, int function_len, int dim, int nodes, int natt,
						   double *x, double *y, double *z, double t, double dt,
						   int name_len, int *nfields, char *field_names, double *field_vals,
						   int *old_nfields, char *old_field_names, double *old_field_vals,
						   int *old_natts, char *old_att_names, int *old_att_dims, double *old_atts,
						   double *result, int *stat)
{
  set_field_from_python_fields(function, function_len, dim, nodes, &natt,
			       x, y, z, t, dt, name_len,
			       nfields, field_names, field_vals,
			       old_nfields, old_field_names, old_field_vals,
			       old_natts, old_att_names, old_att_dims, old_atts,
			       stat, &dim, (void**)&result,
			       set_vector_result_double_array, &natt);
}

void set_tensor_particles_from_python_fields(char *function, int function_len, int dim, int nodes,
					     double *x, double *y, double *z, double t, double dt,
					     int name_len, int *nfields, char *field_names, double *field_vals,
					     int *old_nfields, char *old_field_names, double *old_field_vals,
					     int *old_natts, char *old_att_names, int *old_att_dims, double *old_atts,
					     double *result, int *stat)
{
  int result_dims[] = {dim, dim};
  set_field_from_python_fields(function, function_len, dim, nodes, NULL,
			       x, y, z, t, dt, name_len,
			       nfields, field_names, field_vals,
			       old_nfields, old_field_names, old_field_vals,
			       old_natts, old_att_names, old_att_dims, old_atts,
			       stat, result_dims, (void**)&result,
			       set_tensor_result_double, NULL);
}

void set_tensor_particles_from_python_fields_array(char *function, int function_len, int dim, int nodes, int natt,
						   double *x, double *y, double *z, double t, double dt,
						   int name_len, int *nfields, char *field_names, double *field_vals,
						   int *old_nfields, char *old_field_names, double *old_field_vals,
						   int *old_natts, char *old_att_names, int *old_att_dims, double *old_atts,
						   double *result, int *stat)
{
  int result_dims[] = {dim, dim};
  set_field_from_python_fields(function, function_len, dim, nodes, &natt,
			       x, y, z, t, dt, name_len,
			       nfields, field_names, field_vals,
			       old_nfields, old_field_names, old_field_vals,
			       old_natts, old_att_names, old_att_dims, old_atts,
			       stat, result_dims, (void**)&result,
			       set_tensor_result_double_array, &natt);
}


#define set_detectors_from_python F77_FUNC(set_detectors_from_python, SET_DETECTORS_FROM_PYTHON)
void set_detectors_from_python(char *function, int *function_len, int *dim,
                               int *ndete, double *t,
                               int *result_dim,
                               double result_x[], double result_y[],
                               double result_z[], int* stat)
{
#ifndef HAVE_PYTHON
  int i;
  strncpy(function, "No Python support!\n", (size_t) *function_len);
  for (i=0; i < *function_len; i++)
  {
    if (function[i] == '\0')
      function[i] = ' ';
  }
  *stat=1;
  return;
#else
  PyObject *pMain, *pGlobals, *pLocals, *pFunc, *pCode, *pResult, *pResultItem,
    *pArgs, *px, *pT;
  char *function_c;
  int i;

  // the function string passed down from Fortran needs terminating,
  // so make a copy and fiddle with it (remember to free it)
  function_c = (char *)malloc(*function_len+3);
  memcpy( function_c, function, *function_len );
  function_c[*function_len] = 0;

  // Get a reference to the main module and global dictionary
  pMain = PyImport_AddModule("__main__");
  pGlobals = PyModule_GetDict(pMain);
  // Global and local namespace dictionaries for our code.
  pLocals=PyDict_New();

  // Execute the user's code.
  pCode=PyRun_String(function_c, Py_file_input, pGlobals, pLocals);

  // Extract the function from the code.
  pFunc=PyDict_GetItemString(pLocals, "val");

  // Clean up memory from null termination.
  free(function_c);

  // Check for errors in executing user code.
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  // Python form of time variable.
  pT=PyFloat_FromDouble(*t);

  // Tuple of arguments to function;
  pArgs=PyTuple_New(1);
  PyTuple_SetItem(pArgs, 0, pT);

  // Check for a Python error in the function call
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  pResult=PyObject_CallObject(pFunc, pArgs);

  // Check for a Python error in the function call
   if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }


  for (i = 0; i < *ndete; i++){
    pResultItem = PySequence_GetItem(pResult, i);

    px=PySequence_GetItem(pResultItem, 0);

    result_x[i]=PyFloat_AsDouble(px);
    // Check for a Python error in unpacking tuple.
    if (PyErr_Occurred()){
      PyErr_Print();
      *stat=1;
      return;
    }
    Py_DECREF(px);

    if (*result_dim>1) {
      px=PySequence_GetItem(pResultItem, 1);
      result_y[i]=PyFloat_AsDouble(px);
      // Check for a Python error in unpacking tuple.
      if (PyErr_Occurred()){
         PyErr_Print();
         return;
      }

      Py_DECREF(px);
      if (*result_dim>2) {
        px=PySequence_GetItem(pResultItem, 2);
        result_z[i]=PyFloat_AsDouble(px);
      // Check for a Python error in unpacking tuple.
       if (PyErr_Occurred()){
          PyErr_Print();
          return;
       }
        Py_DECREF(px);
      }
    }

    Py_DECREF(pResultItem);
  }

  Py_DECREF(pResult);

  // Clean up
  Py_DECREF(pArgs);
  Py_DECREF(pLocals);
  Py_DECREF(pCode);

  // Force a garbage collection
  PyGC_Collect();

  *stat=0;
  return;
#endif
}

#define real_from_python F77_FUNC(real_from_python, REAL_FROM_PYTHON)
void real_from_python(char* function, int* function_len,
                        double* t,
                        double* result, int* stat)
{
#ifndef HAVE_PYTHON
  int i;
  strncpy(function, "No Python support!\n", (size_t) *function_len);
  for (i=0; i < *function_len; i++)
  {
    if (function[i] == '\0')
      function[i] = ' ';
  }
  *stat=1;
  return;
#else
  PyObject *pMain, *pGlobals, *pLocals, *pFunc, *pCode, *pResult,
    *pArgs, *pT;

  char *function_c;

  // the function string passed down from Fortran needs terminating,
  // so make a copy and fiddle with it (remember to free it)
  function_c = (char *)malloc(*function_len+3);
  memcpy( function_c, function, *function_len );
  function_c[*function_len] = 0;

  // Get a reference to the main module and global dictionary
  pMain = PyImport_AddModule("__main__");

  pGlobals = PyModule_GetDict(pMain);
  // Global and local namespace dictionaries for our code.
  pLocals=PyDict_New();

  // Execute the user's code.
  pCode=PyRun_String(function_c, Py_file_input, pGlobals, pLocals);

  // Extract the function from the code.
  pFunc=PyDict_GetItemString(pLocals, "val");

  // Clean up memory from null termination.
  free(function_c);

  // Check for errors in executing user code.
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  // Python form of time variable.
  pT=PyFloat_FromDouble(*t);

  // Tuple of arguments to function;
  pArgs=PyTuple_New(1);
  PyTuple_SetItem(pArgs, 0, pT);

  // Check for a Python error in the function call
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  pResult=PyObject_CallObject(pFunc, pArgs);

  // Check for a Python error in the function call
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  *result = PyFloat_AsDouble(pResult);

  // Check for a Python error in result.
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  Py_DECREF(pResult);

  // Clean up
  Py_DECREF(pArgs);
  Py_DECREF(pLocals);
  Py_DECREF(pCode);

  // Force a garbage collection
  PyGC_Collect();


  *stat=0;
  return;
#endif
}

void free_c_vector(void** vector)
{
  free(*vector);
}

void real_vector_from_python(char* function, int* function_len,
                             double* t,
                             void** result,
                             int* result_len,
                             int* stat)
{
 int i;
#ifndef HAVE_PYTHON
  strncpy(function, "No Python support!\n", (size_t) *function_len);
  for (i=0; i < *function_len; i++)
  {
    if (function[i] == '\0')
      function[i] = ' ';
  }
  *stat=1;
  return;
#else
  PyObject *pMain, *pGlobals, *pLocals, *pFunc, *pCode, *pResult,
    *pArgs, *pT, *pResultItem;

  char *function_c;

  // the function string passed down from Fortran needs terminating,
  // so make a copy and fiddle with it (remember to free it)
  function_c = (char *)malloc(*function_len+3);
  memcpy( function_c, function, *function_len );
  function_c[*function_len] = 0;

  // Get a reference to the main module and global dictionary
  pMain = PyImport_AddModule("__main__");

  pGlobals = PyModule_GetDict(pMain);
  // Global and local namespace dictionaries for our code.
  pLocals=PyDict_New();

  // Execute the user's code.
  pCode=PyRun_String(function_c, Py_file_input, pGlobals, pLocals);

  // Extract the function from the code.
  pFunc=PyDict_GetItemString(pLocals, "val");

  // Clean up memory from null termination.
  free(function_c);

  // Check for errors in executing user code.
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  // Python form of time variable.
  pT=PyFloat_FromDouble(*t);

  // Tuple of arguments to function;
  pArgs=PyTuple_New(1);
  PyTuple_SetItem(pArgs, 0, pT);

  // Check for a Python error in the function call
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  pResult=PyObject_CallObject(pFunc, pArgs);

  // Check for a Python error in the function call
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  *result_len = PySequence_Length(pResult);

  // Check for a Python error in result_dim.
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  *result = malloc(*result_len * sizeof(double));

  // Unpack tuple to pointer
  for (i = 0; i < *result_len; i++){
    pResultItem = PySequence_GetItem(pResult, i);
    // Check for a Python error in unpacking tuple.
    if (PyErr_Occurred()){
      PyErr_Print();
      *stat=1;
      return;
    }

    ((double*)*result)[i]=PyFloat_AsDouble(pResultItem);

    // Check we really got a float.
    if (PyErr_Occurred()){
      PyErr_Print();
      *stat=1;
      return;
    }

    Py_DECREF(pResultItem);
  }


  Py_DECREF(pResult);

  // Clean up
  Py_DECREF(pArgs);
  Py_DECREF(pLocals);
  Py_DECREF(pCode);

  // Force a garbage collection
  PyGC_Collect();


  *stat=0;
  return;
#endif
}

void integer_vector_from_python(char* function, int* function_len,
                             double* t,
                             void** result,
                             int* result_len,
                             int* stat)
{
 int i;
#ifndef HAVE_PYTHON
  strncpy(function, "No Python support!\n", (size_t) *function_len);
  for (i=0; i < *function_len; i++)
  {
    if (function[i] == '\0')
      function[i] = ' ';
  }
  *stat=1;
  return;
#else
  PyObject *pMain, *pGlobals, *pLocals, *pFunc, *pCode, *pResult,
    *pArgs, *pT, *pResultItem;

  char *function_c;

  // the function string passed down from Fortran needs terminating,
  // so make a copy and fiddle with it (remember to free it)
  function_c = (char *)malloc(*function_len+3);
  memcpy( function_c, function, *function_len );
  function_c[*function_len] = 0;

  // Get a reference to the main module and global dictionary
  pMain = PyImport_AddModule("__main__");

  pGlobals = PyModule_GetDict(pMain);
  // Global and local namespace dictionaries for our code.
  pLocals=PyDict_New();

  // Execute the user's code.
  pCode=PyRun_String(function_c, Py_file_input, pGlobals, pLocals);

  // Extract the function from the code.
  pFunc=PyDict_GetItemString(pLocals, "val");

  // Clean up memory from null termination.
  free(function_c);

  // Check for errors in executing user code.
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  // Python form of time variable.
  pT=PyFloat_FromDouble(*t);

  // Tuple of arguments to function;
  pArgs=PyTuple_New(1);
  PyTuple_SetItem(pArgs, 0, pT);

  // Check for a Python error in the function call
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  pResult=PyObject_CallObject(pFunc, pArgs);

  // Check for a Python error in the function call
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  *result_len = PySequence_Length(pResult);

  // Check for a Python error in result_dim.
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  *result = malloc(*result_len * sizeof(long));

  // Unpack tuple to pointer
  for (i = 0; i < *result_len; i++){
    pResultItem = PySequence_GetItem(pResult, i);
    // Check for a Python error in unpacking tuple.
    if (PyErr_Occurred()){
      PyErr_Print();
      *stat=1;
      return;
    }

    ((long*)*result)[i]=PyInt_AsLong(pResultItem);

    // Check we really got a float.
    if (PyErr_Occurred()){
      PyErr_Print();
      *stat=1;
      return;
    }

    Py_DECREF(pResultItem);
  }


  Py_DECREF(pResult);

  // Clean up
  Py_DECREF(pArgs);
  Py_DECREF(pLocals);
  Py_DECREF(pCode);

  // Force a garbage collection
  PyGC_Collect();


  *stat=0;
  return;
#endif
}

#define integer_from_python F77_FUNC(integer_from_python, INTEGER_FROM_PYTHON)
void integer_from_python(char* function, int* function_len,
                        double* t,
                        int* result, int* stat)
{
#ifndef HAVE_PYTHON
  int i;
  strncpy(function, "No Python support!\n", (size_t) *function_len);
  for (i=0; i < *function_len; i++)
  {
    if (function[i] == '\0')
      function[i] = ' ';
  }
  *stat=1;
  return;
#else
  PyObject *pMain, *pGlobals, *pLocals, *pFunc, *pCode, *pResult,
    *pArgs, *pT;

  char *function_c;

  // the function string passed down from Fortran needs terminating,
  // so make a copy and fiddle with it (remember to free it)
  function_c = (char *)malloc(*function_len+3);
  memcpy( function_c, function, *function_len );
  function_c[*function_len] = 0;

  // Get a reference to the main module and global dictionary
  pMain = PyImport_AddModule("__main__");

  pGlobals = PyModule_GetDict(pMain);
  // Global and local namespace dictionaries for our code.
  pLocals=PyDict_New();

  // Execute the user's code.
  pCode=PyRun_String(function_c, Py_file_input, pGlobals, pLocals);

  // Extract the function from the code.
  pFunc=PyDict_GetItemString(pLocals, "val");

  // Clean up memory from null termination.
  free(function_c);

  // Check for errors in executing user code.
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  // Python form of time variable.
  pT=PyFloat_FromDouble(*t);

  // Tuple of arguments to function;
  pArgs=PyTuple_New(1);
  PyTuple_SetItem(pArgs, 0, pT);

  // Check for a Python error in the function call
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  pResult=PyObject_CallObject(pFunc, pArgs);

  // Check for a Python error in the function call
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  *result = PyLong_AsLong(pResult);

  // Check for a Python error in result.
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  Py_DECREF(pResult);

  // Clean up
  Py_DECREF(pArgs);
  Py_DECREF(pLocals);
  Py_DECREF(pCode);

  // Force a garbage collection
  PyGC_Collect();


  *stat=0;
  return;
#endif
}

#define string_from_python F77_FUNC(string_from_python, STRING_FROM_PYTHON)
void string_from_python(char* function, int* function_len,
                        int* result_len,
                        double* t,
                        char* result, int* stat)
{
#ifndef HAVE_PYTHON
  int i;
  strncpy(function, "No Python support!\n", (size_t) *function_len);
  for (i=0; i < *function_len; i++)
  {
    if (function[i] == '\0')
      function[i] = ' ';
  }
  *stat=1;
  return;
#else
  PyObject *pMain, *pGlobals, *pLocals, *pFunc, *pCode, *pResult,
    *pArgs, *pT;
  int pResult_len;

  char *function_c;

  // the function string passed down from Fortran needs terminating,
  // so make a copy and fiddle with it (remember to free it)
  function_c = (char *)malloc(*function_len+3);
  memcpy( function_c, function, *function_len );
  function_c[*function_len] = 0;

  // Get a reference to the main module and global dictionary
  pMain = PyImport_AddModule("__main__");

  pGlobals = PyModule_GetDict(pMain);
  // Global and local namespace dictionaries for our code.
  pLocals=PyDict_New();

  // Execute the user's code.
  pCode=PyRun_String(function_c, Py_file_input, pGlobals, pLocals);

  // Extract the function from the code.
  pFunc=PyDict_GetItemString(pLocals, "val");

  // Clean up memory from null termination.
  free(function_c);

  // Check for errors in executing user code.
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  // Python form of time variable.
  pT=PyFloat_FromDouble(*t);

  // Tuple of arguments to function;
  pArgs=PyTuple_New(1);
  PyTuple_SetItem(pArgs, 0, pT);

  // Check for a Python error in the function call
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  pResult=PyObject_CallObject(pFunc, pArgs);

  // Check for a Python error in the function call
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  pResult_len = PyString_Size(pResult);
  if(pResult_len > *result_len){
    fprintf(stderr, "In string_from_python\n");
    fprintf(stderr, "Warning: Truncating returned string\n");
    fflush(stderr);
    memcpy(result, PyString_AsString(pResult), *result_len * sizeof(char));
  }else{
    memcpy(result, PyString_AsString(pResult), pResult_len * sizeof(char));
    *result_len = pResult_len;
  }

  // Check for a Python error in result.
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  Py_DECREF(pResult);

  // Clean up
  Py_DECREF(pArgs);
  Py_DECREF(pLocals);
  Py_DECREF(pCode);

  // Force a garbage collection
  PyGC_Collect();


  *stat=0;
  return;
#endif
}
