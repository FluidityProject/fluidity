/*  Copyright (C) 2006 Imperial College London and others.
⋅⋅⋅⋅
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
#endif
#ifdef HAVE_NUMPY
#include "numpy/arrayobject.h"
#endif
#ifdef HAVE_ADJOINT
#include "libadjoint/libadjoint.h"
#else
#define adj_adjointer int
#endif

void adj_variables_from_python(adj_adjointer* adjointer, char* function, int function_len,
                               double start_time, double end_time, long timestep,
                               void** result,
                               int* result_len,
                               int* stat)
{
#if !defined(HAVE_PYTHON) || !defined(HAVE_ADJOINT)
  int i;
  strncpy(function, "Need both python and libadjoint support.\n", (size_t) function_len);
  for (i = 0; i < function_len; i++)
  {
    if (function[i] == '\0')
      function[i] = ' ';
  }
  *stat=1;
  return;
#else
  PyObject *pMain, *pGlobals, *pLocals, *pFunc, *pCode, *pResult, *pResult2,
    *pArgs, *pStartT, *pEndT, *pTimes, *pTimestep, *pResultItem, *pName;

  char *function_c;

  int i;
  int ierr;
  adj_variable* result_var;
  int auxiliary;

  // the function string passed down from Fortran needs terminating,
  // so make a copy and fiddle with it (remember to free it)
  function_c = (char *)malloc(function_len+3);
  memcpy( function_c, function, function_len );
  function_c[function_len] = 0;

  // Get a reference to the main module and global dictionary
  pMain = PyImport_AddModule("__main__");

  pGlobals = PyModule_GetDict(pMain);
  // Global and local namespace dictionaries for our code.
  pLocals=PyDict_New();

  // Execute the user's code.
  pCode=PyRun_String(function_c, Py_file_input, pGlobals, pLocals);

  // Extract the function from the code.
  pFunc=PyDict_GetItemString(pLocals, "dependencies");

  // Clean up memory from null termination.
  free(function_c);

  // Check for errors in executing user code.
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  // Python form of time variable.
  pStartT=PyFloat_FromDouble(start_time);
  pEndT=PyFloat_FromDouble(end_time);
  pTimes = PyTuple_New(2);
  PyTuple_SetItem(pTimes, 0, pStartT);
  PyTuple_SetItem(pTimes, 1, pEndT);

  // Timestep.
  pTimestep = PyLong_FromLong(timestep);

  // Tuple of arguments to function;
  pArgs=PyTuple_New(2);
  PyTuple_SetItem(pArgs, 0, pTimes);
  PyTuple_SetItem(pArgs, 1, pTimestep);

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

  // OK. Now we need to call fluidity.parse_functional.make_adj_variable on the output.
  PyRun_String("from fluidity.parse_functional import make_adj_variables", Py_file_input, pGlobals, pLocals);
  pFunc=PyDict_GetItemString(pLocals, "make_adj_variables");
  // Check for a Python error.
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  // We need to pack pResult into a tuple so that we can call make_adj_variable on it
  Py_DECREF(pArgs);
  pArgs=PyTuple_New(1);
  PyTuple_SetItem(pArgs, 0, pResult);
  pResult2=PyObject_CallObject(pFunc, pArgs);

  *result_len = PySequence_Length(pResult2);

  // Check for a Python error.
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=1;
    return;
  }

  *result = malloc(*result_len * sizeof(adj_variable));
  result_var = (adj_variable*) *result;

  // Unpack tuple to pointer
  for (i = 0; i < *result_len; i++)
  {
    pResultItem = PySequence_GetItem(pResult2, i);
    // Check for a Python error in unpacking tuple.
    if (PyErr_Occurred()){
      PyErr_Print();
      *stat=1;
      return;
    }

    // pResultItem is a python dictionary, containing two keys: name and timestep.
    if (!PyDict_Check(pResultItem))
    {
      *stat = 1;
      return;
    }

    pName = PyDict_GetItemString(pResultItem, "name");
    pTimestep = PyDict_GetItemString(pResultItem, "timestep");
    // Check we really got our stuff out.
    if (PyErr_Occurred()){
      PyErr_Print();
      *stat=1;
      return;
    }

    /* If we don't have a MaterialPhaseName:: in front, assume it's auxiliary */
    if (strstr(PyString_AsString(pName), "::") == NULL)
      auxiliary = ADJ_TRUE;
    else
      auxiliary = ADJ_FALSE;

    /* If "Coordinate" is in the variable name, assume it's auxiliary */
    if (strstr(PyString_AsString(pName), "Coordinate") != NULL)
      auxiliary = ADJ_TRUE;

    {
      int iteration_count;
      ierr = adj_create_variable(PyString_AsString(pName), (int) PyLong_AsLong(pTimestep), 0, auxiliary, &(result_var[i]));
      adj_chkierr(ierr);
      ierr = adj_iteration_count(adjointer, result_var[i], &iteration_count);
      if (ierr != ADJ_ERR_INVALID_INPUTS)
      {
        ierr = adj_create_variable(PyString_AsString(pName), (int) PyLong_AsLong(pTimestep), iteration_count-1, auxiliary, &(result_var[i]));
        adj_chkierr(ierr);
      }
    }

    Py_DECREF(pResultItem);
  }

  Py_DECREF(pResult2);

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

