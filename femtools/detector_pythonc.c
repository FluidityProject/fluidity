#include "detector_pythonc.h"

// Functions to execute Python over detectors 

void python_run_string_store_locals_c(char *str, int strlen, 
                                     char *dict, int dictlen, 
                                     char *key, int keylen, int *stat){
#ifdef HAVE_PYTHON
  /* Run a python command from Fortran and store local namespace
   * in a global dictionary for later evaluation
   */

  char *c = fix_string(str,strlen);
  int tlen=8+strlen;
  char t[tlen];
  snprintf(t, tlen, "%s\n",c);

  // Get a reference to the main module and global dictionary
  PyObject *pMain = PyImport_AddModule("__main__");
  PyObject *pGlobals = PyModule_GetDict(pMain);

  // Global and local namespace dictionaries for our code.
  PyObject *pLocals = PyDict_New();
  
  // Execute the user's code.
  PyObject *pCode = PyRun_String(t, Py_file_input, pGlobals, pLocals);

  // Check for Python errors
  if(!pCode){
    PyErr_Print();
    *stat=-1;
    return;
  }

  // Create dictionary to hold the locals if it doesn't exist
  char *local_dict = fix_string(dict, dictlen);
  PyObject *pLocalDict= PyDict_GetItemString(pGlobals, local_dict);
  if(!pLocalDict){
    int cmdlen = dictlen + 8;
    char cmd[cmdlen];
    snprintf(cmd, cmdlen, "%s=dict()", local_dict);
    PyRun_SimpleString(cmd);
    pLocalDict= PyDict_GetItemString(pGlobals, local_dict);
  }

  // Now store a copy of our local namespace
  char *local_key = fix_string(key, keylen);
  *stat = PyDict_SetItemString(pLocalDict, local_key, pLocals);

  // Check for Python errors
  if(!pCode){
    PyErr_Print();
    *stat=-1;
    return;
  }

  Py_DECREF(pCode);
  Py_DECREF(pLocals);
  free(c); 
  free(local_dict);
  free(local_key);

#endif
}


void python_run_random_walk_from_locals_c(int ele, int dim, 
                                           double lcoords[], double *dt,
                                           char *dict, int dictlen, 
                                           char *key, int keylen,
                                           double value[], int *stat){
#ifdef HAVE_PYTHON
  /* Evaluate the detector val() function from a previously stored local namespace
   * found in dict under key. The interface is: val(ele, local_coords)
   */

  // Get a reference to the main module and global dictionary
  PyObject *pMain = PyImport_AddModule("__main__");
  PyObject *pGlobals = PyModule_GetDict(pMain);

  // Get the local context from the global dictionary
  char *local_dict = fix_string(dict, dictlen);
  PyObject *pLocalDict= PyDict_GetItemString(pGlobals, local_dict);
  char *local_key = fix_string(key, keylen);
  PyObject *pLocals = PyDict_GetItemString(pLocalDict, local_key);

  // Extract val function from the local dict and its code object
  PyObject *pFunc = PyDict_GetItemString(pLocals, "val");
  PyObject *pFuncCode = PyObject_GetAttrString(pFunc, "func_code");

  // Create ele argument
  PyObject *pEle = PyInt_FromLong( (long)ele );

  // Create local_coords argument
  int i;
  PyObject *pLCoords = PyTuple_New(dim+1);
  for(i=0; i<dim+1; i++){
    PyTuple_SET_ITEM(pLCoords, i, PyFloat_FromDouble(lcoords[i]));
  }

  // Create dt argument
  PyObject *pDt = PyFloat_FromDouble(*dt);

  // Create argument array
  PyObject **pArgs= malloc(sizeof(PyObject*)*3);
  pArgs[0] = pEle;
  pArgs[1] = pLCoords;
  pArgs[2] = pDt;

  // Run val(ele, local_coords)
  PyObject *pResult = PyEval_EvalCodeEx((PyCodeObject *)pFuncCode, pLocals, NULL, pArgs, 3, NULL, 0, NULL, 0, NULL);
 
  // Check for Python errors
  *stat=0;
  if(!pResult){
    PyErr_Print();
    *stat=-1;
    return;
  }

  // Convert the python result
  PyObject *result_ref;
  for(i=0; i<dim; i++){
    // GetItem returns a new reference that needs a DECREF
    result_ref = PySequence_GetItem(pResult, i);
    value[i] = PyFloat_AsDouble( result_ref );
    Py_DECREF(result_ref);
  }

  Py_DECREF(pEle);
  Py_DECREF(pDt);
  Py_DECREF(pLCoords);
  Py_DECREF(pFuncCode);
  Py_DECREF(pResult);
  free(pArgs);
  free(local_dict);
  free(local_key);

#endif
}
