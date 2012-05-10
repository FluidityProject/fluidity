#include "detector_pythonc.h"

// Functions to execute Python over detectors 

void python_run_string_store_locals_c(char *str, int strlen, 
                                     char *dict, int dictlen, 
                                     char *key, int keylen, int *stat){
#ifdef HAVE_PYTHON
  /* Run a python string from Fortran and store local namespace
   * in a global dictionary under a given key for later evaluation
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

void python_run_agent_biology_c(double *dt, char *dict, int dictlen, char *key, int keylen,
                                double biovars[], int n_biovars, double env_vars[],
                                int n_env_vars, int *stat){
#ifdef HAVE_PYTHON
  /* Evaluate the detector biology update function using a previously stored local namespace
   * found in dict under key. The interface is: val(biovars, envvars, dt), 
   * where biovars is a dict of agent variables, envvars a dict of environment field names to 
   * locally evaluated field values, and dt is the timestep.
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

  // Create biology dict argument
  int i;
  PyObject *pPersistent= PyDict_GetItemString(pGlobals, "persistent");
  PyObject *pFGVarNames= PyDict_GetItemString(pPersistent, "fg_var_names");
  PyObject *pVarNames= PyDict_GetItemString(pFGVarNames, local_dict);
  PyObject *pBiology = PyDict_New();
  for(i=0; i<n_biovars; i++){
    PyObject *pVarVal = PyFloat_FromDouble(biovars[i]);
    PyDict_SetItem(pBiology, PyList_GET_ITEM(pVarNames, i), pVarVal);
    Py_DECREF(pVarVal);
  }

  // Create environment dict argument
  PyObject *pFGEnvNames= PyDict_GetItemString(pPersistent, "fg_env_names");
  PyObject *pEnvNames= PyDict_GetItemString(pFGEnvNames, local_dict);
  PyObject *pEnvironment = PyDict_New();
  for(i=0; i<n_env_vars; i++){
    PyObject *pEnvVal = PyFloat_FromDouble(env_vars[i]);
    PyDict_SetItem(pEnvironment, PyList_GET_ITEM(pEnvNames, i), pEnvVal);
    Py_DECREF(pEnvVal);
  }

  // Create dt argument
  PyObject *pDt = PyFloat_FromDouble(*dt);

  // Create argument array
  PyObject **pArgs= malloc(sizeof(PyObject*)*3);
  pArgs[0] = pBiology;
  pArgs[1] = pEnvironment;
  pArgs[2] = pDt;

  // Run val(biovars, envvars, dt)
  PyObject *pResult = PyEval_EvalCodeEx((PyCodeObject *)pFuncCode, pLocals, NULL, pArgs, 3, NULL, 0, NULL, 0, NULL);
 
  // Check for Python errors
  *stat=0;
  if(!pResult){
    PyErr_Print();
    *stat=-1;
    return;
  }

  // Convert the python result
  for(i=0; i<n_biovars; i++){
    biovars[i] = PyFloat_AsDouble( PyDict_GetItem(pResult, PyList_GET_ITEM(pVarNames, i)) );
  }

  Py_DECREF(pDt);
  Py_DECREF(pBiology);
  Py_DECREF(pEnvironment);
  Py_DECREF(pFuncCode);
  Py_DECREF(pResult);
  free(pArgs);
  free(local_dict);
  free(local_key);

#endif
}

void python_get_element_integer_c(int dim, double coords_centre[], 
                                  char *dict, int dictlen, 
                                  char *key, int keylen,
                                  double *result, int *stat){
#ifdef HAVE_PYTHON
  /* 
   * 
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

  // Create coords argument
  int i;
  PyObject *pCoords = PyTuple_New(dim);
  for(i=0; i<dim; i++){
    PyTuple_SET_ITEM(pCoords, i, PyFloat_FromDouble(coords_centre[i]));
  }

  // Create argument array
  PyObject **pArgs= malloc(sizeof(PyObject*));
  pArgs[0] = pCoords;

  // Run val(ele, local_coords)
  PyObject *pResult = PyEval_EvalCodeEx((PyCodeObject *)pFuncCode, pLocals, NULL, pArgs, 1, NULL, 0, NULL, 0, NULL);

  // Check for Python errors
  *stat=0;
  if(!pResult){
    PyErr_Print();
    *stat=-1;
    return;
  }

  // Convert the python result
  *result = (int) PyFloat_AsDouble(pResult);

  Py_DECREF(pCoords);
  Py_DECREF(pFuncCode);
  Py_DECREF(pResult);
  free(pArgs);
  free(local_dict);
  free(local_key);

#endif
}

