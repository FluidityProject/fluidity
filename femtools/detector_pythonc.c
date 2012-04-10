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


void python_run_random_walk_from_locals_c(double position[], int ele, int dim, 
                                           double lcoords[], double *dt,
                                           double vars[], int vars_ind[], int varslen,
                                           char *dict, int dictlen, 
                                           char *key, int keylen,
                                           double value[], int *stat){
#ifdef HAVE_PYTHON
  /* Evaluate the detector random walk function using a previously stored local namespace
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

  // Create local_coords argument
  int i;
  PyObject *pPosition = PyList_New(dim);
  for(i=0; i<dim; i++){
    PyList_SET_ITEM(pPosition, i, PyFloat_FromDouble(position[i]));
  }

  // Create ele argument
  PyObject *pEle = PyInt_FromLong( (long)ele );

  // Create local_coords argument
  PyObject *pLCoords = PyList_New(dim+1);
  for(i=0; i<dim+1; i++){
    PyList_SET_ITEM(pLCoords, i, PyFloat_FromDouble(lcoords[i]));
  }

  // Create dt argument
  PyObject *pDt = PyFloat_FromDouble(*dt);

  // If additional variables are given populate the vars dict
  PyObject *pPersistent, *pFGVarNames, *pVarNames;
  PyObject *pVars = PyDict_New();
  if(varslen > 0){
    // Grab variable names from persistent dict
    pPersistent = PyDict_GetItemString(pGlobals, "persistent");
    pFGVarNames = PyDict_GetItemString(pPersistent, "fg_var_names");
    pVarNames = PyDict_GetItemString(pFGVarNames, local_dict);

    for(i=0; i<varslen; i++){
      PyObject *pVarVal = PyFloat_FromDouble(vars[i]);
      PyDict_SetItem(pVars, PyList_GET_ITEM(pVarNames, vars_ind[i]-1), pVarVal);
      Py_DECREF(pVarVal);
    }
  }

  // Create argument array
  PyObject **pArgs= malloc(sizeof(PyObject*)*4);
  pArgs[0] = pPosition;
  pArgs[1] = pEle;
  pArgs[2] = pLCoords;
  pArgs[3] = pVars;
  pArgs[4] = pDt;

  // Run val(ele, local_coords)
  PyObject *pResult = PyEval_EvalCodeEx((PyCodeObject *)pFuncCode, pLocals, NULL, pArgs, 5, NULL, 0, NULL, 0, NULL);
 
  // Check for Python errors
  *stat=0;
  if(!pResult){
    PyErr_Print();
    *stat=-1;
    return;
  }

  if(varslen > 0){
    for(i=0; i<varslen; i++){
      vars[i] = PyFloat_AsDouble( PyDict_GetItem(pVars, PyList_GET_ITEM(pVarNames, vars_ind[i]-1)) );
    }
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
  Py_DECREF(pPosition);
  Py_DECREF(pVars);
  Py_DECREF(pFuncCode);
  Py_DECREF(pResult);
  free(pArgs);
  free(local_dict);
  free(local_key);

#endif
}

void python_run_agent_biology_init_c(char *function, int function_len,
                                 char *var_list, int var_list_len, 
                                 double biovars[], int n_biovars, 
                                 double *stage_id, int *stat){
#ifdef HAVE_PYTHON
  /* Initialise agent biology variables via a dictionary of name -> value pairs.
   * This will create a dict of all variables stored by the agent, set the stage explicitly,
   * and let the user fill in the blanks in the val(biovars) function.
   */

  // Get a reference to the main module and global dictionary
  PyObject *pMain = PyImport_AddModule("__main__");
  PyObject *pGlobals = PyModule_GetDict(pMain);
  PyObject *pLocals = PyDict_New();

  char *local_key = fix_string(var_list, var_list_len);
  char *py_function = fix_string(function, function_len);

  // Execute the user's code, and get val(...) code
  PyObject *pCode=PyRun_String(py_function, Py_file_input, pGlobals, pLocals);
  PyObject *pFunc=PyDict_GetItemString(pLocals, "val");
  PyObject *pFuncCode = PyObject_GetAttrString(pFunc, "func_code");

  // Grab variable names from persistent dict
  PyObject *pPersistent= PyDict_GetItemString(pGlobals, "persistent");
  PyObject *pFGVarNames= PyDict_GetItemString(pPersistent, "fg_var_names");
  PyObject *pVarNames= PyDict_GetItemString(pFGVarNames, local_key);

  // Create var dictionary; the first entry stage_id is provided, 
  // everything else defaults to 0.0
  int i;
  PyObject *pBiology = PyDict_New();
  PyObject *pStageVal = PyFloat_FromDouble(*stage_id);
  PyDict_SetItem(pBiology, PyList_GET_ITEM(pVarNames, 0), pStageVal);
  Py_DECREF(pStageVal);
  for(i=1; i<n_biovars; i++){
    PyObject *pZeroVal = PyFloat_FromDouble(0.0);
    PyDict_SetItem(pBiology, PyList_GET_ITEM(pVarNames, i), pZeroVal);
    Py_DECREF(pZeroVal);
  }

  // Create argument array
  PyObject **pArgs= malloc(sizeof(PyObject*));
  pArgs[0] = pBiology;

  // Run val(variables)
  PyObject *pResult = PyEval_EvalCodeEx((PyCodeObject *)pFuncCode, pLocals, NULL, pArgs, 1, NULL, 0, NULL, 0, NULL);

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

  Py_DECREF(pLocals);
  Py_DECREF(pBiology);
  Py_DECREF(pFuncCode);
  Py_DECREF(pResult);
  free(pArgs);
  free(local_key);
  free(py_function);
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

