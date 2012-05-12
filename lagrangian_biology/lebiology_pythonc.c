#include "lebiology_pythonc.h"

/**** Initialisation function for the module 'lebiology' ****/
PyMODINIT_FUNC initlebiology(void)
{
#ifdef HAVE_PYTHON
  // Define methods
  (void) Py_InitModule("lebiology", LEBiologyMethods);

  // Create static dicts for variable names, 
  // stage IDs and local namespaces
  pFGLocalsDict = PyDict_New();
  pFGVarNames = PyDict_New();
  pFGEnvNames = PyDict_New();
  pFGStageID = PyDict_New();
#endif
}

/****************************
 * Meta-model utilities
 ****************************/

/* Add a variable name to the array under pFGVarNames['fg'] */
void lebiology_add_fg_varname_c(char *fg, int fglen, 
                                char *var, int varlen, int *stat)
{
#ifdef HAVE_PYTHON
  //Get pFGVarNames['fg'] or create if it doesn't exist
  char *fg_key = fix_string(fg, fglen);
  PyObject *pVarList = PyDict_GetItemString(pFGVarNames, fg_key);
  if (!pVarList) {
     pVarList = PyList_New(0);
     *stat = PyDict_SetItemString(pFGVarNames, fg_key, pVarList);
     if (*stat<0) {
        PyErr_Print();
        return;
     }
     Py_DECREF(pVarList);
  }

  // Append the variable name
  char *varname = fix_string(var, varlen);
  PyObject *pVarName = PyString_FromString(varname);
  *stat = PyList_Append(pVarList, pVarName);

  Py_DECREF(pVarName);
  free(fg_key);
  free(varname); 
  return;
#endif
}

/* Add an environment name to the array under pFGEnvNames['fg'] */
void lebiology_add_fg_envname_c(char *fg, int fglen, 
                                char *env, int envlen, int *stat)
{
#ifdef HAVE_PYTHON
  //Get pFGVarNames['fg'] or create if it doesn't exist
  char *fg_key = fix_string(fg, fglen);
  PyObject *pEnvList = PyDict_GetItemString(pFGEnvNames, fg_key);
  if (!pEnvList) {
     pEnvList = PyList_New(0);
     *stat = PyDict_SetItemString(pFGEnvNames, fg_key, pEnvList);
     if (*stat<0) {
        PyErr_Print();
        return;
     }
     Py_DECREF(pEnvList);
  }

  // Append the environment name
  char *envname = fix_string(env, envlen);
  PyObject *pEnvName = PyString_FromString(envname);
  *stat = PyList_Append(pEnvList, pEnvName);

  Py_DECREF(pEnvName);
  free(fg_key);
  free(envname);
  return;
#endif
}

/* Add a stage ID to pFGStageID['fg'] */
void lebiology_add_fg_stage_id_c(char *fg, int fglen, 
                                char *stage, int stagelen, 
                                double *id, int *stat)
{
#ifdef HAVE_PYTHON
  //Get pFGStageID['fg'] or create if it doesn't exist
  char *fg_key = fix_string(fg, fglen);
  PyObject *pStageID = PyDict_GetItemString(pFGStageID, fg_key);
  if (!pStageID) {
     pStageID = PyDict_New();
     *stat = PyDict_SetItemString(pFGStageID, fg_key, pStageID);
     if (*stat<0) {
        PyErr_Print();
        return;
     }
     Py_DECREF(pStageID);
  }

  // Insert the stage ID
  char *stagename = fix_string(stage, stagelen);
  PyObject *pID = PyFloat_FromDouble(*id);
  *stat = PyDict_SetItemString(pStageID, stagename, pID);

  Py_DECREF(pID);
  free(fg_key);
  free(stagename); 
  return;
#endif
}

/*****************************
 * Embedded Python interfaces
 *****************************/

/* compile_function runs the specified Python function 'def val(...)'
 * from the schema and stores the local namespace under pFGNamespace['fg'].
 * The compiled function object for can be retrieved from the namespace
 * and applied over agents with different arguments. */
void lebiology_compile_function_c(char *fg, int fglen, 
                                  char *key, int keylen, 
                                  char *func, int funclen, 
                                  int *stat)
{
#ifdef HAVE_PYTHON
  char *c = fix_string(func,funclen);
  int tlen=8+funclen;
  char t[tlen];
  snprintf(t, tlen, "%s\n",c);

  // Runc function with global and local namespace
  PyObject *pMain = PyImport_AddModule("__main__");
  PyObject *pGlobals = PyModule_GetDict(pMain);
  PyObject *pLocals = PyDict_New();
  PyObject *pCode = PyRun_String(t, Py_file_input, pGlobals, pLocals);

  if(!pCode){
    PyErr_Print();
    *stat=-1;
    return;
  }

  // Get pFGLocalsDict['fg'] or create if it doesn't exist
  char *fg_key = fix_string(fg, fglen);
  PyObject *pLocalsDict = PyDict_GetItemString(pFGLocalsDict, fg_key);
  if (!pLocalsDict) {
     pLocalsDict = PyDict_New();
     *stat = PyDict_SetItemString(pFGLocalsDict, fg_key, pLocalsDict);
     if (*stat<0) {
        PyErr_Print();
        return;
     }
     Py_DECREF(pLocalsDict);
  }

  // Add local namespace under 'key'
  char *keystr = fix_string(key,keylen);
  PyDict_SetItemString(pLocalsDict, keystr, pLocals);

  // Check for exceptions
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=-1;
    return;
  }

  Py_DECREF(pCode);
  Py_DECREF(pLocals);
  free(c);
  free(fg_key);
  free(keystr);
  return;
#endif
}

void lebiology_agent_init_c(char *fg, int fglen, 
                            char *key, int keylen, 
                            double vars[], int n_vars, 
                            int *stat)
{
#ifdef HAVE_PYTHON
  // Get variable names and local namespace
  char *fg_key = fix_string(fg, fglen);
  PyObject *pLocalsDict = PyDict_GetItemString(pFGLocalsDict, fg_key);
  PyObject *pVarNames = PyDict_GetItemString(pFGVarNames, fg_key);

  // Get compiled code object from local namespace
  char *keystr = fix_string(key,keylen);
  PyObject *pLocals = PyDict_GetItemString(pLocalsDict, keystr);
  PyObject *pFunc = PyDict_GetItemString(pLocals, "val");
  PyObject *pFuncCode = PyObject_GetAttrString(pFunc, "func_code");

  // Create agent dict
  int i;
  PyObject *pAgent = PyDict_New();
  for(i=0; i<n_vars; i++){
    PyObject *pZeroVal = PyFloat_FromDouble(0.0);
    PyDict_SetItem(pAgent, PyList_GET_ITEM(pVarNames, i), pZeroVal);
    Py_DECREF(pZeroVal);
  }

  // Create args and execute kernel function 'def val(agent):'
  PyObject **pArgs= malloc(sizeof(PyObject*));
  pArgs[0] = pAgent;

  PyObject *pResult = PyEval_EvalCodeEx((PyCodeObject *)pFuncCode, pLocals, NULL, pArgs, 1, NULL, 0, NULL, 0, NULL);

  // Check for correct return
  *stat=0;
  if(!pResult){
    PyErr_Print();
    *stat=-1;
    return;
  }

  // Convert the python result
  for(i=0; i<n_vars; i++){
    vars[i] = PyFloat_AsDouble( PyDict_GetItem(pResult, PyList_GET_ITEM(pVarNames, i)) );
  }

  // Check for exceptions
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=-1;
    return;
  }

  Py_DECREF(pAgent);
  Py_DECREF(pFuncCode);
  Py_DECREF(pResult);
  free(pArgs);
  free(fg_key);
  free(keystr);
#endif
}

void lebiology_agent_update_c(char *fg, int fglen, 
                              char *key, int keylen, 
                              double vars[], int n_vars, 
                              double envvals[], int n_envvals, 
                              double *dt, int *stat)
{
#ifdef HAVE_PYTHON
  // Get variable names and local namespace
  char *fg_key = fix_string(fg, fglen);
  PyObject *pLocalsDict = PyDict_GetItemString(pFGLocalsDict, fg_key);
  PyObject *pVarNames = PyDict_GetItemString(pFGVarNames, fg_key);
  PyObject *pEnvNames = PyDict_GetItemString(pFGEnvNames, fg_key);

  // Get compiled code object from local namespace
  char *keystr = fix_string(key,keylen);
  PyObject *pLocals = PyDict_GetItemString(pLocalsDict, keystr);
  PyObject *pFunc = PyDict_GetItemString(pLocals, "val");
  PyObject *pFuncCode = PyObject_GetAttrString(pFunc, "func_code");

  // Create agent dict
  int i;
  PyObject *pAgent = PyDict_New();
  for(i=0; i<n_vars; i++){
    PyObject *pVarVal = PyFloat_FromDouble(vars[i]);
    PyDict_SetItem(pAgent, PyList_GET_ITEM(pVarNames, i), pVarVal);
    Py_DECREF(pVarVal);
  }

  // Create environment dict
  PyObject *pEnvironment = PyDict_New();
  for(i=0; i<n_envvals; i++){
    PyObject *pEnvVal = PyFloat_FromDouble(envvals[i]);
    PyDict_SetItem(pEnvironment, PyList_GET_ITEM(pEnvNames, i), pEnvVal);
    Py_DECREF(pEnvVal);
  }

  // Create dt argument
  PyObject *pDt = PyFloat_FromDouble(*dt);

  // Create args and execute kernel function 'def val(agent, env, dt):'
  PyObject **pArgs= malloc(sizeof(PyObject*)*3);
  pArgs[0] = pAgent;
  pArgs[1] = pEnvironment;
  pArgs[2] = pDt;

  PyObject *pResult = PyEval_EvalCodeEx((PyCodeObject *)pFuncCode, pLocals, NULL, pArgs, 3, NULL, 0, NULL, 0, NULL);
 
  // Check for Python errors
  *stat=0;
  if(!pResult){
    PyErr_Print();
    *stat=-1;
    return;
  }

  // Convert the python result
  for(i=0; i<n_vars; i++){
    vars[i] = PyFloat_AsDouble( PyDict_GetItem(pResult, PyList_GET_ITEM(pVarNames, i)) );
  }

  // Check for exceptions
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=-1;
    return;
  }

  Py_DECREF(pAgent);
  Py_DECREF(pEnvironment);
  Py_DECREF(pDt);
  Py_DECREF(pFuncCode);
  Py_DECREF(pResult);
  free(pArgs);
  free(fg_key);
  free(keystr);
#endif
}

// Note: Disbaled until dependencies are sorted out
/*
void lebiology_agent_move_c(char *fg, int fglen, 
                            char *key, int keylen, 
                            double pos[], int n_pos, 
                            double vars[], int n_vars, int var_inds[],
                            double *dt, double vector[], int *stat)
{
#ifdef HAVE_PYTHON
  // Get variable names and local namespace
  char *fg_key = fix_string(fg, fglen);
  PyObject *pLocalsDict = PyDict_GetItemString(pFGLocalsDict, fg_key);

  // Get compiled code object from local namespace
  char *keystr = fix_string(key,keylen);
  PyObject *pLocals = PyDict_GetItemString(pLocalsDict, keystr);
  PyObject *pFunc = PyDict_GetItemString(pLocals, "val");
  PyObject *pFuncCode = PyObject_GetAttrString(pFunc, "func_code");

  // Create position argument
  int i;
  PyObject *pPosition = PyList_New(n_pos);
  for(i=0; i<n_pos; i++){
    PyList_SET_ITEM(pPosition, i, PyFloat_FromDouble(pos[i]));
  }

  // If additional variables are given populate the vars dict
  PyObject *pVarNames, *pVarVal;
  PyObject *pVariables = PyDict_New();
  if(n_vars > 0){
    // Grab variable names from persistent dict
    pVarNames = PyDict_GetItemString(pFGVarNames, fg_key);

    for(i=0; i<n_vars; i++){
      pVarVal = PyFloat_FromDouble(vars[i]);
      PyDict_SetItem(pVariables, PyList_GET_ITEM(pVarNames, var_inds[i]-1), pVarVal);
      Py_DECREF(pVarVal);
    }
  }

  // Create dt argument
  PyObject *pDt = PyFloat_FromDouble(*dt);

  // Create args and execute kernel function 'def val(agent, env, dt):'
  PyObject **pArgs= malloc(sizeof(PyObject*)*3);
  pArgs[0] = pPosition;
  pArgs[1] = pVariables;
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

  if(n_vars > 0){
    for(i=0; i<n_vars; i++){
      vars[i] = PyFloat_AsDouble( PyDict_GetItem(pVariables, PyList_GET_ITEM(pVarNames, var_inds[i]-1)) );
    }
  }

  // Convert the python result
  PyObject *result_ref;
  for(i=0; i<n_pos; i++){
    result_ref = PySequence_GetItem(pResult, i);
    vector[i] = PyFloat_AsDouble( result_ref );
    Py_DECREF(result_ref);
  }

  // Check for exceptions
  if (PyErr_Occurred()){
    PyErr_Print();
    *stat=-1;
    return;
  }

  Py_DECREF(pPosition);
  Py_DECREF(pVariables);
  Py_DECREF(pDt);
  Py_DECREF(pFuncCode);
  Py_DECREF(pResult);
  free(pArgs);
  free(fg_key);
  free(keystr);
#endif
}
*/

/*****************************************
 * Callbacks via Python module "lebiology"
 *****************************************/

static PyObject *lebiology_stage_id(PyObject *self, PyObject *args)
{
#ifdef HAVE_PYTHON
  const char *fg, *stage;
  if (!PyArg_ParseTuple(args, "ss", &fg, &stage)) {
     return NULL;
  }

  // Check if given FGroup exists
  PyObject *pStageIDDict = PyDict_GetItemString(pFGStageID, fg);
  if (!pStageIDDict) {
     PyErr_SetString(PyExc_KeyError, "FGroup not found in lebiology.stage_id");
     Py_INCREF(PyExc_KeyError);
     return NULL;
  }

  // Check if given stage exists
  PyObject *pStageID = PyDict_GetItemString(pStageIDDict, stage);
  if (!pStageID) {
     PyErr_SetString(PyExc_KeyError, "Stage not found in lebiology.stage_id");
     Py_INCREF(PyExc_KeyError);
     return NULL;
  }

  Py_INCREF(pStageID);
  return pStageID;
#endif
}

static PyObject *lebiology_add_agent(PyObject *self, PyObject *args)
{
#ifdef HAVE_PYTHON
  PyObject *pAgentDict, *pPositionList, *pAgentVar;
  const char *fg;
  if (!PyArg_ParseTuple(args, "sOO", &fg, &pAgentDict, &pPositionList)) {
    return NULL;
  }

  // Check for correctly parsed agent argument
  if (!pAgentDict) {
     PyErr_SetString(PyExc_TypeError, "Error parsing agent dictionary in lebiology.add_agent");
     Py_INCREF(PyExc_TypeError);
     return NULL;
  }

  // Check for correctly parsed agent argument
  if (!pPositionList) {
     PyErr_SetString(PyExc_TypeError, "Error parsing position list in lebiology.add_agent");
     Py_INCREF(PyExc_TypeError);
     return NULL;
  }

  // Check if given FGroup exists
  PyObject *pVarNames = PyDict_GetItemString(pFGVarNames, fg);
  if (!pVarNames) {
     PyErr_SetString(PyExc_KeyError, "FGroup not found in lebiology.add_agent");
     Py_INCREF(PyExc_KeyError);
     return NULL;
  }

  // Convert the agent dict
  int *n_vars = malloc(sizeof(int));
  *n_vars = PyList_Size(pVarNames);
  double *vars = malloc(*n_vars * sizeof(double));
  int i;
  for (i=0; i<*n_vars; i++) {
    pAgentVar = PyDict_GetItem(pAgentDict, PyList_GET_ITEM(pVarNames, i));
    if (pAgentVar) {
      vars[i] = PyFloat_AsDouble(pAgentVar);
    } else {
      vars[i] = 0.0;
    }
  }

  //Convert the position list
  int *n_pos = malloc(sizeof(int));
  *n_pos = PyList_Size(pPositionList);
  double *pos = malloc(*n_pos * sizeof(double));
  for (i=0; i<*n_pos; i++) {
    pos[i] = PyFloat_AsDouble( PyList_GET_ITEM(pPositionList, i) );
  }

  // Pass the new agent on to the Fortran
  fl_add_agent_c(vars, n_vars, pos, n_pos);

  free(vars);
  free(n_vars);
  free(pos);
  free(n_pos);

  Py_INCREF(Py_None);
  return Py_None;
#endif
}
