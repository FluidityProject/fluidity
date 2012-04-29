#include "lebiology_pythonc.h"

/**** Initialisation function for the module 'lebiology' ****/
PyMODINIT_FUNC initlebiology(void){
  // Define methods
  (void) Py_InitModule("lebiology", LEBiologyMethods);

  // Create static dicts for variable names, 
  // stage IDs and local namespaces
  pFGLocalsDict = PyDict_New();
  pFGVarNames = PyDict_New();
  pFGEnvNames = PyDict_New();
  pFGStageID = PyDict_New();
  
}

/****************************
 * Meta-model utilities
 ****************************/

/* Add a variable name to the array under pFGVarNames['fg'] */
void lebiology_add_fg_varname_c(char *fg, int fglen, 
                                char *var, int varlen, int *stat)
{
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
}

/* Add an environment name to the array under pFGEnvNames['fg'] */
void lebiology_add_fg_envname_c(char *fg, int fglen, 
                                char *env, int envlen, int *stat)
{
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
}

/* Add a stage ID to pFGStageID['fg'] */
void lebiology_add_fg_stage_id_c(char *fg, int fglen, 
                                char *stage, int stagelen, 
                                double *id, int *stat)
{
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
}

void lebiology_agent_init_c(char *fg, int fglen, 
                            char *key, int keylen, 
                            double vars[], int n_vars, 
                            int *stat)
{
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
}

void lebiology_agent_update_c(char *fg, int fglen, 
                              char *key, int keylen, 
                              double vars[], int n_vars, 
                              double envvals[], int n_envvals, 
                              double *dt, int *stat)
{
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

  //   // Create args and execute kernel function 'def val(agent, env, dt):'
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
}

/*****************************************
 * Callbacks via Python module "lebiology"
 *****************************************/

static PyObject *lebiology_stage_id(PyObject *self, PyObject *args)
{
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
}

static PyObject *lebiology_add_agent(PyObject *self, PyObject *args)
{
  PyObject *pAgentDict;
  int ok;

  if (!PyArg_ParseTuple(args, "O", pAgentDict)) {
    return NULL;
  }

  PyErr_SetString(PyExc_Exception, "lebiology.add_agent not implemented yet!");
  Py_INCREF(PyExc_Exception);
  return NULL;
}
