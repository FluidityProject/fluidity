#include "lebiology_pythonc.h"

/**** Initialisation function for the module 'lebiology' ****/
PyMODINIT_FUNC initlebiology(void){
  // Define methods
  (void) Py_InitModule("lebiology", LEBiologyMethods);

  // Create static dicts for variable names and stage IDs
  pFGVarNames = PyDict_New();
  pFGStageID = PyDict_New();
  
}

/****************************
 * Meta-model utilities
 ****************************/

/* Add a variable name to the array under pFGVarNames['fg'] */
void lebiology_add_fg_varname_c(char *fg, int fglen, 
                                char *var, int varlen, int *stat){
  char *fg_key = fix_string(fg, fglen);
  PyObject *pVarList = PyDict_GetItemString(pFGVarNames, fg_key);

  // Create list for 'fg' if it doesn't exist
  if (!pVarList) {
     pVarList = PyList_New(0);
     *stat = PyDict_SetItemString(pFGVarNames, fg_key, pVarList);
     if (*stat<0) {
        PyErr_Print();
        return;
     }
  }

  // Append the variable name
  char *varname = fix_string(var, varlen);
  PyObject *pVarName = PyString_FromString(varname);
  *stat = PyList_Append(pVarList, pVarName);
  return;
}

/* Add a stage ID to pFGStageID['fg'] */
void lebiology_add_fg_stage_id_c(char *fg, int fglen, 
                                char *stage, int stagelen, 
                                double *id, int *stat){
  char *fg_key = fix_string(fg, fglen);
  PyObject *pStageID = PyDict_GetItemString(pFGStageID, fg_key);

  // Create dict for 'fg' if it doesn't exist
  if (!pStageID) {
     pStageID = PyDict_New();
     *stat = PyDict_SetItemString(pFGStageID, fg_key, pStageID);
     if (*stat<0) {
        PyErr_Print();
        return;
     }
  }

  // Insert the stage ID
  char *stagename = fix_string(stage, stagelen);
  PyObject *pID = PyFloat_FromDouble(*id);
  *stat = PyDict_SetItemString(pStageID, stagename, pID);
  return;
}

/*****************************
 * Embedded Python interfaces
 *****************************/




/*****************************************
 * Callbacks via Python module "lebiology"
 *****************************************/

static PyObject *lebiology_addagent(PyObject *self, PyObject *args){
  PyObject *pAgentDict;
  int ok;

  if (!PyArg_ParseTuple(args, "O", pAgentDict)) {
    return NULL;
  }
  return Py_BuildValue("i", 0);
}
