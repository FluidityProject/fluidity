#include "confdefs.h"
#include "python_statec.h"

#ifdef HAVE_PYTHON
#include "Python.h"
#endif

/* Initialisation function for the 'lebiology' module */
PyMODINIT_FUNC initlebiology(void);

/* Static dictionaries that hold meta-data for access at runtime.
 * The first index is always the FGroup name, eg. pFGVarNames['FGroup']
 */
static PyObject *pFGVarNames;  // List of variable names
static PyObject *pFGStageID;   // Dict of stage->ID mappings

/* Add a variable name to the array under pFGVarNames['fg'] */
void lebiology_add_fg_varname_c(char *fg, int fglen, char *var, int varlen, int *stat);

/* Add a stage ID to pFGStageID['fg'] */
void lebiology_add_fg_stage_id_c(char *fg, int fglen, char *stage, int stagelen, double *id, int *stat);

/* Callback function to add new agents to the system from inside the embedded Python biology update */
static PyObject *lebiology_addagent(PyObject *self, PyObject *args);

static PyMethodDef LEBiologyMethods[] = {
  {"addagent",  lebiology_addagent, METH_VARARGS,
   "Add a new agent to the system"},
  {NULL, NULL, 0, NULL}        /* Sentinel */
};
