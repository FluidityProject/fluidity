#include "confdefs.h"
#include "python_statec.h"

#ifdef HAVE_PYTHON
#include "Python.h"
#endif

/* Initialisation function for the 'lebiology' module */
PyMODINIT_FUNC initlebiology(int dim);

static int lebiology_dim;

/* Static dictionaries that hold meta-data for access at runtime.
 * The first index is always the FGroup name, eg. pFGVarNames['FGroup'] 
 */
static PyObject *pFGLocalsDict; // Dict of local namespaces with compiled function objects
static PyObject *pFGVarNames;   // List of agent variable names
static PyObject *pFGEnvNames;   // List of environment field names
static PyObject *pFGFoodNames;  // Dict of foodname pointing to list of variety names
static PyObject *pFGStageID;    // Dict of stage->ID mappings

/* Add a variable name to the array under pFGVarNames['fg'] */
void lebiology_add_fg_varname_c(char *fg, int fglen, char *var, int varlen, int *stat);

/* Add an environment name to the array under pFGEnvNames['fg'] */
void lebiology_add_fg_envname_c(char *fg, int fglen, char *env, int envlen, int *stat);

/* Add a food variety name to the array under pFGFoodNames['fg']['food'] */
void lebiology_add_fg_foodname_c(char *fg, int fglen, char *food, int foodlen, char *variety, int varietylen, int *stat);

/* Add a stage ID to pFGStageID['fg'] */
void lebiology_add_fg_stage_id_c(char *fg, int fglen, char *stage, int stagelen, double *id, int *stat);

/* The main compile function runs the outer Python code and stores the local namespace, 
 * including the compiled function object for 'val', in pFGLocalsDict 
 */
void lebiology_compile_function_c(char *fg, int fglen, char *key, int keylen, char *func, int funclen, int *stat);

/* Initialise agent biology from a Python function
 * Usage: def val(agent):
 *          agent['var'] = value
 *          return agent
 */
void lebiology_agent_init_c(char *fg, int fglen, char *key, int keylen, double vars[], int n_vars, int *stat);

/* Update agent biology from a Python function
 * Usage: def val(agent, env, dt):
 *          agent['var'] = f( env['field'], dt )
 *          return agent
 */
void lebiology_agent_update_c(char *fg, int fglen, char *key, int keylen, char *food, int foodlen, 
                              double vars[], int n_vars, double envvals[], int n_envvals, 
                              double fvariety[], double frequest[], double fthreshold[], double fingest[], int n_fvariety, 
                              double *dt, int *stat);

/* Agent motion interface:
 * Usage: def val(position, vars, dt):
 *          vector = f(position, vars, dt)
 *          return vector
 */
// Note: Disabled until dependencies are sorted
void lebiology_agent_move_c(char *fg, int fglen, char *key, int keylen, double pos[], int n_pos, double vars[], int n_vars, int var_inds[], double *dt, double vector[], int *stat);

/* Callback function to add new agents to the system from inside the embedded Python biology update */
static PyObject *lebiology_add_agent(PyObject *self, PyObject *args);

/* Wrapper for Fortran function to add converted new agents to the system */
void fl_add_agent_c(double vars[], int *n_vars, double pos[], int *n_pos);

/* Callback function to translate stage names to internal IDs */
static PyObject *lebiology_stage_id(PyObject *self, PyObject *args);

/* Method definitions for lebilogy module */
static PyMethodDef LEBiologyMethods[] = {
  {"add_agent",  lebiology_add_agent, METH_VARARGS,
   "Add a new agent to the system"},
  {"stage_id",  lebiology_stage_id, METH_VARARGS,
   "Add a new agent to the system"},
  {NULL, NULL, 0, NULL}        /* Sentinel */
};
