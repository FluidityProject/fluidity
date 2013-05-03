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
static PyObject *pFGKernelFunc; // Dict of kernel function objects
static PyObject *pFGParamDicts ;// Dict of parameter->value mappings
static PyObject *pFGVarNames;   // List of agent variable names
static PyObject *pFGEnvNames;   // List of environment field names
static PyObject *pFGFoodNames;  // Dict of foodname pointing to list of variety names
static PyObject *pFGStageID;    // Dict of stage->ID mappings

static PyObject *pPersistent;   // Pointer to the current "persistent" dictionary

static PyObject *pJobServer;    // PP job server for threaded agent processing
static PyObject *pJobDict;      // Dict holding the task objects associated with each agent

static bool dropout_agent = false;

extern "C" {
/* Load the kernel fucntion and associated paramter set from the module provided and store in pFGKernelFunc and pFGParamDicts */
void lebiology_fg_kernel_load_c(char *fg, char *key, char *module, char *kernel, char *param, int *stat);

/* Reload the static referecen to the persisten dictionary */
void lebiology_reload_persistent_c();

/* Add a variable name to the array under pFGVarNames['fg'] */
void lebiology_add_fg_varname_c(char *fg, char *var, int *stat);

/* Add an environment name to the array under pFGEnvNames['fg'] */
void lebiology_add_fg_envname_c(char *fg, char *env, int *stat);

/* Add a food variety name to the array under pFGFoodNames['fg']['food'] */
void lebiology_add_fg_foodname_c(char *fg, char *food, char *variety, int *stat);

/* Add a stage ID to pFGStageID['fg'] */
void lebiology_add_fg_stage_id_c(char *fg, char *stage, double *id, int *stat);

/* The main compile function runs the outer Python code and stores the local namespace, 
 * including the compiled function object for 'val', in pFGLocalsDict 
 */
void lebiology_compile_function_c(char *fg, char *key, char *func, int *stat);

/* Initialise agent biology from a Python function
 * Usage: def val(agent):
 *          agent['var'] = value
 *          return agent
 */
void lebiology_agent_init_c(char *fg, char *key, double vars[], int n_vars, int *stat);

void lebiology_kernel_update_c(char *fg, char *key, char *food, double vars[], int n_vars, double envvals[], int n_envvals, 
			       double fvariety[], double frequest[], double fthreshold[], double fingest[], int n_fvariety, 
			       double *dt, int persistent, int *stat);

void lebiology_parallel_prepare_c(char *fg, char *key, char *food, double vars[], int n_vars, 
				  double envvals[], int n_envvals, double fvariety[], double fingest[], int n_fvariety, int agent_id, 
				  double *dt, int persistent, int *stat);

void lebiology_parallel_finish_c(char *fg, char *food, double vars[], int n_vars, double frequest[], 
				 double fthreshold[], int n_fvariety, int agent_id, int *stat);

/* Update agent biology from a Python function
 * Usage: def val(agent, env, dt):
 *          agent['var'] = f( env['field'], dt )
 *          return agent
 */
void lebiology_agent_update_c(char *fg, char *key, char *food, double vars[], int n_vars, double envvals[], int n_envvals, 
                              double fvariety[], double frequest[], double fthreshold[], double fingest[], int n_fvariety, 
                              double *dt, int *dropout, int *stat);

/* Agent motion interface:
 * Usage: def val(position, vars, dt):
 *          vector = f(position, vars, dt)
 *          return vector
 */
// Note: Disabled until dependencies are sorted
void lebiology_agent_move_c(char *fg, char *key, double pos[], int n_pos, double vars[], int n_vars, int var_inds[], double *dt, double vector[], int *stat);

/* Wrapper for Fortran function to add converted new agents to the system */
void fl_add_agent_c(double vars[], int *n_vars, double pos[], int *n_pos);
}

/* Callback function to add new agents to the system from inside the embedded Python biology update */
static PyObject *lebiology_add_agent(PyObject *self, PyObject *args);

/* Callback function to mark the current agent as a dropout to be removed at the end of the timestep */
static PyObject *lebiology_dropout_agent(PyObject *self, PyObject *args);

/* Callback function to translate stage names to internal IDs */
static PyObject *lebiology_stage_id(PyObject *self, PyObject *args);

/* Method definitions for lebilogy module */
static PyMethodDef LEBiologyMethods[] = {
  {"add_agent",  lebiology_add_agent, METH_VARARGS,
   "Add a new agent to the system"},
  {"dropout_agent",  lebiology_dropout_agent, METH_VARARGS,
   "Mark the current agent to be removed from the system"},
  {"stage_id",  lebiology_stage_id, METH_VARARGS,
   "Add a new agent to the system"},
  {NULL, NULL, 0, NULL}        /* Sentinel */
};
