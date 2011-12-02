#include "confdefs.h"
#include "python_statec.h"

#ifdef HAVE_PYTHON
#include "Python.h"
#endif

// Run a python string and store local namespace in a global dictionary for later evaluation
void python_run_string_store_locals_c(char *str, int strlen, char *dict, int dictlen, char *key, int keylen, int *stat);

// Evaluate the detector random walk via the val() function from a previously stored local namespace found in dict under key. The interface is: val(ele, local_coords)
void python_run_detector_val_from_locals_c(int ele, int dim, double lcoords[], double *dt, char *dict, int dictlen, char *key, int keylen, double value[], int *stat);

/* Initialise agent biology variables via a dictionary of name -> value pairs. This will create a dict of all variables stored by the agent, set the stage explicitly,
 * and let the user fill in the blanks in the val(biovars) function. */
void python_run_agent_biology_init_c(char *function, int function_len, char *var_list, int var_list_len, double biovars[], int n_biovars, double *stage_id, int *stat);

/* Evaluate the detector biology update function using a previously stored local namespace found in dict under key. 
 * The interface is: val(biovars, envvars, dt), where biovars is a dict of agent variables, 
 * envvars a dict of environment field names to locally evaluated field values, and dt is the timestep. */
void python_run_agent_biology_c(double *dt, char *dict, int dictlen, char *key, int keylen, double biovars[], int n_biovars, double env_vars[], int n_env_vars, int *stat);
