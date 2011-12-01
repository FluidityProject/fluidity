#ifdef HAVE_PYTHON
#include "Python.h"
#endif

// Run a python string and store local namespace in a global dictionary for later evaluation
void python_run_string_store_locals_c(char *str, int strlen, char *dict, int dictlen, char *key, int keylen, int *stat);

// Evaluate the detector random walk via the val() function from a previously stored local namespace found in dict under key. The interface is: val(ele, local_coords)
void python_run_detector_val_from_locals_c(int ele, int dim, double lcoords[], double *dt, char *dict, int dictlen, char *key, int keylen, double value[], int *stat);
