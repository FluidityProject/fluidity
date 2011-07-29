#include <Python.h>
#include <string.h>


//#define PyPickers_NUM 0

static PyObject *PickersError;
static void *cNodeOwnerFinderResetFunction;
static void *cNodeOwnerFinderSetInputFunction;

//void cNodeOwnerFinderReset(const int* id)

//void cNodeOwnerFinderSetInput(int* id, const double* positions, const int* enlist, const int* dim, const int* loc, const int* nnodes, const int* nelements)

static PyObject *
pickers_cNodeOwnerFinderReset(PyObject *self, PyObject *args)
{
  printf("cNodeOwnerFinderReset entered\n");

  void (*cNodeOwnerFinderReset)(const int* id) = NULL;
  int id;
  
  cNodeOwnerFinderReset = cNodeOwnerFinderResetFunction;

  if (!PyArg_ParseTuple(args, "i", &id))
    return NULL;
  
  (*cNodeOwnerFinderReset)(&id);
  
  Py_RETURN_NONE;
}

static PyObject *
pickers_cNodeOwnerFinderSetInput(PyObject *self, PyObject *args)
{
  printf("cNodeOwnerFinderSetInput entered\n");

  void (*cNodeOwnerFinderSetInput)(int* id, const double* positions, const int* enlist, const int* dim, 
          const int* loc, const int* nnodes, const int* nelements) = NULL;
  int id;
  double positions;
  int enlist;
  int dim;
  int loc;
  int nnodes;
  int nelements;
  
  cNodeOwnerFinderSetInput = cNodeOwnerFinderSetInputFunction;

  if (!PyArg_ParseTuple(args, "ifiiiii", &id, &positions, &enlist, &dim, &loc, &nnodes, &nelements))
    return NULL;
  
  (*cNodeOwnerFinderSetInput)(&id, &positions, &enlist, &dim, &loc, &nnodes, &nelements);
  
  Py_RETURN_NONE;
}




static PyMethodDef pickersMethods[] = {
  {"cNodeOwnerFinderReset",  pickers_cNodeOwnerFinderReset, METH_VARARGS,
   "cNodeOwnerFinderReset in Node_Owner_Finder.cpp."},
  {"cNodeOwnerFinderSetInput",  pickers_cNodeOwnerFinderSetInput, METH_VARARGS,
   "cNodeOwnerFinderSetInput in Node_Owner_Finder.cpp."},
  {NULL, NULL, 0, NULL},
            /* Sentinel */
};


PyMODINIT_FUNC
initpickers(void)
{
  PyObject *m;

  m = Py_InitModule("pickers", pickersMethods);
  if (m == NULL)
    return;
        
  cNodeOwnerFinderResetFunction = (void *)PyCapsule_Import("fluidity_api._cNodeOwnerFinderReset", 0);
  if (cNodeOwnerFinderResetFunction == NULL){
    printf("cNodeOwnerFinderResetFunction is NULL\n");
    return NULL;
  }
  cNodeOwnerFinderSetInputFunction = (void *)PyCapsule_Import("fluidity_api._cNodeOwnerFinderSetInput", 0);
  if (cNodeOwnerFinderSetInputFunction == NULL){
    printf("cNodeOwnerFinderSetInputFunction is NULL\n");
    return NULL;
  }
    
  PickersError = PyErr_NewException("Pickers.error", NULL, NULL);

  Py_INCREF(PickersError);

  PyModule_AddObject(m, "PickersError", PickersError);

}







 //void cNodeOwnerFinderFind(const int* id, const double* position, const int* dim)
 
 //void cNodeOwnerFinderQueryOutput(const int* id, int* nelms)
 
//  void cNodeOwnerFinderGetOutput(const int* id, int* ele_id, const int* index)
