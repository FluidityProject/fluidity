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
  double *positions;
  int *enlist;
  int dim;
  int loc;
  int nnodes;
  int nelements;
  PyObject* pyid;
  PyObject* pypositions;
  PyObject* pyenlist;
  PyObject* pydim;
  PyObject* pyloc;
  PyObject* pynnodes;
  PyObject* pynelements;
  int pypositionsSize;
  int enlistSize;
  
  pyid = PyTuple_GetItem(args, 0);
  pypositions = PyTuple_GetItem(args, 1);
  pyenlist = PyTuple_GetItem(args, 2);
  pydim = PyTuple_GetItem(args, 3);
  pyloc = PyTuple_GetItem(args, 4);
  pynnodes = PyTuple_GetItem(args, 5);
  pynelements = PyTuple_GetItem(args, 6);
  
  PyArg_Parse(pyid, "i", &id);
  
  pypositionsSize = PyList_Size(pypositions);
  positions = malloc(pypositionsSize*sizeof(double));
  int i;
  for (i = 0; i < pypositionsSize; i++){
    double element;
    PyObject* pylistElement = PyList_GetItem(pypositions, i);
    PyArg_Parse(pylistElement, "f", &element);
    positions[i] = element;
  }
  
  enlistSize = PyList_Size(pyenlist);
  enlist = malloc(enlistSize*sizeof(double));
  int j;
  for (j = 0; j < enlistSize; j++){
    int element;
    PyObject* pylistElement = PyList_GetItem(pyenlist, j);
    PyArg_Parse(pylistElement, "i", &element);
    enlist[j] = element;
  }
  
  PyArg_Parse(pydim, "i", &dim);
  PyArg_Parse(pyloc, "i", &loc);
  PyArg_Parse(pynnodes, "i", &nnodes);
  PyArg_Parse(pynelements, "i", &nelements);
  
  cNodeOwnerFinderSetInput = cNodeOwnerFinderSetInputFunction;

  (*cNodeOwnerFinderSetInput)(&id, positions, enlist, &dim, &loc, &nnodes, &nelements);
  
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
