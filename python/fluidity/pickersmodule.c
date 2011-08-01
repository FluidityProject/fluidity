#include <Python.h>
#include <string.h>




static PyObject *PickersError;
static PyObject *PickersCapsuleImportError;


static void *cNodeOwnerFinderResetFunction;
static void *cNodeOwnerFinderSetInputFunction;
static void *cNodeOwnerFinderQueryOutputFunction;
static void *cNodeOwnerFinderGetOutputFunction;
static void *cNodeOwnerFinderFindFunction;

//void cNodeOwnerFinderReset(const int* id)

//void cNodeOwnerFinderSetInput(int* id, const double* positions, const int* enlist, const int* dim, const int* loc, const int* nnodes, const int* nelements)
//(const int* id, int* nelms)

//  void cNodeOwnerFinderGetOutput(const int* id, int* ele_id, const int* index)
// //void cNodeOwnerFinderFind(const int* id, const double* position, const int* dim)

static PyObject *
pickers_cNodeOwnerFinderFind(PyObject *self, PyObject *args)
{
  printf("cNodeOwnerFinderFind entered\n");

  void (*cNodeOwnerFinderFind)(const int* id, const double* position, const int* dim) = NULL;
  int id;
  int dim;
  
  PyObject* pyid;
  PyObject* pyposition;
  PyObject* pydim;
  
  pyid = PyTuple_GetItem(args, 0);
  pyposition = PyTuple_GetItem(args, 1);
  int check = PyList_Check(pyposition);
  printf("%i\n", check);
  pydim = PyTuple_GetItem(args, 2);
  
  PyArg_Parse(pyid, "i", &id);
  PyArg_Parse(pydim, "i", &dim);
    
  double positions [dim];
  int i;
  for (i = 0; i < dim; i++){
    double element = -1.0;
    PyObject* pylistElement = PyList_GetItem(pyposition, i);
    
    PyArg_Parse(pylistElement, "d", &element);
    positions[i] = element;
    printf("what is element? %f\n", element);
  }
  //positions[dim] = '\0';
  
  printf("cNodeOwnerFinderFind before assignment\n");
  cNodeOwnerFinderFind = cNodeOwnerFinderFindFunction;
  if (cNodeOwnerFinderFind == NULL){
    printf("cNodeOwnerFinderFind is NULL in function\n");
  }
  printf("what is id? %i\n", id);
  printf("what is dim? %i\n", dim);
  printf("what is positions[0]? %f\n", positions[0]);
  
  (*cNodeOwnerFinderFind)(&id, positions, &dim);
  
  printf("after function pointer executed\n");
  Py_RETURN_NONE;
}

static PyObject *
pickers_cNodeOwnerFinderQueryOutput(PyObject *self, PyObject *args)
{
  printf("cNodeOwnerFinderQueryOutput entered\n");
  int id;
  int nelms;
  PyObject *pyid;
  void (*cNodeOwnerFinderQueryOutput)(const int* id, int* nelms) = NULL;
  
  pyid = PyTuple_GetItem(args, 0);
  PyArg_Parse(pyid, "i", &id);
  
  cNodeOwnerFinderQueryOutput = cNodeOwnerFinderQueryOutputFunction;

  (*cNodeOwnerFinderQueryOutput)(&id, &nelms);
  
  return Py_BuildValue("i", nelms);
}

static PyObject *
pickers_cNodeOwnerFinderGetOutput(PyObject *self, PyObject *args)
{
  printf("cNodeOwnerFinderGetOutput entered\n");
  int id;
  int ele_id;
  int index;
  PyObject *pyid;
  PyObject *pyindex;
  void (*cNodeOwnerFinderGetOutput)(const int* id, int* ele_id, const int* index) = NULL;
  
  pyid = PyTuple_GetItem(args, 0);
  pyindex = PyTuple_GetItem(args, 1);

  PyArg_Parse(pyid, "i", &id);
  PyArg_Parse(pyindex, "i", &index);
  
  cNodeOwnerFinderGetOutput = cNodeOwnerFinderGetOutputFunction;

  (*cNodeOwnerFinderGetOutput)(&id, &ele_id, &index);
  
  return Py_BuildValue("i", ele_id);
}



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
  double positions [pypositionsSize];
  int i;
  for (i = 0; i < pypositionsSize; i++){
    double element;
    PyObject* pylistElement = PyList_GetItem(pypositions, i);
    PyArg_Parse(pylistElement, "d", &element);
    positions[i] = element;
  }
  
  enlistSize = PyList_Size(pyenlist);
  int enlist[enlistSize];
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
  {"cNodeOwnerFinderQueryOutput",  pickers_cNodeOwnerFinderQueryOutput, METH_VARARGS,
   "cNodeOwnerFinderQueryOutput in Node_Owner_Finder.cpp."},
  {"cNodeOwnerFinderGetOutput",  pickers_cNodeOwnerFinderGetOutput, METH_VARARGS,
   "cNodeOwnerFinderGetOutput in Node_Owner_Finder.cpp."},
  {"cNodeOwnerFinderFind",  pickers_cNodeOwnerFinderFind, METH_VARARGS,
   "cNodeOwnerFinderFind in Node_Owner_Finder.cpp."},
  {NULL, NULL, 0, NULL},
            /* Sentinel */
};


PyMODINIT_FUNC
initpickers(void)
{
  printf("initialisation entered\n");
  PyObject *m;

  m = Py_InitModule("pickers", pickersMethods);
  if (m == NULL)
    return;
        
  cNodeOwnerFinderResetFunction = (void *)PyCapsule_Import("fluidity_api._cNodeOwnerFinderReset", 0);
  if (cNodeOwnerFinderResetFunction == NULL){
    printf("cNodeOwnerFinderResetFunction is NULL in initialisation\n");
    PyErr_SetString(PickersCapsuleImportError, "cNodeOwnerFinderResetFunction is NULL");
    return;
  }
  else{
    printf("cNodeOwnerFinderResetFunction is not NULL in initialisation\n");
  }
  cNodeOwnerFinderSetInputFunction = (void *)PyCapsule_Import("fluidity_api._cNodeOwnerFinderSetInput", 0);
  if (cNodeOwnerFinderSetInputFunction == NULL){
    printf("cNodeOwnerFinderSetInputFunction is NULL in initialisation\n");
    PyErr_SetString(PickersCapsuleImportError, "cNodeOwnerFinderSetInputFunction is NULL");
    return;
  }
  else{
    printf("cNodeOwnerFinderSetInputFunction is not NULL in initialisation\n");
  }
  cNodeOwnerFinderQueryOutputFunction = (void *)PyCapsule_Import("fluidity_api._cNodeOwnerFinderQueryOutput", 0);
  if (cNodeOwnerFinderQueryOutputFunction == NULL){
    printf("cNodeOwnerFinderQueryOutputFunction is NULL in initialisation\n");
    PyErr_SetString(PickersCapsuleImportError, "cNodeOwnerFinderQueryOutputFunction is NULL");
    return;
  }
  else{
    printf("cNodeOwnerFinderQueryOutputFunction is not NULL in initialisation\n");
  }
  cNodeOwnerFinderGetOutputFunction = (void *)PyCapsule_Import("fluidity_api._cNodeOwnerFinderGetOutput", 0);
  if (cNodeOwnerFinderGetOutputFunction == NULL){
    printf("cNodeOwnerFinderGetOutputFunction is NULL in initialisation\n");
    PyErr_SetString(PickersCapsuleImportError, "cNodeOwnerFinderGetOutputFunction is NULL");
    return;
  }
  else{
    printf("cNodeOwnerFinderGetOutputFunction is not NULL in initialisation\n");
  }
  cNodeOwnerFinderFindFunction = (void *)PyCapsule_Import("fluidity_api._cNodeOwnerFinderFind", 0);
  if (cNodeOwnerFinderFindFunction == NULL){
    printf("cNodeOwnerFinderFindFunction is NULL in initialisation\n");
    PyErr_SetString(PickersCapsuleImportError, "cNodeOwnerFinderFindFunction is NULL");
    return;
  }
  else{
    printf("cNodeOwnerFinderFindFunction is not NULL in initialisation\n");
  }
    
  PickersError = PyErr_NewException("Pickers.error", NULL, NULL);

  Py_INCREF(PickersError);

  PyModule_AddObject(m, "PickersError", PickersError);

}








 

 

