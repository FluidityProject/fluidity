#define ALLOW_IMPORT_ARRAY
#include "python_statec.h"

void python_init_(void){
#ifdef HAVE_PYTHON
  // Initialize the Python interpreter
  Py_Initialize();
  PyRun_SimpleString("import string");  
#endif
#ifdef HAVE_NUMPY
  // Enable use of NumPy arrays in C
  import_array();

  // Import the NumPy module in our Python interpreter
  if(PyRun_SimpleString("import numpy") == -1)
    fprintf(stderr,"Error: Importing the NumPy module failed.\n");

  // Initialize a persistent dictionary
  PyRun_SimpleString("persistent = {}");

  // Add the working directory to the module search path.
  PyRun_SimpleString("import sys");
  PyRun_SimpleString("sys.path.append('.')");
  
  init_vars();
#endif
}


void init_vars(void){
#ifdef HAVE_PYTHON
  // Import the types and setup the states dictionary
  // This is called every time the data is reset
  if(PyRun_SimpleString("from fluidity.state_types import *") == -1){
    fprintf(stderr, "Warning: The 'state_types.py' module could not be loaded. Make sure the PYTHONPATH environment variable is set.\n");
    fprintf(stderr, "This is a problem if you have ocean biology or Python diagnostic fields.\n");
    fprintf(stderr, "It will not otherwise affect the running of Fluidity.\n");
    PyErr_Clear();
  }
  PyRun_SimpleString("states = dict()");
#endif
}


void python_reset_(void){
#ifdef HAVE_PYTHON
  if(Py_IsInitialized()){
    // Create a list of items to be kept
    PyRun_SimpleString("keep = [];keep.append('keep');keep.append('rem');keep.append('__builtins__');keep.append('__name__'); keep.append('__doc__');keep.append('string');keep.append('numpy');keep.append('persistent')");

    // Create a list of items to  be removed
    PyRun_SimpleString("rem = []");
    PyRun_SimpleString("for i in globals().keys():\n if(not (i in keep)): rem.append(i)");

    // Delete every item except the ones we want to keep
    PyRun_SimpleString("for i in rem: del globals()[i]");
    PyRun_SimpleString("del globals()['keep'];del globals()['rem'];del globals()['i']");

    // Reinitialize the variables
    init_vars();

    // And run a garbage collection
    PyGC_Collect();
  }
#endif
}


void python_end_(void){
#ifdef HAVE_PYTHON
  if(Py_IsInitialized()){
    // Garbage collection
    PyGC_Collect();
    // Finalize the Python interpreter
    Py_Finalize();
  }
#endif
}


void python_run_stringc_(char *s,int *slen, int *stat){
#ifdef HAVE_PYTHON
  // Run a python command from Fortran
  char *c = fix_string(s,*slen);
  char t[8+*slen];
  sprintf(t,"%s\n",c);
  *stat = PyRun_SimpleString(t);
  if(*stat != 0){
    PyErr_Print();
  }
  free(c); 
#endif
}


void python_run_filec_(char *f,int *flen, int *stat){
#ifdef HAVE_PYTHON
  // Run a python file from Fortran
  char *filename = fix_string(f,*flen);
  FILE *pyfile;
  if ((pyfile = fopen(filename, "r")) == NULL){
    fprintf(stderr, "Error: cannot open '%s'. \n", filename);
    *stat = 1;
  }else {
    *stat = PyRun_SimpleFileExFlags(pyfile,filename,1,NULL);
    if(*stat != 0){
      PyErr_Print();
    }
  }
  free(filename);
#endif
}


char* fix_string(char *s,int len){
  char *ns = (char *)malloc(len+3);
  memcpy( ns, s, len );
  ns[len] = 0;
  return ns;
}








// Functions to add a state and fields: scalar, vector, tensor, mesh, quadrature, polynomial


void python_add_statec_(char *name,int *len){
#ifdef HAVE_PYTHON
  // Add a new state object to the Python environment
  char *n = fix_string(name,*len);
  char t[20+*len+*len];
  // 'state' in Python will always be the last state added while the 'states' dictionary 
  // includes all added states
  sprintf(t,"states[\"%s\"] = State(\"%s\")",n,n);
  PyRun_SimpleString(t);
  sprintf(t,"state = states[\"%s\"]",n);
  PyRun_SimpleString(t);

  free(n); 
#endif
}


void python_add_scalar_(int *sx,double x[],char *name,int *nlen, int *field_type, 
  char *option_path, int *oplen, char *state,int *slen,
  char *mesh_name, int *mesh_name_len){
#ifdef HAVE_NUMPY
  // Add the Fortran scalar field to the dictionary of the Python interpreter
  PyObject *pMain = PyImport_AddModule("__main__");
  PyObject *pDict = PyModule_GetDict(pMain);
  // Fix the Fortran strings for C and Python
  char *namec = fix_string(name,*nlen);
  char *opc = fix_string(option_path,*oplen);
  char *meshc = fix_string(mesh_name,*mesh_name_len);

  // Create the array
  python_add_array_double_1d(x,sx,"s");

  PyObject *pname = PyString_FromString(namec);
  PyDict_SetItemString(pDict,"n",pname); 
  PyObject *poptionp = PyString_FromString(opc);
  PyDict_SetItemString(pDict,"op",poptionp);  
  PyObject *pft = PyInt_FromLong(*field_type);
  PyDict_SetItemString(pDict,"ft",pft);  

  PyRun_SimpleString("n = string.strip(n)");
  PyRun_SimpleString("op = string.strip(op)");

  char *n = fix_string(state,*slen);
  char t[150+*slen+*mesh_name_len];
  sprintf(t,"field = ScalarField(n,s,ft,op); states['%s'].scalar_fields['%s'] = field",n,namec);
  PyRun_SimpleString(t);

  // Set the mesh for this field
  sprintf(t,"field.set_mesh(states['%s'].meshes['%s'])",n,meshc);
  PyRun_SimpleString(t);

  // Clean up
  PyRun_SimpleString("del n; del op; del ft; del s; del field");
  free(namec); 
  free(opc);
  free(n);
  free(meshc);
  Py_DECREF(pname);
  Py_DECREF(poptionp);
  Py_DECREF(pft);
#endif
}


void python_add_vector_(int *s, int *num_dim, 
  double x[], double y[], double z[],
  char *name,int *nlen, int *field_type, char *option_path, int *oplen, char *state,int *slen,
  char *mesh_name, int *mesh_name_len){
#ifdef HAVE_NUMPY
  // Make the Fortran vector field availabe to the Python interpreter
  PyObject *pMain = PyImport_AddModule("__main__");
  PyObject *pDict = PyModule_GetDict(pMain);

  if(*num_dim>0)
    python_add_array_double_1d(x,s,"vector1");
  else
    PyRun_SimpleString("vector1 = numpy.matrix([])");

  if(*num_dim>1)
    python_add_array_double_1d(y,s,"vector2");
  else
    PyRun_SimpleString("vector2 = numpy.matrix([])");

  if(*num_dim>2)
    python_add_array_double_1d(z,s,"vector3");
  else
    PyRun_SimpleString("vector3 = numpy.matrix([])");

  // Fix the Fortran strings for C and Python
  char *namec = fix_string(name,*nlen);
  char *opc = fix_string(option_path,*oplen);
  char *meshc = fix_string(mesh_name,*mesh_name_len);

  PyObject *pname = PyString_FromString(namec);
  PyDict_SetItemString(pDict,"n",pname); 
  PyObject *poptionp = PyString_FromString(opc);
  PyDict_SetItemString(pDict,"op",poptionp);  
  PyObject *pft = PyInt_FromLong(*field_type);
  PyDict_SetItemString(pDict,"ft",pft);  
  PyObject *pnd = PyInt_FromLong(*num_dim);
  PyDict_SetItemString(pDict,"nd",pnd);  

  PyRun_SimpleString("n = string.strip(n)");
  PyRun_SimpleString("op = string.strip(op)");

  char *n = fix_string(state,*slen);
  char t[150+*slen+*mesh_name_len];
  sprintf(t,"field = VectorField(n,vector1,vector2,vector3,ft,op,nd); states[\"%s\"].vector_fields['%s'] = field",n,namec);
  PyRun_SimpleString(t);

  // Set the mesh for this field
  sprintf(t,"field.set_mesh(states['%s'].meshes['%s'])",n,meshc);
  PyRun_SimpleString(t);

  // Clean up
  PyRun_SimpleString("del n; del op; del ft; del nd; del vector1; del vector2; del vector3; del field");
  free(n);
  free(namec); 
  free(opc);
  free(meshc);

  Py_DECREF(pname);
  Py_DECREF(poptionp);
  Py_DECREF(pft);
  Py_DECREF(pnd);
#endif
}


void python_add_tensor_(int *sx,int *sy,int *sz, double *x, int *num_dim,
  char *name,int *nlen, int *field_type, char *option_path, int *oplen, char *state,int *slen,
  char *mesh_name, int *mesh_name_len){
#ifdef HAVE_NUMPY
  // Expose a Fortran tensor field to the Python interpreter
  PyObject *pMain = PyImport_AddModule("__main__");
  PyObject *pDict = PyModule_GetDict(pMain);

  python_add_array_double_3d(x,sx,sy,sz,"val");
  PyRun_SimpleString("val = val.transpose(2, 0, 1)");

  // Fix the Fortran strings for C and Python

  char *namec = fix_string(name,*nlen);
  char *opc = fix_string(option_path,*oplen);
  char *meshc = fix_string(mesh_name,*mesh_name_len);

  PyObject *pname = PyString_FromString(namec);
  PyDict_SetItemString(pDict,"n",pname); 
  PyObject *poptionp = PyString_FromString(opc);
  PyDict_SetItemString(pDict,"op",poptionp);  
  PyObject *pft = PyInt_FromLong(*field_type);
  PyDict_SetItemString(pDict,"ft",pft);  
  PyObject *pnd = PyInt_FromLong(*num_dim);
  PyDict_SetItemString(pDict,"nd",pnd);  

  PyRun_SimpleString("n = string.strip(n)");
  PyRun_SimpleString("op = string.strip(op)");

  char *n = fix_string(state,*slen);
  char t[150+*slen+*mesh_name_len];
  sprintf(t,"field = TensorField(n,val,ft,op,nd); states[\"%s\"].tensor_fields['%s'] = field",n,namec);
  PyRun_SimpleString(t);  

  // Set the mesh for this field
  sprintf(t,"field.set_mesh(states['%s'].meshes['%s'])",n,meshc);
  PyRun_SimpleString(t);

  // Clean up
  PyRun_SimpleString("del n; del op; del val; del field");
  free(n);
  free(namec); 
  free(opc); 
  free(meshc);

  Py_DECREF(pname);
  Py_DECREF(poptionp);
  Py_DECREF(pft);
  Py_DECREF(pnd);
#endif
}


void python_add_mesh_(int ndglno[],int *sndglno, int *elements, int *nodes, 
  char *name,int *nlen, char *option_path, int *oplen, 
  int *continuity, int region_ids[], int *sregion_ids,
  char *state_name, int *state_name_len){
#ifdef HAVE_NUMPY
  // Add the Mesh to the interpreter
  PyObject *pMain = PyImport_AddModule("__main__");
  PyObject *pDict = PyModule_GetDict(pMain);

  // Fix the Fortran strings for C and Python
  char *namec = fix_string(name,*nlen);
  char *opc = fix_string(option_path,*oplen);

  python_add_array_integer_1d(ndglno, sndglno,"mesh_array");

  python_add_array_integer_1d(region_ids, sregion_ids,"region_ids");

  PyObject *pname = PyString_FromString(namec);
  PyDict_SetItemString(pDict,"n",pname);
  PyObject *poptionp = PyString_FromString(opc);
  PyDict_SetItemString(pDict,"op",poptionp);

  PyRun_SimpleString("n = string.strip(n)");
  PyRun_SimpleString("op = string.strip(op)");

  char *n = fix_string(state_name,*state_name_len);
  char t[150+*state_name_len];

  sprintf(t,"states[\"%s\"].meshes[\"%s\"] = Mesh(mesh_array,%d,%d,%d,n,op,region_ids)",n,namec,*elements,*nodes,*continuity);
  PyRun_SimpleString(t);
  sprintf(t,"m = states[\"%s\"].meshes[\"%s\"]",n,namec);
  PyRun_SimpleString(t);

  // Clean up
  PyRun_SimpleString("del n; del op; del m");
  free(namec); 
  free(n); 
  free(opc);

  Py_DECREF(pname);
  Py_DECREF(poptionp);
#endif
}


void python_add_element_(int *dim, int *loc, int *ngi, int *degree,
  char *state_name, int *state_name_len, char *mesh_name, int *mesh_name_len,
  double *n,int *nx, int *ny, double *dn, int *dnx, int *dny, int *dnz,
  int *size_spoly_x,int *size_spoly_y,int *size_dspoly_x,int *size_dspoly_y ){
#ifdef HAVE_NUMPY
  // Fix the Fortran strings for C and Python
  char *meshc = fix_string(mesh_name,*mesh_name_len);
  char *statec = fix_string(state_name,*state_name_len);
  char t[80+*mesh_name_len+*state_name_len];

  // Set n
  python_add_array_double_2d(n,nx,ny,"n_array");
  // Set dn
  python_add_array_double_3d(dn,dnx,dny,dnz,"dn_array");

  // Add the element to the interpreter and make the element variable available
  // so the other attributes in Fortran can be passed in
  sprintf(t,"element = Element(%d,%d,%d,%d,n_array,dn_array,%d,%d,%d,%d); states['%s'].meshes['%s'].shape = element",
    *dim,*loc,*ngi,*degree,
    *size_spoly_x,*size_spoly_y,*size_dspoly_x,*size_dspoly_y,
    statec, meshc
  );
  PyRun_SimpleString(t);

  PyRun_SimpleString("del n_array; del dn_array");
#endif
}


void python_add_quadrature_(int *dim,int *degree,int *loc, int *ngi, 
  double *weight, int *weight_size, double *locations, int *l_size, int *is_surfacequadr){
  // Only being called right after an element has been added
#ifdef HAVE_NUMPY
  // Set weights
  python_add_array_double_1d(weight,weight_size,"weight");

  // Set locations
  python_add_array_double_1d(locations,l_size,"locations");

  char c[150];
  if(*is_surfacequadr == 1)
    sprintf(c,"element.set_surface_quadrature(Quadrature(weight,locations,%d,%d,%d,%d))",
      *dim,*degree,*loc,*ngi);
  else if(*is_surfacequadr == 0)
    sprintf(c,"element.set_quadrature(Quadrature(weight,locations,%d,%d,%d,%d))",
      *dim,*degree,*loc,*ngi);
  PyRun_SimpleString(c);

  PyRun_SimpleString("del weight; del locations");
#endif
}


void python_add_polynomial_(double *coefs,int *size,int *degree, int *x,int *y, int *spoly){
#ifdef HAVE_NUMPY
  // Add a polynomial to the latest element
  // Set the coefs array
  python_add_array_double_1d(coefs, size, "coefs");

  // Set the polynomial to the element
  char c[120];
  if(*spoly == 1)
    sprintf(c, "element.set_polynomial_s(Polynomial(coefs,%d),%d,%d)",*degree,*x,*y);
  else if(*spoly == 0)
    sprintf(c, "element.set_polynomial_ds(Polynomial(coefs,%d),%d,%d)",*degree,*x,*y);
  PyRun_SimpleString(c);
  PyRun_SimpleString("del coefs");
//   Py_DECREF(arr); 
#endif
}



// Interface for adding arrays

void python_add_array_double_1d(double *arr, int *size, char *name){
#ifdef HAVE_NUMPY
  // Add an array in Python which will be availabe under the variable name 'name'
  PyObject *pMain = PyImport_AddModule("__main__");
  PyObject *pDict = PyModule_GetDict(pMain);

  // Create our NumPy matrix struct:
  // Arguments are: 
  //  number of dimensions (int),
  //  size of each dimension (int[]),
  //  data type to determine the width of each element in memory (int)
  //  the actual data as a byte array(char*)

  // Set the array
  PyObject *a = PyArray_SimpleNewFromData(1, (npy_intp[]){*size}, PyArray_DOUBLE, (char*)arr);
  PyDict_SetItemString(pDict,name,a);
  Py_DECREF(a);
#endif
}

void python_add_array_double_2d(double *arr, int *sizex, int *sizey, char *name){
#ifdef HAVE_NUMPY
  // Add an array in Python which will be availabe under the variable name 'name'
  PyObject *pMain = PyImport_AddModule("__main__");
  PyObject *pDict = PyModule_GetDict(pMain);

  // Set the array
  npy_intp dims[] = {*sizey,*sizex};
  PyObject *a = PyArray_SimpleNewFromData(2, dims, PyArray_DOUBLE, (char*)arr);
  PyDict_SetItemString(pDict,name,a);
  char c[200];
  sprintf(c,"%s = numpy.transpose(%s,(1,0))",name,name);
  PyRun_SimpleString(c);
  Py_DECREF(a);
#endif
}

void python_add_array_double_3d(double *arr, int *sizex, int *sizey, int *sizez, char *name){
#ifdef HAVE_NUMPY
  // Add an array in Python which will be availabe under the variable name 'name'
  PyObject *pMain = PyImport_AddModule("__main__");
  PyObject *pDict = PyModule_GetDict(pMain);

  // Set the array
  npy_intp dims[] = {*sizez,*sizey,*sizex};
  PyObject *a = PyArray_SimpleNewFromData(3, dims, PyArray_DOUBLE, (char*)arr);
  PyDict_SetItemString(pDict,name,a);
  char c[200];
  sprintf(c,"%s = numpy.transpose(%s,(2,1,0))",name,name);
  PyRun_SimpleString(c);
  Py_DECREF(a);
#endif
}

void python_add_array_integer_1d(int *arr, int *size, char *name){
#ifdef HAVE_NUMPY
  // Add an array in Python which will be availabe under the variable name 'name'
  PyObject *pMain = PyImport_AddModule("__main__");
  PyObject *pDict = PyModule_GetDict(pMain);

  // Set the array
  PyObject *a = PyArray_SimpleNewFromData(1, (npy_intp[]){*size}, PyArray_INT, (char*)arr);
  PyDict_SetItemString(pDict,name,a);
  Py_DECREF(a);
#endif
}

void python_add_array_integer_2d(int *arr, int *sizex, int *sizey, char *name){
#ifdef HAVE_NUMPY
  // Add an array in Python which will be availabe under the variable name 'name'
  PyObject *pMain = PyImport_AddModule("__main__");
  PyObject *pDict = PyModule_GetDict(pMain);

  // Set the array
  npy_intp dims[] = {*sizey,*sizex};
  PyObject *a = PyArray_SimpleNewFromData(2, dims, PyArray_INT, (char*)arr);
  PyDict_SetItemString(pDict,name,a);
  char c[200];
  sprintf(c,"%s = numpy.transpose(%s,(1,0))",name,name);
  PyRun_SimpleString(c);
  Py_DECREF(a);
#endif
}

void python_add_array_integer_3d(int *arr, int *sizex, int *sizey, int *sizez, char *name){
#ifdef HAVE_NUMPY
  // Add an array in Python which will be availabe under the variable name 'name'
  PyObject *pMain = PyImport_AddModule("__main__");
  PyObject *pDict = PyModule_GetDict(pMain);

  // Set the array
  npy_intp dims[] = {*sizez,*sizey,*sizex};
  PyObject *a = PyArray_SimpleNewFromData(3, dims, PyArray_INT, (char*)arr);
  PyDict_SetItemString(pDict,name,a);
  char c[200];
  sprintf(c,"%s = numpy.transpose(%s,(2,1,0))",name,name);
  PyRun_SimpleString(c);
  Py_DECREF(a);
#endif
}




// Wrapper functions

void python_add_array_double_1d_(double *arr, int *size, char *name, int *name_len){
  // Called from Fortran
  char *namec = fix_string(name,*name_len);
  python_add_array_double_1d(arr, size, namec);
  free(namec);
}

void python_add_array_double_2d_(double *arr, int *sizex, int *sizey, char *name, int *name_len){
  // Called from Fortran
  char *namec = fix_string(name,*name_len);
  python_add_array_double_2d(arr, sizex,sizey, namec);
  free(namec);
}

void python_add_array_double_3d_(double *arr, int *sizex, int *sizey, int *sizez, char *name, int *name_len){
  // Called from Fortran
  char *namec = fix_string(name,*name_len);
  python_add_array_double_3d(arr, sizex,sizey,sizez, namec);
  free(namec);
}
void python_add_array_integer_1d_(int *arr, int *size, char *name, int *name_len){
  // Called from Fortran
  char *namec = fix_string(name,*name_len);
  python_add_array_integer_1d(arr, size, namec);
  free(namec);
}

void python_add_array_integer_2d_(int *arr, int *sizex, int *sizey, char *name, int *name_len){
  // Called from Fortran
  char *namec = fix_string(name,*name_len);
  python_add_array_integer_2d(arr, sizex,sizey, namec);
  free(namec);
}

void python_add_array_integer_3d_(int *arr, int *sizex, int *sizey, int *sizez, char *name, int *name_len){
  // Called from Fortran
  char *namec = fix_string(name,*name_len);
  python_add_array_integer_3d(arr, sizex,sizey,sizez, namec);
  free(namec);
}
