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
//cNodeOwnerFinderQueryOutput(const int* id, int* nelms)

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


//cNodeOwnerFinderQueryOutput(const int* id, int* nelms)

 //void cNodeOwnerFinderFind(const int* id, const double* position, const int* dim)
//integer, dimension(size(positions, 2)), intent(out) :: ele_ids

static PyObject *
pickers_findserial(PyObject *self, PyObject *args)
{
  int closest_ele_id;
  int id;
  double out_of_bounds_tolerance = 0.1
  int i; 
  int j; 
  int r;
  int c;
  int nele_ids; 
  int possible_ele_id;
  double closest_miss; 
  double miss;
  int ele_ids = -1 ;
  PyObject *pyid;
  PyObject *pypositions;
  PyObject *pysublist;
  int positionsSizeDim1;
  int positionsSizeDim2;
  int ele_ids [positionsSizeDim2];
  void (*cNodeOwnerFinderFind)(const int* id, const double* position, const int* dim) = NULL;
  void (*cNodeOwnerFinderQueryOutput)(const int* id, int* nelms) = NULL;
  void (*cNodeOwnerFinderGetOutput)(const int* id, int* ele_id, const int* index) = NULL;
  cNodeOwnerFinderFind = cNodeOwnerFinderFindFunction;
  cNodeOwnerFinderQueryOutput = cNodeOwnerFinderQueryOutputFunction;
  cNodeOwnerFinderGetOutput = cNodeOwnerFinderGetOutputFunction;
  
  pyid = PyTuple_GetItem(args, 0);
  PyArg_Parse(pyid, "i", &id);
  
  positions = PyTuple_GetItem(args, 1);
  positionsSizeDim1 = PyList_Size(pypositions);
  pysublist = PyList_GetItem(pypositions,0);
  
  if (PyList_Check(pysublist)){
    positionsSizeDim2 = PyList_Size(pysublist);
  }
  else{
    positionsSizeDim2 = 1;
  }
  
  for (i = 0; i < positionsSizeDim2; i++){
  
    double positions[positionsSizeDim1];
    for (r = 0; r < positionsSizeDim1; r++){
      double element;
      pysublist = PyList_GetItem(pypositions,r);    
      PyObject *pysublistelement = PyList_GetItem(pysublist,i);
      PyArg_Parse(pysublistelement, "d", &element);
      positions[r] = element;
    }
    
    (*cNodeOwnerFinderFind)(&id, positions, &positionsSizeDim1);
    (*cNodeOwnerFinderQueryOutput)(&id, &nele_ids);
    closest_ele_id = -1;
    closest_miss = out_of_bounds_tolerance;
    for (j = 0; j < nele_ids; j++){
      (*cNodeOwnerFinderGetOutput)(&id, &possible_ele_id, &j);
        
      if(1){
        closest_ele_id = possible_ele_id;
        break;
      }
      else if (miss < closest_miss){
        closest_ele_id = possible_ele_id;
        closest_miss = miss;        
      }
    }
    
  ele_ids[i] = closest_ele_id;  
  }
    
}
  
  
  
  
  
  /*
  
 
  positions_loop: do i = 1, size(positions, 2)//how do i know the size of positions dim 2
  call cnode_owner_finder_find(id, positions(:, i), size(positions, 1))
  call cnode_owner_finder_query_output(id, nele_ids)

  closest_ele_id = -1
  ! We don't tolerate very large ownership failures
        closest_miss = out_of_bounds_tolerance
        do j = 1, nele_ids
          call cnode_owner_finder_get_output(id, possible_ele_id, j)
          ! Zero tolerance - we're not using an "epsilon-ball" approach here
          if(ownership_predicate(positions_a, possible_ele_id, positions(:, i), 0.0, miss = miss)) then
            ele_ids(i) = possible_ele_id
            ! We've found an owner - no need to worry about the closest miss
            cycle positions_loop
          else if(miss < closest_miss) then
            ! We didn't find an owner, but did find the closest miss so far
            closest_ele_id = possible_ele_id
            closest_miss = miss
          end if
        end do

        ! We didn't find an owner, so choose the element with the closest miss
        ele_ids(i) = closest_ele_id
          
      end do positions_loop

    end subroutine find_serial


}





subroutine node_owner_finder_find_multiple_positions(id, positions_a, positions, ele_ids, global)
    !!< For the node owner finder with ID id corresponding to positions
    !!< positions_a, find the element IDs owning the given positions.
    !!< This does not use ownership tolerances - instead, it determines the
    !!< "best" owning elements (those that are the smallest distance in ideal
    !!< space from test nodes).
    
    integer, intent(in) :: id
    type(vector_field), intent(in) :: positions_a
    real, dimension(:, :), intent(in) :: positions
    integer, dimension(size(positions, 2)), intent(out) :: ele_ids
    !! If present and .false., do not perform a global ownership test across all
    !! processes
    logical, optional, intent(in) :: global
   
    if(.not. present_and_false(global) .and. isparallel()) then
      call find_parallel()
    else
      call find_serial()
    end if

subroutine find_serial()
      integer :: closest_ele_id, i, j, nele_ids, possible_ele_id
      real :: closest_miss, miss
      
      ele_ids = -1
      positions_loop: do i = 1, size(positions, 2)
        call cnode_owner_finder_find(id, positions(:, i), size(positions, 1))
        call cnode_owner_finder_query_output(id, nele_ids)

        closest_ele_id = -1
        ! We don't tolerate very large ownership failures
        closest_miss = out_of_bounds_tolerance
        do j = 1, nele_ids
          call cnode_owner_finder_get_output(id, possible_ele_id, j)
          ! Zero tolerance - we're not using an "epsilon-ball" approach here
          if(ownership_predicate(positions_a, possible_ele_id, positions(:, i), 0.0, miss = miss)) then
            ele_ids(i) = possible_ele_id
            ! We've found an owner - no need to worry about the closest miss
            cycle positions_loop
          else if(miss < closest_miss) then
            ! We didn't find an owner, but did find the closest miss so far
            closest_ele_id = possible_ele_id
            closest_miss = miss
          end if
        end do

        ! We didn't find an owner, so choose the element with the closest miss
        ele_ids(i) = closest_ele_id
          
      end do positions_loop

    end subroutine find_serial
    
  function ownership_predicate_position(positions_a, ele_a, position, ownership_tolerance, miss, l_coords) result(owned)
    !!< Node ownership predicate. Returns .true. if the given position is
    !!< contained within element ele_a of positions_a to within tolerance
    !!< ownership_tolerance.
  
    type(vector_field), intent(in) :: positions_a
    integer, intent(in) :: ele_a
    real, dimension(positions_a%dim), intent(in) :: position
    real, intent(in) :: ownership_tolerance
    !!< Return the "miss" - the distance (in ideal space) of the test position
    !!< from the test element
    real, optional, intent(out) :: miss
    !!< Return the coordinate (in ideal space) of the test position
    !!< in the test element
    real, dimension(positions_a%dim + 1), optional, intent(out) :: l_coords
    
    logical :: owned
    
    real :: lmiss
    real, dimension(positions_a%dim + 1) :: ll_coords
   
    assert(ownership_tolerance >= 0.0)

    ll_coords = local_coords(positions_a, ele_a, position)
    
    assert(ele_numbering_family(positions_a, ele_a) == FAMILY_SIMPLEX)
    if(any(ll_coords < 0.0)) then
      lmiss = -minval(ll_coords)
      if(lmiss < ownership_tolerance) then
        owned = .true.
      else
        owned = .false.
      end if
      if(present(miss)) miss = lmiss
    else
      owned = .true.
      if(present(miss)) miss = 0.0
    end if

    if(present(l_coords)) l_coords = ll_coords
    
  end function ownership_predicate_position

  function ownership_predicate_node(positions_a, positions_b, ele_a, node_b, ownership_tolerance, miss, l_coords) result(owned)
    !!< Node ownership predicate. Returns .true. if the given node in
    !!< positions_b is contained within element ele_a of positions_a to within
    !!< tolerance ownership_tolerance.
    
    type(vector_field), intent(in) :: positions_a
    type(vector_field), intent(in) :: positions_b
    integer, intent(in) :: ele_a
    integer, intent(in) :: node_b
    real, intent(in) :: ownership_tolerance
    !!< Return the "miss" - the distance (in ideal space) of the test position
    !!< from the test element
    real, optional, intent(out) :: miss
    !!< Return the coordinate (in ideal space) of the test position
    !!< in the test element
    real, dimension(positions_a%dim + 1), optional, intent(out) :: l_coords
    
    logical :: owned
    
    owned = ownership_predicate(positions_a, ele_a, node_val(positions_b, node_b), ownership_tolerance, &
      & miss = miss, l_coords = l_coords)
    
  end function ownership_predicate_node
*/

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
  {"findserial",  pickers_findserial, METH_VARARGS,
   "findserial in Node_Owner_Finder_Fortran.F90."},
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








 

 

