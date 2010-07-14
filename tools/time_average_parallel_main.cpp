// wrapper for time_average_parallel main
#include "confdefs.h"

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <cstdlib>
#include <iostream>
#include <string.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#ifdef HAVE_PETSC
#include <petsc.h>
#endif

using std::cerr;
using std::endl;
using std::string;


extern "C"{
void time_average_parallel_();
#ifdef HAVE_PYTHON
#include "python_statec.h"
#endif
}

int main(int argc, char **argv){

#ifdef HAVE_MPI
  // This must be called before we process any arguments
  MPI::Init(argc,argv);

  // Undo some MPI init shenanigans
  chdir(getenv("PWD"));
#endif

#ifdef HAVE_PETSC
  static char help[] = "Use -help to see the help.\n\n";
  PetscErrorCode ierr = PetscInitialize(&argc, &argv, NULL, help);
  // PetscInitializeFortran needs to be called 
  // when initialising PETSc from C, but calling it from Fortran
  // This sets all kinds of objects such as PETSC_NULL_OBJECT, 
  // PETSC_COMM_WORLD, etc., etc.
  ierr = PetscInitializeFortran();
  
#ifdef HAVE_PYTHON
  // Initialize the Python Interpreter
  python_init_();
#endif
  
  time_average_parallel_();
  PetscFinalize();
  return 0;
#else
  cerr << "ERROR: Not configured with PETSc, so petsc_readnsolve is not gonna work!" << endl; 
  return 1;
#endif

}
