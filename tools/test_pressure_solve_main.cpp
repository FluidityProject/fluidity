// wrapper for test_pressure_solve main
#include "confdefs.h"

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <cstdlib>
#include <iostream>
using std::cerr;
using std::endl;

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#ifdef HAVE_PETSC
#include <petsc.h>
#endif

extern "C"{
#ifdef HAVE_PYTHON
#include "python_statec.h"
#endif
void test_pressure_solve_();
}

int main(int argc, char **argv){

#ifdef HAVE_MPI
  // This must be called before we process any arguments
  MPI_Init(&argc, &argv);

  // Undo some MPI init shenanigans
  int cderr = chdir(getenv("PWD"));
  if (cderr == -1) {
        cerr << "Unable to switch to directory " << getenv("PWD");
        abort();
  }
#endif

#ifdef HAVE_PYTHON
  // Initialize the Python Interpreter
  python_init_();
#endif

#ifdef HAVE_PETSC
  static char help[] = "Use -help to see the help.\n\n";
  PetscErrorCode ierr = PetscInitialize(&argc, &argv, NULL, help);
  // PetscInitializeFortran needs to be called when initialising PETSc from C, but calling it from Fortran
  ierr = PetscInitializeFortran();
  
  test_pressure_solve_();
  PetscFinalize();
#ifdef HAVE_PYTHON
  // Finalize the Python Interpreter
  python_end_();
#endif
  return 0;
#else
  cerr << "ERROR: Not configured with PETSc, so test_pressure_solve is not gonna work!" << endl; 
  return 1;
#endif

}
