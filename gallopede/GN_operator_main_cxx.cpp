#ifdef HAVE_MPI
#include <mpi.h>
#endif

#ifdef HAVE_PETSC
#include <petsc.h>
#include "../include/c++misc.h"
#endif

#include <map>
#include <iostream>
#include <string>

#include "fmangle.h"

using namespace std;

int main(int argc, char **argv){

  int level = 3;

  PetscErrorCode ierr;
  static char help[] = "No help yet.\n\n";
  ierr = PetscInitialize(&argc,&argv,NULL,help);CHKERRQ(ierr);
  set_global_debug_level_fc(&level);
  gn_operator_main_fc();

  PetscFinalize();
  exit(0);
}
