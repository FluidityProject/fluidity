// wrapper for petsc_readnsolve main
#include "confdefs.h"

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <cstdlib>
#include <iostream>
#include <string.h>
#include <sstream>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#ifdef HAVE_PETSC
#include <petsc.h>
#endif

using std::cerr;
using std::endl;
using std::string;
using std::ostringstream;


extern "C"{
void petsc_readnsolve_();
#ifdef HAVE_PYTHON
#include "python_statec.h"
#endif
}

void usage(int argc, char **argv){
// Users are allowed to specify the flml and field name as 2 consecutive
// *non-optional* arguments. As PETSc doesn't really provide a way to
// access those, and we can't access the command-line from fortran
// We read those here and stick them in the PETSc options database to be
// read from fortran
  
  char flml_extension[]=".flml";
  char *flml_file=NULL;
  PetscErrorCode ierr;
  PetscTruth     flg;
  
  // if it's already specified as a PETSc option, we do nothing:
  ierr = PetscOptionsHasName("prns_","-flml",&flg);
  if (flg) {
    return;
  }
  // search for any argument ending in .flml
  int i;
  for (i=0; i<argc; i++) {
    int l=strlen(argv[i]);
    if( (l>5) && !strcmp( &(argv[i][l-5]), flml_extension)) {
      flml_file=argv[i];
      break;
    }
  }
  if (flml_file) {
    string my_PETSc_options="-prns_flml ";
    my_PETSc_options+= flml_file;
    // see if next argument is a valid fieldname
    // but only if not already in the PETSc options database
    ierr = PetscOptionsHasName("prns_","-field",&flg);
    if( !flg && (i+1<argc) && (argv[i+1][0]!='-') ) {
      my_PETSc_options+= " -prns_field " + string(argv[i+1]);
    }
    ierr = PetscOptionsInsertString( my_PETSc_options.c_str() );
  }
  
  // -l option needs to be dealt with in c++ already
  ierr = PetscOptionsHasName("","-l",&flg);
  if (flg) {
    int rank = 0;
#ifdef HAVE_MPI
    if(MPI::Is_initialized()){
      rank = MPI::COMM_WORLD.Get_rank();
    }
#endif
    ostringstream buffer;
    buffer << "petsc_readnsolve.log-" << rank;
    if(freopen(buffer.str().c_str(), "w", stdout) == NULL){
      cerr << "Failed to redirect stdout" << endl;
      exit(-1);
    }
    buffer.str("");
    buffer << "petsc_readnsolve.err-" << rank;
    if(freopen(buffer.str().c_str(), "w", stderr) == NULL){
      cerr << "Failed to redirect stderr" << endl;
      exit(-1);
    }
    buffer.str("");
  }
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
  // PetscInitializeFortran needs to be called when initialising PETSc from C, but calling it from Fortran
  // This sets all kinds of objects such as PETSC_NULL_OBJECT, PETSC_COMM_WORLD, etc., etc.
  ierr = PetscInitializeFortran();
  
  usage(argc, argv);
  
#ifdef HAVE_PYTHON
  // Initialize the Python Interpreter
  python_init_();
#endif
  
  petsc_readnsolve_();
  PetscFinalize();
  return 0;
#else
  cerr << "ERROR: Not configured with PETSc, so petsc_readnsolve is not gonna work!" << endl; 
  return 1;
#endif

}
