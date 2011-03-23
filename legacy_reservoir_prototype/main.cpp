/* A cpp main for calling MultiphasePrototype */

#include "Usage.h"
#include "fmangle.h"
//#ifdef HAVE_PYTHON
//#include "Python.h"
//#endif

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <unistd.h>
#include <string.h>

#define F77_FUNC(name,NAME) name ## _
#define F77_FUNC_(name,NAME) name ## _

extern "C" {
#define multiphase_prototype_wrapper_fc F77_FUNC(multiphase_prototype_wrapper, MULTIPHASE_PROTOTYPE_WRAPPER)
  void multiphase_prototype_wrapper_fc(const char *, const int *);
#ifdef USING_GFORTRAN
  /* gfortran hack to ensure 4-byte record marker for unformatted files */
  void _gfortran_set_record_marker(int);
#endif
}


int main(int argc, char **argv){
  
#ifdef USING_GFORTRAN
  /* gfortran hack to ensure 4-byte record marker for unformatted files */
  _gfortran_set_record_marker(4);
#endif
     
  // Get any command line arguments.
  ParseArguments(argc, argv);

  // Start fortran main
  if(fl_command_line_options.count("simulation_name")){
    int filenamelen = fl_command_line_options["simulation_name"].size();
    multiphase_prototype_wrapper_(fl_command_line_options["simulation_name"].c_str(), &filenamelen);    
  }else{
    usage(argv[0]);
    exit(-1);
  }

  exit(0);
}









