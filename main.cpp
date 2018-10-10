/*  Copyright (C) 2006 Imperial College London and others.
    
    Please see the AUTHORS file in the main source directory for a full list
    of copyright holders.

    Prof. C Pain
    Applied Modelling and Computation Group
    Department of Earth Science and Engineering
    Imperial College London

    amcgsoftware@imperial.ac.uk
    
    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation,
    version 2.1 of the License.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
    USA
*/

#include "Usage.h"
#include "fmangle.h"
//#ifdef HAVE_PYTHON
//#include "Python.h"
//#endif

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <unistd.h>

extern "C" {
#ifdef HAVE_PYTHON
#include "python_statec.h"
#endif
  void mainfl();

/* we probably can get rid of this hack nowadays, cos the default in gfortran
(since 4.2, I think) is to use four-byte record markers */
#ifdef USING_GFORTRAN
  /* gfortran hack to ensure 4-byte record marker for unformatted files */
  void _gfortran_set_record_marker(int);
#endif
}

#include "Profiler.h"

int main(int argc, char **argv){
  
#ifdef HAVE_MPI
  // This must be called before we process any arguments
  MPI_Init(&argc,&argv);

  // Undo some MPI init shenanigans
  int ierr = chdir(getenv("PWD"));
  if (ierr == -1) {
	std::cerr << "Unable to switch to directory " << getenv("PWD");
	abort();
  }
  
#endif

  flprofiler.tic("/fluidity");
  
#ifdef USING_GFORTRAN
  /* gfortran hack to ensure 4-byte record marker for unformatted files */
  _gfortran_set_record_marker(4);
#endif

  // Get any command line arguments.
  ParseArguments(argc, argv);

  if(atoi(fl_command_line_options["verbose"].c_str()) >= 2){
    print_version(std::cout);
  }    
  
  // Initialise PETSc (this also parses PETSc command line arguments)
  PetscInit(argc, argv);

#ifdef HAVE_PYTHON
  // Initialize the Python Interpreter
  python_init_();
#endif

  // Start fortran main
  if(fl_command_line_options.count("simulation_name")){
    mainfl();    
  }else{
    usage(argv[0]);
    exit(-1);
  }

#ifdef HAVE_PYTHON
  // Finalize the Python Interpreter
  python_end_();
#endif
  
#ifdef HAVE_PETSC
  PetscFinalize();
#endif
  
  flprofiler.toc("/fluidity");
  if ( fl_command_line_options.count("profile") ) {
    flprofiler.print();
  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return(0);
}
