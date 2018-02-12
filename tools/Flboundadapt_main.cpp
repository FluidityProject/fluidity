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

#include <iostream>
#include <map>
#include <stdlib.h>

#include "confdefs.h"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include "Usage.h"
#include "c++debug.h"

using namespace std;

extern "C"{
  void flboundadapt(const char *, const char *, const char *, const char *);
#ifdef HAVE_PYTHON
#include "python_statec.h"
#endif
}

void Usage(){
  cerr << "Usage: flboundadapt INPUT_FLML INPUT_GEO_FILE OUTPUT_MESH\n"
       << "or\n"
       << "flboundadapt -i INPUT_VTU METRIC_NAME INPUT_GEO_FILE OUTPUT_MESH\n"
       << "\n"
       << "Performs mesh generation based on the input options file (which may be a\n"
       << "checkpoint). Outputs the resulting mesh.\n"
       << "\n"
       << "Options:\n"
       << "\n"
       << "-h\t\tDisplay this help\n"
       << "-i\t\t read from VTU file\n"
       << "-v\t\tVerbose mode" << endl;
}

int main(int argc, char** argv){
#ifdef HAVE_MPI
  MPI::Init(argc, argv);
  // Undo some MPI init shenanigans
  chdir(getenv("PWD"));
#endif
  PetscInit(argc, argv);
#ifdef HAVE_PYTHON
  // Initialize the Python Interpreter
  python_init_();
#endif

  // Modified version of flredecomp argument parsing
  // Get any command line arguments
  // Reset optarg so we can detect changes
  optarg = NULL;  
  char c;
  map<char, string> args;
  bool vtk_mode = false;
  while((c = getopt(argc, argv, "hiv")) != -1){
    if (c != '?'){
      if (c == 'i') vtk_mode=true;
      if(optarg == NULL){
        args[c] = "true";
      }else{
        args[c] = optarg;
      }
    }else{
      if(isprint(optopt)){
        cerr << "Unknown option " << optopt << endl;
      }else{
        cerr << "Unknown option " << hex << optopt << endl;
      }
      Usage();
      exit(-1);
    }
  }

  // Help
  if(args.count('h')){
    Usage();
    exit(0);
  }
  
  // Verbosity
  int verbosity = 0;
  if(args.count('v') > 0){
    verbosity = 3;
  }
  set_global_debug_level_fc(&verbosity);
  
  // Input and output base names
  int I = vtk_mode?1:0;
  string input_name, input_geometryname, output_meshname, metric_name="";
  if(argc > optind + 3 +I){
    input_name = argv[optind + 1];
    if (vtk_mode) metric_name = argv[optind + 2];
    input_geometryname = argv[optind + 2 + I];
    output_meshname = argv[optind + 3 + I];
  }else if(argc == optind + 3+(vtk_mode?1:0) ){
    input_name = argv[optind];
    if (vtk_mode) metric_name = argv[optind + 1];
    input_geometryname = argv[optind + 1 +I];
    output_meshname = argv[optind + 2 + I];
  }else{
    Usage();
    exit(-1);
  }
      


  flboundadapt(input_name.c_str(),
	       input_geometryname.c_str(),
	       output_meshname.c_str(),
	       metric_name.c_str());

#ifdef HAVE_PYTHON
  // Finalize the Python Interpreter
  python_end_();
#endif
  
#ifdef HAVE_PETSC
  PetscFinalize();
#endif
  
#ifdef HAVE_MPI
  MPI::Finalize();
#endif
  return 0;
}
