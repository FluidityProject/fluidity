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
  void vertical_integration(const char* target_basename, size_t target_basename_len,
                            const char* integrated_filename, size_t integrated_filename_len,
                            const char* output_basename, size_t output_basename_len,
                            double top, double bottom, double sizing, int result_continuity, int result_degree);

#ifdef HAVE_PYTHON
#include "python_statec.h"
#endif
}

void Usage(){
  cerr << "Usage: vertical_intergation [OPTIONS] ... TARGET INTEGRATED OUTPUT\n"
       << "\n"
       << "Performs a vertical integration of a vtu via the use of supermeshing a\n"
       << "vertical extrusion.\n"
       << "\n"
       << "Options:\n"
       << "\n"
       << "-b\t\tBottom\n"
       << "-d\t\tDiscontinuous result (default false for degree > 0, true for degree == 0)\n"
       << "-h\t\tDisplay this help\n"
       << "-p\t\tResult degree (default 1)\n"
       << "-s\t\tSizing\n"
       << "-t\t\tTop (default 0.0)\n"
       << "-v\t\tVerbose mode" << endl;
}

int main(int argc, char** argv){
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  // Undo some MPI init shenanigans
  int ierr = chdir(getenv("PWD"));
  if (ierr == -1) {
        cerr << "Unable to switch to directory " << getenv("PWD");
        abort();
  }
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
  while((c = getopt(argc, argv, "b:dhp:s:t:v")) != -1){
    if (c != '?'){
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
  
  // Options
  double bottom, sizing, top = 0.0;
  if(args.count('b') > 0){
    bottom = atof(args['b'].c_str());
  }else{
    cerr << "Bottom required" << endl;
    Usage();
    exit(-1);
  }
  if(args.count('s') > 0){
    sizing = atof(args['s'].c_str());
  }else{
    cerr << "Sizing required" << endl;
    Usage();
    exit(-1);
  }  
  if(args.count('t') > 0){
    top = atof(args['t'].c_str());
  }
  int result_degree = 1;
  if(args.count('p') > 0){
    result_degree = atoi(args['p'].c_str());
  }
  int result_continuity = args.count('d') > 0 ? -1 : (result_degree == 0 ? -1 : 0);
  
  // Input / output
  string target_basename, integrated_filename, integrated_fieldname, output_basename;
  if(argc > optind + 3){
    target_basename = argv[optind + 1];
    integrated_filename = argv[optind + 2];
    output_basename = argv[optind + 3];
  }else if(argc == optind + 3){
    target_basename = argv[optind];
    integrated_filename = argv[optind + 1];
    output_basename = argv[optind + 2];
  }else{
    Usage();
    exit(-1);
  }
      
  size_t target_basename_len = target_basename.size();
  size_t integrated_filename_len = integrated_filename.size();
  size_t output_basename_len = output_basename.size();
  vertical_integration(target_basename.c_str(), target_basename_len,
                       integrated_filename.c_str(), integrated_filename_len,
                       output_basename.c_str(), output_basename_len,
                       top, bottom, sizing, result_continuity, result_degree);

#ifdef HAVE_PYTHON
  // Finalize the Python Interpreter
  python_end_();
#endif
  
#ifdef HAVE_PETSC
  PetscFinalize();
#endif
  
#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}
