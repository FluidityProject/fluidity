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

#include <cstring>
#include <iostream>
#include <stdlib.h>

#include "confdefs.h"

#include "c++debug.h"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

using namespace std;

extern "C"{
#define pseudo_mesh F77_FUNC(pseudo_mesh, PSEUDO_MESH)
  void pseudo_mesh(const char* filename, const int* filename_len, const int* target_elements);
}

void Usage(){
  cout << "Usage: pseudo_mesh TRIANGLE_BASENAME\n"
       << "\n"
       << "Generates a new mesh with approximately the same local nodal density as the\n"
       << "input mesh\n"
       << "\n"
       << "Options:\n"
       << "\n"
       << "-h\t\tDisplay this help\n"
       << "-t\t\tLimit the adaptivity metric to target the same number of\n"
       <<   "\t\telements in the pseudo mesh as in the target mesh\n"
       << "-v\t\tVerbose mode" << endl;
}

int main(int argc, char** argv){
#ifdef HAVE_MPI
  MPI::Init(argc, argv);
  // Undo some MPI init shenanigans
  chdir(getenv("PWD"));
#endif
    
  if(argc < 2){
    Usage();
    return -1;
  }

  // Modified version of flredecomp argument parsing
  // Get any command line arguments
  // Reset optarg so we can detect changes
  optarg = NULL;  
  char c;
  map<char, string> args;
  while((c = getopt(argc, argv, "htv")) != -1){
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
  
  if (optind != argc - 1){
    cerr << "Need exactly one non-option argument" << endl;
    Usage();
    exit(-1);
  }

  // Verbosity
  int verbosity = 0;
  if(args.count('v') > 0){
    verbosity = 3;
  }
  set_global_debug_level_fc(&verbosity);
  
  int target_elements;
  target_elements = args.count('t') > 0 ? 1 : 0;

  int filename_len = strlen(argv[optind]);
  pseudo_mesh(argv[optind], &filename_len, &target_elements);

#ifdef HAVE_MPI
  MPI::Finalize();
#endif
  return 0;
}
