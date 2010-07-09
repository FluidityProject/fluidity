/*  Copyright (C) 2006 Imperial College London and others.
    
    Please see the AUTHORS file in the main source directory for a full list
    of copyright holders.

    Prof. C Pain
    Applied Modelling and Computation Group
    Department of Earth Science and Engineering
    Imperial College London

    C.Pain@Imperial.ac.uk
    
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
#include "fmangle.h"
#include "Usage.h"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

extern "C" {
#define project_vtu_fc F77_FUNC(project_vtu, PROJECT_VTU)
  void project_vtu_fc(const char*, int*, const char*, int*, const char*, int*, const char*, int*);
}

#ifdef _AIX
#include <unistd.h>
#else
#include <getopt.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include <map>
#include <string>
#include <iostream>

using namespace std; 

void project_vtu_usage(char *binary){
  cerr<<"Usage: "<<binary<<" [OPTIONS] input_filename donor_basename target_basename output_filename\n"
      <<"Project an input vtu onto a triangle mesh\n"
      <<"\t-h\t\tPrints out this message\n"
      <<"\t-v\t\tVerbose mode\n";
}

int main(int argc, char **argv){
 #ifdef HAVE_MPI
  // This must be called before we process any arguments
  MPI::Init(argc,argv);

  // Undo some MPI init shenanigans
  chdir(getenv("PWD"));
#endif
  
  // Initialise PETSc (this also parses PETSc command line arguments)
  PetscInit(argc, argv);
 
  // Get any command line arguments
  // reset optarg so we can detect changes
  optarg = NULL;  
  char c;
  map<char, string> args;
  while((c = getopt(argc, argv, "hv")) != -1){
    if (c != '?'){
      if (optarg == NULL){
        args[c] = "true";
      }else{
        args[c] = optarg;
      };
    }else{
      if (isprint(optopt)){
        cerr << "Unknown option " << optopt << endl;
      }else{
        cerr << "Unknown option " << hex << optopt << endl;
      }
      project_vtu_usage(argv[0]);
      exit(-1);
    }
  }

  // Help?
  if(args.count('h')){
    project_vtu_usage(argv[0]);
    exit(-1);
  }
  
  if (optind != argc - 4){
    cerr << "Need exactly four non-option arguments" << endl;
    project_vtu_usage(argv[0]);
    exit(-1);
  }

  // What to do with stdout?
  int val=3;
  if(args.count('v')==0)
    val = 0;
  set_global_debug_level_fc(&val);

  string input_filename = argv[optind];
  int input_filename_len = input_filename.length();  
  
  string donor_basename = argv[optind + 1];
  int donor_basename_len = donor_basename.length();
  
  string target_basename = argv[optind + 2];
  int target_basename_len = target_basename.length();
  
  string output_filename = argv[optind + 3];
  int output_filename_len = output_filename.length();  

  project_vtu_fc(input_filename.c_str(), &input_filename_len, donor_basename.c_str(), &donor_basename_len, target_basename.c_str(), &target_basename_len, output_filename.c_str(), &output_filename_len);
    
#ifdef HAVE_PETSC
  PetscFinalize();
#endif

#ifdef HAVE_MPI
  MPI::Finalize();
#endif
  return(0);
}
