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
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include "confdefs.h"

#include "c++debug.h"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

using namespace std;

extern "C"{
  void flredecomp(const char* input_basename, size_t input_basename_len, const char* output_basename, size_t output_basename_len,
                  int input_nprocs, int target_nprocs);
}

void Usage(){
  cerr << "Usage: flredecomp [OPTIONS] ... INPUT OUTPUT\n"
       << "\n"
       << "Performs a re-decomposition of a Fluidity checkpoint. Must be run on\n"
       << "max(input no. procs, target no. procs) processors.\n"
       << "\n"
       << "Options:\n"
       << "\n"
       << "-h\t\tDisplay this help\n"
       << "-i\t\tInput number of processors\n"
       << "-l\t\tWrite output to log files\n"
       << "-o\t\tTarget number of processors\n"
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

  // Modified version of fldecomp argument parsing
  // Get any command line arguments
  // Reset optarg so we can detect changes
  optarg = NULL;  
  char c;
  map<char, string> args;
  while((c = getopt(argc, argv, "i:o:hlv")) != -1){
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
  
  // Logging
  if(args.count('l')){
    int rank = 0;
#ifdef HAVE_MPI
    int init_flag;
    MPI_Initialized(&init_flag);
    if(init_flag){
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    }
#endif
    ostringstream buffer;
    buffer << "flredecomp.log-" << rank;
    if(freopen(buffer.str().c_str(), "w", stdout) == NULL){
      cerr << "Failed to redirect stdout" << endl;
      exit(-1);
    }
    buffer.str("");
    buffer << "flredecomp.err-" << rank;
    if(freopen(buffer.str().c_str(), "w", stderr) == NULL){
      cerr << "Failed to redirect stderr" << endl;
      exit(-1);
    }
    buffer.str("");
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
  string input_basename, output_basename;
  if(argc > optind + 2){
    input_basename = argv[optind + 1];
    output_basename = argv[optind + 2];
  }else if(argc == optind + 2){
    input_basename = argv[optind];
    output_basename = argv[optind + 1];
  }else{
    Usage();
    exit(-1);
  }
  
  // Input number of processors
  if(args.count('i') == 0){
    Usage();
    exit(-1);
  }
  int input_nprocs = atoi(args['i'].c_str());
  
  // Target number of processors
  if(args.count('o') == 0){
    Usage();
    exit(-1);
  }
  int target_nprocs = atoi(args['o'].c_str());
    
  size_t input_basename_len = input_basename.size();
  size_t output_basename_len = output_basename.size();
  flredecomp(input_basename.c_str(), input_basename_len, output_basename.c_str(), output_basename_len,
             input_nprocs, target_nprocs);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}
