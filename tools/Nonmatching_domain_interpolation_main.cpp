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
  void nonmatching_domain_interpolation(const char* input_basename_1, size_t input_basename_1_len, const char* input_basename_2, 
                  size_t input_basename_2_len, const char* input_mesh_format, size_t input_mesh_format_len, int ifinder);
}

void Usage(){
  cerr << "Usage: nonmatching_domain_interpolation [OPTIONS] ... mesh1 mesh2 \n"
       << "with 'mesh1' and 'mesh2' being the basenames of the input mesh files\n"
       << "\n"
       << "Options:\n"
       << "\n"
       << "-h\t\tDisplay this help\n"
       << "-m\t\tInput mesh format (both mesh files must be of the same format!)\n"
       << "-f\t\tIntersection finder method to be used.\n"
       << "-l\t\tWrite output to log files\n"
       << "-v\t\tVerbose mode\n"
       << "\n"
       << "With the option 'm' being one of the following:\n"
       << " * triangle\n"
       << " * gmsh\n"
       << " * exodusii,\n"
       << "and option 'f' being one of the following:\n"
       << " * 0 (rtree intersection finder).\n"
       << " * 1 (advancing front for matching domains)\n"
       << " * 2 (advancing front for non-matching domains)\n"
       << " * 3 (advancing front for non-matching domains with brute force search when no clues available)\n"
       << " * 4 (advancing front for non-matching domains with boundary bboxes search).\n"
       << endl;
}

int main(int argc, char** argv){

#ifdef HAVE_MPI
  MPI::Init(argc, argv);
  // Undo some MPI init shenanigans
  chdir(getenv("PWD"));
#endif

  // Modified version of fldecomp argument parsing
  // Get any command line arguments
  // Reset optarg so we can detect changes
  optarg = NULL;
  char c;
  map<char, string> args;
  while((c = getopt(argc, argv, "m:f:hlv")) != -1){
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
    if(MPI::Is_initialized()){
      rank = MPI::COMM_WORLD.Get_rank();
    }
#endif

    ostringstream buffer;
    buffer << "nonmatching_domain_interpolation.log-" << rank;
    if(freopen(buffer.str().c_str(), "w", stdout) == NULL){
      cerr << "Failed to redirect stdout" << endl;
      exit(-1);
    }
    buffer.str("");
    buffer << "nonmatching_domain_interpolation.err-" << rank;
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

  // Input base names
  string input_basename_1, input_basename_2;
  if(argc == optind + 2){
    input_basename_1 = argv[optind];
    input_basename_2 = argv[optind + 1];
  }else if(argc == optind + 3){
    input_basename_1 = argv[optind + 1];
    input_basename_2 = argv[optind + 2];
  }else{
    Usage();
    exit(-1);
  }


  // Input mesh format:
  string input_mesh_format;
  if(args.count('m') == 0){
    Usage();
    exit(-1);
  }
  input_mesh_format = args['m'].c_str();

  // Intersection finder algorithm to be used:
  int ifinder = 0;
  if(args.count('f') == 0){
    Usage();
    exit(-1);
  }
  ifinder = atoi(args['f'].c_str());
  if (ifinder < 0 or ifinder >=5){
    Usage();
    exit(-1);
  }


  size_t input_basename_1_len = input_basename_1.size();
  size_t input_basename_2_len = input_basename_2.size();
  size_t input_mesh_format_len = input_mesh_format.size();
  nonmatching_domain_interpolation(input_basename_1.c_str(), input_basename_1_len, input_basename_2.c_str(), input_basename_2_len, input_mesh_format.c_str(), input_mesh_format_len, ifinder);


#ifdef HAVE_MPI
  MPI::Finalize();
#endif

  return 0;
}
