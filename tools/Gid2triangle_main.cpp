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

#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <map>

#include "confdefs.h"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

using namespace std;

extern "C"{
#define gid2triangle F77_FUNC(gid2triangle, GID2TRIANGLE)
  void gid2triangle(const char* input_filename, const int* input_filename_len,
                    const char* output_basename, const int* output_basename_len,
                    const int* reorder_nodes = NULL);
}

void Usage(){
  cout << "Usage: gid2triangle [OPTIONS] INPUT\n"
       << "\n"
       << "Converts a GiD mesh to a triangle mesh. Does not process halo data.\n"
       << "\n"
       << "Options:\n"
       << "\n"
       << "-h  Display this message\n"
       << "-n  Disable GiD to Fluidity node reordering" << endl;
}

int main(int argc, char** argv){
#ifdef HAVE_MPI
  MPI::Init(argc, argv);

  // Undo some MPI init shenanigans
  chdir(getenv("PWD"));
#endif

  // Based on fldecomp.cpp argument handling
  // Get any command line arguments
  // reset optarg so we can detect changes
  optarg = NULL;  
  char c;
  map<string, string> args;
  while((c = getopt(argc, argv, "hn")) != -1){
    if (c != '?'){
      if (optarg == NULL){
        args[string(&c, 1)] = "true";
      }else{
        args[string(&c, 1)] = optarg;
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
  if(argc > optind + 1){
    args["input"] = argv[optind + 1];
  }else if(argc == optind + 1){
    args["input"] = argv[optind];
  }
  args["input"] = string(getenv("PWD")) + string("/") + args["input"];
  if(args.count("h") > 0){
    Usage();
    exit(0);
  }
  
  string inputFilename = args["input"];
  if(inputFilename.size() < 4 or inputFilename.substr(inputFilename.size() - 4, 4) != ".dat"){
    cerr << "Input GiD mesh must have .dat extension" << endl;
    exit(-1);
  }
  string outputBasename = inputFilename.substr(0, inputFilename.size() - 4);
  
  int inputFilename_len = inputFilename.size(), outputBasename_len = outputBasename.size();
  int reorder_nodes_handle = args.count("n") == 0 ? 1 : 0;
  gid2triangle(inputFilename.c_str(), &inputFilename_len, outputBasename.c_str(), &outputBasename_len, &reorder_nodes_handle);

#ifdef HAVE_MPI
  MPI::Finalize();
#endif
  return 0;
}
