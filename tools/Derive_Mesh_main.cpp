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

#ifdef HAVE_MPI
#include <mpi.h>
#endif

extern "C" {
#define derive_mesh_fc F77_FUNC(derive_mesh, DERIVE_MESH)
  void derive_mesh_fc(const char*, int*, const char*, int*, int*, int*);
  void set_global_debug_level_fc(int* val);
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

void usage(char *binary){
  cerr<<"Usage: "<<binary<<" [OPTIONS] input_trianglename output_trianglename\n"
      <<"\t-h\t\tPrints out this message\n"
      <<"\t-d\t\tDerive a discontinuous mesh\n"
      <<"\t-p DEGREE\tDerive a mesh of degree DEGREE\n"
      <<"\t-v\t\tVerbose mode\n";
}

int main(int argc, char **argv){
 #ifdef HAVE_MPI
  // This must be called before we process any arguments
  MPI::Init(argc,argv);

  // Undo some MPI init shenanigans
  chdir(getenv("PWD"));
#endif
 
  // Get any command line arguments
  // reset optarg so we can detect changes
  optarg = NULL;  
  char c;
  map<char, string> args;
  while((c = getopt(argc, argv, "dhp:v")) != -1){
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
      usage(argv[0]);
      exit(-1);
    }
  }
  
  if (optind != argc - 2){
    cerr << "Need exactly two non-option arguments" << endl;
    usage(argv[0]);
    exit(-1);
  }

  // Help?
  if(args.count('h')){
    usage(argv[0]);
    exit(-1);
  }

  // What to do with stdout?
  int val=3;
  if(args.count('v')==0)
    val = 0;
  set_global_debug_level_fc(&val);

  int degree;
  if(args.count('p') > 0){
    degree = atoi(args['p'].c_str());
  }else{
    cerr << "Mesh degree not specified" << endl;
    exit(-1);
  }
  
  int continuity = 0;
  if(args.count('d') == 1){
    continuity = 1;
  }

  string input_trianglefile = argv[optind];
  int input_trianglelen = input_trianglefile.length();  
  
  string output_trianglefile = argv[optind + 1];
  int output_trianglelen = output_trianglefile.length();  

  derive_mesh_fc(input_trianglefile.c_str(), &input_trianglelen, output_trianglefile.c_str(), &output_trianglelen, &degree, &continuity);
  
#ifdef HAVE_MPI
  MPI::Finalize();
#endif
  return(0);
}
