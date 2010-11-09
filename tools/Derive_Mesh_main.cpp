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
#include "fmangle.h"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

extern "C" {
#define derive_mesh_fc F77_FUNC(derive_mesh, DERIVE_MESH)
  void derive_mesh_fc(const char*, int*, const char*, int*, int*, int*, int*);
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
  cerr<<"Usage: "<<binary<<" [OPTIONS] input_filename output_filename\n"
      <<"Derive a mesh from a supplied triangle mesh file\n"
      <<"\t-h\t\tPrints out this message\n"
      <<"\t-d\t\tDerive a discontinuous mesh\n"
      <<"\t-k\t\tOutput to vtu rather than a triangle file\n"
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
  while((c = getopt(argc, argv, "dhkp:v")) != -1){
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

  // Help?
  if(args.count('h')){
    usage(argv[0]);
    exit(-1);
  }
  
  if (optind != argc - 2){
    cerr << "Need exactly two non-option arguments" << endl;
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
  
  int vtu = 0;
  if(args.count('k') > 0){
    vtu = 1;
  }
  
  int continuity = 0;
  if(args.count('d') == 1){
    continuity = -1;
  }

  string input_filename = argv[optind];
  int input_filename_len = input_filename.length();  
  
  string output_filename = argv[optind + 1];
  int output_filename_len = output_filename.length();  

  derive_mesh_fc(input_filename.c_str(), &input_filename_len, output_filename.c_str(), &output_filename_len, &degree, &continuity, &vtu);
  
#ifdef HAVE_MPI
  MPI::Finalize();
#endif
  return(0);
}
