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

#include <stdlib.h>

extern "C" {
  void project_to_continuous(const char *, size_t, const char *, size_t);
}

#include <unistd.h>

#ifndef _AIX
#include <getopt.h>
#endif

#include <stdio.h>
#include <errno.h>

#include <map>
#include <string>
#include <iostream>
#include "c++debug.h"
using namespace std; 

void usage(char *binary){
  cerr<<"Usage: "<<binary<<" [OPTIONS] vtufile meshname.\n"
      <<"\twhere 'vtufile' is the name of the discontinuous vtu\n"
      <<"\tand 'meshname' is the base name of the mesh file.\n"
      <<"\n"
      <<"\t-h\t\tPrints out this message\n"
      <<"\t-v\t\tVerbose mode\n";
}

int main(int argc, char **argv){
 #ifdef HAVE_MPI
  // This must be called before we process any arguments
  MPI_Init(&argc, &argv);

  // Undo some MPI init shenanigans
  int ierr = chdir(getenv("PWD"));
  if (ierr == -1) {
        cerr << "Unable to switch to directory " << getenv("PWD");
        abort();
  }
#endif
 
  // Get any command line arguments
  // reset optarg so we can detect changes
  optarg = NULL;  
  char c;
  int no_option_args=0;
  map<char, string> flArgs;
  while((c = getopt(argc, argv, "f:hn:v")) != -1){
    if (c != '?'){
      if (optarg == NULL){
        flArgs[c] = "true";
      }else{
        flArgs[c] = optarg;
      };
      no_option_args++;
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
  if(flArgs.count('h')){
    usage(argv[0]);
    exit(-1);
  }
  
  if (argc-no_option_args != 3){
    cerr << "Need exactly two non-option arguments" << endl;
    usage(argv[0]);
    exit(-1);
  }

  // What to do with stdout?
  int val=3;
  if(flArgs.count('v')==0)
    val = 0;
  set_global_debug_level_fc(&val);

  string vtufile;
  vtufile.append(argv[no_option_args+1]);
  size_t vtulen = vtufile.size();

  string meshfile;
  meshfile.append(argv[no_option_args+2]);
  size_t meshlen = meshfile.size();

  project_to_continuous(vtufile.c_str(),vtulen,meshfile.c_str(),meshlen);
  
#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return(0);
}
