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

#include "Usage.h"
#include "Precision.h"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

extern "C" {
  void vtu_bins(const char*, size_t, const char*, size_t, double*, size_t);
}

#include <unistd.h>

#ifndef _AIX
#include <getopt.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include <map>
#include <string>
#include <iostream>

using namespace std; 

void vtu_bins_usage(char *binary){
  cerr<<"Usage: "<<binary<<" [OPTIONS] input_filename input_fieldname BOUND1 [BOUND2 BOUND3 ...]\n"
      <<"Return the volumes between field values.\n"
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
      vtu_bins_usage(argv[0]);
      exit(-1);
    }
  }

  // Help?
  if(args.count('h')){
    vtu_bins_usage(argv[0]);
    exit(-1);
  }
  
  if(optind > argc - 3){
    cerr << "Need an input vtu, input field name and at least one boundary value" << endl;
    vtu_bins_usage(argv[0]);
    exit(-1);
  }

  // What to do with stdout?
  int val=3;
  if(args.count('v')==0)
    val = 0;
  set_global_debug_level_fc(&val);

  string input_filename = argv[optind];
  size_t input_filename_len = input_filename.size();  
  
  string input_fieldname;
  input_fieldname = argv[optind + 1];  
  int input_fieldname_len = input_fieldname.size();
  
  size_t nbounds =  argc - optind - 2;
  double bounds[nbounds];
  for(int i = 0;i < nbounds;i++){
    bounds[i] = atof(argv[optind + 2 + i]);
  }

  vtu_bins(input_filename.c_str(), input_filename_len, input_fieldname.c_str(), input_fieldname_len, bounds, nbounds);
    
#ifdef HAVE_PETSC
  PetscFinalize();
#endif

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return(0);
}
