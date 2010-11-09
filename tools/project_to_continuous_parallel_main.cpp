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

#include "confdefs.h"
#include "c++debug.h"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#ifdef HAVE_PETSC
#include <petsc.h>
#endif

extern "C" {
#define project_to_continuous_parallel_fc F77_FUNC(project_to_continuous_parallel, PROJECT_TO_CONTINUOUS_PARALLEL_FC)
  void project_to_continuous_parallel_fc(const char *, const int *, const char *, const int *);
}
void PetscInit(int argc, char** argv);

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

void project_to_continuous_parallel_usage(char *binary){
  cerr<<"Usage: "<<binary<<" [OPTIONS] vtufile trianglename.\n"
      <<"\twhere 'vtufile' is the name of the discontinuous vtu (without extension)\n"
      <<"\tand 'trianglename' is the base name of the triangle files.\n"
      <<"\n"
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
      project_to_continuous_parallel_usage(argv[0]);
      exit(-1);
    }
  }

  // Help?
  if(flArgs.count('h')){
    project_to_continuous_parallel_usage(argv[0]);
    exit(-1);
  }
  
  if (argc-no_option_args != 3){
    cerr << "Need exactly two non-option arguments" << endl;
    project_to_continuous_parallel_usage(argv[0]);
    exit(-1);
  }

  // What to do with stdout?
  int val=3;
  if(flArgs.count('v')==0)
    val = 0;
  set_global_debug_level_fc(&val);

  string vtufile;
  if(argv[no_option_args+1][0]!='/')
    vtufile = getenv("PWD");
  vtufile.append("/");
  vtufile.append(argv[no_option_args+1]);
  vtufile.append(" "); // needed by fluidity to find end of file
  int vtulen = vtufile.length();  

  string trianglefile;
  if(argv[no_option_args+1][0]!='/')
    trianglefile = getenv("PWD");
  trianglefile.append("/");
  trianglefile.append(argv[no_option_args+2]);
  trianglefile.append(" "); // needed by fluidity to find end of file
  int trianglelen = trianglefile.length();  

  project_to_continuous_parallel_fc(vtufile.c_str(),&vtulen,trianglefile.c_str(),&trianglelen);
  
#ifdef HAVE_PETSC
  PetscFinalize();
#endif
  
#ifdef HAVE_MPI
  MPI::Finalize();
#endif
  return(0);
}
