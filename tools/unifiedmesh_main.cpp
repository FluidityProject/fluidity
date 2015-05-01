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
#include <unistd.h>

#include "confdefs.h"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

using namespace std;

extern "C"{
#define unifiedmesh F77_FUNC(unifiedmesh, UNIFIEDMESH)
  void unifiedmesh_(const char* filename1, const int* filename1_len,
                   const char* filename2, const int* filename2_len,
                   const char* output, const int* output_len);
}

void usage(){
  cerr<<"usage:\n"
      <<"unifiedmesh <mesh-file-1> <mesh-file-2> <output-file-name>\n"
      <<"Dumps the supermesh constructed from the two input meshes\n"
      <<"to discontinuous gmsh mesh files and a vtu file,\n"
      <<"both called output-file-name.\n"
      <<"All file names should not have extensions."<<endl;
}

int main(int argc, char** argv){
#ifdef HAVE_MPI
  MPI::Init(argc, argv);
  // Undo some MPI init shenanigans
  chdir(getenv("PWD"));
#endif
    
  if(argc<3){
    usage();
    return -1;
  }

  int mesh_file_1_len=strlen(argv[1]);
  int mesh_file_2_len=strlen(argv[2]);
  int output_file_name_len=strlen(argv[3]);
  unifiedmesh_(argv[1], &mesh_file_1_len,
              argv[2], &mesh_file_2_len,
              argv[3], &output_file_name_len);

#ifdef HAVE_MPI
  MPI::Finalize();
#endif
  return 0;
}
