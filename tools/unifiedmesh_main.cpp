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
  void unifiedmesh(const char* filename1, size_t filename1_len,
                   const char* filename2, size_t filename2_len,
                   const char* output, size_t output_len);
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
  MPI_Init(&argc, &argv);
  // Undo some MPI init shenanigans
  int ierr = chdir(getenv("PWD"));
  if (ierr == -1) {
        cerr << "Unable to switch to directory " << getenv("PWD");
        abort();
  }
#endif
    
  if(argc<3){
    usage();
    return -1;
  }

  string filename1, filename2, outputfilename;
  filename1.append(argv[1]);
  filename2.append(argv[2]);
  outputfilename.append(argv[3]);
  size_t mesh_file_1_len = filename1.size();
  size_t mesh_file_2_len = filename2.size();
  size_t output_file_name_len = outputfilename.size();
  unifiedmesh(filename1.c_str(), mesh_file_1_len,
              filename2.c_str(), mesh_file_2_len,
              outputfilename.c_str(), output_file_name_len);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}
