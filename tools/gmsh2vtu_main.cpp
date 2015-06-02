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
  void gmsh2vtu(const char* filename, size_t filename_len);
}

void usage(){
  cerr<<"usage: gmsh2vtu <gmsh_file_name>\n"
      <<"The gmsh file name should be without the .msh suffix."<<endl;
}

int main(int argc, char** argv){
#ifdef HAVE_MPI
  MPI::Init(argc, argv);
  // Undo some MPI init shenanigans
  chdir(getenv("PWD"));
#endif
    
  if(argc<2){
    usage();
    return -1;
  }

  string filename;
  filename.append(argv[1]);
  size_t filename_len= filename.size();
  gmsh2vtu(filename.c_str(), filename_len);

#ifdef HAVE_MPI
  MPI::Finalize();
#endif
  return 0;
}
