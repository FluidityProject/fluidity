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

#include <cstring>
#include <iostream>
#include <stdlib.h>

#include "confdefs.h"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

using namespace std;

extern "C"{
#define vtu2triangle F77_FUNC(vtu2triangle, TRIANGLE2VTU)
  void vtu2triangle(const char* filename, const int* filename_len);
}

void usage(){
  cerr<<"usage: vtu2triangle <vtu_file_name>\n"
      <<"The vtu file name should be without the .vtu suffix."<<endl;
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

  int filename_len=strlen(argv[1]);
  vtu2triangle(argv[1], &filename_len);

#ifdef HAVE_MPI
  MPI::Finalize();
#endif
  return 0;
}
