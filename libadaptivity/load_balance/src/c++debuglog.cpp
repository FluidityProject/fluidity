/* Copyright (C) 2006 Imperial College London and others.

 Please see the AUTHORS file in the main source directory for a full list
 of copyright holders.

 Dr Gerard J Gorman
 Applied Modelling and Computation Group
 Department of Earth Science and Engineering
 Imperial College London

 g.gorman@imperial.ac.uk

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
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
#include <mpi.h>
#include <string.h>
#include "c++debug.h"

#include <fstream>

using namespace std;

ofstream debugfile;

void openlog(const bool all2null){
  
  int MyRank;
  char filename[20];

  if(all2null){
    strcpy(filename, "/dev/null");
  }else{
    if( MPI::Is_initialized() ){
      MPI_Comm_rank(MPI_COMM_WORLD, &MyRank);
      sprintf(filename, "c++debuglog%06d\n", MyRank);
    }else{
      sprintf(filename, "c++debuglog\n");
    }
  }

  debugfile.open( filename, ofstream::app);
}

void closelog(){
  debugfile.close();
  return;
}
