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
#include <string>
#include <deque>
#include <vector>
#include <cassert>

#include <string.h>

#include "Mesh.h"
#include "giovtk.h"
#include "c++debug.h"

using namespace std;

// Currently this only supports VTK_TRIANGLE and VTK_TETRA
int Mesh::writeVTK(const string vtkfilename){
  vtkug vtkfile;
  string title = "SAM VTK dump";

  vtkfile.openw(vtkfilename, title, false);
  unsigned nlen = node_list.size();
  unsigned elen = element_list.size();

  // Write out points
  for(unsigned i=0; i<nlen; i++){
    float x = node_list[i].get_x();
    float y = node_list[i].get_y();
    float z = node_list[i].get_z();
    vtkfile.putPoint(x,y,z); 
  }
 
  // Write out elements
  for(unsigned i=0; i<elen; i++){
    unsigned char flag = element_list[i].get_flags();
    int type;
    if(flag&ELM_SURFACE){
      type = VTK_TRIANGLE;
    }else if(flag&ELM_VOLUME){
      type = VTK_TETRA;
    }else{
      ERROR("unreconised element type in writeVTK()");
      //      exit(-1);
    }
    vtkfile.putCell(element_list[i].get_enlist(), type);
  }

  // Write out the fields
  for(unsigned i=0; i<nlen; i++){

    const vector<samfloat_t>& flds = node_list[i].get_fields();
    const unsigned flen = flds.size();
    for(unsigned j=0; j<flen; j++){
      char str[80];
      sprintf(str, "field_%d\n", j);
      string tag(str, strlen( str ));

      vtkfile.putScalarPointData(tag, (float)flds[j]);
    }

  }

  // Finish
  vtkfile.close();

  return(0);
}
