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
#include <assert.h>
#include <vector>
#include <deque>
#include <string>
#include <map>
#include <set>

using namespace std;

#include "Mesh.h"
#include "c++debug.h"

// Invent a pressure mesh.
void Mesh::invent_pressure_mesh(){
  ECHO("Inventing a pressure mesh.");
  
  MFnode_list.clear();
        
  // Generate a pressure mesh based on the velocity.
  MFnode_list.expand( node_list.size() );      
  for(unsigned i=0; i<node_list.size(); i++){
    MFnode_list[i].set_unn( node_list[i].get_unn() );
    MFnode_list[i].set_gnn( node_list[i].get_gnn() );
    assert(node_list[i].get_gnn() == i);
    MFnode_list[i].set_flags( node_list[i].get_flags() );
    MFnode_list[i].set_owner( node_list[i].get_current_owner() );
  }
  
  for(ElementVector<Element>::iterator it = element_list.begin(); it != element_list.end(); ++it){
    if((*it).get_flags() & ELM_VOLUME)
      (*it).set_MFenlist( (*it).get_enlist() );
  }
      
  shared_pnodes = shared_nodes;
  halo_pnodes   = halo_nodes;
  
  assert( mesh_consistant() );
  
  return;
} // finished inventing pressure mesh
