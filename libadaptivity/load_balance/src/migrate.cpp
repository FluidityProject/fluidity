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

#include "Mesh.h"
#include "MI5.h"
#include "packing.h"
#include "c++debug.h"

using namespace std;

//
// The main migrate routine.
//
void Mesh::migrate(const vector<int> noddom){
  unsigned paranoid_element = 0;

  //
  // write in new node decomposition  
  //
  {
    unsigned ncnt = num_nodes("total");
    assert(ncnt == noddom.size());       // if this assertion is wrong, there is a bug
    
    for(unsigned i = 0; i<ncnt; i++)
      node_list[i].set_future_owner( noddom[i] );
  }
  
  //
  // Based on a new domain decomposition; 
  // + calculate the submesh at remains on this processor
  // + calculate the submeshs that have to be sent to other processors
  //
  MI5 intelligance_report(num_nodes("total"), num_elements("total")); // Initalize.
  intelligance_report.spy( *this ); // Do calculation based on *this mesh
  
  //
  // Using the result of the above calculation;
  // + Pack each of the submeshs that have to be sent, 
  //   and send them to their respective destinations.
  // + Wait and receive any submeshs that may have been 
  //   sent to this processor. The received graphs are appended 
  //   to the origional
  // 
  { // This scope contains packing, sending, receiving and unpacking.
    packing meshpack(intelligance_report, *this);    // Initalize packing.
    meshpack.pack(intelligance_report,    *this);    // Pack submeshes
    
    
    { // In the interest of memory conservation, the node_list is contracted now.
      
      deque<unsigned> free_space;
      unsigned stay_cnt = intelligance_report.new_owned_nodes.size();
      unsigned pos      = 0;

      // We never need free spaces beyond the number of nodes that will be 
      // actually staying.
      for(unsigned i = 0; i<stay_cnt; i++){

        if( i == (unsigned)( intelligance_report.new_owned_nodes[pos] )){
          pos++;
        }else{
          free_space.push_back(i);
        }

      }
     
      unsigned cnt=0;
      unsigned free_cnt = free_space.size();
      for(unsigned i = stay_cnt-1; cnt<free_cnt; i--){
        node_list[ free_space[cnt] ] = node_list[ intelligance_report.new_owned_nodes[i] ];
        cnt++;
      }
      node_list.shrink( stay_cnt );

    }
    
    { // In the interest of memory conservation, the element_list is contracted now.
      
      deque<unsigned> free_space;
      unsigned stay_cnt = intelligance_report.new_elements.size();
      unsigned pos      = 0;
      
      // We never need free spaces beyond the number of elements that will be 
      // actually staying.
      for(unsigned i = 0; i<stay_cnt; i++){

        if(i == (unsigned)( intelligance_report.new_elements[pos] )){
          pos++;
        }else{
          free_space.push_back(i);
        }

      }

      unsigned cnt=0;
      unsigned free_cnt = free_space.size();

      for(unsigned i = stay_cnt-1; cnt<free_cnt; i--){
        element_list[ free_space[cnt] ] = element_list[ intelligance_report.new_elements[i] ];
        cnt++;
     } 
      element_list.shrink( stay_cnt );
      paranoid_element = stay_cnt;
    }
    
    // Now everything is prepared to scommunicate all the data.
    meshpack.send();
    
    // Unpack all the data received
    ECHO("Unpacking submeshes.");
    meshpack.unpack(intelligance_report, *this);
    ECHO("Submeshes unpacked.");

    meshpack.clear();
    
    // At this point, all the submeshes have been redistributed.
  }
  
  //
  // Re-arrange element list
  //
  fixate_elements(paranoid_element, intelligance_report.node_cache, intelligance_report.new_halo_nodes);

  //
  // Clean-up node list
  //
  fixate_nodes(intelligance_report.new_halo_nodes);
  
  //
  // For mixed formulation
  //
#ifdef UNTESTED
  if( mixed_formulation() && (!MFnode_list.empty()) ){
    ECHO("Fixating migrated pressure mesh...");
    // Clean-up pressure-node list
    fixate_pressure();
    
    // Obtain aditional nodes required for halo II
    formHalo2();
    ECHO("...fixated pressure mesh.");
  }
#endif

  ECHO("Finished main migration routine.");
  
  return;
} // end















