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
#include "comTools.h"

using namespace std;

// 
// This file primarily contains post-migration clean-up functions.
//

// Re-arrange element list
void Mesh::fixate_elements(const unsigned paranoid_element, 
			   map<unsigned, Node>& node_cache, 
			   vector< map<unsigned, Node> >& new_halo_nodes){
  
  unsigned ecnt=0;
  for(deque<Element>::iterator it = element_list.begin(); it != element_list.end(); ++it){
    (*it).set_eid(ecnt);
    
    // In the case of newly received elements: We may have previously
    // excluded an old halo node that this new element may need. Thus 
    // we have to check the unn's of this element against the nodes in the cache.
    if(ecnt >= paranoid_element){

      // Check all of the nodes in this element.
      const vector<unn_t>& nodes = (*it).get_enlist();
      for(vector<unn_t>::const_iterator in = nodes.begin(); in != nodes.end(); ++in){
	
	// Check if this node was taged to be deleted.
	map<unn_t, Node>::iterator pos = node_cache.find( *in );
	if( pos != node_cache.end()){
	  
	  unsigned owner = (*pos).second.get_future_owner();
	  
	  assert( new_halo_nodes[ owner ].find( *in ) == 
		  new_halo_nodes[ owner ].end() );
	  ECHO("Found worst case! Using a node from cache (unn = "<< *in << ")");

	  // Put the node into the halo.
	  new_halo_nodes[ owner ][ *in ] = (*pos).second;
	  
	  // Remove the entry from the cache. This will also shorten the search.
	  node_cache.erase( pos );
	  
	}
      }
    }
      
    ecnt++;
  }
  
#ifndef NDEBUG
  // Verify that there is only one instance of each element.
  for(deque<Element>::iterator i=element_list.begin(); i != element_list.end(); ++i){
    deque<Element>::iterator j=i; ++j;
    for(;j!=element_list.end(); ++j)
      if( i != j )
	if( (*i)==(*j) ){
	  ECHO("Two Elements have been found to be equal! See below:");
	  ECHO(*i);
	  ECHO(*j);
	  MPI::COMM_WORLD.Abort(-1);
	} 
  }
#endif

  return;
}

// Clean-up node list
void Mesh::fixate_nodes(vector< map<unsigned, Node> >& new_halo_nodes){
  
  { // Write in the new global numbering for the nodes.
    int ncnt=0;
    
    // Write new id's and ownerships...
    for(deque<Node>::iterator it = node_list.begin(); it != node_list.end(); ++it){
      
      (*it).set_gnn(ncnt++);
      (*it).set_current_owner(MyRank);
  
    }
    
    // Add in halo nodes...
    for(int p=0; p<NProcs; p++){
      for(map<unn_t, Node>::iterator it=new_halo_nodes[p].begin(); it != new_halo_nodes[p].end(); ++it){
	(*it).second.set_gnn(ncnt++);
	(*it).second.set_current_owner(p);
	node_list.push_back( (*it).second );
	assert( (*it).second.get_future_owner() == p);
      }
    }
    
    do_node_headcount();
    
  } // end the building of new node list and halo_nodes
  
  //
  // Build new shared_nodes/halo_nodes vectors.
  //
  {
    ECHO("find new halo/shared nodes ");
    shared_nodes.clear();
    shared_nodes.resize( NProcs );
    halo_nodes.clear();
    halo_nodes.resize(NProcs);

    for(deque<Element>::iterator ie = element_list.begin(); ie != element_list.end(); ++ie){
      
      vector<unn_t> nodes = (*ie).get_enlist();
      set< unsigned > neigh_procs;
      
      // What domains know this element
      for(vector<unn_t>::iterator in = nodes.begin(); in != nodes.end(); ++in)
	neigh_procs.insert( node_list[ unn2gnn(*in) ].get_current_owner() );
      
      assert( neigh_procs.find(MyRank) != neigh_procs.end() );
      
      if( neigh_procs.size() == 1 )
	continue;
      
      neigh_procs.erase( MyRank );

      // Write in the halo/shared nodes
      for(vector<unn_t>::iterator in = nodes.begin(); in != nodes.end(); ++in){
	unsigned owner = node_list[ unn2gnn(*in) ].get_current_owner();
	
	if( owner == (unsigned)MyRank ){ // it's ours to share	
	  for(set<unsigned>::iterator ip = neigh_procs.begin(); ip != neigh_procs.end(); ++ip){
	    shared_nodes[ *ip ].insert( *in );
	  }
	}else{ // we're getting it from...
	  halo_nodes[ owner ].insert( *in );
	}
	
      }
    }
  }
}

// NOTE: The owner of a pressure node is defined to be the minimum node owner
// of all the elements that surround a pressure node.
void Mesh::fixate_pressure(){
  
  // Through local information, find the set of processors that should know each pressure node.
  map<unsigned, set<unsigned> > new_pnodes;
  for(deque<Element>::const_iterator elm=element_list.begin(); elm!=element_list.end(); ++elm){
    
    // Find the set of processors that this element is distributed across.
    const vector<unn_t>& enlist = (*elm).get_enlist();
    set<unsigned> known_by;
    for(vector<unn_t>::const_iterator nod = enlist.begin(); nod != enlist.end(); ++nod){
      const Node& __node__ = node_list.unn(*nod);
      known_by.insert( __node__.get_current_owner() );
    }
    
    // Add this set to the set of processors that know each of the pressure nodes.
    const vector<unn_t>& penlist = (*elm).get_MFenlist();
    for(vector<unn_t>::const_iterator nod = penlist.begin(); nod != penlist.end(); ++nod)
      new_pnodes[ *nod ].insert(known_by.begin(), known_by.end()); 
    
  }
  
  //
  // Enter in negiotions about the ownership of the pressure nodes
  // of elements on the halo.
  //
  
  vector< vector<unsigned> > ownership(NProcs);
  // Go through all pressure nodes that are being kept
  for(map<unsigned, set<unsigned> >::const_iterator it=new_pnodes.begin(); it!=new_pnodes.end(); ++it){
    
    // Negiotions are required if more than one guy may lay claim on a node.
    if((*it).second.size() > 1){
      
      // A unn/minimum-owner pair needs to be sent to all interested parties.
      for(set<unsigned>::const_iterator proc=(*it).second.begin(); proc!=(*it).second.end(); ++proc){
	if(*proc == (unsigned)MyRank) 
	  continue;
	
	ownership[*proc].push_back( (*it).first );
	ownership[*proc].push_back( *( (*it).second.begin() ) );
      }
    }
  }
  
  // send/receive everything
  vector< vector<unsigned> > FinalWord(NProcs);
  allSendRecv(ownership, FinalWord);
  
  // Update owner data
  for(unsigned i=0; i<(unsigned)NProcs; i++)
    for(unsigned j=0; j<FinalWord.size(); j+=2)
      new_pnodes[ FinalWord[i][j] ].insert( FinalWord[i][j+1] );
  
  // Clean-up.
  ownership.clear();
  FinalWord.clear() ;
  
  // Now we know enough to write consistant ownership data and a new MFnode_list.
  deque<PressureNode> halo;
  unsigned end=0;
  
  for(unsigned n=0; n<MFnode_list.size(); n++){
    const unsigned unn = MFnode_list[n].get_unn();
    
    // Check to see if this node is required any longer.
    map<unsigned, set<unsigned> >::const_iterator micky = new_pnodes.find(unn);
    if(micky == new_pnodes.end())
      continue;
    
    // Who is the Minimum-Node-Owner of this pressure node.
    unsigned mno = *( (*micky).second.begin() );
    MFnode_list[n].set_owner(mno);
    
    if(mno==(unsigned)MyRank){
      // Add into list straight away.
      if(end != n){
	MFnode_list[end] = MFnode_list[n];
      }
      end++;
    }else{
      // This node will be passed in at the end.
      halo.push_back( MFnode_list[n] );
    }
  }
  
  // Using this as a short-cut later...
  unsigned private_top = end;
  
  // Append all the halo nodes.
  for(deque<PressureNode>::const_iterator it=halo.begin(); it!=halo.end(); ++it)
    MFnode_list[end++] = *it;
  
  new_pnodes.clear();  
  MFnode_list.shrink(end);
 
  // Form halo
  // std::vector< set<unn> > halo_pnodes;
  halo_pnodes.clear();
  halo_pnodes.resize(NProcs);
  for(unsigned n=private_top; n<MFnode_list.size(); n++)
    halo_pnodes[ MFnode_list[n].get_owner() ].insert(MFnode_list[n].get_unn());
  
  allSendRecv(halo_pnodes, shared_pnodes);
  
  // Finished ;-)

  return;
}









