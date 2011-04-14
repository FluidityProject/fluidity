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
#include <vector>
#include <deque>
#include <string>
#include <set>
#include <mpi.h>

using std::set;
using std::vector;
using std::deque;
using std::string;
using std::endl;

#include "Mesh.h"
#include "c++debug.h"

void Mesh::calculate_submeshes(vector<Mesh>& submeshes){
  
  unsigned len_node_list = node_list.size();
  unsigned len_elem_list = element_list.size();
  
  // Given a node number, keep a list of partitions that 
  // already know its particulars.
  deque< set<unsigned short> > nodeKnowers( len_node_list );
  
  // Given a node number, who will need to know this node in the 
  // next round.
  deque< set<unsigned short> > nodeKnowers_next( len_node_list );
  
  // Given an element number, keep a list of partitions that 
  // already know its particulars.
  deque< set<unsigned short> > elementKnowers( len_elem_list );
  
  // Given an element number, keep a list of partitions that 
  // will need to know this node in the next round.
  deque< set<unsigned short> > elementKnowers_next( len_elem_list );
  
    // List of nodes to be sent each rank.
  vector< deque<int> > nodes2send( NProcs );

  // List of elements to be sent each rank.
  vector< deque<int> > elems2send( NProcs );

  // List of pressure nodes to be sent each rank.
  vector< deque<int> > pnodes2send( NProcs );

   // List of nodes that we are keeping.
  deque<int> new_owned_nodes;
  
  // New halo nodes...They are being copied to that they
  // don't get overwritten!
  vector< map<unn_t, Node> > new_halo_nodes;
  
  // List of elements that we are keeping.
  deque<int> new_elements;

  // In the case or halo nodes, there are some which may be
  // required, thus they are stored away in a cache.
  map<unn_t, Node> node_cache;
    
  // Consider each element.
  unsigned num_managed_elements = 0;
  unsigned eid = 0;
  for(deque<Element>::iterator ie = element_list.begin(); ie != element_list.end(); ++ie){
    vector<bool> PPresent(NProcs, false);
    vector<bool> PFuture(NProcs, false);
    
    // Consider all the nodes that are in this element, and note the
    //  current and future owners of these nodes.
    vector<unn_t> nodes = (*ie).get_enlist();      
    set< unsigned > doms;
    
    // What domains know this element
    for(vector<unn_t>::iterator in = nodes.begin(); in != nodes.end(); ++in)
      doms.insert( node_list[ unn2gnn(*in) ].get_current_owner() );
    unsigned min_node_owner = *(doms.begin()); 
    
    for(set<unsigned>::iterator it = doms.begin(); it != doms.end(); ++it){
      elementKnownBy_insert(eid, *it);
      PPresent[ *it ] = true;
    }

    if(doms.size() == 2){
      ECHO("halo ");
      for(vector<unn_t>::iterator in = nodes.begin(); in != nodes.end(); ++in)
	ECHO(*in);
    }

    doms.clear();

    // What domains will know this element
    for(vector<unn_t>::iterator in = nodes.begin(); in != nodes.end(); ++in)
      doms.insert( node_list[ unn2gnn(*in) ].get_future_owner() );
    for(set<unsigned>::iterator it = doms.begin(); it != doms.end(); ++it){
      element2bKnownBy_insert(eid, *it);
      PFuture[ *it ] = true;
    }

    // A nice list containing elements
    // which will be used in the next round.
    if( PFuture[MyRank] )
      new_elements.push_back( eid );
    
    if( min_node_owner == MyRank){ // ie. The element is our responcibility.
      num_managed_elements++;
      
      // Another list containing elements that we have to send.
      for(unsigned p = 0; p<NProcs; p++){
        if( (PPresent[p] == false) && (PFuture[p]==true) ){
	  elems2send[ p ].push_back(eid);
        }
      }
    }

    // If a element is to be on a processor, then all the
    // nodes of that element must also be known to that processor.
    vector<unn_t> enl = get_enlist(eid);
    for(vector<unn_t>::iterator it = enl.begin(); it != enl.end(); ++it){
      gnn_t node = unn2gnn( *it );
      
      for(unsigned p=0; p<NProcs; p++){
	if( PPresent[ p ] )
	  nodeKnownBy_insert(node,   p);
	if( PFuture[ p ] )
	  node2bKnownBy_insert(node, p);
      }
    }

    eid++;
  }
  
  
  // Determine the ranks that need to know a node
  // but don't already know it. Also, make a count of the 
  // number of nodes to be sent to each rank.
  vector< set<unsigned> > halo(NProcs);
  
  for(unsigned n=0; n<len_node_list; n++){
    
    vector<bool> PPresent(NProcs, false);
    vector<bool> PFuture( NProcs, false);
    
    // Build the two bit maps.
    for( set<unsigned short>::iterator it = nodeKnowers[n].begin();      it != nodeKnowers[n].end();      ++it)
      PPresent[ *it ] = true;
    for( set<unsigned short>::iterator it = nodeKnowers_next[n].begin(); it != nodeKnowers_next[n].end(); ++it)
      PFuture[ *it ] = true;
    
    unsigned owner = get_node_current_owner(n);
    if( owner == MyRank )     // MyRank is responcible for sending 
      for(unsigned p=0; p<NProcs; p++)
	if( (PPresent[p]==false) && (PFuture[p]==true) )
	  nodes2send[ p ].push_back( n );
    
    // Storing node for future?
    if( PFuture[ MyRank ]==true ){
      unsigned future_owner = get_node_future_owner(n);
      
      if( future_owner == MyRank ){
	// Note its position
	new_owned_nodes.push_back( n );
      }else{
	// Put the new halo nodes somewhere save.
	unn_t u = get_node_unn(n);
        ECHO("MI5:: " << u << " to be known through a halo with " << future_owner);
	new_halo_nodes[ future_owner ][u] = get_node(n);
      }
    }else{
      
      // There is the possibility that some halo nodes will be appear 
      // to be unnecessary in the future when infact they will be again 
      // required as halo nodes, being part of previously unknown elements. 
      // For this reason, these nodes are backed up so that they won't get
      // deleted in the clean-up phase.
      
      if( owner != MyRank ){ // a halo
	node_cache[ get_node_unn(n) ] = get_node(n);
      }
      
    }

  }

  // Fixed formulation.
  if( using_mixed_formulation() && (!MFnode_list.empty()) ){
    pnodes2send.resize(NProcs);
    
    //
    // If an element, e, has to be sent to a process, d, then the
    // pressure nodes must also be sent.
    //
    for(unsigned d=0; d<NProcs; d++){
      
      // Get the set of pressure nodes to be sent to d.
      set<unsigned> pnodes;
      for(deque<int>::const_iterator it = elems2send[d].begin(); it != elems2send[d].end(); ++it){
	const Element& Elem = get_element( *it );
	const vector<unn_t>& pns = Elem.get_MFenlist();
	for(vector<unn_t>::const_iterator p = pns.begin(); p != pns.end(); ++p)
	  pnodes.insert( *p );
      }
      
      // Put set into store.
      for(set<unsigned>::const_iterator pn = pnodes.begin(); pn != pnodes.end(); pn++){
	pnodes2send[d].push_back( *pn );
      }
      
    }
  }
  
  ECHO("All necessary calculations for migration has been collected.");  
}















