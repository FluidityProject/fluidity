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
#include "packing.h"
#include "c++debug.h"
#include "comTools.h"

//#define CREATE_VELOCITY_HALO2 1

void Mesh::formHalo2(){
#ifndef NDEBUG
  unsigned len45 = element_list.size();
  for(unsigned i=0; i<len45; i++)
    for(unsigned j=0; j<len45; j++)
      if(i!=j)
	assert(element_list[i] != element_list[j]);
  #endif

  // This is what we are looking for.
  vector< deque<unsigned> > halo2Elements(NProcs);
  set<unsigned> dangerElements;
  {  
    // Wizz through all elements.
    unsigned pos = 0;
    for(deque<Element>::const_iterator elm = element_list.begin(); elm != element_list.end(); ++elm, ++pos){

      const vector<unn_t>& nodes = (*elm).get_enlist();
      vector<bool> halo1(NProcs, false);
      bool interesting=false;

      for(int p=0; p<NProcs; p++){
	if(p==MyRank) continue;
	
	// Record if elm has a halo node with p.
	for( vector<unn_t>::const_iterator jt = nodes.begin(); jt!=nodes.end(); ++jt){
	  halo1[p] = (halo_nodes[p].find(*jt)!= halo_nodes[p].end());
	  
	  // Early exit...
	  if(halo1[p])
	    break;
	}
	
	if(!halo1[p]){
	  // Record if elm shares a node with p, but is not itself a halo1 element for p.
	  for(vector<unn_t>::const_iterator jt=nodes.begin(); jt!=nodes.end(); ++jt){
	    if(shared_nodes[p].find(*jt)!=shared_nodes[p].end()){
	      interesting = true;
	      halo2Elements[p].push_back(pos);
	      break;
	    }
	  }
	}

      }
      
      if(interesting){
	// A halo2 element may be multiably sent if it contains halo1
	// nodes for any domain.
	for(unsigned p=0;p<(unsigned)NProcs;p++){
	  if(halo1[p]){
	    dangerElements.insert( pos );
	    break;
	  }
	}	  
      }
    }
  } // Have halo2 elements and a set of dangerElements.
  

  //
  // Identify all the nodes that must be sent.
  //
  vector< set<unsigned> > sendhalo2nodes(NProcs);
  vector< set<unsigned> > sendhalo2pnodes(NProcs);
  for(unsigned p=0; p<(unsigned)NProcs; p++){

    for(deque<unsigned>::const_iterator elm=halo2Elements[p].begin(); elm!=halo2Elements[p].end(); ++elm){
      // It's enough to just check the volume elements
      if(element_list[*elm].get_flags() & ELM_SURFACE)
        continue;
      
      {// Regular nodes	
	const vector<unn_t>& nodes = element_list[*elm].get_enlist();
	for(vector<unn_t>::const_iterator nod=nodes.begin(); nod!=nodes.end(); ++nod){
	  sendhalo2nodes[p].insert( *nod );
	}
      }
      {// Pressure nodes
	const vector<unn_t>& nodes = element_list[*elm].get_MFenlist();
	for(vector<unn_t>::const_iterator nod=nodes.begin(); nod!=nodes.end(); ++nod){
	  sendhalo2pnodes[p].insert( *nod );
	}
      }
    }
    
    {// Remove nodes that p should already know through halo1.
      set<unsigned> toDel;
      for(set<unsigned>::const_iterator it=sendhalo2nodes[p].begin(); it != sendhalo2nodes[p].end(); ++it){
	if(shared_nodes[p].find( *it ) != shared_nodes[p].end()){
	  toDel.insert( *it );
	}
      }
      for(set<unsigned>::const_iterator it=toDel.begin(); it != toDel.end(); ++it){
	sendhalo2nodes[p].erase( *it );
      }
    }
    
    {// Remove pressure nodes that p should already know through halo1.
      set<unsigned> toDel;
      for(set<unsigned>::const_iterator it=sendhalo2pnodes[p].begin(); it != sendhalo2pnodes[p].end(); ++it){
	if( shared_pnodes[p].find( *it ) != shared_pnodes[p].end() ){
	  toDel.insert( *it );
	}
      }   
      for(set<unsigned>::const_iterator it=toDel.begin(); it != toDel.end(); ++it){
	sendhalo2pnodes[p].erase( *it );
      }
    }
  }

  //
  // At this point we have identified all the information which we
  // want to communicate: 
  // vector< deque<unsigned> > elems2send( NProcs );
  // vector< set<unsigned> > nodes2send( NProcs );
  //
  
  // Make the send-packs
  vector< vector<char> > SendRecvBuffer(NProcs);
  
  { // Allocate space for buffers
    unsigned max_node_size      = max_nodepack_size();
    unsigned max_pnode_size     = max_pressurepack_size();
    unsigned max_element_size   = max_elementpack_size();
    unsigned space_for_unsigned = MPI::UNSIGNED.Pack_size(1, MPI::COMM_WORLD);
    
    for(int i=0; i<NProcs; i++){
      unsigned nbytes = space_for_unsigned          +
	space_for_unsigned*halo2Elements[i].size() +
	space_for_unsigned                          +
	max_element_size*halo2Elements[i].size()    +
	space_for_unsigned                          +
	max_node_size*sendhalo2nodes[i].size()      +
	space_for_unsigned                          +
	max_pnode_size*sendhalo2pnodes[i].size();
      
      SendRecvBuffer[i].resize( (unsigned)(1.1*nbytes) );
    }
  }
  
  vector<int> offsets(NProcs, 0); // int because of mpi calls
  
  // Pack.
  for(int i=0; i<NProcs; i++){
    int len = SendRecvBuffer[i].size();
    
    if( (i == MyRank)||(len == 0) )
      continue;
    
    char *buff = &(SendRecvBuffer[i][0]);
    
    // Elements
    unsigned cnt = halo2Elements[i].size();
    MPI::UNSIGNED.Pack(&cnt, 1, buff, len, offsets[i], MPI::COMM_WORLD);
    
    ECHO("Packing "<<cnt<<" halo2 elements for "<<i<<".");
    for(unsigned j=0; j<cnt; j++){

      // Dangerious?
      unsigned danger = 0;
      if( dangerElements.find( halo2Elements[i][j] ) != dangerElements.end() )
	danger = 1;
      MPI::UNSIGNED.Pack(&danger, 1, buff, len, offsets[i], MPI::COMM_WORLD);

      // Pack element
      element_list[ halo2Elements[i][j] ].pack(buff, len, offsets[i]);
    }
    // Nodes
    cnt = sendhalo2nodes[i].size();
    MPI::UNSIGNED.Pack(&cnt, 1, buff, len, offsets[i], MPI::COMM_WORLD);

    ECHO("Packing "<<cnt<<" halo2 nodes for "<<i<<".");    
    for(set<unsigned>::const_iterator it=sendhalo2nodes[i].begin(); it!=sendhalo2nodes[i].end(); ++it){
      const Node& node = node_list.unn( *it );
      node.pack(buff, len, offsets[i]);
    }
    
    // Pressure nodes
    cnt = sendhalo2pnodes[i].size();
    MPI::UNSIGNED.Pack(&cnt, 1, buff, len, offsets[i], MPI::COMM_WORLD);
    ECHO("Packing "<<cnt<<" halo2 pressure nodes for "<<i<<".");

    for(set<unsigned>::const_iterator it=sendhalo2pnodes[i].begin(); it!=sendhalo2pnodes[i].end(); ++it){
      const PressureNode& node = MFnode_list.unn( *it );
      node.pack(buff, len, offsets[i]);
    }
    
    assert(offsets[i] <= (int)SendRecvBuffer[i].size());
  }
  
  // Clean-up.
  dangerElements.clear();
  halo2Elements.clear();
  sendhalo2nodes.clear();
  sendhalo2pnodes.clear();
  
  // Send/recieve everything.
  allSendRecv(SendRecvBuffer);

  { // Unpacking.
    ECHO("Starting unpacking.");
    
    for(int p=0; p<NProcs; p++){
      offsets[p] = 0;
    }

    set<Element>  DangerElements;
    set<unsigned> ReceivedNodes;
    set<unsigned> ReceivedPressureNodes;
    vector< set<unsigned> > extraHalo(NProcs);
    
    // Paranoid
    assert(SendRecvBuffer.size() == (unsigned)NProcs);

#ifndef NDEBUG
    CHECK( element_list.size() );	
    do_element_headcount();
    CHECK( num_elements("total") );
    CHECK( num_elements("volume") );
    CHECK( num_elements("surface") );
#endif
   
    for(int p=0; p<NProcs; p++){
      int nbytes = SendRecvBuffer[p].size();
      if( (p == MyRank)||(nbytes == 0) )
	continue;
      
      char *buffer = &(SendRecvBuffer[p][0]);
      { // elements
	unsigned cnt;
	MPI::UNSIGNED.Unpack(buffer, nbytes, &cnt, 1, offsets[p], MPI::COMM_WORLD);
	
	ECHO("Unpacking "<<cnt<<" elements from "<<p<<".");
	
	for(unsigned j=0; j<cnt; j++){
	  Element element;
	  
	  ECHO("Unpacking "<<j<<"...");
	  
	  // Unpack danger flag.
	  unsigned danger;
	  MPI::UNSIGNED.Unpack(buffer, nbytes, &danger, 1, offsets[p], MPI::COMM_WORLD);
	  CHECK(danger);
	  
	  // Unpack element...      
	  element.unpack(buffer, nbytes, offsets[p]);
	  ECHO("unpacked.");
	  
	  if(danger){
	    ECHO("Danger element...taking evasive action."); 
	    DangerElements.insert( element );	      
	  }else{
	    element_list.push_back(element);
	  }  
	}
      } // finished unpacking elements.
      
      { // nodes
	unsigned cnt;
	MPI::UNSIGNED.Unpack(buffer, nbytes, &cnt, 1, offsets[p], MPI::COMM_WORLD);
	ECHO("Unpacking "<<cnt<<" nodes from "<<p<<".");
	
	for(unsigned j=0; j<cnt; j++){
	  Node node;
	  node.unpack(buffer, nbytes, offsets[p]);
	  unsigned unn   = node.get_unn();
	  unsigned owner = node.get_owner();

	  if(halo_nodes[owner].find( unn ) == halo_nodes[owner].end())
	    if(ReceivedNodes.find( unn ) == ReceivedNodes.end()){
	      assert( !node_list.contains_unn(unn) );
	      node_list.push_back( node );
	      ReceivedNodes.insert( unn );
	    }
	}
	
      } // finished unpacking nodes
      
      { // pressure nodes
	unsigned cnt;
	MPI::UNSIGNED.Unpack(buffer, nbytes, &cnt, 1, offsets[p], MPI::COMM_WORLD);
	ECHO("Unpacking "<<cnt<<" pressure nodes from "<<p<<".");
	
	for(unsigned j=0; j<cnt; j++){
	  PressureNode node;
	  node.unpack(buffer, nbytes, offsets[p]);
	  unsigned unn = node.get_unn();
	  unsigned owner = node.get_owner();

	  if(halo_pnodes[owner].find( unn ) == halo_pnodes[owner].end())
	    if(ReceivedPressureNodes.find( node.get_unn() ) == ReceivedPressureNodes.end()){
	      assert( !MFnode_list.contains_unn(unn) );
	      MFnode_list.push_back( node );
	      ReceivedPressureNodes.insert( node.get_unn() );
	      
	      // This next line is tricky as "node.get_owner()" is used
	      // rather than p. This is because the owner of this node is
	      // not necessarly the same as the processor that sent
	      // it.
	      extraHalo[owner].insert( node.get_unn() );
	    }
	}
	
      } // finished unpacking pressure nodes
      
    }
  
    //
    // Finally, add in the danger elements
    //
    if( !DangerElements.empty() ){
      ECHO("Add in danger elements.");
      
      for(set<Element>::const_iterator it=DangerElements.begin(); it!=DangerElements.end(); ++it)
	element_list.push_back( *it ); 
      DangerElements.clear();
    }
    
#ifndef NDEBUG
    CHECK( element_list.size() );	
    do_element_headcount();
    CHECK( num_elements("total") );
    CHECK( num_elements("volume") );
    CHECK( num_elements("surface") );
#endif

    //
    // Communicate extra halo stuff.
    //
    vector< set<unsigned> > extraShared(NProcs);
    ECHO("Communicating the extended halo...");
    assert(extraHalo[MyRank].size() == 0);    
    allSendRecv(extraHalo, extraShared);
    assert(extraShared[MyRank].size() == 0);
    ECHO("...done.");
    
    // 
    // Add in extra halo values...
    //
    for(int p=0; p<NProcs; p++){
      if(p==MyRank)
	continue;
      
      // add in extra halo values
      ECHO("Adding halo nodes for "<<p);
      halo_pnodes[p].insert(extraHalo[p].begin(), extraHalo[p].end());
      extraHalo[p].clear();
      
      ECHO("Adding shared nodes for "<<p);
      shared_pnodes[p].insert(extraShared[p].begin(), extraShared[p].end());
      extraShared[p].clear();
      
    }
    extraHalo.clear();
    extraShared.clear();
    
  } // finished unpacking

#ifndef NDEBUG
  len45 = element_list.size();
  for(unsigned i=0; i<len45; i++)
    for(unsigned j=0; j<len45; j++)
      if(i!=j)
	assert(element_list[i] != element_list[j]);
#endif
  
  ECHO("Halo2 created.");
  return; 
 
}
  











