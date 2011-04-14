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
//
// This function is intend to be called from fluidity. It basically
// strips out all the elements and nodes that are not part of a
// regular mesh - as in we remove anything related to halo-2 elements.
//
#include <confdefs.h>

#include <mpi.h>

#include <cassert>
#include <deque>
#include <map>
#include <set>
#include <vector>

#include "comTools.h"
#include "c++debug.h"
#include "samtypes.h"

using namespace std;

// Function used to find all the halo1 and halo2 elements. 
int getH1andH2_elems(deque<unsigned>& halo1_elems, deque<unsigned>& halo2_elems,
		     const unsigned nNodes, const unsigned nPrivateNodes, 
		     const unsigned nElements, const unsigned nloc,
		     const int * const enlist){
  
  halo1_elems.clear();
  halo2_elems.clear();
  
  // If there is no nodes then there are no elements.
  if( (nNodes == 0)||(nElements == 0) )
    return(0);
  
  // Identify all the halo1 and halo2 elements.
  for(unsigned i=0; i<nElements; i++){
    
    // For each element, start by assuming that everything is a halo2
    // element and that no element might be halo1.
    bool halo2       = true;
    bool halo1_maybe = false;
    
    // Test these assumptions by examining each of the nodes in that
    // element.
    const int * const element = enlist + nloc*i;
    for(unsigned j=0; j<nloc; j++){
      
#ifndef NDEBUG
      if( element[j] < 1 ){
	ERROR("Corrupted element-node list. Found node id "
	      <<element[j]<<" in element "<<i<<".");
	return(-1);
      }
      if( element[j] > (int)nNodes ){
	ERROR("Found a node id, "<<element[j]<< ", greater "
	      << "than the total number of nodes, "<<nNodes<<".");
	return(-2);
      }
#endif
      
      if(element[j] <= (int)nPrivateNodes){
	// If a element contains a node which is priviatly owned then
	// it is not a halo2 element.
	halo2 = false;
      }else{
	// If an element contains any node that is not privatly owned
	// then it is either a halo1 or halo2 element.
	halo1_maybe = true;
      }
      
    }
    
    // Record any halo1 elements and halo2 elements
    if(halo2){
      // Record as a halo2 element.
      halo2_elems.push_back( i );
    }else{
      if(halo1_maybe){
	// If an element is not halo2 but is in the halo, then it's
	// obviously a halo1 element which we wish to store.
	halo1_elems.push_back( i );
      }
    }
    
  }
  
  return(0);
}

// Create the set of all nodes at are in halo2 but not halo1.
int mk_h2node_set(set<unsigned>& halo2_nodes, 
		  const deque<unsigned>& halo1_elems, const deque<unsigned>& halo2_elems,
		  const unsigned nloc, const int * const enlist){
  
  if(halo2_elems.empty())
    return(0);

  // First identify all the
  // nodes in the halo2 elements.
  for(deque<unsigned>::const_iterator it=halo2_elems.begin(); it!=halo2_elems.end(); ++it){
    const int * const element = enlist + nloc*(*it);
    
    for(unsigned i=0; i<nloc; i++)
      halo2_nodes.insert( element[i] );
  }
  
  // Identify all nodes in the halo1 elements.
  set<unsigned> halo1_nodes;
  for(deque<unsigned>::const_iterator it=halo1_elems.begin(); it!=halo1_elems.end(); ++it){
    const int * const element = enlist + nloc*(*it);
    
    for(unsigned i=0; i<nloc; i++)
      halo1_nodes.insert( element[i] );
  }
  
  // If a node is found in both sets, halo1_nodes and halo2_nodes, when
  // remove node from the halo2_nodes set.
  for(set<unsigned>::const_iterator it = halo1_nodes.begin(); it != halo1_nodes.end(); ++it){
    set<unsigned>::iterator pos = halo2_nodes.find( *it );
    
    if( pos != halo2_nodes.end() )
      halo2_nodes.erase( pos );
  }
  
  return(0);
}

// returns the new number of elements
int cleanup_enlist(deque<unsigned>& halo2_elems, int *enlist, int nelems, unsigned nloc){ 
  
  unsigned del_cnt = halo2_elems.size();  
  if(del_cnt == 0) return(nelems);
  
  // Start deleting
  do{
      
    deque<unsigned>::reverse_iterator pos_back = halo2_elems.rbegin();
    deque<unsigned>::iterator pos_front = halo2_elems.begin();
    
    if(*pos_back == (unsigned)(nelems - 1)){
      
      // If the last element is to be deleted then delete from back.
      halo2_elems.pop_back();

    }else{
      
      // Copy element forward from back.
      int *element1 = enlist + nloc*(*pos_front);
      int *element2 = enlist + nloc*(nelems - 1);

      for(unsigned i=0; i<nloc; i++)
	element1[i] = element2[i];

      halo2_elems.pop_front();

    }
    
    // One less element to worry about.
    nelems--;

  }while( !halo2_elems.empty() );

  return(nelems);
}

#define flstriph2_fc F77_FUNC(flstriph2, FLSTRIPH2)
extern "C" {
  void flstriph2_fc(int *nNodes, const int *nPrivateNodes, const int *nprocs,
		    int VolumeENlist[],   int *nVolumeElems,   const int *nloc, 
		    int SurfaceENlist[],  int SurfaceIDs[], int *nSurfaceElems,  const int *snloc, 
		    samfloat_t X[], samfloat_t Y[], samfloat_t Z[],
		    samfloat_t fields[], int *nfields,  int *fstride,
		    samfloat_t metric[],
		    int scatter[],  int *nscatter){
  
    { //
      // Identify volume elements to be deleted, and delete.
      //
      unsigned last_element = 0;
      for(size_t i=0; i<(size_t)*nVolumeElems; i++){
	bool keep = false;
	int *elem = VolumeENlist + i*(*nloc);
	for(size_t j=0; j<(size_t)*nloc; j++){
	  if(elem[j] <= *nPrivateNodes){
	    keep = true;
	    break;
	  }
	}
      
	if(keep){
	  if(last_element != i){
	    int *last = VolumeENlist + last_element*(*nloc);
	    for(size_t j=0; j<(size_t)*nloc; j++){
	      last[j] = elem[j];
	    }
	  }
	  last_element++;
	}
      }
    
      // Update the number of volume elements
      assert(*nVolumeElems >= last_element);
      *nVolumeElems = last_element;
    } // unnecessary volume elements deleted

    {  //
      // Identify surface elements to be deleted, and delete.
      //
      unsigned last_element = 0;
      for(size_t i=0; i<(size_t)*nSurfaceElems; i++){
	bool keep = false;
	int *elem = SurfaceENlist + i*(*snloc);
	for(size_t j=0; j<(size_t)*snloc; j++){
	  if(elem[j] <= *nPrivateNodes){
	    keep = true;
	    break;
	  }
	}
      
	if(keep){
	  if(last_element != i){
	    int *last = SurfaceENlist + last_element*(*snloc);
	    for(size_t j=0; j<(size_t)*snloc; j++){
	      last[j] = elem[j];
	    }
	    SurfaceIDs[last_element] = SurfaceIDs[i];
	  }
	  last_element++;
	}
      }

      // Update the number of surface elements
      assert(*nSurfaceElems >= last_element);
      *nSurfaceElems = last_element;
    } // unnecessary surface elements deleted

    // Find all halo1 nodes  
    set<unsigned> halo1Nodes;
    for(size_t i=0; i<(size_t)*nscatter; i++)
      halo1Nodes.insert( scatter[i] );
  
    // Establish the renumbering of the halo nodes
    int last_node = *nPrivateNodes;
    map<int, int> renumbering;
    for(set<unsigned>::iterator it=halo1Nodes.begin(); it!=halo1Nodes.end(); ++it){
      last_node++;
      renumbering[ *it ] = last_node;
    }
  
    // renumber nodes in elements-node list
    for(size_t i=0; i<(size_t)*nVolumeElems; i++){
      int *elem = VolumeENlist + i*(*nloc);
      for(size_t j=0; j<(size_t)*nloc; j++)
	if(elem[j] > *nPrivateNodes){
	  assert(renumbering.find(elem[j]) != renumbering.end());
	  elem[j] = renumbering[ elem[j] ];
	}
    }
    for(size_t i=0; i<(size_t)*nSurfaceElems; i++){
      int *elem = SurfaceENlist + i*(*snloc);
      for(size_t j=0; j<(size_t)*snloc; j++)
	if(elem[j] > *nPrivateNodes){
	  assert(renumbering.find(elem[j]) != renumbering.end());
	  elem[j] = renumbering[ elem[j] ];
	}
    }

    // Sort out coordinates and field values
    for(size_t i=*nPrivateNodes; i<(size_t)*nNodes; i++)
      if(renumbering.find(i+1) != renumbering.end()){
	if(renumbering[i+1]<1){
	  cerr<<"ERROR "<<__FILE__<<", "<<__LINE__<<": renumbering foobar.\n";
	  exit(-1);
	}
	size_t gid = renumbering[i+1] - 1;
	
	if(gid==i) 
	  continue;
      
	X[gid] = X[i];
	Y[gid] = Y[i];
	Z[gid] = Z[i];
      
	for(size_t j=0; j<(size_t)*nfields; j++)
	  fields[j*(*fstride) + gid] = fields[j*(*fstride) + i];
      
	for(size_t j=0; j<9; j++)
	  metric[gid*9 + j] = metric[i*9 + j];
      }
  
    // Sort out the halo numbering for the halo
    for(size_t i=0; i<(size_t)*nscatter; i++){
      assert(renumbering.find(scatter[i]) != renumbering.end());
#ifndef NDEBUG
      cout<<i<<" "<<scatter[i]<<" --> "<<renumbering[ scatter[i] ]<<endl;
#endif
      scatter[i] = renumbering[ scatter[i] ];
    }
  
    *nNodes = last_node;
  
    return;
  }
}
