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

// CSR - Compressed Storage Row format

#include <mpi.h>
#include "c++debug.h"

#include <vector>
#include <set>
#include <deque>
#include <map>
#include <algorithm>
#include <cassert>
#include <cmath>

#include <string.h>
#include "Graph.h"
#include "samtypes.h"
#include <string.h>

using namespace std;
using csr::Graph; 

//
// Stuff for a csr Node-Node list
//
Graph::Graph(){}

Graph::Graph(const map<unsigned, set<unsigned> >& _edges){
  unsigned nverts         = _edges.size();
  
  bptr.resize( nverts+1 );
  bptr[0] = 0;
  
  unsigned ecnt = 0;
  unsigned vcnt = 0;
  for(map<unsigned, set<unsigned> >::const_iterator vit = _edges.begin(); vit != _edges.end(); ++vit){
    
    const set<unsigned> &nodes = (*vit).second;
    
    // Base pointer
    bptr[vcnt+1] = bptr[vcnt] + nodes.size(); 
    
    for(set<unsigned>::const_iterator eit = nodes.begin(); eit != nodes.end(); ++eit){      
      edges.push_back( *eit );
      ecnt++;
    }
    
    vcnt++;
  }
  
  nnodes = vcnt;
  nedges = ecnt;
  return;
}


Graph::~Graph(){
  bptr.clear();
  edges.clear();
  eweight.clear();
}

void Graph::build(const map<unsigned, set<unsigned> >& _edges, const deque<int>& _eweights){
  unsigned nverts         = _edges.size();
  bool using_edge_weights = !(_eweights.empty());
  
  bptr.resize( nverts+1 );
  bptr[0] = 0;
  
  unsigned ecnt = 0;
  unsigned vcnt = 0;
  for(map<unsigned, set<unsigned> >::const_iterator vit = _edges.begin(); vit != _edges.end(); ++vit){
    
    const set<unsigned> &nodes = (*vit).second;

    // Base pointer
    bptr[vcnt+1] = bptr[vcnt] + nodes.size(); 
    
    for(set<unsigned>::const_iterator eit = nodes.begin(); eit != nodes.end(); ++eit){      
      edges.push_back( *eit );
      
      if( using_edge_weights ){
	assert(_eweights.begin()+ecnt != _eweights.end());
	eweight.push_back( _eweights[ecnt] );
      }
      ecnt++;
    }

    vcnt++;
  }
  
  nnodes = vcnt;
  return;
}

void Graph::clear(){
  bptr.clear();
  edges.clear();
  eweight.clear();
}
int* Graph::get_cptr_bptr(){
  return &(bptr[0]);
}

int* Graph::get_cptr_edges(){
  return &(edges[0]);
}

int Graph::getNnodes(){
  return nnodes ;
}

int Graph::getNedges(){
  return nedges;
}

// You are probably going to have to see the documentation to see the
// full reasoning behind this function. Basically it identifies all
// the elements about each element (note that this is consistant for
// parallel work) and choses the edge weight to be the maximum of the
// functionals of these functionals.
void Graph::buildEdgeWeights_adapt(const deque< set<unsigned> >& nelist, 
				   const map<unsigned, unsigned>& unn2gnn, 
				   const deque<samfloat_t>& elem_fxnl,
				   samfloat_t functional_tolerence){
#ifndef NDEBUG
  unsigned max_eid    = elem_fxnl.size();
#endif
  unsigned max_edge   = edges.size();
  unsigned num_nodes  = bptr.size() -1;  

  // Calculate a edge weight which is equal to the maximum element
  // functional surrounding an edge divided by the functional
  // tolerence used by adaptivity.
  vector<float> weight(max_edge, 0.0);

  for(unsigned n=0; n<num_nodes; n++){
    
    for(int j = bptr[n]; j<bptr[n+1]; j++){
      assert((unsigned)j<max_edge);

      // get the gnn
      unsigned gnn;
      {
	const unsigned unn = edges[j];	
	
	map<unsigned, unsigned>::const_iterator ignn = unn2gnn.find( unn );
	assert( ignn != unn2gnn.end() );
	gnn = (*ignn).second;
	
	// the edge in question is defined by (n, gnn)
      }
      
      // Identify the elements about the edge
      set<unsigned> elems;
      {
	for(set<unsigned>::const_iterator ie = nelist[n].begin(); ie != nelist[n].end(); ++ie){
	  if( nelist[ gnn ].find( *ie ) != nelist[ gnn ].end() )
	    elems.insert( *ie );
	}
      }
      
      // Calculate the L-infinity norm.      
      {
	for(set<unsigned>::iterator it = elems.begin(); it != elems.end(); it++){
	  unsigned eid = *it;
	  assert(eid<max_eid);
	  // weight[j] = max((double)weight[j], (double)elem_fxnl[ eid ]/functional_tolerence);
	  weight[j] = max(weight[j], (float)elem_fxnl[ eid ]);
	}
      }
    }
  }
  
  float max_w = weight[0];
  // float min_w = weight[0];
  for(size_t i=1;i<(size_t)max_edge;i++){
    max_w = max(max_w, weight[i]);
    // min_w = max(min_w, weight[i]);
  }

  // Edges that have a value less or equal to 1.0 and get a weight of
  // 1. Greater values indicate regions that need to be adapted and so
  // are given a nominal weight of 20.
  eweight.resize(max_edge);
  
  for(unsigned i=0;i<max_edge;i++){
    if(weight[i]<=1.0)
      eweight[i] = 1;
    else
      eweight[i] = (int)(2*ceil(weight[i]));
    
    // Make the top 10% of edges uncutable
    if(0.9*max_w > 1.0)
      if(weight[i]>(0.9*max_w))
	eweight[i] = max_edge;
  }
  
  return;
}

void Graph::setVertexWeights(const std::deque<float> in){
  unsigned len = in.size();
  nweight.resize(len);
  for(unsigned i=0;i<len;i++){
    nweight[i] = max(1, (int)(in[i] + 0.5));
    assert(nweight[i] >= 1);
  }
}















