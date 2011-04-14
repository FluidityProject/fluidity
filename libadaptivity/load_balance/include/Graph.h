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
#ifndef H_GRAPH
#define H_GRAPH

#include <confdefs.h>
#include "samtypes.h"

#include <deque>
#include <map>
#include <set>
#include <vector>

namespace csr{

class Graph{
 public:  
  Graph();
  Graph(const std::map<unsigned, std::set<unsigned> >&);
  ~Graph();
  
  void build(const std::map<unsigned, std::set<unsigned> >& _edges, 
	     const std::deque<int>& _eweights);
  void setVertexWeights(const std::deque<float>);
  void buildEdgeWeights_adapt(const std::deque< std::set<unsigned> >&, 
			      const std::map<unsigned, unsigned>&, 
			      const std::deque<samfloat_t>&,
			      samfloat_t);
  

  void clear();
  int *get_cptr_bptr();
  int *get_cptr_edges();
  int getNnodes();
  int getNedges();

  //
  // Provides interface for METIS
  //

  // Use for small numbers of domains (>=8)
  void PartGraphRecursive(int nparts, std::vector<int>& noddom);
  // Otherwise use...
  void PartGraphKway(int nparts, std::vector<int>& noddom);
  // Minimize communication volume rather than edge-cut
  void PartGraphVKway(int nparts, std::vector<int>& noddom);
  
  // Multi-constraint partitioning
  void mCPartGraphRecursive(int nparts, int nconstraints, std::vector<int>& noddom);
  void mCPartGraphKway(int nparts, int nconstraints, std::vector<samfloat_t>& vcon, 
		       std::vector<int>& noddom);

  // Partitions of different sizes.
  void WPartGraphRecursive(int nparts, std::vector<samfloat_t>& pweights, std::vector<int>& noddom);
  void WPartGraphKway(int nparts, const std::vector<samfloat_t>& pweights, std::vector<int>& noddom);
  void WPartGraphVKway(int nparts, std::vector<samfloat_t>& pweights, std::vector<int>& noddom);

  // Sparse matrix reording Routines
  void NodeND(std::vector<int>& perm, std::vector<int>& iperm);

  //
  // Provides interface for parMETIS
  //
  void Repart(std::vector<int>& vtxdist, std::vector<int>& noddom, int nparts);
  void RepartRemap(std::vector<int>& vtxdist,   std::vector<int>& noddom);

 private:
  unsigned nnodes;
  std::vector<int> nweight; // Node weights

  unsigned nedges;
  std::vector<int> bptr;    // Base pointer into edges
  std::vector<int> edges;   // Edges
  std::vector<int> eweight; // Edge weights
};

} // end of namespace
#endif
