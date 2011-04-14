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
#ifndef MI5_H
#define MI5_H

#include <confdefs.h>

#include <vector>
#include <map>
#include <set>
#include <string>
#include <deque>

#include "samtypes.h"
#include "Node.h"
#include "PressureNode.h"

// This class deals with 
// - Who knows what nodes and elements
// - Who will need to know what nodes and elements
// - What nodes and elements are we responsible for
//   sending and to whom.

class MI5{
 private:
  unsigned __nnodes;
  unsigned __nelems;
  unsigned MyRank;
  unsigned NProcs;
  
  // Given a node number, keep a list of partitions that 
  // already know its particulars.
  std::vector< std::set<unsigned short> > nodeKnowers;
  
  // Given a node number, who will need to know this node in the 
  // next round.
  std::vector< std::set<unsigned short> > nodeKnowers_next;
  
  // Given an element number, keep a list of partitions that 
  // already know its particulars.
  std::vector< std::set<unsigned short> > elementKnowers;
  
  // Given an element number, keep a list of partitions that 
  // will need to know this node in the next round.
  std::vector< std::set<unsigned short> > elementKnowers_next;

  void nodeKnownBy_insert(unsigned node, unsigned short knower);
  void node2bKnownBy_insert(unsigned node, unsigned short knower);
  void elementKnownBy_insert(unsigned elem, unsigned short knower);
  void element2bKnownBy_insert(unsigned elem, unsigned short knower);

  bool nodeKnownBy(unsigned node, unsigned short wiseguy);
  bool rankNeedsNode(unsigned node, unsigned short wiseguy);
  bool elementKnownBy(unsigned elem, unsigned short wiseguy);
  bool rankNeedsElement(unsigned elem, unsigned short wiseguy);
  
 public:  
  void clear();

  // List of nodes to be sent each rank.
  std::vector< std::deque<int> > nodes2send;

  // List of elements to be sent each rank.
  std::vector< std::deque<int> > elems2send;

  // List of pressure nodes to be sent each rank.
  std::vector< std::deque<int> > pnodes2send;

  // List of nodes that we are keeping.
  std::deque<int> new_owned_nodes;
  
  // New halo nodes...They are being copied to that they
  // don't get overwritten!
  std::vector< std::map<unsigned, Node> > new_halo_nodes;
  
  // List of elements that we are keeping.
  std::deque<int> new_elements;

  // In the case or halo nodes, there are some which may be
  // required, thus they are stored away in a cache.
  std::map<unsigned, Node> node_cache;

  MI5(int nnodes, int nelems);
  ~MI5();
  
  // This function compiles all the information to determine who
  // knows what and what I need to tell others.
  void spy(Mesh& );
  int get_nprocs(){
    return( NProcs);
  }
  int nnodes() const{
    return( __nnodes );
  }
  int nelems() const{
    return( __nelems );
  }
  int get_myrank() const{
    return( MyRank );
  }
};

#endif










