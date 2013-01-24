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
#ifndef NODEVECTOR_H
#define NODEVECTOR_H

#include "confdefs.h"
#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <vector>
#include <deque>
#include <map>
#include <set>
#include <string>
#include <cmath>
#include <iostream>
#include <assert.h>

#include "samtypes.h"

// This describes the extra bits ans pieces that are required for the 
// mixed formulation stuff. In particular, issues related to the 
// pressure nodes with in mixed formulation.

template<class NodeType>
class NodeVector: public std::deque<NodeType>{
 public:
  NodeVector();
  NodeVector(const NodeVector&);
  ~NodeVector();
  
  // Overloaded operators.
  NodeVector &operator=(const NodeVector<NodeType> &);
#ifndef NDEBUG
  const NodeType &operator[](const unsigned i) const;
  NodeType &operator[](const unsigned i);
#endif

  // Other member functions.
  NodeType& unn(const unsigned);
  void refresh_unn2gnn();
  void refresh_pcnt();
  bool contains_unn(const unsigned);

  void made_edit();
  void shrink(const unsigned);
  void expand(const unsigned);

  void push_back(const NodeType&);
  unsigned size() const;
  unsigned psize();
  void clear();

 private:
  int MyRank;
  
  bool complete_unn2gnn;
  bool complete_pcnt;

  unsigned len, plen;
  std::map<unsigned, unsigned> unn2gnn;
};

template<class NodeType>
NodeVector<NodeType>::NodeVector(){
#ifdef HAVE_MPI
  if( MPI::Is_initialized() ){
    MyRank = MPI::COMM_WORLD.Get_rank();
  }else{
    MyRank = 0;
  }
#else
  MyRank = 0;
#endif
  plen = 0;
  len = 0;
  complete_pcnt    = true;
  complete_unn2gnn = true;
}

template<class NodeType>
NodeVector<NodeType>::NodeVector(const NodeVector& newkid){
  *this = newkid;
}

template<class NodeType>
NodeVector<NodeType>::~NodeVector(){
  clear();
}

// Overloaded operators.
template<class NodeType>
NodeVector<NodeType>& NodeVector<NodeType>::operator=(const NodeVector<NodeType>& in){
  len = in.len;
  plen = in.plen;
  this->insert( (*this).end(), in.begin(), in.end() );
  complete_unn2gnn = in.complete_unn2gnn;
  complete_pcnt    = in.complete_pcnt;
  
  if(complete_unn2gnn){
    unn2gnn = in.unn2gnn;
  }
  
  return *this;
}

#ifndef NDEBUG
template<class NodeType>
NodeType& NodeVector<NodeType>::operator[](const unsigned i){  
  // don't forgive out of bounds referencing
  assert(i<len); 
  return *(this->begin() + i);
}

template<class NodeType>
const NodeType& NodeVector<NodeType>::operator[](const unsigned i) const{  
  // don't forgive out of bounds referencing
  assert(i<len); 
  
  return *(this->begin() + i);
}
#endif

// Other member functions
template<class NodeType>
void NodeVector<NodeType>::push_back(const NodeType& node){
#ifdef HAVE_MPI
  complete_unn2gnn = false;
  
  // Add node.
  this->insert(this->end(), node);
  
  // Update count
  len++;
  if( MPI::Is_initialized() ){
    if( ((unsigned)node.get_owner()) == (unsigned)(MyRank) )
      plen++;
  }
#endif
  return;
}

template<class NodeType>
void NodeVector<NodeType>::shrink(const unsigned newlen){
  made_edit();

  assert(newlen <= len);

  this->resize( newlen );
  len = newlen;

}

template<class NodeType>
void NodeVector<NodeType>::expand(const unsigned newlen){
  made_edit();

  assert(newlen >= len);
  
  this->resize( newlen );
  len = newlen;
  
}

template<class NodeType>
NodeType& NodeVector<NodeType>::unn(const unsigned in){
  
  // If the LUT is incomplete, generate a new one before continuing.
  if(!complete_unn2gnn)
    refresh_unn2gnn();
  
  // Locate the entry
  std::map<unsigned, unsigned>::const_iterator n = unn2gnn.find(in);
  
  // If we cannot get a reference, try building a new map
  if(n == unn2gnn.end()){
    refresh_unn2gnn();
    n = unn2gnn.find(in);
    
    if(n == unn2gnn.end()){
      assert(false);
    }
  }
  
  // Return the node reference.
  return *(this->begin() + (*n).second);
  
}

template<class NodeType>
void NodeVector<NodeType>::made_edit(){
  complete_unn2gnn = false;
  complete_pcnt    = false;
}

template<class NodeType>
unsigned NodeVector<NodeType>::size() const{
  return len;
}

template<class NodeType>
unsigned NodeVector<NodeType>::psize(){
  if(!complete_pcnt)
    refresh_pcnt();
  
  return plen;
}

template<class NodeType>
void NodeVector<NodeType>::clear(){
  made_edit();
  len = 0;
  unn2gnn.clear();
  this->erase(this->begin(), this->end());
}

template<class NodeType>
bool NodeVector<NodeType>::contains_unn(const unsigned unn){
  if(!complete_unn2gnn)
    refresh_unn2gnn();
  
  return( unn2gnn.find(unn) != unn2gnn.end() );
}

template<class NodeType>
void NodeVector<NodeType>::refresh_unn2gnn(){
  
  unn2gnn.clear();
  
  unsigned gnn = 0;
  for(typename NodeVector<NodeType>::iterator it = this->begin(); it != this->end(); ++it){
    
    // Get the 
    unsigned _unn = (*it).get_unn();
    
    // Re-set the gnn
    (*it).set_gnn( gnn );

    // Once again, we don't forgive mistakes;
    // duplicate entries are unacceptable.
    assert(unn2gnn.find(_unn) == unn2gnn.end());
    
    // Set the map entry
    unn2gnn[ _unn ] = gnn;
    
    gnn++;
  }
  
  complete_unn2gnn = true;
  
  return;
}

template<class NodeType>
void NodeVector<NodeType>::refresh_pcnt(){
#ifdef HAVE_MPI
  plen = 0;
  if( MPI::Is_initialized() ){
    for(typename NodeVector<NodeType>::iterator it = this->begin(); it != this->end(); ++it){
      if( (*it).get_owner() == (unsigned)(MyRank))
	plen++;
    }
  }
#endif
}
#endif


