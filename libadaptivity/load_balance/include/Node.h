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
#ifndef H_NODES
#define H_NODES

#include "confdefs.h"
#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <vector>
#include <deque>
#include <map>
#include <set>
#include <string>
#include <iostream>

#include <cmath>

#include "samtypes.h"

/* **********************************
   NODES:***************************
   Node class
*/
// Node flags {0 2 4 8 10 20 40 80}
#define NODE_DEFAULT 0
#define NODE_DELETE  0x2
#define NODE_CHECKED 0x4
#define NODE_HALO    0x8

class Node{
 public:
  Node();                                   // "No argument" constructor
  Node(const Node& node);                   // Copy constructor
  ~Node();                                  // destructer
  
  // Overloaded operators.
  friend std::ostream &operator<<(std::ostream& out, const Node& n);
  const Node& operator=(const Node &in); 
  bool operator==(const Node& in) const;
  bool operator!=(const Node& in) const;
  bool operator<(const Node& in) const;

  // Member functions.
  void append_field(const std::vector<int>&);
  void append_field(const std::vector<samfloat_t>&);

  void append_field(const int*, const size_t);
  void append_field(const samfloat_t*, const size_t);

  void set_unn(const unn_t&);
  unn_t get_unn() const;
  
  void set_gnn(const gnn_t&);
  gnn_t get_gnn() const;
  
  void set_flags(const unsigned char);
  unsigned char get_flags() const;
  
  void set_current_owner(const unsigned short);
  unsigned short get_current_owner() const;
  
  void set_owner(const unsigned short);
  unsigned short get_owner() const;
   
  void set_future_owner(const unsigned short);
  unsigned short get_future_owner() const;
  
  void set_CE(const std::vector<eid_t>&);
  const std::vector<eid_t>& get_CE() const;
  const eid_t* get_cptr_CE() const;
  unsigned get_size_CE() const;

  void set_fields(const std::vector<samfloat_t>&);
  void set_fields(const samfloat_t* f, const unsigned flen);

  const std::vector<int>& get_ifields() const;
  const std::vector<samfloat_t>& get_fields() const;
  
  int pop_ifield(void);
  samfloat_t pop_field(void);

  const int* get_cptr_ifields() const;
  const samfloat_t* get_cptr_fields() const;

  size_t get_size_ifields() const;
  size_t get_size_fields() const;

  void set_coord(const std::vector<samfloat_t>&);
  void set_coord(const samfloat_t, const samfloat_t);
  void set_coord(const samfloat_t, const samfloat_t, const samfloat_t);
  const std::vector<samfloat_t>& get_coord() const;
  const samfloat_t* get_cptr_coord() const;
  samfloat_t get_x() const;
  samfloat_t get_y() const;
  samfloat_t get_z() const;
  size_t get_size_x() const;

  void set_metric(const std::vector<samfloat_t>&);
  void set_metric(const samfloat_t* m, const unsigned mlen);
  const std::vector<samfloat_t>& get_metric() const;
  const samfloat_t* get_cptr_metric() const;
  unsigned get_size_metric() const;
  
  void pack(char *buffer,   int& bsize, int& offset) const;
  void unpack(char *buffer, int& bsize, int& offset);
  unsigned pack_size() const;
  
 private:
  unn_t unn;                   // Universal Node Number.        
  gnn_t gnn;                   // Partition-wide node number.   
  unsigned char  flags;        // flags! 
  unsigned short owner[2];     // owner[0] == current owner, owner[1] == future owner 
  std::vector<eid_t> CE;       // Connected elements
  
  // Nodewise integer field values.
  std::vector<int> ifields;

  // Nodewise field values.
  std::vector<samfloat_t> fields;
  
  // Position.
  std::vector<samfloat_t> x;
  
  // Nodewise metric used for adaptivity
  std::vector<samfloat_t> metric;
};

#endif

