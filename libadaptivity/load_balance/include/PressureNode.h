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
#ifndef H_PRESSURENODES
#define H_PRESSURENODES

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

class PressureNode{  
 public:
  PressureNode();                                   // "No argument" constructor
  PressureNode(const PressureNode& node);           // Copy constructor
  ~PressureNode();                                  // destructer
  
  // Overloaded operators.
  friend std::ostream &operator<<(std::ostream& out, const PressureNode& in);
  PressureNode &operator=(const PressureNode &in); 
  bool operator==(const PressureNode& in) const;
  bool operator!=(const PressureNode& in) const;
  bool operator<(const PressureNode& in) const;

  // Member functions.
  void            pack(char *buffer,   int& bsize, int& offset) const;
  void            unpack(char *buffer, int& bsize, int& offset);
  unsigned        pack_size() const;

  void            set_unn(const unsigned);
  unsigned        get_unn() const;
  
  void            set_gnn(const unsigned);
  unsigned        get_gnn() const;
  
  void            set_flags(const unsigned char);
  unsigned char   get_flags() const;

  void            set_owner(const unsigned);
  unsigned        get_owner() const;

  void            set_pressure(const samfloat_t);
  void            set_pressure(const std::vector<samfloat_t>);
  samfloat_t      get_pressure(const int) const;
  std::vector<samfloat_t> get_pressure() const;
  
 private:
  unsigned unn;                             // Universal Node Number.        
  unsigned gnn;                             // Partition-wide node number.   
  unsigned char  flags;                     // flags! 
  unsigned owner;
  // For multiphase calculations, pressure can have n-degrees of freedom
  std::vector<samfloat_t> pressure;
};

#endif









