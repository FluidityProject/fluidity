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
#ifndef H_ELEMENTS
#define H_ELEMENTS

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

#include "samtypes.h"

// **********************************
// ELEMENTS:*************************
// Element class
#define ELM_DEFAULT 0
#define ELM_HALO    0x2
#define ELM_SURFACE 0x4
#define ELM_VOLUME  0x8
#define ELM_MNO     0x10 // minimum node owner (ie. managed element)

class Element{
 public:
  Element();
  Element(const Element& element);                   // Copy constructor
  ~Element();
  
  // Overloaded operators.
  friend std::ostream &operator<<(std::ostream& out, const Element& in);
  Element &operator=(const Element& in);
  bool operator==(const Element& in) const;
  bool operator!=(const Element& in) const;
  bool operator<(const Element& in) const;

  // Member functions.
  void append_field(const std::vector<int>&);
  void append_field(const std::vector<samfloat_t>&);

  void append_field(const int*, const size_t);
  void append_field(const samfloat_t*, const size_t);

  void  set_eid(const eid_t);
  eid_t get_eid() const;
  
  void set_flags(const unsigned char flags);
  unsigned char get_flags() const;
  
  void set_enlist(const std::vector<unn_t>&);
  void set_MFenlist(const std::vector<unn_t>&);

  const std::vector<unn_t>& get_enlist() const;
  const std::vector<unn_t>& get_MFenlist() const;
  
  unsigned get_size_enlist() const;
  unsigned get_size_MFenlist() const;
  
  void set_ifields(const std::vector<int>&);
  void set_fields(const std::vector<samfloat_t>&);
  
  const std::vector<int>& get_ifields() const;  
  const std::vector<samfloat_t>& get_fields() const;
  
  size_t get_size_ifields() const;
  size_t get_size_fields() const;

  void pack(char *buffer, int& bsize, int& offset) const;
  void unpack(char *buffer, int& bsize, int& offset);
  unsigned pack_size() const;

 private:
  eid_t eid;                   // element id number.
  unsigned char flags;         // states: {ELM_DEFAULT, ELM_CHECK, ELM_DELETE} */
  std::vector<unn_t> nodes;    // the element-node list. 
  std::vector<unn_t> MFnodes;  // Stores the varies mixed-formulation elements aligned with this element.
  
  // All integer variables associated with this element. This will
  // include region id information for example.
  std::vector<int> ifields;
  
  // All real variables associated with this element.  NOTE: in the
  // case of discontinious methods, information regarding field values
  // would be included here.
  std::vector<samfloat_t> fields;
};


#endif
