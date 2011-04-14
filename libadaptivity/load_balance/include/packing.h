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
#ifndef PACKING_H
#define PACKING_H

#include <confdefs.h>

#include <algorithm>
#include <map>
#include <vector>

#include "Mesh.h"
#include "MI5.h"

class packing{
 private:
  unsigned int NProcs;
  unsigned int MyRank;
  unsigned int __nnodes;
  unsigned int __nelems;

  void pack_nodes(const MI5 &,   Mesh &);
  void pack_pnodes(const MI5 &,   Mesh &);
  void pack_elems(const MI5 &,   Mesh &);
  void unpack_nodes(MI5 &, Mesh &);
  void unpack_pnodes(MI5 &, Mesh &);
  void unpack_elems(MI5 &, Mesh &);

 public:
  packing(const MI5&, const Mesh&);
  ~packing();
  
  unsigned int packCnt2Bytes(unsigned int);

  void pack(const MI5& intelligance_report,   Mesh &mesh);
  void send();
  void unpack(MI5& intelligance_report, Mesh &mesh);
  void clear();
  
  std::vector<int> offsets; // int because of mpi calls
  std::vector< std::vector<char> > SendRecvBuffer;
};

#endif














