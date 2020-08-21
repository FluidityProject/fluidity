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
#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <set>
#include <assert.h>

#include "sam_mpi.h"
#include "c++debug.h"
#include "Mesh.h"
#include "PressureNode.h"

using namespace std;

PressureNode::PressureNode(){
  unn   = 0;
  gnn   = 0;
  flags = 0;
  owner = 0;
}
PressureNode::PressureNode(const PressureNode& node){
  *this = node;
}
PressureNode::~PressureNode(){}

// Assignment operator for PressureNode.
PressureNode &PressureNode::operator=( const PressureNode &n ){
  unn      = n.unn;
  gnn      = n.gnn;
  flags    = n.flags;
  owner    = n.owner;
  pressure = n.pressure;
  
  return( *this );
}

// Equivalence operator for pressure node
bool PressureNode::operator==(const PressureNode& node) const{
  if(unn != node.unn)
    return false;
  
  return(true);
}

bool PressureNode::operator!=(const PressureNode& node) const{
  return( !(*this == node) );
}

bool PressureNode::operator<(const PressureNode& node) const{
  return(unn<node.unn);
}

// UNN stuff
void PressureNode::set_unn(const unn_t in){ 
  unn = in; 
}
unn_t PressureNode::get_unn() const{
  return unn; 
}

// GNN stuff
void PressureNode::set_gnn(const unn_t in){ 
  gnn = in; 
}
unn_t PressureNode::get_gnn() const{ 
  return gnn; 
}

// Flags??
void PressureNode::set_flags(const unsigned char in){ 
  flags = in; 
}
unsigned char PressureNode::get_flags() const{ 
  return flags; 
}

// Owner stuff.
void PressureNode::set_owner(const unsigned in){ 
  owner = in; 
}
unsigned PressureNode::get_owner() const{ 
  return owner; 
}

// Pressure value
void PressureNode::set_pressure(const samfloat_t p){ 
  pressure.clear(); 
  pressure.push_back(p); 
}
void PressureNode::set_pressure(const vector<samfloat_t> p){ 
  pressure = p; 
}
samfloat_t PressureNode::get_pressure(const int i) const{ 
  assert( i < (int)pressure.size() );
  return pressure[i]; 
}
vector<samfloat_t> PressureNode::get_pressure() const{ 
  return pressure; 
}

void PressureNode::pack(char *buffer, int& bsize, int& offset) const{
  MPI_Pack(&unn,           1, UNN_T, buffer, bsize, &offset, MPI_COMM_WORLD);
  MPI_Pack(&gnn,           1, GNN_T, buffer, bsize, &offset, MPI_COMM_WORLD);
  MPI_Pack(&flags, 1, MPI_UNSIGNED_CHAR, buffer, bsize, &offset, MPI_COMM_WORLD);
  MPI_Pack(&owner,     1, MPI_UNSIGNED, buffer, bsize, &offset, MPI_COMM_WORLD);
  
  unsigned nDofF = pressure.size();
  MPI_Pack(&nDofF, 1, MPI_UNSIGNED, buffer, bsize, &offset, MPI_COMM_WORLD);
  if( nDofF>0 )
    MPI_Pack(&(pressure[0]), nDofF, SAMFLOAT, buffer, bsize, &offset, MPI_COMM_WORLD);

  return;
}

void PressureNode::unpack(char *buffer, int& bsize, int& offset){
  MPI_Unpack(buffer, bsize, &offset, &unn,  1, UNN_T, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bsize, &offset, &gnn,  1, GNN_T, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bsize, &offset, &flags,1, MPI_UNSIGNED_CHAR, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bsize, &offset, &owner,1, MPI_UNSIGNED, MPI_COMM_WORLD);
  
  unsigned nDofF;
  MPI_Unpack(buffer, bsize, &offset, &nDofF,1, MPI_UNSIGNED, MPI_COMM_WORLD);
  pressure.resize( nDofF );
  if( nDofF>0 )
    MPI_Unpack(buffer, bsize, &offset, &(pressure[0]), nDofF, SAMFLOAT, MPI_COMM_WORLD);
  
  return;
}

unsigned PressureNode::pack_size() const{
  int total=0,s1,s2,s3,s4;
  
  MPI_Pack_size(1, UNN_T, MPI_COMM_WORLD, &s1);
  MPI_Pack_size(1, GNN_T, MPI_COMM_WORLD, &s2);
  MPI_Pack_size(1, MPI_UNSIGNED_CHAR, MPI_COMM_WORLD, &s3);
  MPI_Pack_size(1, MPI_UNSIGNED, MPI_COMM_WORLD, &s4);
  total = s1 + s2 + s3 + s4;
  {
    unsigned nfields = pressure.size();      
    total += MPI_Pack_size(1, MPI_UNSIGNED, MPI_COMM_WORLD, &s1);
    total += MPI_Pack_size(nfields, SAMFLOAT, MPI_COMM_WORLD, &s2);
    total += (s1 + s2);
  }

  return total;
}

// printing operator for PressureNode
ostream& operator<<(ostream& s, const PressureNode& n){
  s << "PressureNode unn number\t\t"   << n.unn                  << endl;
  s << "PressureNode gnn number\t\t"   << n.gnn                  << endl;
  s << "PressureNode flags\t\t\t"      << (unsigned)(n.flags)    << endl;
  s << "PressureNode pressure\t\t";

  for(vector<samfloat_t>::const_iterator it=n.pressure.begin(); it != n.pressure.end(); ++it)
    s << *it << " ";
  s << endl;
   
  return(s);
}
