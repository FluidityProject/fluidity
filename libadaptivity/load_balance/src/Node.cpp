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
#include <cassert>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "sam_mpi.h"
#include "c++debug.h"
#include "Node.h"

using namespace std;

// Node "No argument" constructor
Node::Node(){
  unn   = 0;
  gnn   = 0;
  flags = NODE_DEFAULT;
  owner[0] = 0;
  owner[1] = 0;
} 

// Node copy constructor
Node::Node(const Node& node){
  *this = node;
}

// Node destructer
Node::~Node(){}

// Assignment operator for Node.
const Node &Node::operator=( const Node &n ){  
  unn       = n.unn;
  gnn       = n.gnn;
  flags     = n.flags;
  owner[0]  = n.owner[0];
  owner[1]  = n.owner[1];
  CE        = n.CE;
  ifields   = n.ifields;
  fields    = n.fields;
  x         = n.x;
  metric    = n.metric;
  
  return( *this );
}

// Equivalance operator for Node. Two Nodes will only be consider 
//  equal if they have the same uid
bool Node::operator==(const Node& in) const{
  return(unn == in.unn);
}
bool Node::operator!=(const Node& in) const{
  return( !(*this == in) );
}
bool Node::operator<(const Node& in) const{
  return(unn<in.unn);
}

// This is the printing operator for node
ostream &operator<<(ostream& s, const Node& n){
  s << "Node profile:" << endl;
  s << "UNN\t\t"       << n.unn               << endl;
  s << "GNN\t\t"       << n.gnn               << endl;
  s << "Flag code\t\t\t"         << (unsigned)(n.flags) << endl;
  s << "Current node owner\t\t"  << n.owner[0]          << endl;
  s << "Future node owner\t\t"   << n.owner[1]          << endl;
  { // element info
    s << "Number of connected elements\t" << n.CE.size() << endl;
    for(vector<eid_t>::const_iterator it = n.CE.begin(); it != n.CE.end(); ++it)
      s << *it << endl;
  }
  { // Field info
    int cnt = 0;
    s << "Number of floating point field components for this node\t" << n.fields.size() << endl;
    for(vector<samfloat_t>::const_iterator it = n.fields.begin(); it != n.fields.end(); ++it)
      s << " Value of field " << cnt++ << "\t\t" << *it << endl;

    cnt = 0;
    s << "Number of integer field components for this node\t" << n.ifields.size() << endl;
    for(vector<int>::const_iterator it = n.ifields.begin(); it != n.ifields.end(); ++it)
      s << " Value of field " << cnt++ << "\t\t" << *it << endl;
  }
  { // Position
    s << "Position\t\t\t(";
    for(vector<samfloat_t>::const_iterator it = n.x.begin(); it != n.x.end(); ++it)
      s << *it << " ";
    s << ")" << endl;
  }
  { // metric
    int cnt = n.metric.size();
    int sqrtcnt = (int)( sqrt( (double)cnt ) + 0.5 );
    assert( sqrtcnt*sqrtcnt == cnt ); // expecting a S X S matrix.
    s << "Metric: \t";
    cnt = 0;
    for(vector<samfloat_t>::const_iterator it = n.metric.begin(); it != n.metric.end(); ++it){
      s << *it << "\t" ;
      cnt++;
      if(cnt%sqrtcnt == 0){
	if(it+1 != n.metric.end())
	  s << endl << "\t\t";
      }
    }
  }

  return(s);
}

void Node::append_field(const std::vector<int>& _field){
  ifields.insert(ifields.end(), _field.begin(), _field.end());
  
}

void Node::append_field(const std::vector<samfloat_t>& _field){
  fields.insert(fields.end(), _field.begin(), _field.end());
  
}

void Node::append_field(const int* _field, const size_t block_size){
  for(size_t i=0;i<block_size;i++)
    ifields.push_back(_field[i]);
}

void Node::append_field(const samfloat_t* _field, const size_t block_size){
  for(size_t i=0;i<block_size;i++)
    fields.push_back(_field[i]);
}

void  Node::set_unn(const unn_t& _unn){ 
  unn = _unn; 
}
unn_t Node::get_unn() const{ 
  return unn; 
}

void  Node::set_gnn(const gnn_t& _gnn){ 
  gnn = _gnn; 
}
gnn_t Node::get_gnn() const{ 
  return gnn; 
}

void          Node::set_flags(const unsigned char _flags){ 
  flags = _flags; 
}
unsigned char Node::get_flags() const{ 
  return flags; 
}

void           Node::set_current_owner(const unsigned short in){ 
  owner[0] = in; 
}
unsigned short Node::get_current_owner() const{ 
  return owner[0]; 
}
void           Node::set_owner(const unsigned short in){ 
  set_current_owner(in); 
}
unsigned short Node::get_owner() const{ 
  return( get_current_owner() ); 
}
void           Node::set_future_owner(const unsigned short _new){ 
  owner[1] = _new; 
}
unsigned short Node::get_future_owner() const{ 
  return owner[1]; 
}

void                      Node::set_CE(const vector<eid_t>& _CE){ CE = _CE; }
const vector<eid_t>& Node::get_CE() const{ 
  return CE; 
}
const eid_t*              Node::get_cptr_CE() const{ 
  return &(CE[0]); 
}
unsigned                  Node::get_size_CE() const{ 
  return CE.size(); 
}

void Node::set_fields(const vector<samfloat_t>& _fields){ 
  fields = _fields; 
}

void Node::set_fields(const samfloat_t* f, const unsigned flen){
  fields.resize( flen );
  for(unsigned i=0; i<flen; i++)
    fields[i] = f[i];
}

const vector<samfloat_t>& Node::get_fields() const{ 
  return fields; 
}

const vector<int>& Node::get_ifields() const{ 
  return ifields; 
}

int Node::pop_ifield(void)
{
  int i;
  i = *(ifields.rbegin());
  ifields.pop_back();
  return i;
}

samfloat_t Node::pop_field(void)
{
  samfloat_t r;
  r = *(fields.rbegin());
  fields.pop_back();
  return r;
}

const samfloat_t* Node::get_cptr_fields() const{ 
  return &(fields[0]); 
}

const int* Node::get_cptr_ifields() const{ 
  return &(ifields[0]); 
}

size_t Node::get_size_fields() const{ 
  return fields.size(); 
}

size_t Node::get_size_ifields() const{ 
  return ifields.size(); 
}

void Node::set_coord(const vector<samfloat_t>& _pos){
  x = _pos; 
}

void Node::set_coord(const samfloat_t _x, const samfloat_t _y){ 
  x.resize(2); 
  x[0] = _x; 
  x[1] = _y;
}

void Node::set_coord(const samfloat_t _x, const samfloat_t _y, const samfloat_t _z){
  x.resize(3);
  x[0] = _x;
  x[1] = _y;
  x[2] = _z;
}

const vector<samfloat_t>& Node::get_coord() const{ 
  return x; 
}

const samfloat_t* Node::get_cptr_coord() const{ 
  return &(x[0]); 
}

samfloat_t Node::get_x() const{ 
  return x[0]; 
}

samfloat_t Node::get_y() const{ 
  return x[1]; 
}

samfloat_t Node::get_z() const{ 
  return x[2]; 
}

size_t Node::get_size_x() const{ 
  return x.size(); 
}

void Node::set_metric(const vector<samfloat_t>& _metric){
  metric = _metric;
}

void Node::set_metric(const samfloat_t* m, const unsigned mlen){
  metric.resize(mlen);
  for(unsigned i = 0; i<mlen; i++)
    metric[i] = m[i];
}
const vector<samfloat_t>& Node::get_metric() const{ 
  return metric; 
}
const samfloat_t*              Node::get_cptr_metric() const{ 
  return &(metric[0]); 
}
unsigned                     Node::get_size_metric() const{ 
  return metric.size(); 
}

void Node::pack(char *buffer, int& bsize, int& offset) const{
  MPI_Pack(&unn, 1, UNN_T, buffer, bsize, &offset, MPI_COMM_WORLD);
  MPI_Pack(&gnn, 1, GNN_T, buffer, bsize, &offset, MPI_COMM_WORLD);
  MPI_Pack(&flags,  1, MPI_UNSIGNED_CHAR, buffer, bsize, &offset, MPI_COMM_WORLD);
  MPI_Pack(&owner,  2, MPI_UNSIGNED_SHORT, buffer, bsize, &offset, MPI_COMM_WORLD);
  
  { // Pack ifields.
    unsigned nifields = ifields.size();      
    MPI_Pack(&nifields, 1, MPI_UNSIGNED, buffer, bsize, &offset, MPI_COMM_WORLD);
    if(nifields>0)
      MPI_Pack(&(ifields[0]), nifields,MPI_INT, buffer, bsize, &offset, MPI_COMM_WORLD);

  }

  { // Pack fields.
    unsigned nfields = fields.size();      
    MPI_Pack(&nfields, 1, MPI_UNSIGNED, buffer, bsize, &offset, MPI_COMM_WORLD);
    if(nfields>0)
      MPI_Pack(&(fields[0]), nfields, SAMFLOAT, buffer, bsize, &offset, MPI_COMM_WORLD);
  }
  
  { // Pack coordinates.
    unsigned ndim = x.size();
    MPI_Pack(&ndim, 1, MPI_UNSIGNED, buffer, bsize, &offset, MPI_COMM_WORLD);
    if(ndim>0)
      MPI_Pack(&(x[0]), ndim,SAMFLOAT, buffer, bsize, &offset, MPI_COMM_WORLD);
  }

  { // Pack metric.
    unsigned len = metric.size();
    MPI_Pack(&len, 1, MPI_UNSIGNED, buffer, bsize, &offset, MPI_COMM_WORLD);
    if(len>0)
      MPI_Pack(&(metric[0]), len, SAMFLOAT, buffer, bsize, &offset, MPI_COMM_WORLD);
  }

}

void Node::unpack(char *buffer, int& bsize, int& offset){
  MPI_Unpack(buffer, bsize, &offset, &unn, 1, UNN_T, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bsize, &offset, &gnn, 1, GNN_T, MPI_COMM_WORLD );
  MPI_Unpack(buffer, bsize, &offset, &flags, 1, MPI_UNSIGNED, MPI_COMM_WORLD );
  MPI_Unpack(buffer, bsize, &offset, owner, 2, MPI_UNSIGNED_SHORT, MPI_COMM_WORLD );
  
  { // Unpack ifields.
    unsigned nifields;
    MPI_Unpack(buffer, bsize, &offset, &nifields, 1, MPI_UNSIGNED, MPI_COMM_WORLD);
    if(nifields>0){
      ifields.resize(nifields);  
      MPI_Unpack(buffer, bsize, &offset, &(ifields[0]), nifields, MPI_INT, MPI_COMM_WORLD);
    }
  }

  { // Unpack fields.
    unsigned nfields;
    MPI_Unpack(buffer, bsize, &offset, &nfields, 1, MPI_UNSIGNED, MPI_COMM_WORLD);
    if(nfields>0){
      fields.resize( nfields );  
      MPI_Unpack(buffer, bsize, &offset, &(fields[0]), nfields, SAMFLOAT, MPI_COMM_WORLD);
    }
  }

  { // Unpack coordinates.
    unsigned ndim;
    MPI_Unpack(buffer, bsize, &offset, &ndim, 1, MPI_UNSIGNED, MPI_COMM_WORLD);
    if(ndim>0){
      x.resize(ndim);
      MPI_Unpack(buffer, bsize, &offset, &(x[0]), ndim, SAMFLOAT, MPI_COMM_WORLD);
    }
  }

  { // Unpack metric.
    unsigned len;
    MPI_Unpack(buffer, bsize, &offset, &len, 1, MPI_UNSIGNED, MPI_COMM_WORLD);
    if(len>0){
      metric.resize(len);
      MPI_Unpack(buffer, bsize, &offset, &(metric[0]), len, SAMFLOAT, MPI_COMM_WORLD);
    }
  }
     
}

// Return an estimate of the number of bytes required to pack this guy.
unsigned Node::pack_size() const{
  int total=0,s1,s2,s3,s4;
  
  MPI_Pack_size(1, UNN_T, MPI_COMM_WORLD, &s1);
  MPI_Pack_size(1, GNN_T, MPI_COMM_WORLD, &s2);         
  MPI_Pack_size(1, MPI_UNSIGNED_CHAR, MPI_COMM_WORLD, &s3); 
  MPI_Pack_size(2, MPI_UNSIGNED_SHORT, MPI_COMM_WORLD, &s4);
  total = s1 + s2 + s3 + s4;
  {
    unsigned nifields = ifields.size();      
    MPI_Pack_size(1, MPI_UNSIGNED, MPI_COMM_WORLD, &s1);      
    MPI_Pack_size(nifields, MPI_INT, MPI_COMM_WORLD, &s2); 
    total += (s1 + s2);
  }
  {
    unsigned nfields = fields.size();      
    total += MPI_Pack_size(1, MPI_UNSIGNED, MPI_COMM_WORLD, &s1);      
    total += MPI_Pack_size(nfields, SAMFLOAT, MPI_COMM_WORLD, &s2);
    total += (s1 + s2); 
  }
  {
    unsigned ndim = x.size();
    MPI_Pack_size(1, MPI_UNSIGNED, MPI_COMM_WORLD, &s1);      
    MPI_Pack_size(ndim, SAMFLOAT, MPI_COMM_WORLD, &s2);       
    total += (s1 + s2);
  }
  {
    unsigned len = metric.size();
    MPI_Pack_size(1, MPI_UNSIGNED, MPI_COMM_WORLD, &s1);      
    MPI_Pack_size(len, SAMFLOAT, MPI_COMM_WORLD, &s2);
    total += (s1 + s2);    
  }

  return total;
}

